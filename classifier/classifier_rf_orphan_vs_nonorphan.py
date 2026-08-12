import numpy as np
import pandas as pd
import sys
import os
import json
import pickle
from joblib import dump

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.metrics import f1_score, precision_score, recall_score, roc_auc_score, confusion_matrix, classification_report

# Retrieve the dataframe and prepare it 

df = pd.read_csv("feature_table.tsv", sep="\t")
df = df[
    ((df["Orphan_status"] == 1) & (df["Orthogroup"].notna())) |
    (df["Orphan_status"] == 0)
].copy()
df = df.drop(columns=[
    "Species", "Emergence", "InterPro", "SignalP", "Number_of_homologs", 
    "Orthogroup", "Root", "Transcription", "Proteomics"
])
df = df.set_index("Sequence_id")

X = df.drop(columns=["Orphan_status"])
y = df["Orphan_status"].values

# 60/20/20 train-validation-test split

X_80, X_test, y_80, y_test = train_test_split(
    X, y, test_size=0.20, stratify=y, random_state=42
)

X_train, X_val, y_train, y_val = train_test_split(
    X_80, y_80, test_size=0.25, stratify=y_80, random_state=42
)

# GridSearch

param_grid = {
    "n_estimators": [200, 400, 600],
    "max_features": ["sqrt", None],
    "min_samples_leaf": [1, 2],
    "min_samples_split": [2, 10],
    "max_depth": [None, 20]
}

rf = RandomForestClassifier(
    class_weight="balanced",
    n_jobs=32,
    random_state=42
)

gs = GridSearchCV(
    rf,
    param_grid=param_grid,
    scoring="f1",
    cv=StratifiedKFold(n_splits=3, shuffle=True, random_state=42),
    n_jobs=1,
    verbose=1,
    refit=True
)

gs.fit(X_train, y_train)

print("Best params:", gs.best_params_)
print("Best CV F1:", gs.best_score_)

y_val_pred = gs.best_estimator_.predict(X_val)
print("Validation F1:", f1_score(y_val, y_val_pred))

# Final training with best params

X_final = pd.concat([X_train, X_val])
y_final = np.concatenate([y_train, y_val])

final_model = RandomForestClassifier(
    **gs.best_params_,
    class_weight="balanced",
    n_jobs=-1,
    random_state=42
)

final_model.fit(X_final, y_final)

# Test performance

y_test_pred = final_model.predict(X_test)

test_f1 = f1_score(y_test, y_test_pred)
test_prec = precision_score(y_test, y_test_pred)
test_rec = recall_score(y_test, y_test_pred)

print("Scores on test")
print(f"F1:        {test_f1:.4f}")
print(f"Precision: {test_prec:.4f}")
print(f"Recall:    {test_rec:.4f}")

print("\nConfusion matrix:")
print(confusion_matrix(y_test, y_test_pred))
print("\nReport:\n", classification_report(y_test, y_test_pred, digits=4))

# Save model and outputs

pd.DataFrame(gs.cv_results_).to_csv("gridsearch_cv_results.csv", index=False)

with open("best_params.json", "w") as f:
    json.dump(gs.best_params_, f, indent=2)

with open("rf_orphan_best.pkl", "wb") as f:
    pickle.dump(final_model, f)

dump(final_model, "rf_orphan_best.joblib")

print("\nModel saved as rf_orphan_best.pkl and rf_orphan_best.joblib")


