import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu, chi2_contingency, fisher_exact, shapiro, anderson
import matplotlib.pyplot as plt
import seaborn as sns

# Retrieve all information

# Orphan lists
with open("/path/to/orphans_8melo_verified.txt") as f:
    verified_orphans = set(line.strip() for line in f)
    
with open("/path/to/orphans_8melo_nonverified.txt") as f:
    nonverified_orphans = set(line.strip() for line in f)
    
# Diverged or de novo
with open("/path/to/diverged_orphans.txt") as f:
    diverged_orphans = set(line.strip() for line in f)
    
with open("/path/to/denovo_orphans.txt") as f:
    denovo_orphans = set(line.strip() for line in f)
    
# Orthogroup info
og_info = pd.read_csv("/path/to/sequence_to_orthogroup.tsv", sep="\t", names=["Sequence_id", "Orthogroup"])
genecount_df = pd.read_csv("/path/to/Orthogroup.GeneCount.tsv", sep="\t")

# Age info

root_df = pd.read_csv("/path/to/orthogroup_age.tsv", sep="\t")

# Homolog count
hom_info = pd.read_csv("/path/to/homologcount.tsv", sep="\t", names=["Sequence_id", "Number_of_homologs"])

# Transcription
with open("/path/to/expressed_contigs.txt") as f:
    expressed_contigs = set(line.strip() for line in f)

# Translation 

with open("/path/to/translated_contigs.txt") as f:
    translated_contigs = set(line.strip() for line in f)
    

# Physicochemical features
physicochemprop = pd.read_csv("/path/to/physicochemprop_all.tsv", sep="\t")

# Localization
deeploc = pd.read_csv("/path/to/deeploc_results.csv")

# Peptide signals
signalp = pd.read_csv("/path/to/signalp_results.txt", sep="\t")
signalp["ID"] = signalp["ID"].str.strip().str.split().str[0]
signalp['SignalP'] = signalp['Prediction'].apply(lambda x: 1 if x == 'SP' else 0)

# InterPro domains
with open("/path/to/interproscan_ipr_domains.txt") as f:
    interproscan_sequences = set(line.strip() for line in f)
# Add InterPro column with 1 if present in file, else 0
physicochemprop['InterPro'] = physicochemprop['seq_id'].apply(lambda x: 1 if x in interproscan_sequences else 0)

main_df = pd.DataFrame()
main_df["Sequence_id"] = physicochemprop["seq_id"]

main_df["Orphan_status"] = main_df["Sequence_id"].apply(
    lambda x: 1 if x in verified_orphans or x in nonverified_orphans else 0
)


orphan_mask = main_df["Orphan_status"] == 1

verified_mask = main_df["Sequence_id"].isin(verified_orphans)
main_df.loc[orphan_mask & verified_mask, "Transcription"] = 1

nonverified_expressed_mask = main_df["Sequence_id"].isin(nonverified_orphans & expressed_contigs)
main_df.loc[orphan_mask & nonverified_expressed_mask, "Transcription"] = 1

nonverified_unexpressed_mask = main_df["Sequence_id"].isin(nonverified_orphans - expressed_contigs)
main_df.loc[orphan_mask & nonverified_unexpressed_mask, "Transcription"] = 0


translated_mask = main_df["Sequence_id"].isin(translated_contigs)
untranslated_mask = ~main_df["Sequence_id"].isin(translated_contigs)
main_df.loc[orphan_mask & translated_mask, "Translation"] = 1
main_df.loc[orphan_mask & untranslated_mask, "Translation"] = 0

main_df = main_df.merge(og_info, on="Sequence_id", how="left")
main_df = main_df.merge(hom_info, on="Sequence_id", how="left")

main_df = main_df.merge(signalp[["ID", "SignalP"]], left_on="Sequence_id", right_on="ID", how="left")
main_df.drop(columns=["ID"], inplace=True)

main_df = main_df.merge(physicochemprop, left_on="Sequence_id", right_on="seq_id", how="left")
main_df.drop(columns=["seq_id"], inplace=True)

main_df = main_df.merge(deeploc, left_on="Sequence_id", right_on="Protein_ID", how="left")
main_df.drop(columns=["Protein_ID"], inplace=True)

main_df = main_df.merge(root_df, left_on="Orthogroup", right_on="OG", how="left")
main_df.drop(columns=["OG"], inplace=True)

species_columns = genecount_df.columns[1:] 
genecount_df[species_columns] = genecount_df[species_columns].apply(pd.to_numeric, errors='coerce')
all_melo_ogs = genecount_df[genecount_df[species_columns].gt(0).all(axis=1)]["OG"].tolist()
main_df.loc[main_df["Orthogroup"].isin(all_melo_ogs) & main_df["Root"].isna(), "Root"] = "test_all_melo"

only_1_mask = main_df["Root"] == "only_1"
main_df.loc[only_1_mask, ["Root", "Orthogroup"]] = np.nan

conditions = [
    main_df["Sequence_id"].isin(denovo_orphans),
    main_df["Sequence_id"].isin(diverged_orphans)
]
choices = ["denovo", "diverged"]

main_df["Emergence"] = np.select(conditions, choices, default=np.nan)
main_df["Emergence"] = main_df["Emergence"].replace("nan", np.nan)

main_df = main_df.drop_duplicates(subset="Sequence_id", keep="first")
main_df = main_df.drop(columns=["Localizations", "Signals"])

main_df = main_df[['Sequence_id', 'Orphan_status', 'Orthogroup', 'Number_of_homologs', 'Root', 'Emergence',
 'Transcription', 'Translation', 'SignalP', 'InterPro',
 'pp_seq_len', 'pp_mol_wt', 'pp_instab_idx', 'pp_gravy',
 'pp_isoelec_point', 'pp_charge_at_pH',
 'pp_A_number', 'pp_C_number', 'pp_D_number', 'pp_E_number', 'pp_F_number',
 'pp_G_number', 'pp_H_number', 'pp_I_number', 'pp_K_number', 'pp_L_number',
 'pp_M_number', 'pp_N_number', 'pp_P_number', 'pp_Q_number', 'pp_R_number',
 'pp_S_number', 'pp_T_number', 'pp_V_number', 'pp_W_number', 'pp_Y_number',
 'pp_tiny', 'pp_small', 'pp_aliphatic', 'pp_aromatic', 'pp_non_polar',
 'pp_polar', 'pp_charged', 'pp_basic', 'pp_acidic',
 'pp_mean_vihinen_flex', 'pp_mean_kd_hydro', 'pp_mean_rose_hydro',
 'pp_helix', 'pp_turn', 'pp_sheet',
 'Cytoplasm', 'Nucleus', 'Extracellular', 'Cell membrane', 'Mitochondrion',
 'Plastid', 'Endoplasmic reticulum', 'Lysosome/Vacuole',
 'Golgi apparatus', 'Peroxisome']
]


main_df.to_csv("feature_table.tsv", sep="\t")


