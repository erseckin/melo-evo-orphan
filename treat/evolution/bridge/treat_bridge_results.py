import os
import pandas as pd

"""Treat bridge output files to get the number of emerged orthogroups at each node"""

def map_root(filename, root):
    if 'mgra' in filename or 'mchi' in filename:
        mappings = {1: 'only_1', 2: 'mgra_mchi', 3: 'all_melo'}
    elif 'mhap' in filename:
        mappings = {1: 'only_1', 2: 'mhap_clade1', 3: 'all_melo'}
    elif 'ment' in filename:
        mappings = {1: 'only_1', 2: 'clade1', 3: 'mhap_clade1', 4: 'all_melo'}
    elif 'mare' in filename:
        mappings = {1: 'only_1', 2: 'clade1_wo_ment', 3: 'clade1', 4: 'mhap_clade1', 5: 'all_melo'}
    elif 'mjav' in filename:
        mappings = {1: 'only_1', 2: 'clade1_wo_ment_mare', 3: 'clade1_wo_ment', 4: 'clade1', 5: 'mhap_clade1', 6: 'all_melo'}
    elif 'minc' in filename or 'mluc' in filename:
        mappings = {1: 'only_1', 2: 'clade1_wo_ment_mare_mjav', 3: 'clade1_wo_ment_mare', 4: 'clade1_wo_ment', 5: 'clade1', 6: 'mhap_clade1', 7: 'all_melo'}
    else:
        # Default mapping if no specific case matches
        mappings = {1: 'only_1'}
    return mappings.get(root, 'unknown')

def to_align_mapping(root):
    mapping = {
        'only_1': 'not_interested',
        'clade1_wo_ment_mare_mjav': 'mjav',
        'clade1_wo_ment_mare': 'mare',
        'clade1_wo_ment': 'ment',
        'clade1': 'mhap',
        'mhap_clade1': 'mgra_mchi',
        'mgra_mchi' : 'ppen',
        'all_melo': 'ppen'
    }
    return mapping.get(root, 'unknown')

def process_files():
    combined_df = pd.DataFrame()
    for file in os.listdir():
        if file.endswith(".tsv"):
            df = pd.read_csv(file, sep="\t", engine='python')
            df['Root'] = df['Root'].apply(lambda x: map_root(file, x))
            df['To_align'] = df['Root'].apply(to_align_mapping)
            df = df[['OG', 'Root', 'To_align']]
            combined_df = pd.concat([combined_df, df])

    # The order of preference for 'To_align'
    preference_order = ['ppen', 'mgra_mchi', 'mhap', 'ment', 'mare', 'mjav', 'not_interested']
    combined_df['Preference'] = combined_df['To_align'].apply(lambda x: preference_order.index(x))
    
    # Keeping the highest preference entry for each OG
    combined_df.sort_values(by=['OG', 'Preference'], ascending=[True, True], inplace=True)
    combined_df.drop_duplicates(subset=['OG'], keep='first', inplace=True)
    combined_df.drop(columns=['Preference'], inplace=True)  # Clean up the DataFrame
    
    # Save the combined file
    combined_df.to_csv("root_ident_bridge.txt", index=False, sep='\t')
    print("Combined and processed file saved as: root_ident_bridge.txt")

process_files()

