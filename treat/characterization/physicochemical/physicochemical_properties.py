import os
import numpy as np
import pandas as pd  
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import Bio
import matplotlib.pyplot as plt
import seaborn as sns 

# Adapted from code written by Djampa KL Kozlowski.

def read_fasta_file(fpath, case='upper'):
    """
    Read a fasta file and return a dict where keys correspond to the sequences'
    headers and the value are the sequences.
    
    Parameters :
    ------------
    fpath (str) -- input fasta file path
    case (str) -- convert every characters to the desired case : 'upper' for 
    uppercase, 'lower' for lowercase. Else, sequence will not be formated.
    (default : 'upper')
    
    Returns :
    ---------
    (dict) -- a dictionary where keys correspond to the sequences' headers and 
    the value are the sequences.
    """
    if case == 'upper':
        dct_seq = {rec.id:str(rec.seq).upper() for rec in SeqIO.parse(fpath, "fasta")}
    elif case == 'lower':
        dct_seq = {rec.id:str(rec.seq).lower() for rec in SeqIO.parse(fpath, "fasta")}
    else:
        dct_seq = {rec.id:str(rec.seq) for rec in SeqIO.parse(fpath, "fasta")}
    return dct_seq

def write_fasta_file(fpath, sequences, n=60):
    """
    Create a fasta file from a dictionary of sequences
    
    Parameters:
    -----------
    fpath (str) --  output fasta file path
    sequences (dict) -- sequences dictionary where values are the formated 
                        sequences and keys are the headers.
    n (int) --  number of characters per line (default : 60)
    """
    
    with  open(fpath, 'w') as f:
        for header in sequences:
            f.write(
                f">{header}\n{self.split_sequence(sequences[header], n)}\n"
                )

def write_fasta(dictionary, filename):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when caling the function
    """

    import textwrap
    with open(filename, "w") as outfile:
        for key, value in dictionary.items():
            outfile.write(">" + key + "\n")
            outfile.write("\n".join(textwrap.wrap(value, 60)))
            outfile.write("\n")



class PhysicochemicalProperties:
    """
    Extract features/measurments from protein sequences. 
    Rely on Bio.SeqUtils.ProtParam package.
    see https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html. for 
    more informations.
    
    NB : The 'X' (undifined) amino acid are removed from the sequence to compute 
    the properties but are counted in the sequence's length.
    
    Todo : add he hydrophobic dipole moment (see 
    https://gist.github.com/JoaoRodrigues/568c845915aea3efa3578babfd72423c)
    
    
    
    The computed features (columns of the output dataframe) are: 

    - 'id' (str) :                  the sequence id 
    - 'seq_len' (int) :             the sequence length 
    - 'mol_wt' (float) :            the total molecular weight
    - 'instab_idx' (float) :        "Calculate the instability index according 
                                    to Guruprasad et al 1990.Implementation of 
                                    the method of Guruprasad et al. 1990 to test 
                                    a protein for stability Any value above 40 
                                    means the protein is unstable (has a short 
                                    half life). See: Guruprasad K., 
                                    Reddy B.V.B., Pandit M.W. Protein 
                                    Engineering 4:155-161(1990)." 
                                    (Bio.SeqUtils.ProtParam documentation).
    - 'mean_vihinen_flex' (float) : mean value of the flexibility index 
                                    according to Vihinen, 1994 (window's size : 
                                    9)."Protein average flexibility indices are 
                                    inversely correlated to protein stability" 
                                    (https://doi.org/10.1002/prot.340190207)
    - 'mean_kd_hydro' (float) :     mean value of the hydrophobicity profile 
                                    calculated according to the Kyte-Doolittle 
                                    scale (window's size : 6). "The 
                                    Kyte-Doolittle scale is widely used for 
                                    detecting hydrophobic regions in proteins. 
                                    Regions with a positive value are 
                                    hydrophobic. This scale can be used for 
                                    identifying both surface-exposed regions as 
                                    well as transmembrane regions, depending on
                                    the window size used. Short window sizes of 
                                    5-7 generally work well for predicting 
                                    putative surface-exposed regions. Large 
                                    window sizes of 19-21 are well suited for 
                                    finding transmembrane domains if the values 
                                    calculated are above 1.6 [Kyte and 
                                    Doolittle, 1982]. These values should be 
                                    used as a rule of thumb and deviations from 
                                    the rule may occur".
                                    (https://resources.qiagenbioinformatics \
                                    .com/manuals/clcgenomicsworkbench/650/ \
                                    Hydrophobicity_scales.html)
    - 'mean_rose_hydro' (float) :   mean value of the hydrophobicity profile 
                                    calculated according to the Rose scale 
                                    (window's size : 6). "The hydrophobicity 
                                    scale by Rose et al. is correlated to the 
                                    average area of buried amino acids in 
                                    globular proteins [Rose et al., 1985]. 
                                    This results in a scale which is not showing 
                                    the helices of a protein, but rather the 
                                    surface accessibility."
                                    (https://resources.qiagenbioinformatics \
                                    .com/manuals/clcgenomicsworkbench/650/ \
                                    Hydrophobicity_scales.html)
    - 'gravy' (float) :             (grand average of hydropathy) calculated by 
                                    adding the hydropathy value for each residue 
                                    and dividing by the length of the sequence 
                                    (Kyte and Doolittle; 1982) 
                                    (Bio.SeqUtils.ProtParam documentation).
                                    NB : probably redundant with 'mean_kd_hydro' 
                                    and one should be remooved in the futur.
    - 'isoelec_point' (float) :     calculate the isoelectric point. 
                                    (Bio.SeqUtils.ProtParam documentation)
    - 'charge_at_pH' (float) :      calculate the charge of a protein at given 
                                    pH (default : 7.5).
                                    (Bio.SeqUtils.ProtParam documentation) 
                                    "The intracellular pH of living cells is 
                                    strictly controlled in each compartment. 
                                    Under normal conditions, the cytoplasmic pH 
                                    (pHc) and the vacuolar pH (pHv) of typical 
                                    plant cells are maintained at slightly 
                                    alkaline (typically 7.5) and acidic 
                                    (typically 5.5) values, respectively 
                                    (https://doi.org/10.1007/978-3-7091-1254 \
                                    -0_4).
    - 'molar_ext_coef' (float) :    calculates the molar extinction coefficient 
                                    assuming cysteines (reduced) and cystines 
                                    residues (Cys-Cys-bond)
                                    (Bio.SeqUtils.ProtParam documentation)
    - 'tiny' (float) :              cumulative percentage of the smallest amino 
                                    acids (e.g 'A', 'C', 'G', 'S', 'T')
    - 'small' (float) :             cumulative percentage of the small amino 
                                    acids (e.g A', 'C', 'F', 'G', 'I', 'L', 'M', 
                                    'P', 'V', 'W', 'Y').
    - 'aliphatic' (float) :         cumulative percentage of the aliphatic amino 
                                    acids (e.g 'A', 'I', 'L', 'V')
    - 'aromatic' (float) :          cumulative percentage of the aromatic amino 
                                    acids (e.g 'F', 'H', 'W', 'Y')
    - 'non_polar' (float) :         cumulative percentage of the non-polar amino 
                                    acids (e.g 'A', 'C', 'F', 'G', 'I', 'L', 
                                    'M', 'P', 'V', 'W', 'Y')
    - 'polar' (float) :             cumulative percentage of the polar amino 
                                    acids (e.g 'D','E', 'H', 'K', 'N', 'Q', 'R', 
                                    'S', 'T', 'Z')
    - 'charged' (float) :           cumulative percentage of the charged amino 
                                    acids (e.g 'B', 'D', 'E', 'H', 'K', 'R', 
                                    'Z')
    - 'basic' (float) :             cumulative percentage of the basic amino 
                                    acids (e.g 'H', 'K', 'R')
    - 'acidic' (float) :            cumulative percentage of the acidic amino 
                                    acids (e.g 'B', 'D', 'E', 'Z')
    - 'helix' (float) :             cumulative percentage of amino acids which 
                                    tend to be in Helix (e.g. 'V', 'I', 'Y', 
                                    'F', 'W', 'L')
    - 'turn' (float) :              cumulative percentage of amino acids which 
                                    tend to be in turn (e.g 'N', 'P', 'G', 'S')
    - 'sheet' (float) :             cumulative percentage of amino acids which 
                                    tend to be in turn (e.g  'E', 'M', 'A', 'L')
    - 'X_number' (float) :          cumulative percentage of amino acid X 
                                    (for all 20 amino acids)
    """
    
    def __init__(self, fpath, ph=7.5):
        self.lst_prot_prop = []
        self.dict_prop = {
            'tiny':['A', 'C', 'G', 'S', 'T'],
            'small':['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'],
            'aliphatic':['A', 'I', 'L', 'V'],
            'aromatic':['F', 'H', 'W', 'Y'],
            'non_polar':['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'],
            'polar':['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z'],
            'charged':['B', 'D', 'E', 'H', 'K', 'R', 'Z'],
            'basic':['H', 'K', 'R'],
            'acidic':['B', 'D', 'E', 'Z'],
            'A_number':['A'],
            'C_number':['C'],
            'D_number':['D'],
            'E_number':['E'],
            'F_number':['F'],
            'G_number':['G'],
            'H_number':['H'],
            'I_number':['I'],
            'K_number':['K'],
            'L_number':['L'],
            'M_number':['M'],
            'N_number':['N'],
            'P_number':['P'],
            'Q_number':['Q'],
            'R_number':['R'],
            'S_number':['S'],
            'T_number':['T'],
            'V_number':['V'],
            'W_number':['W'],
            'Y_number':['Y'],
            'kyte_doolittle':{'A': 1.8,'C': 2.5,'D': -3.5,'E': -3.5,
                              'F': 2.8,'G': -0.4,'H': -3.2,'I': 4.5,
                              'K': -3.9,'L': 3.8,'M': 1.9,'N': -3.5,
                              'P': -1.6,'Q': -3.5,'R': -4.5,'S': -0.8,
                              'T': -0.7,'V': 4.2,'W': -0.9,'Y': -1.3},
            'rose':{'A': 0.74,'C': 0.91,'D': 0.62,'E': 0.62,'F': 0.88,
                    'G': 0.72,'H': 0.78,'I': 0.88,'K': 0.52,'L': 0.85,
                    'M': 0.85,'N': 0.63,'P': 0.64,'Q': 0.62,'R': 0.64,
                    'S': 0.66,'T': 0.70,'V': 0.86,'W': 0.85,'Y': 0.76}}
        self.fpath = fpath
        self.ph = ph 
    
    def extract_properties(self):
        """
        Read a fasta file and iterate through the sequence to compute 
        several protein properties
        """
        
        self.dct_seqs = read_fasta_file(self.fpath)
        for self.id in self.dct_seqs:
            self.compute_single_seq_properties()
            
            
    def compute_single_seq_properties(self):
        """
        Compute protein properties from the sequence.
        The computed properties are the one described in the class 
        description.
        """
        # remove potential '\t' and '*' char at the end of the sequence      
        seq = str(self.dct_seqs[self.id]).strip().translate({ord('*'): None})
        # compute seq len
        seq_len = len(seq)
        # remove XOU aa from the sequence as they are not necessarly suported 
        # by scales...
        seq = seq.translate({ord(i): None for i in 'XOU'})
        prot_analysis = ProteinAnalysis(seq)
        aa_percent = prot_analysis.get_amino_acids_percent()
        sec_struct_fraction = prot_analysis.secondary_structure_fraction()
        self.lst_prot_prop.append({
            'seq_id':self.id,
            'pp_seq_len':seq_len,
            'pp_mol_wt':prot_analysis.molecular_weight(),
            'pp_instab_idx':prot_analysis.instability_index(),
            'pp_mean_vihinen_flex':np.mean(prot_analysis.flexibility()),
            'pp_mean_kd_hydro':np.mean(prot_analysis.protein_scale( \
                self.dict_prop['kyte_doolittle'], 6, edge=1.0)),
            'pp_mean_rose_hydro':np.mean(prot_analysis.protein_scale( \
                self.dict_prop['rose'], 6, edge=1.0)),
            'pp_gravy':prot_analysis.gravy(),
            'pp_isoelec_point':prot_analysis.isoelectric_point(),
            'pp_charge_at_pH':prot_analysis.charge_at_pH(self.ph),
            'pp_molar_ext_coef':prot_analysis.molar_extinction_coefficient()[1],
            'pp_tiny':np.sum([aa_percent[aa] for aa in self.dict_prop['tiny'] if \
                aa in aa_percent]),
            'pp_small':np.sum([aa_percent[aa] for aa in self.dict_prop['small'] \
                if aa in aa_percent]),
            'pp_aliphatic':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['aliphatic'] if aa in aa_percent]),
            'pp_aromatic':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['aromatic'] if aa in aa_percent]),
            'pp_non_polar':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['non_polar'] if aa in aa_percent]),
            'pp_polar':np.sum([aa_percent[aa] for aa in self.dict_prop['polar'] \
                if aa in aa_percent]),
            'pp_charged':np.sum([aa_percent[aa] for aa in \
                self.dict_prop['charged'] if aa in aa_percent]),
            'pp_basic':np.sum([aa_percent[aa] for aa in self.dict_prop['basic'] \
                if aa in aa_percent]),
            'pp_acidic':np.sum([aa_percent[aa] for aa in self.dict_prop['acidic'] \
                if aa in aa_percent]),
             'pp_A_number':np.sum([aa_percent[aa] for aa in self.dict_prop['A_number'] \
                if aa in aa_percent]),
            'pp_C_number':np.sum([aa_percent[aa] for aa in self.dict_prop['C_number'] \
                if aa in aa_percent]),
            'pp_D_number':np.sum([aa_percent[aa] for aa in self.dict_prop['D_number'] \
                if aa in aa_percent]),
            'pp_E_number':np.sum([aa_percent[aa] for aa in self.dict_prop['E_number'] \
                if aa in aa_percent]),
            'pp_F_number':np.sum([aa_percent[aa] for aa in self.dict_prop['F_number'] \
                if aa in aa_percent]),
            'pp_G_number':np.sum([aa_percent[aa] for aa in self.dict_prop['G_number'] \
                if aa in aa_percent]),
            'pp_H_number':np.sum([aa_percent[aa] for aa in self.dict_prop['H_number'] \
                if aa in aa_percent]),
            'pp_I_number':np.sum([aa_percent[aa] for aa in self.dict_prop['I_number'] \
                if aa in aa_percent]),
            'pp_K_number':np.sum([aa_percent[aa] for aa in self.dict_prop['K_number'] \
                if aa in aa_percent]),
            'pp_L_number':np.sum([aa_percent[aa] for aa in self.dict_prop['L_number'] \
                if aa in aa_percent]),
            'pp_M_number':np.sum([aa_percent[aa] for aa in self.dict_prop['M_number'] \
                if aa in aa_percent]),
            'pp_N_number':np.sum([aa_percent[aa] for aa in self.dict_prop['N_number'] \
                if aa in aa_percent]),
            'pp_P_number':np.sum([aa_percent[aa] for aa in self.dict_prop['P_number'] \
                if aa in aa_percent]),
            'pp_Q_number':np.sum([aa_percent[aa] for aa in self.dict_prop['Q_number'] \
                if aa in aa_percent]),
            'pp_R_number':np.sum([aa_percent[aa] for aa in self.dict_prop['R_number'] \
                if aa in aa_percent]),
            'pp_S_number':np.sum([aa_percent[aa] for aa in self.dict_prop['S_number'] \
                if aa in aa_percent]),
            'pp_T_number':np.sum([aa_percent[aa] for aa in self.dict_prop['T_number'] \
                if aa in aa_percent]),
            'pp_V_number':np.sum([aa_percent[aa] for aa in self.dict_prop['V_number'] \
                if aa in aa_percent]),
            'pp_W_number':np.sum([aa_percent[aa] for aa in self.dict_prop['W_number'] \
                if aa in aa_percent]),
            'pp_Y_number':np.sum([aa_percent[aa] for aa in self.dict_prop['Y_number'] \
                if aa in aa_percent]),
            'pp_KEN_number':np.sum([aa_percent[aa] for aa in self.dict_prop['KEN_number'] \
                if aa in aa_percent]),
            'pp_helix':sec_struct_fraction[0],
            'pp_turn':sec_struct_fraction[1],
            'pp_sheet':sec_struct_fraction[2]
        })
        

    def export_as_dataframe(self):
        """
        Return a pandas dataframe from a list of dict.
        """
        
        return pd.DataFrame(self.lst_prot_prop)
    
    
def compute_physicochemical_properties(fpath, outdir, ph=7.5):
    """
    Execute the PhysicochemicalProperties class, compute peptide sequence 
    properties and store the results in a file. 
    
    Parameters :
    ------------
    fpath (str)   --  input fasta file
    outdir (str)    --  output directory path
    ph (float)  --  pH value to compute protein charge
    
    Returns :
    ---------
    (str)   --  output file path
    """
    output_path = os.path.join(outdir, f'physicochemprop_all.tsv')
    # init the class 
    physicochem_prop = PhysicochemicalProperties(fpath, ph=ph)
    # compute protein sequences properties for each sequence in the file. 
    physicochem_prop.extract_properties()
    # store the information (dataframe) in a file
    physicochem_prop.export_as_dataframe()\
        .to_csv(
            output_path, 
            sep='\t', 
            index=False, 
            header=True)
    return output_path


def plot_feature_importance(importance,names,model_type):

    #Create arrays from feature importance and feature names
    feature_importance = np.array(importance)
    feature_names = np.array(names)

    #Create a DataFrame using a Dictionary
    data={'feature_names':feature_names,'feature_importance':feature_importance}
    fi_df = pd.DataFrame(data)

    #Sort the DataFrame in order decreasing feature importance
    fi_df.sort_values(by=['feature_importance'], ascending=False,inplace=True)
    
    #Define size of bar plot
    plt.figure(figsize=(15,10))
    #Plot Seaborn bar chart
    sns.barplot(x=fi_df['feature_importance'], y=fi_df['feature_names'])
    #Add chart labels
    plt.title(model_type + 'FEATURE IMPORTANCE')
    plt.xlabel('FEATURE IMPORTANCE')
    plt.ylabel('FEATURE NAMES')

if __name__ == "__main__":
    compute_physicochemical_properties("/path/to/fasta", "/path/to/outdir")
