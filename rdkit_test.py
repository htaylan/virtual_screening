import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import DataStructs
import seaborn as sns
from rdkit.Chem import AllChem

def get_targets(target_name):
    """
    :param target_name: The name of the inhibitor target
    :return: This function gets targets from ChemBL database
    and return as a dataframe
    """
    target = new_client.target
    target_query = target.search('kras')
    targets = pd.DataFrame.from_dict(target_query)
    return targets

def select_target(n):
    """
    :return: Selects the target
    """
    selected_target = targets.target_chembl_id[0]
    return selected_target

def get_inhibitors(selected_target):
    """
    :param selected_target: ChemBl Id of selected target
    :return: The inhibitors for selected target
    from ChemBl database as a dataframe
    """
    activity = new_client.activity
    res = activity.filter(target_chembl_id=selected_target)
    df = pd.DataFrame.from_dict(res)
    return df

def select_potency_type(df,pot_type):
    """
    :parameter pot_type: type of potency
    :return dataframe based on selected potency type and drop NaN :
    """
    valid_types = 'IC50', 'EC50', 'Kd'
    if pot_type not in valid_types:
        raise ValueError("Type should be one of valid types")

    df.drop(df[(df.type != pot_type)].index, axis=0, inplace=True)
    df = df[df.standard_value.notna()]
    df = df.reset_index(drop=True)
    return df

def select_nec_columns(df,liste):
    """
    :return Select only necessary columns:
    """
    df = df[liste]
    return df

def calculate_bioactivity(df,in_limit,ac_limit):
    """
    :return Calculate bioactivity (ac_limit), inactivity (in_limit)
     and intermediate activity
     based on standard_values. Add those standard values to the dataframe
     as a new column called bioactivity_class :
    """

    df.standard_value = df.standard_value.astype("float")
    bioactivity_class = []
    for i in df['standard_value']:
        if float(i) >= in_limit:
            bioactivity_class.append('inactive')
        elif float(i) <= ac_limit:
            bioactivity_class.append('active')
        else:
            bioactivity_class.append('intermediate')

    print(len(bioactivity_class))
    df = pd.concat([df, pd.Series(bioactivity_class)], axis=1)
    df.rename({0: "bioactivity_class"}, axis=1, inplace=True)
    return df


def lipinski(smiles, verbose=False):
    """
    :param smiles: The smiles strings of the molecules
    :param verbose: To not explicitly write descriptors to the screen
    :return: Molecular Weight, LogP, H donors, and H acceptors descriptors
    of each molecule as np.array then generates a dataframe
    """
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = []

    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Descriptors.NumHDonors(mol)
        desc_NumHAcceptor = Descriptors.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptor])

        baseData.append(row)

    columnnames = ['MW', 'LogP', 'NumHDonors', 'NumHAcceptors']
    descriptors = pd.DataFrame(data=baseData, columns=columnnames)

    return descriptors

def combine_dataframes(df1,df2):
    """
    :param df1: First dataframe
    :param df2: Second dataframe
    :return: Combine two dataframe column-wise
    """
    df_combined = pd.concat([df1, df2], axis=1)
    return df_combined

def subset_df_for_lipinski(df):
    """
    :param df: The dataframe to subset
    :return: Subset data frame based on Lipinski rules
    and bioactivity
    """
    df_sub_lipinski = df[
        (df['MW'] <= 500) &
        (df['LogP'] <= 5) &
        (df['NumHDonors'] <= 5) &
        (df['NumHAcceptors'] <= 10) &
        (df['bioactivity_class'] == 'active')]

    df_sub_lipinski = df_sub_lipinski.reset_index(drop=True)
    return df_sub_lipinski


def get_scaffold(smiles, verbose=False):
    """
    :param smiles: The smiles string of the molecules
    :param verbose: Keep function quiet
    :return: The Murckoscaffold of a molecule
    """
    scaffold = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        core = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold.append(Chem.MolToSmiles(core, isomericSmiles=True))

    return scaffold

def get_min_std_val_idx(df):
    """
    :param df: Name of the dataframe
    :return: the ID of the molecule with the lowest standard value
    """
    min_idx = df['standard_value'].idxmin()
    return min_idx

def get_min_std_val_smiles(df,min_idx):
    """
    :param df: Name of the dataframe
    :param min_idx: the ID of the molecule with the lowest standard value
    :return: The smiles of the molecule with the lowest standard value
    """
    min_sml = df['standard_value'][min_idx]
    return min_sml

def calc_Tanimoto_sim(df,min_sml):
    """
    :param df:  Name of the dataframe
    :param min_sml: The smiles of the molecule with the lowest standard value
    :return: Return the list of Tanimoto similarity values compared to
    reference molecule
    """
    ref = Chem.MolFromSmiles(min_sml)
    fp_ref = Chem.RDKFingerprint(ref)

    Tan_list = []
    for mol in df['canonical_smiles']:
        mol_comp = Chem.MolFromSmiles(mol)
        fp_comp = Chem.RDKFingerprint(mol_comp)
        Tanimoto = DataStructs.TanimotoSimilarity(fp_ref, fp_comp)
        Tan_list.append(Tanimoto)

    return Tan_list

def get_n_mol_with_smallest_std_val(df, n):
    """
    :param df: Name of the dataframe
    :param n: Select n smallest standard value
    :return: The smiles of the n molecule with smallest standard value
    """
    smallest_n = list(df.nsmallest(n, 'standard_value').index)
    smiles_n = []
    for i in smallest_n:
        mol = df['canonical_smiles'][i]
        smiles_n.append(mol)
    return smiles_n

def add_H_to_smiles(smile_list):
    """
    :param smile_list: The list that contains the smiles of the molecules
    :return: Return the objects of the molecules with H
    """
    mols = [Chem.MolFromSmiles(smi) for smi in smile_list]
    mols_addH = [Chem.AddHs(mol) for mol in mols]
    return mols_addH


targets = get_targets('kras')
selected_target = select_target(0)
print(selected_target)
df = get_inhibitors(selected_target)
df_ic50 = select_potency_type(df,"IC50")
print(df.shape)
df_nec = select_nec_columns(df_ic50,['molecule_chembl_id','canonical_smiles','standard_value'])
df_bio = calculate_bioactivity(df_nec,10000,1000)
df_lipinski = lipinski(df_bio.canonical_smiles)
df_combined = combine_dataframes(df_bio,df_lipinski)
df_sub_lipinski = subset_df_for_lipinski(df_combined)
print(df.shape)
#
# min_idx = get_min_std_val_idx(df_sub_lipinski)
# min_sml = get_min_std_val_smiles(df_sub_lipinski,min_idx)
# print(min_sml)