"""
This example code allow users to screen KRAS targets and their potential inhibitors
Then screen the inhibitors based on IC50 values and bioactivities.
df_sub_lipinski dataframe contains 33 potential candidates, you can play with them
with the functions provided (Scaffold search, similarity scores etc.) in similarity.py as you wish.
"""

from get_targets import *
from screening import *
from similarity import *


def screen_kras_inhibitors(verbose=False):
    chembl_client = ChemBLClient()
    targets = chembl_client.get_targets('kras')
    selected_target = chembl_client.select_target(0)
    inhibitors = chembl_client.get_inhibitors()

    bio = BioactivityData()
    df = bio.select_potency_type(inhibitors, 'IC50')
    df = bio.select_nec_columns(df, ['molecule_chembl_id', 'canonical_smiles', 'standard_value'])
    df = bio.calculate_bioactivity(df, 10000, 1000)
    descriptors = bio.lipinski(df['canonical_smiles'].tolist())
    df_combined = bio.combine_dataframes(df, descriptors)
    df_sub_lipinski = bio.subset_df_for_lipinski(df_combined)
    print(df_sub_lipinski)

    similarity_score = SimilarityScore(df_sub_lipinski, df_sub_lipinski['canonical_smiles'], 5, verbose)
    scaffolds = similarity_score.get_scaffold()
    min_idx = similarity_score.get_min_std_val_idx()
    min_smiles = similarity_score.get_min_std_val_smiles(min_idx)
    smallest_n_smiles = similarity_score.get_n_mol_with_smallest_std_val()
    mols_with_H = similarity_score.add_H_to_smiles(smallest_n_smiles)

    df_sub_lipinski['scaffold'] = scaffolds

    if verbose:
        print("Selected target:", selected_target)
        print("Potential inhibitors data frame (subset based on Lipinski's rule):\n", df_sub_lipinski)

    return df_sub_lipinski, mols_with_H


if __name__ == '__main__':
    df_sub_lipinski, mols_with_H = screen_kras_inhibitors(verbose=True)
