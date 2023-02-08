import get_targets
import screening
import similiarity


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
