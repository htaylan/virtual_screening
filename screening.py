import pandas as pd
import numpy as np
import rdkit.Chem as Chem
from rdkit.Chem import Descriptors

class BioactivityData:
    def __init__(self):
        self.valid_types = ('IC50', 'EC50', 'Kd')
        self.bioactivity_class = []

    def select_potency_type(self, df, pot_type):
        if pot_type not in self.valid_types:
            raise ValueError("Type should be one of valid types")
        df.drop(df[(df.type != pot_type)].index, axis=0, inplace=True)
        df = df[df.standard_value.notna()]
        df = df.reset_index(drop=True)
        return df

    def select_nec_columns(self, df, liste):
        df = df[liste]
        return df

    def calculate_bioactivity(self, df, in_limit, ac_limit):
        df.standard_value = df.standard_value.astype("float")
        for i in df['standard_value']:
            if float(i) >= in_limit:
                self.bioactivity_class.append('inactive')
            elif float(i) <= ac_limit:
                self.bioactivity_class.append('active')
            else:
                self.bioactivity_class.append('intermediate')
        df = pd.concat([df, pd.Series(self.bioactivity_class)], axis=1)
        df.rename({0: "bioactivity_class"}, axis=1, inplace=True)
        return df

    def lipinski(self, smiles, verbose=False):
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

    def combine_dataframes(self, df1, df2):
        df_combined = pd.concat([df1, df2], axis=1)
        return df_combined

    def subset_df_for_lipinski(self, df):
        df_sub_lipinski = df[
                (df['MW'] <= 500) &
                (df['LogP'] <= 5) &
                (df['NumHDonors'] <= 5) &
                (df['NumHAcceptors'] <= 10) &
                (df['bioactivity_class'] == 'active')]
        df_sub_lipinski = df_sub_lipinski.reset_index(drop=True)
        return df_sub_lipinski


