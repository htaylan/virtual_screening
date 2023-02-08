import pandas as pd
import numpy as np
import rdkit.Chem as Chem
from rdkit import Descriptors

class DataProcessor:
    def __init__(self, df, pot_type, nec_columns, in_limit, ac_limit, smiles):
        self.df = df
        self.pot_type = pot_type
        self.nec_columns = nec_columns
        self.in_limit = in_limit
        self.ac_limit = ac_limit
        self.smiles = smiles

    def select_potency_type(self):
        """
        :return dataframe based on selected potency type and drop NaN :
        """
        valid_types = 'IC50', 'EC50', 'Kd'
        if self.pot_type not in valid_types:
            raise ValueError("Type should be one of valid types")

        self.df.drop(self.df[(self.df.type != self.pot_type)].index, axis=0, inplace=True)
        self.df = self.df[self.df.standard_value.notna()]
        self.df = self.df.reset_index(drop=True)
        return self.df

    def select_nec_columns(self):
        """
        :return Select only necessary columns:
        """
        self.df = self.df[self.nec_columns]
        return self.df

    def calculate_bioactivity(self):
        """
        :return Calculate bioactivity (ac_limit), inactivity (in_limit)
        and intermediate activity
        based on standard_values. Add those standard values to the dataframe
        as a new column called bioactivity_class :
        """
        self.df.standard_value = self.df.standard_value.astype("float")
        bioactivity_class = []
        for i in self.df['standard_value']:
            if float(i) >= self.in_limit:
                bioactivity_class.append('inactive')
            elif float(i) <= self.ac_limit:
                bioactivity_class.append('active')
            else:
                bioactivity_class.append('intermediate')

        self.df = pd.concat([self.df, pd.Series(bioactivity_class)], axis=1)
        self.df.rename({0: "bioactivity_class"}, axis=1, inplace=True)
        return self.df
    
    def lipinski(self):
        for elem in self.smiles:
            mol = Chem.MolFromSmiles(elem)
            self.moldata.append(mol)
        
        for mol in self.moldata:
            desc_MolWt = Descriptors.MolWt(mol)
            desc_MolLogP = Descriptors.MolLogP(mol)
            desc_NumHDonors = Descriptors.NumHDonors(mol)
            desc_NumHAcceptor = Descriptors.NumHAcceptors(mol)

            row = np.array([desc_MolWt,
                            desc_MolLogP,
                            desc_NumHDonors,
                            desc_NumHAcceptor])

            self.baseData.append(row)
        self.descriptors = pd.DataFrame(data=self.baseData, columns=self.columnnames)
        
    def combine_dataframes(self, df2):
        df_combined = pd.concat([self.descriptors, df2], axis=1)
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
