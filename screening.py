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
