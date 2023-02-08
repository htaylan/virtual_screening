import rdkit.Chem as Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdmolops import SanitizeFlags
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs

class SimilarityScore:
    def __init__(self, df, smiles, n, verbose=False):
        self.df = df
        self.smiles = smiles
        self.n = n
        self.verbose = verbose
    
    def get_scaffold(self):
        """
        :return: The Murckoscaffold of a molecule
        """
        scaffold = []
        for smile in self.smiles:
            mol = Chem.MolFromSmiles(smile)
            core = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold.append(Chem.MolToSmiles(core, isomericSmiles=True))

        return scaffold

    def get_min_std_val_idx(self):
        """
        :return: the ID of the molecule with the lowest standard value
        """
        min_idx = self.df['standard_value'].idxmin()
        return min_idx

    def get_min_std_val_smiles(self,min_idx):
        """
        :param min_idx: the ID of the molecule with the lowest standard value
        :return: The smiles of the molecule with the lowest standard value
        """
        min_sml = self.df['standard_value'][min_idx]
        return min_sml

    def calc_Tanimoto_sim(self, min_sml):
        """
        :return: Return the list of Tanimoto similarity values compared to
        reference molecule
        """
        ref = Chem.MolFromSmiles(min_sml)
        fp_ref = Chem.RDKFingerprint(ref)

        Tan_list = []
        for mol in self.df['canonical_smiles']:
            mol_comp = Chem.MolFromSmiles(mol)
            fp_comp = Chem.RDKFingerprint(mol_comp)
            Tanimoto = DataStructs.TanimotoSimilarity(fp_ref, fp_comp)
            Tan_list.append(Tanimoto)

        return Tan_list

    def get_n_mol_with_smallest_std_val(self):
        """
        :return: The smiles of the n molecule with smallest standard value
        """
        smallest_n = list(self.df.nsmallest(self.n, 'standard_value').index)
        smiles_n = []
        for i in smallest_n:
            mol = self.df['canonical_smiles'][i]
            smiles_n.append(mol)
        return smiles_n

    def add_H_to_smiles(self, smile_list):
        """
        :return: Return the objects of the molecules with H
        """
        mols = [Chem.MolFromSmiles(smi) for smi in smile_list]
        mols_addH = [Chem.AddHs(mol) for mol in mols]
        return mols_addH

