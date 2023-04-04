import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator
from typing import List, Tuple


def read_molecule_csv(file_path: str) -> pd.DataFrame:
    """Read the CSV file containing the molecule data and return a pandas dataframe"""
    df = pd.read_csv(file_path)
    return df


def create_molecule_array(df: pd.DataFrame) -> List[Tuple[Chem.rdchem.Mol, str]]:
    """Create molecules from SMILES and store them in a list"""
    molecules = [(Chem.MolFromSmiles(smiles), chembl_id) for _, chembl_id, smiles in
                 df[["molecule_chembl_id", "canonical_smiles"]].itertuples()]
    return molecules


def create_fingerprints(molecules: List[Tuple[Chem.rdchem.Mol, str]]) -> List[DataStructs.cDataStructs.ExplicitBitVect]:
    """Create fingerprints for all molecules"""
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
    fingerprints = [rdkit_gen.GetFingerprint(mol) for mol, _ in molecules]
    return fingerprints


def calculate_tanimoto_distance(fp_list: List[DataStructs.cDataStructs.ExplicitBitVect]) -> List[float]:
    """Calculate distance matrix for a list of fingerprints using the Tanimoto coefficient"""
    distance_matrix = []
    for i in range(1, len(fp_list)):
        # Compare the current fingerprint against all the previous ones in the list
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        distance_matrix.extend([1 - x for x in similarities])
    return distance_matrix


def cluster_fingerprints(fingerprints: List[DataStructs.cDataStructs.ExplicitBitVect], cutoff: float = 0.2) -> List[List[int]]:
    """Cluster fingerprints using the Butina algorithm"""
    distance_matrix = calculate_tanimoto_distance(fingerprints)
    clusters = Butina.ClusterData(distance_matrix, len(fingerprints), cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters


def calculate_intra_similarity(fps_clusters: List[List[DataStructs.cDataStructs.ExplicitBitVect]]) -> List[List[float]]:
    """Calculate the Tanimoto similarity for all pairs of fingerprints in each cluster"""
    intra_similarity = []
    for cluster in fps_clusters:
        similarities = [1 - x for x in calculate_tanimoto_distance(cluster)]
        intra_similarity.append(similarities)
    return intra_similarity


def plot_intra_cluster_similarity(intra_sim: List[List[float]], cutoff: float = 0.4) -> None:
    """Plot the intra-cluster Tanimoto similarity"""
    fig, ax = plt.subplots(figsize=(10, 5))
    indices = list(range(1, len(intra_sim) + 1))
    ax.set_xlabel("Cluster index")
    ax.set_ylabel("Similarity")
    ax.set_xticks(indices)
    ax.set_xticklabels(indices)
    ax.set_yticks(np.arange(0.4, 1.0, 0.1))
    ax.set_title("Intra-cluster Tanimoto similarity", fontsize=13)
    r = ax.violinplot(intra_sim, indices, showmeans=True, showmedians=True, showextrema=False)
    r["cmeans"].set_color("red")
    plt.show()
   
def visualize_clusters(compounds, clusters):
    """Draw molecules in clusters"""
    Draw.MolsToGridImage(
        [compounds[i][0] for i in clusters[1][:10]],
        legends=[compounds[i][1] for i in clusters[1][:10]],
        
    
    
    
    
    
    
    
    
    
    
    
    
    
import time
import random
from pathlib import Path

import pandas as pd
import numpy
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdFingerprintGenerator

compound_df = pd.read_csv('clean_molecule_library.csv')

# Create molecules from SMILES and store in array
compounds = []

for _, chembl_id, smiles in compound_df[["molecule_chembl_id", "canonical_smiles"]].itertuples():
    compounds.append((Chem.MolFromSmiles(smiles), chembl_id))
compounds[:5]

# Create fingerprints for all molecules
rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5)
fingerprints = [rdkit_gen.GetFingerprint(mol) for mol, idx in compounds]

print("Number of compounds converted:", len(fingerprints))
print("Fingerprint length per compound:", len(fingerprints[0]))

def tanimoto_distance_matrix(fp_list):
    """Calculate distance matrix for fingerprint list"""
    dissimilarity_matrix = []
    
    for i in range(1, len(fp_list)):
        # Compare the current fingerprint against all the previous ones in the list
        similarities = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        dissimilarity_matrix.extend([1 - x for x in similarities])
    return dissimilarity_matrix
  
def cluster_fingerprints(fingerprints, cutoff=0.2):
    """Cluster fingerprints
    Parameters:
        fingerprints
        cutoff: threshold for the clustering
    """
    # Calculate Tanimoto distance matrix
    distance_matrix = tanimoto_distance_matrix(fingerprints)
    # Now cluster the data with the implemented Butina algorithm:
    clusters = Butina.ClusterData(distance_matrix, len(fingerprints), cutoff, isDistData=True)
    clusters = sorted(clusters, key=len, reverse=True)
    return clusters
  
  
def intra_tanimoto(fps_clusters):
    """Function to compute Tanimoto similarity for all pairs of fingerprints in each cluster"""
    intra_similarity = []
    # Calculate intra similarity per cluster
    for cluster in fps_clusters:
        # Tanimoto distance matrix function converted to similarity matrix (1-distance)
        intra_similarity.append([1 - x for x in tanimoto_distance_matrix(cluster)])
    return intra_similarity
  
  # Recompute fingerprints for 10 first clusters
mol_fps_per_cluster = []
for cluster in clusters[:10]:
    mol_fps_per_cluster.append([rdkit_gen.GetFingerprint(compounds[i][0]) for i in cluster])

# Compute intra-cluster similarity
intra_sim = intra_tanimoto(mol_fps_per_cluster)

fig, ax = plt.subplots(figsize=(10, 5))
indices = list(range(1,11))
ax.set_xlabel("Cluster index")
ax.set_ylabel("Similarity")
ax.set_xticks(indices)
ax.set_xticklabels(indices)
ax.set_yticks(numpy.arange(0.4, 1.0, 0.1))
ax.set_title("Intra-cluster Tanimoto similarity", fontsize=13)
r = ax.violinplot(intra_sim, indices, showmeans=True, showmedians=True, showextrema=False)
r["cmeans"].set_color("red")
# mean=red, median=blue

 #Visualize it 
cutoff = 0.4
clusters = cluster_fingerprints(fingerprints, cutoff=cutoff)

# Plot the size of the clusters - save plot
fig, ax = plt.subplots(figsize=(15, 4))
ax.set_xlabel("Cluster index")
ax.set_ylabel("# molecules")
ax.bar(range(1, len(clusters) + 1), [len(c) for c in clusters])
ax.set_title(f"Threshold: {cutoff:3.1f}")
fig.savefig(
    f"cluster_dist_cutoff_{cutoff:4.2f}.png",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)

print(
    f"Number of clusters: {len(clusters)} from {len(compounds)} molecules at distance cut-off {cutoff:.2f}"
)
print("Number of molecules in largest cluster:", len(clusters[0]))
print(
    f"Similarity between two random points in same cluster: {DataStructs.TanimotoSimilarity(fingerprints[clusters[0][0]], fingerprints[clusters[0][1]]):.2f}"
)
print(
    f"Similarity between two random points in different cluster: {DataStructs.TanimotoSimilarity(fingerprints[clusters[0][0]], fingerprints[clusters[1][0]]):.2f}"
)

# Draw molecules
Draw.MolsToGridImage(
    [compounds[i][0] for i in clusters[1][:10]],
    legends=[compounds[i][1] for i in clusters[1][:10]],
    molsPerRow=5,
)

Draw.MolsToGridImage(
    [compounds[clusters[i][0]][0] for i in range(10)],
    legends=[compounds[clusters[i][0]][1] for i in range(10)],
    molsPerRow=5,
)



