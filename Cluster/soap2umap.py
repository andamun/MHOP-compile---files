import os
import pandas as pd
import numpy as np
from dscribe.descriptors import SOAP
import matplotlib.pyplot as plt
import umap
import hdbscan
from ase.io import read,write

# input
FRAMES_DIR = "frames"
CSV_FILE = "mh_energies_labeled.csv"


# SOAP Parameters
SOAP_PARAMS = {
        "species": ["Pd", "C", "H", "N", "O", "Cl"],
    "r_cut": 6.0,
    "n_max": 6,
    "l_max": 8,
    "sigma": 0.5,
    "periodic": True,
    "sparse": False,
    "average": "off"
}
SLAB_ELEMENT = "Pd"
# UMAP Parameters
UMAP_PARAMS = {
    "n_neighbors": 15,
    "min_dist": 0.1,
    "n_components": 2,  # reducing to 2D for plotting
    "metric": "euclidean",
    "random_state": 42
}

# HDBSCAN Parameters
HDBSCAN_PARAMS = {
    "min_cluster_size": 30,
    "min_samples": 5,
    "metric": "euclidean",
    "cluster_selection_method": "eom" # 'eom' (Excess of Mass) prefers big clusters
}

def load_data(csv_path, frames_dir):
    """Loads CSV and reads existing .extxyz files."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)
    required = ["hop", "energy_eV", "label"]
    if not all(col in df.columns for col in required):
        raise ValueError(f"CSV missing columns. Needs: {required}")

    structures = []
    valid_indices = []

    print(f"Loading frames from {frames_dir}...")
    for idx, row in df.iterrows():
        fname = f"frame_{int(row['hop']):03d}.extxyz"
        fpath = os.path.join(frames_dir, fname)

        if os.path.exists(fpath):
            structures.append(read(fpath))
            valid_indices.append(idx)

    return df.loc[valid_indices].reset_index(drop=True), structures

def compute_soap_vectors(structures, soap_params, ignore_element):
    """Computes SOAP and averages only the molecule atoms."""
    print("Computing SOAP descriptors...")
    soap = SOAP(**soap_params)
    vectors = []

    for atoms in structures:
        features = soap.create(atoms)
        symbols = atoms.get_chemical_symbols()
        # Find indices of atoms that are NOT the slab (ignore_element)
        mol_idx = [i for i, s in enumerate(symbols) if s != ignore_element]

        if mol_idx:
            vectors.append(np.mean(features[mol_idx], axis=0))
        else:
            vectors.append(np.mean(features, axis=0))

    return np.array(vectors)
def plot_clusters_title_only(df, n_clusters, filename="umap_clusters_count.png"):
    """Creates the cluster plot with the count in the Title."""
    plt.figure(figsize=(10, 8))

    # 1. Plot Noise (-1) in grey
    noise = df[df["cluster"] == -1]
    if not noise.empty:
        plt.scatter(noise["umap_x"], noise["umap_y"], c="lightgrey", s=20, label="Noise")

    # 2. Plot Clusters
    clustered = df[df["cluster"] != -1]
    if not clustered.empty:
        scatter = plt.scatter(
            clustered["umap_x"], clustered["umap_y"],
            c=clustered["cluster"], cmap='tab10', s=60, alpha=0.9, edgecolor='k', linewidth=0.5
        )
        plt.colorbar(scatter, label="Cluster ID")

    # [MODIFIED] Title now displays the exact count
    plt.title(f"UMAP\n{n_clusters} Distinct Groups", fontsize=16)
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"[OK] Saved cluster plot to '{filename}'")

def plot_energy(df, filename="umap_energy.png"):
    """Creates the plot colored by Energy (eV)"""
    plt.figure(figsize=(10, 8))

    scatter = plt.scatter(
        df["umap_x"], df["umap_y"],
        c=df["energy_eV"], cmap='plasma_r',
        s=60, alpha=0.9, edgecolor='none'
    )

    cbar = plt.colorbar(scatter)
    cbar.set_label("Energy (eV)")

    plt.title(f"UMAP Energy Landscape")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"[OK] Saved energy plot to '{filename}'")

def main():
    # 1. Load Data
    df, structures = load_data(CSV_FILE, FRAMES_DIR)
    if len(df) == 0: return

    # 2. Compute SOAP
    X = compute_soap_vectors(structures, SOAP_PARAMS, SLAB_ELEMENT)

    # 3. Run UMAP
    print("Running UMAP...")
    reducer = umap.UMAP(**UMAP_PARAMS)
    embedding = reducer.fit_transform(X)
    df["umap_x"] = embedding[:, 0]
    df["umap_y"] = embedding[:, 1]

    # 4. Run HDBSCAN
    print("Running HDBSCAN...")
    clusterer = hdbscan.HDBSCAN(**HDBSCAN_PARAMS)
    df["cluster"] = clusterer.fit_predict(embedding)

    # Calculate number of clusters (ignoring noise -1)
    unique_labels = set(df["cluster"])
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)
    print(f"Algorithm detected {n_clusters} clusters.")

    # 5. Generate Plots
    print("-" * 30)
    # Pass the n_clusters count to the plotting function
    plot_clusters_title_only(df, n_clusters)
    plot_energy(df)
    print("-" * 30)
    
    # Save Data
    df.to_csv("final_results.csv", index=False)

if __name__ == "__main__":
    main()
