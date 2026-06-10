from mace.calculators import mace_mp
from ase import units
from ase.io.trajectory import Trajectory
from ase.md.langevin import Langevin #NVT ensemble for equilibration
from ase.md.verlet import VelocityVerlet #NVE for running
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import MDLogger
from ase.io import read,write
from ase.visualize import view
import warnings
from ase.neighborlist import neighbor_list,natural_cutoffs
from ase.constraints import FixAtoms, Hookean

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import numpy as np

# Suppress all warnings
warnings.filterwarnings("ignore")

# def for detecting the cluster of multiple pcnb attach on teh surface
def get_molecular_clusters(atoms, adsorbate_indices):
    """
    Separates adsorbate atoms into distinct molecules using graph connectivity.
    Returns a list of lists: [[indices of mol A], [indices of mol B], ...]
    """
    # 1. Build a connectivity matrix only for adsorbate atoms
    cutoffs = natural_cutoffs(atoms, mult=1.05)
    i_list, j_list = neighbor_list('ij', atoms, cutoffs)

    # Filter bonds to only include those within the adsorbate list
    bonds = []
    map_idx = {original: i for i, original in enumerate(adsorbate_indices)}

    for k in range(len(i_list)):
        a, b = i_list[k], j_list[k]
        if a in adsorbate_indices and b in adsorbate_indices:
            bonds.append([map_idx[a], map_idx[b]])

    if not bonds:
        return [[i] for i in adsorbate_indices] # No bonds, treat all as separate

    # 2. Use Scipy to find connected components (clusters)
    data = np.ones(len(bonds), dtype=int)
    rows = [b[0] for b in bonds]
    cols = [b[1] for b in bonds]
    N = len(adsorbate_indices)
    graph = csr_matrix((data, (rows, cols)), shape=(N, N))

    n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)

    # 3. Group original indices by label
    clusters = [[] for _ in range(n_components)]
    for local_idx, label in enumerate(labels):
        original_idx = adsorbate_indices[local_idx]
        clusters[label].append(original_idx)

    return clusters

# Define the Hookean constraint among the molecule
def apply_pcnb_constraints_universal(atoms,slab_atom,length_atom_slab):
    """
    Applies constraints to MULTIPLE PCNB molecules automatically.
    """
    # --- 1. Define PCNB Parameters ---
    bond_params = {
        ('C', 'Cl'): {'rt': 2.50, 'k': 3.0},
        ('C', 'N'):  {'rt': 2.15, 'k': 5.0},
        ('N', 'O'):  {'rt': 1.75, 'k': 10.0},
        ('C', 'C'):  {'rt': 2.00, 'k': 5.0},
        ('C', 'H'):  {'rt': 1.59, 'k': 7.0}
    }

    # Identify the surface height (Preparing to add the constraints )
    pd_z_values = [atom.position[2] for atom in atoms if atom.symbol == slab_atom]
    top_pd_z = max(pd_z_values)

    # Identify all non-surface atoms
    # (Assuming surface is Pd, change if using Cu)
    adsorbate_indices = [a.index for a in atoms if a.symbol != slab_atom]

    # Separate into distinct molecules (e.g., Mol A and Mol B)
    molecular_clusters = get_molecular_clusters(atoms, adsorbate_indices)

    constraints = []

    print(f"Detected {len(molecular_clusters)} distinct molecule(s).")

    # --- Loop over EACH molecule found ---
    for mol_id, cluster_indices in enumerate(molecular_clusters):
        print(f"  Processing Molecule {mol_id+1} (Indices: {cluster_indices})")

        # A. Apply Internal Bond Constraints (The Springs)
        # We only check bonds within this specific cluster
        subset = atoms[cluster_indices]
        cutoffs = natural_cutoffs(subset, mult=1.05)
        i_sub, j_sub = neighbor_list('ij', subset, cutoffs)

        processed_pairs = set()

        for k in range(len(i_sub)):
            # Map subset indices back to global indices
            idx_i = cluster_indices[i_sub[k]]
            idx_j = cluster_indices[j_sub[k]]

            pair_id = tuple(sorted((idx_i, idx_j)))
            if pair_id in processed_pairs: continue
            processed_pairs.add(pair_id)

            sym_i = atoms[idx_i].symbol
            sym_j = atoms[idx_j].symbol
            bond_type = tuple(sorted((sym_i, sym_j)))

            if bond_type in bond_params:
                p = bond_params[bond_type]
                constraints.append(Hookean(a1=int(idx_i), a2=int(idx_j), rt=p['rt'], k=p['k']))

        # B. Apply Anchor Constraint (The Leash to Surface)
        # Goal: Find a Ring Carbon. In PCNB, all carbons are ring carbons.
        # We pick the one closest to the surface (lowest Z) for stability.

        cluster_atoms = [(idx, atoms[idx]) for idx in cluster_indices]
        carbons = [item for item in cluster_atoms if item[1].symbol == 'C']

        if carbons:
            # Sort carbons by Z height and pick the lowest one
            carbons.sort(key=lambda x: x[1].position[2])
            anchor_idx = carbons[0][0]
            anchor_pos = atoms[anchor_idx].position

            # This prevents THIS specific molecule from flying away
            constraints.append(Hookean(a1=int(anchor_idx), a2=(0.,0.,1.,-(top_pd_z+length_atom_slab)), k=10.0))
            print(f"    -> Anchored Molecule {mol_id+1} via Carbon {anchor_idx} (z={anchor_pos[2]:.2f})")
        else:
            print(f"    ! Warning: Molecule {mol_id+1} has no Carbons! Cannot anchor.")

    return constraints

def verify_constraints(atoms, constraints):
    print("\n--- CONSTRAINT VERIFICATION REPORT ---")

    # 1. Re-run cluster detection to see what the computer thinks exist
    adsorbate_indices = [a.index for a in atoms if a.symbol != 'Pd']
    clusters = get_molecular_clusters(atoms, adsorbate_indices)

    # Create a map: Atom Index -> Molecule ID
    atom_to_mol_id = {}
    for mol_id, indices in enumerate(clusters):
        print(f"Molecule {mol_id+1} found: {len(indices)} atoms. Indices: {indices}")
        for idx in indices:
            atom_to_mol_id[idx] = mol_id + 1

    print("-" * 30)

    # 2. Check every constraint
    errors = 0
    for k, constr in enumerate(constraints):
        # We only care about Hookean constraints for this check
        if not isinstance(constr, Hookean):
            continue

        # Hookean constraints store indices in .indices (list) or .index (single)
        # For a bond (2 atoms), it is usually in 'constr.indices' or 'constr.a1/a2' depending on version
        # ASE Hookean typically uses a1 and a2 as attributes

        try:
            # Check if it is a "Two Atom" constraint (Bond)
            if hasattr(constr, 'a1') and hasattr(constr, 'a2'):
                # Some versions of ASE store a2 as a vector (for plane constraints)
                # We need to check if a2 is an integer (atom index) or a tuple (coordinate)

                if isinstance(constr.a2, (int, np.integer)):
                    idx1 = constr.a1
                    idx2 = constr.a2

                    mol1 = atom_to_mol_id.get(idx1, "Surface/Unknown")
                    mol2 = atom_to_mol_id.get(idx2, "Surface/Unknown")

                    if mol1 != mol2:
                        print(f" [CRITICAL ERROR] Constraint {k} connects Molecule {mol1} (Atom {idx1}) to Molecule {mol2} (Atom {idx2})!")
                        errors += 1
                    # else:
                    #    print(f" [OK] Constraint {k}: Mol {mol1} internal bond.")

                elif isinstance(constr.a2, (tuple, list, np.ndarray)):
                    # This is likely the Anchor constraint (Atom to Point)
                    mol1 = atom_to_mol_id.get(constr.a1, "Surface/Unknown")
                    print(f" [INFO] Anchor Constraint detected on Molecule {mol1} (Atom {constr.a1}).")

        except Exception as e:
            print(f"Could not parse constraint {k}: {e}")

    if errors == 0:
        print("\nSUCCESS: All constraints are chemically valid. No cross-linking detected.")
    else:
        print(f"\nFAILURE: Found {errors} cross-linked constraints!")

#--- main ----#
# import the adsorbed system 
slab = read('Thermalized_struct.traj')
#slab = read('opt_struct.traj')

# Identify the slab_atom and length_atom_slab
slab_atom='Pd'
length_atom_slab = 7 #Angstorm


# Adding constraints

# Create constraint to make sure the molecule not fly out from the surface
# Define the constraint for fixing the bottom two layers
cu_indices = [atom.index for atom in slab if atom.symbol == slab_atom]
mask = [i for i in range(0,int(len(cu_indices)/2),1)]
fixlayers = FixAtoms(mask=mask)

# Create the list for containing the fixlayers first
constraints = [fixlayers]

# Generate the Hookean constrains making sure that the molecule will not explode
glycerol_constraints = apply_pcnb_constraints_universal(slab,slab_atom,length_atom_slab)
constraints.extend(glycerol_constraints)

# Apply all constaints
slab.set_constraint(constraints)

# --- RUN THE CHECK ---
verify_constraints(slab, slab.constraints)

# Identify model
calc = mace_mp(model='medium',
               dispersion=True,
               damping='zero',
               default_dtype="float32")
slab.calc = calc

required_temp = 4000.0 # Kelvin
MaxwellBoltzmannDistribution(atoms=slab,temp=required_temp*units.kB,force_temp=True) # To ensure that the temperature correctly distributed

# Identify thermostat
dyn = VelocityVerlet(atoms=slab,timestep=0.2*units.fs)
dyn.attach(MDLogger(dyn,slab,'md.log' #log file name
                    ,header=True,stress=False,
                    peratom=False, # make it print out as total energy
                    mode='w')
                    ,interval = 10) # interval timestep

# save the position of atom every 10 interval
traj = Trajectory('VelocityVerlet_md.traj', 'w', slab)
dyn.attach(traj.write, interval=50)

dyn.run(25000) # running for 25 PS
