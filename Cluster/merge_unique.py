from ase import io
import os

# --- Configuration ---
# specific files you want to merge
files_to_merge = ['unique_default.traj', 'unique_eclass.traj']
output_combined = 'combined_minima.traj'
output_unique = 'unique_minima.traj'

# 1. Combine them
all_atoms = []
total_count = 0

print(f"Starting merge of: {files_to_merge}")

for f in files_to_merge:
    try:
        if os.path.exists(f):
            # Read all frames from the file
            traj = io.read(f, index=':')
            count = len(traj)
            print(f"  Reading {f}: found {count} structures")
            all_atoms.extend(traj)
            total_count += count
        else:
            print(f"  Warning: File '{f}' not found.")
            
    except Exception as e:
        print(f"  Error reading {f}: {e}")

if not all_atoms:
    print("No structures loaded. Exiting.")
    exit()

# 2. Save the combined file (All trajectories)
io.write(output_combined, all_atoms)

print("-" * 40)
print(f"Successfully merged {total_count} structures into '{output_combined}'")
print("-" * 40)

# 3. Create Unique File (Filter based on Potential Energy)
print("Filtering for unique energies...")

unique_atoms = []
seen_energies = set()

for atoms in all_atoms:
    # Round to 4 decimal places to avoid floating point differences
    e = round(atoms.get_potential_energy(), 4) 
    
    if e not in seen_energies:
        unique_atoms.append(atoms)
        seen_energies.add(e)

# 4. Save the unique file
io.write(output_unique, unique_atoms)

print(f"Filtered down to {len(unique_atoms)} unique structures.")
print(f"Saved unique structures to '{output_unique}'")
print("-" * 40)
