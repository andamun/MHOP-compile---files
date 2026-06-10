from ase import io
import os
import glob

# 1. Find all trajectory files in run subdirectories
# This looks for run0/minima.traj, run1/minima.traj, etc.
files = glob.glob('run*/minima.traj')
files.sort()

print(f"Found {len(files)} trajectory files to merge.")

# 2. Combine them
all_atoms = []
total_count = 0

for f in files:
    try:
        # Read all frames from the file
        traj = io.read(f, index=':')
        count = len(traj)
        print(f"  Reading {f}: found {count} structures")
        all_atoms.extend(traj)
        total_count += count
    except Exception as e:
        print(f"  Error reading {f}: {e}")

# 3. Save the combined file
output_file = 'combined_minima.traj'
io.write(output_file, all_atoms)

print("-" * 40)
print(f"Successfully merged {total_count} structures into '{output_file}'")
print("-" * 40)

# Optional: Deduplicate (remove identical energies)
# This is a simple check based on potential energy
unique_atoms = []
seen_energies = set()

for atoms in all_atoms:
    e = round(atoms.get_potential_energy(), 4)  # Round to 4 decimal places
    if e not in seen_energies:
        unique_atoms.append(atoms)
        seen_energies.add(e)

io.write('unique_minima.traj', unique_atoms)
print(f"Filtered down to {len(unique_atoms)} unique structures based on energy.")
print(f"Saved unique structures to 'unique_minima.traj'")
