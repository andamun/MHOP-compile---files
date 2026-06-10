import os
from ase.io import read, write

# Load all frames
frames = read('unique_minima.traj', index=':')

base_dir = "frames"

os.makedirs(base_dir,exist_ok=True)

for i, atoms in enumerate(frames):
    filename=f"frame_{i:03d}.extxyz"
    output_path=os.path.join(base_dir,filename)
    write(output_path, atoms,format="extxyz")
