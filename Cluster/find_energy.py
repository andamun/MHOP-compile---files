import pandas as pd
from ase.io import read
from mace.calculators import mace_mp

# The traj file
TRAJ_FILE = "unique_minima.traj"  # <--- CHANGE THIS to your actual filename

# 2. Output CSV name (must match what your SOAP script expects)
OUTPUT_CSV = "mh_energies_labeled.csv"

def main():
    traj = read(TRAJ_FILE, index=':')
    print("Loading the model")
    mace_calc=mace_mp(model='medium'
                ,dispersion=True
                ,default_dtype='float32')

    # 3. Compute Energies
    results = []
    print(f"\n{'Frame':<6} | {'Total Energy (eV)':<18}")
    print("-" * 30)

    for i, atoms in enumerate(traj):
        try:
            atoms.calc = mace_calc
            total_energy = atoms.get_potential_energy()

            # Store the data
            results.append({
                "hop": i,               # Frame number
                "energy_eV": total_energy,
                "label": "selected"    # Default label for your SOAP script
                })
            # Print progress every 10 frames
            if i % 10 == 0:
                print(f"{i:<6} | {total_energy:.5f}")
        except Exception as e:
            print(f"Error on frame {i}: {e}")

    # 4. Save to CSV
        if results:
            df = pd.DataFrame(results)
            df.to_csv(OUTPUT_CSV, index=False)
            print("\n" + "="*50)
            print(f"[OK] Saved energies for {len(df)} frames to '{OUTPUT_CSV}'")
            print(f"You can now run your SOAP/UMAP script.")
        else:
            print("No results calculated.")

if __name__ == "__main__":
    main()
