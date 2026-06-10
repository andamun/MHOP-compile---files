import os
import glob
import pandas as pd
import shutil
import re

def get_cluster_number(folder_name):
    """Extracts the integer from the folder name for proper numerical sorting."""
    match = re.search(r'\d+', folder_name)
    return int(match.group()) if match else -1

def main():
    # Find all directories that start with "cluster"
    cluster_dirs = glob.glob("cluster*")
    
    if not cluster_dirs:
        print("No cluster directories found in the current folder.")
        return

    # Sort the directories numerically (0, 1, 2, ..., 10, 11)
    cluster_dirs = sorted(cluster_dirs, key=get_cluster_number)

    # List to store the results for the final CSV
    summary_data = []

    for cluster_dir in cluster_dirs:
        if not os.path.isdir(cluster_dir):
            continue
            
        print(f"--- Processing {cluster_dir} ---")

        # 1. Find and read the CSV file inside the cluster folder
        csv_files = glob.glob(os.path.join(cluster_dir, "final_results_*.csv"))
        if not csv_files:
            print(f"[WARN] No CSV file found in {cluster_dir}. Skipping.")
            continue
            
        csv_path = csv_files[0]
        df = pd.read_csv(csv_path)
        
        if df.empty or "energy_eV" not in df.columns or "hop" not in df.columns:
            print(f"[WARN] Missing required columns or empty CSV in {csv_path}. Skipping.")
            continue

        # 2. Find the average energy and the closest frame
        mean_energy = df["energy_eV"].mean()
        df["diff_from_mean"] = (df["energy_eV"] - mean_energy).abs()
        
        closest_idx = df["diff_from_mean"].idxmin()
        closest_row = df.loc[closest_idx]
        
        best_hop = int(closest_row["hop"])
        best_energy = closest_row["energy_eV"]
        
        print(f"Cluster Average Energy : {mean_energy:.5f} eV")
        print(f"Closest Frame Selected : frame_{best_hop:03d}.extxyz (Energy: {best_energy:.5f} eV)")

        # Record the data for the summary CSV
        summary_data.append({
            "cluster_dir": cluster_dir,
            "hop_avg_selected": best_hop,
            "eV": best_energy
        })

        # 3. Create the 'select' subdirectory
        select_dir = os.path.join(cluster_dir, "select")
        os.makedirs(select_dir, exist_ok=True)
        
        # 4. Copy the frame into the 'select' folder
        frame_filename = f"frame_{best_hop:03d}.extxyz"
        src_path = os.path.join(cluster_dir, frame_filename)
        dst_path = os.path.join(select_dir, frame_filename)
        
        if os.path.exists(src_path):
            shutil.copy2(src_path, dst_path)
            print(f"[OK] Copied {frame_filename} to {select_dir}/")
        else:
            print(f"[ERROR] Source file not found: {src_path}")
            
    # 5. Generate the Summary CSV
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        
        # Ensure the columns are ordered exactly how you want them
        summary_df = summary_df[["cluster_dir", "hop_avg_selected", "eV"]]
        
        # Save to file
        output_csv = "summary_average_frames.csv"
        summary_df.to_csv(output_csv, index=False)
        print(f"\n[SUCCESS] Wrote sorted summary of all selected frames to '{output_csv}'")
    
    print("\nAll clusters processed successfully.")

if __name__ == "__main__":
    main()


