import os
import pandas as pd
import shutil

# ==========================================
# 1. CONFIGURATION
# ==========================================
INPUT_CSV = "final_results.csv"
SOURCE_FRAMES_DIR = "frames"
OUTPUT_BASE_DIR = "."  # Where to create the cluster folders (current dir)

# Set to True to MOVE files (deletes from source). 
# Set to False to COPY files (safer, keeps original frames).
MOVE_FILES = False 

def main():
    # 1. Load Data
    if not os.path.exists(INPUT_CSV):
        print(f"Error: {INPUT_CSV} not found.")
        return
        
    df = pd.read_csv(INPUT_CSV)
    
    # Ensure necessary columns exist
    if "cluster" not in df.columns or "energy_eV" not in df.columns:
        print("Error: CSV must contain 'cluster' and 'energy_eV' columns.")
        return

    # 2. Process each cluster group
    # Get unique cluster IDs (sorted)
    unique_clusters = sorted(df["cluster"].unique())
    
    print(f"Found {len(unique_clusters)} clusters (including noise if present). Starting organization...")
    print("-" * 50)

    for cluster_id in unique_clusters:
        # Skip noise if you don't want to organize it (usually cluster -1)
        # if cluster_id == -1: continue 

        # Filter data for this cluster
        cluster_df = df[df["cluster"] == cluster_id].copy()
        
        # 3. Sort by Energy (Lowest/Most Stable first)
        cluster_df = cluster_df.sort_values(by="energy_eV", ascending=True)
        
        # Define Directory Name (e.g., "cluster0", "cluster1", "cluster_noise")
        if cluster_id == -1:
            dir_name = "cluster_noise"
        else:
            dir_name = f"cluster{cluster_id}"
            
        target_dir = os.path.join(OUTPUT_BASE_DIR, dir_name)
        
        # Create Directory
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
            
        # 4. Save the sorted CSV for this cluster
        csv_name = f"final_results_cluster{cluster_id}.csv"
        csv_path = os.path.join(target_dir, csv_name)
        cluster_df.to_csv(csv_path, index=False)
        
        # 5. Move/Copy Frames
        files_processed = 0
        for _, row in cluster_df.iterrows():
            # Assuming 'hop' or 'frame_index' holds the file number
            # Adjust column name if you used 'frame_index' in the previous step
            if "frame_index" in row:
                idx = int(row["frame_index"])
            elif "hop" in row:
                idx = int(row["hop"])
            else:
                print("Error: Could not find 'hop' or 'frame_index' to identify files.")
                return

            filename = f"frame_{idx:03d}.extxyz"
            src_path = os.path.join(SOURCE_FRAMES_DIR, filename)
            dst_path = os.path.join(target_dir, filename)
            
            if os.path.exists(src_path):
                if MOVE_FILES:
                    shutil.move(src_path, dst_path)
                else:
                    shutil.copy2(src_path, dst_path)
                files_processed += 1
            else:
                print(f"[WARN] Frame not found: {src_path}")

        action = "Moved" if MOVE_FILES else "Copied"
        print(f"[{dir_name}] Saved CSV with {len(cluster_df)} rows. {action} {files_processed} frame files.")

    print("-" * 50)
    print("Organization complete.")

if __name__ == "__main__":
    main()
