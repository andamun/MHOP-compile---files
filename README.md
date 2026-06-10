# MHOP-compile---files

licensed and written by -- Anchalee CHAICHANAVONGSAROJ

For easy implication, please feel free to write your own automation script. I didnot uplode mine, in case you may want to edit some procedures depends on your style.

To select the presentative structure of sampling using MHOP, I divide the procedures into 3 main steps.

1. Find the limited temperature you gonna set for running MHOP in MD - check T folder.

  For this part, you should check and set the constraints of the system including:
    (1) 'lashing' adsorbate to not let it fly out of the surface
    (2) bonding each atom in adsorbate for not make it break to each other
    (3) fix two bottom layers of the slab
  For (1) and (2), hookean constraints are used. Please make sure it appropriate with your system. In this folder, I place the one appropriate for PCNB on Pd.

  Run the MD with High Temp (normally 5000 K) for 25 PS, using NVE ensemble - for this case I use velocity verlet. 
  
  Then plot the Temp, check whether the system will drop and run around at which Temperature use this at Tmax for MHOP running.

2. Performing MHOP using code in MHOP folder. What I have done are
  (1) Create several of sub directories e.g., run0, run1, .., run10
  (2) copy initial structure file, both of python code in MHOP folder, in each dir
  (3) sub the job in each dir for 6,12,or 24 hours (up to the DOF of the system if less DOF should increase the time), if running it on Wildfly using - ncpu = 4 cores or 8  cores in each dir is fine, but have to keep monitor the structure to prevent it running on the bad node I normally use count*.sh files in these folder to do so. For NSCC I used 16 core.

 3. For compile data and clustering in Cluster folder.
  (1) I merge the data in one .traj file using merge_minima.py - this screen the replicate structure that could obtain from running in several dir, pick the unique_minima.traj one to go into next step
  (2) separate into frame using separated_frame.py, and fine the energy using find_energy.py
  (3) run soap2umap.py
  (4) run organize_clusters.py
  (5) pick representative each group using select_avg.py

And then, DFT optimize each structures obtain the energy and pray for finding more stable configuration.
