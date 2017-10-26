#############################
#  Run in text mode
#  vmd -dispdev text -eofexit -e this_script.tcl topology_PDB_file trajectory_file
#############################

#beg=1 because the PDB file adds one extra conformation
mol new sh4.pdb
animate read xtc traj_comp0_noPBC.xtc beg 0 end -1 waitfor all
animate write pdb traj_comp0_noPBC.pdb beg 1 end -1 waitfor all top
quit

