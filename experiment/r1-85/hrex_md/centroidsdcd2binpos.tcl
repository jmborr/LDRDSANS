#############################
#  Run in text mode
#  vmd -dispdev text -eofexit -e this_script.tcl topology_PDB_file trajectory_file
#############################

#beg=1 because the PDB file adds one extra conformation
animate write binpos centroids.binpos beg 1 end -1 waitfor all top
quit


