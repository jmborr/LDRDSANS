#############################
#  Run in text mode
#  vmd -dispdev text -eofexit -e this_script.tcl topology_file trajectory_file
#############################

set HA [atomselect top "protein and not element H"]
animate write pdb HA.pdb beg 0 end 0 sel $HA top
animate write binpos HA.binpos beg 0 end -1 sel $HA waitfor all top
quit


