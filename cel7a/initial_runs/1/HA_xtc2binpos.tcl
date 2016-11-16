#############################
#  Run in text mode
#  vmd -dispdev text -eofexit -e this_script.tcl
#############################

# File names
set inpdb "HA.pdb"
set inxtc   "HA.xtc"
set outbinpos  "HA.binpos"

# INPUT PDB
mol new $inpdb

# INPUT trajectory
set begframe 0
set endframe -1
animate read xtc $inxtc beg $begframe end $endframe waitfor all


# OUTPUT trajectory
set begframe 1
set endframe -1
animate write binpos $outbinpos beg $begframe end $endframe waitfor all

quit