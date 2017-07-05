#!/bin/bash

topology=$1
trajectory=$2

################################
#  Auxiliary scripts:
#    partition_trajectory.py
#    kmeans_cluster.py
#    centroidsdcd2binpos.tcl
#    fpch2scph.py
#    plot_linkage_matrix.py
################################
source /etc/profile.d/modules.sh
module load vmd
module load fast_protein_cluster

#Partition trajectory in 10 chunks retaining only the protein
nchunks=10
maxChunkIndex=`echo $[$nchunks-1]`
prefixroot="p"
python partition_trajectory.py $topology $trajectory $nchunks --prefix $prefixroot


# Kmeans cluster of each chunk. Save cluster centroids to DCD file
# (Parallel execution)
centroidTrajs=""  # store file names for the centroid trajectories
centroidNames=""  # store the names of the centroids
for i in `seq 0 $maxChunkIndex`;do
    prefix=`printf "%s%03d" $prefixroot $i`
    echo "Kmeans clustering chunk $i. Prefix is $prefix"
    python kmeans_cluster.py ./$prefix $prefix --nclusters 100 --centroids #&
    centroidTrajs+=" ./$prefix/${prefix}.cluster.centroids.dcd"
    centroidNames+=" ./$prefix/${prefix}.cluster.centroids.names"
done
wait  # wait for all child processes
echo "Finished Kmeans clustering"
sleep 2s

# Gather the centroids into a single trajectory and names file
catdcd -o centroids.dcd $centroidTrajs  # gather all centroid DCD trajectories
cat $centroidNames > centroids.names    # gather all centroid names

# create .pdb and .binpos files
/bin/cp ./$prefix/${prefix}.pdb centroids.pdb  # representative PDB of the protein
vmd -dispdev text -eofexit -e centroidsdcd2binpos.tcl centroids.pdb centroids.dcd  # create centroids.binpos

# hierarchical clustering
fast_protein_cluster -i centroids -o centroids --binary_coords --cluster_rmsd --nclusters 2 --nthreads 8 --hcomplete --sse3

# Create linkage matrix from clustering output
python fpch2scph.py centroids.agglomeration.history centroids.linkage_matrix

# Move all centroids files to centroids directory
mkdir -p ./centroids && /bin/mv centroids.* ./centroids/

# Plot a dendogram with the last 10 merged clusters
python plot_linkage_matrix.py ./centroids/centroids.linkage_matrix --p 10 &



