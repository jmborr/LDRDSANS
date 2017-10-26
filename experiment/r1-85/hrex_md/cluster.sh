#!/bin/bash

topology=$1  # Best is to use a PDB file of only the protein
trajectory=$2  # Input MD trajectory file. Frames should contain only the protein

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
nclusters=100  # number of clusters for each trajectory chunk
for i in `seq 0 $maxChunkIndex`;do
    prefix=`printf "%s%03d" $prefixroot $i`
    echo "Kmeans clustering chunk $i. Prefix is $prefix"
    python kmeans_cluster.py ./$prefix $prefix --nclusters $nclusters --centroids #&
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

# Inform of the output
maxprefix=`printf "%s%03d" $prefixroot $maxChunkIndex`
echo " SUMMARY OF RESULTS
* Trajectory was divided into $nchunks for processing. Directories
  p000 to $maxprefix were created for this.
* $nclusters mini-clusters were obtained from each chunk with a k-means
  clustering algorithm using fast_protein_cluster. For each mini-cluster,
  a representative conformation is selected as cluster-centroid.
* All cluster representatives ($chunks x $nclusters) were gathered into
  file centroids.dcd. File centroids.binpos contains the same info but
  in a format amenable to fast_protein_cluster.
* Mini-cluster representatives are clustered with a hiearchical algorithm
  running fast_protein_cluster. These representatives become the leafs of
  the resulting hiearchical cluster. Results are in directory centroids/.
* centroids.pdb: one of the representatives. It's used as a topology file.
* centroids.agglomeration.history: how the leafs where clustered one
  after the other.
* centroids.names: each leaf is one particular conformation of the original
  input trajectory file. Here we log its associated frame number (the first
  conformation has frame number 1).
* centroids.linkage_matrix: hierarchical cluster in a format that can be
  loaded with python package scipy.cluster.hierarchy.
* To visualize the topmost 10 clusters, run
    python plot_linkage_matrix.py ./centroids/centroids.linkage_matrix --p 10
"

