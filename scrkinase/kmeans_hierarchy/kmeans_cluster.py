# -*- coding: utf-8 -*-

"""Run kmeans algorithm of fast_protein_cluster.
"""

import os
import sys
import re
import subprocess
import argparse
from partition_trajectory import inputHAfiles
import mdtraj
import numpy as np
from pdb import set_trace as tr

def get_centroids(cluster_stats_file):
    """Scan {{prefix}}.cluster.stats and fetch the centroid name
    :returns names of the centroids
    """
    pattern = "Cluster:  \d+ of size \d+ with centroid (\w+) with cluster density score"
    with open(cluster_stats_file, 'r') as handle:
        data = handle.read()
        return re.findall(pattern, data)

def save_centroids(centroid_names, chunkset, extension="dcd"):
    """Extract the centroid conformations and save to file. Also save names file
    Assumes {{prefix}}.binpos and {{prefix}}.pdb, and {{prefix}}.names exist
    """
    # Find frame indices for the centroids
    frame_names_file = os.path.join(chunkset.directory, "{}.names".format(chunkset.prefix))
    with open(frame_names_file,'r') as handle:
        frame_names = handle.read().split('\n')
    centroid_indices = np.where(np.in1d(frame_names, centroid_names))  # Note: order of centroid_names is not preseved
    # Find names associated to the centroid_indices list
    reordered_centroid_names = (np.array(frame_names)[centroid_indices]).tolist()
    # Load binpos trajectory and extract frames for the centroids
    topology_file = os.path.join(chunkset.directory, "{}.pdb".format(chunkset.prefix))
    trajectory_file = os.path.join(chunkset.directory, "{}.binpos".format(chunkset.prefix))
    trajectory = mdtraj.load(trajectory_file, top=topology_file)
    trajectory_file_name = os.path.join(chunkset.directory, "{}.cluster.centroids.{}".format(chunkset.prefix, extension))
    trajectory[centroid_indices].save(trajectory_file_name)
    names_file_name = os.path.join(chunkset.directory, "{}.cluster.centroids.names".format(chunkset.prefix))
    with open(names_file_name, 'w') as handle:
        handle.write('\n'.join(reordered_centroid_names)+'\n')
    
def run_kmeans(chunksets, cluster_options={}):
    """
    Produces files {{prefix}}.clusters, and {{prefix}}.cluster.stats 
    :chunksets : list of named tuple with attributes 'directory' and 'prefix'
    
    :returns name of the centroid for each cluster
    """
    kmeans_options_dict={ "--nclusters" : 10,        # number of output clusters
                          "--min_cluster_size" : 1,  # minimum cluster size
                          "--max_iterations" : 9,    # max iterations for each random starting partition to converge to a final partition
                          "--total_seeds" : 9,       # number of different starting partitions to tr
                          "--nthreads": 4
                        }
    #update kmeans_options_dict with argument cluster_options
    kmeans_options_dict.update({k:v for k, v in cluster_options.items() if k in kmeans_options_dict})
    kmeans_options = " ".join([str(item) for pair in kmeans_options_dict.items() for item in pair])
    fixed_options = " --cluster_rmsd --fine_parallel --sse3 --binary_coords"
    for chunkset in chunksets:
        io_options = " -i {0} -o {0}".format(chunkset.prefix)
        commands = """cd {}
fast_protein_cluster {} {} {}
""".format(chunkset.directory, fixed_options, kmeans_options, io_options)
        subprocess.check_output(commands.replace('\n',';'), shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Kmeans clustering on heavy atoms binpos set")
    parser.add_argument("directory", type=str, help="Storing .binpos, .pdb, and .names files")
    parser.add_argument("basename", type=str, help="prefix for the files, e.g., HA001 for file HA001.pdb")
    parser.add_argument("--nclusters", type=int, default=10, help="number of clusters. Default is 10 clusters")
    parser.add_argument("--centroids", action='store_true', help="Save centroids to DCD file")
    args = parser.parse_args()

    chunkset = inputHAfiles(args.directory, args.basename)
    run_kmeans([chunkset,], {"--nclusters":args.nclusters})
    centroid_names = get_centroids(os.path.join(args.directory, "{}.cluster.stats".format(args.basename)))
    if(args.centroids):
        save_centroids(centroid_names, chunkset)

    sys.exit(0)
