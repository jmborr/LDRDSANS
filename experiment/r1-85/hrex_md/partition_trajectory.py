# -*- coding: utf-8 -*-

"""Given an input trajectory, partition it into N chunks and save"""

import collections
import argparse
import subprocess
import mdtraj
import os
from pdb import set_trace as tr

inputHAfiles = collections.namedtuple("inputHAfiles", "directory prefix".split())

def partition_into_HA_chunks(topology, trajectory, n, prefix="p", atom_selection="protein"):
    """
    
    """

    atom_selections = {"protein": "protein",
                      "protein heavy atoms": "protein and not element H" }
    if atom_selection not in atom_selections.keys():
        atom_selection = "protein"
    cwd = os.getcwd()  # current workding directory
    print "Loading trajectory {}...".format(args.trajectory)
    traj = mdtraj.load(args.trajectory, top=args.topology)
    
    atom_indices = traj.topology.select(atom_selection)  # protein heavy atoms
    chunk_size = int(traj.n_frames/args.n)
    chunksets = list()
    for index in range(args.n):
        print "Processing chunk number {}. Left {} chunks...".format(index, args.n-index)
        dirname = "{}{:03d}".format(prefix, index)
        subprocess.call(["mkdir", "-p", dirname])
        basename = dirname
        begin = index*chunk_size
        chunk = traj[begin:begin+chunk_size].atom_slice(atom_indices)
        chunk.save_binpos("{}/{}.binpos".format(dirname, basename))  # save to binpos format
        chunk[0].save_pdb("{}/{}.pdb".format(dirname, basename))  # save first frame to PDB
        names = '\n'.join(["{:05d}".format(i) for i in range(begin, begin+chunk_size)])
        open("{}/{}.names".format(dirname, basename), 'w').write(names)
        chunksets.append(inputHAfiles(os.path.join(cwd,dirname), basename))
    return chunksets
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Partition trajectory in N chunks")
    parser.add_argument("topology", type=str, help="topology file, PDB's are OK")
    parser.add_argument("trajectory", type=str, help="trajectory file")
    parser.add_argument("--prefix", type=str, default="p", help="name prefix for the chunks")
    parser.add_argument("n", type=int, help="number of chunks to divide the trajectory")
    args = parser.parse_args()

    chunksets = partition_into_HA_chunks(args.topology, args.trajectory, args.n, prefix=args.prefix)
    print chunksets
    
        
    


    
