# -*- coding: utf-8 -*-

"""Fast protein cluster hierarchical to scypi cluster hierarchy

This module reads the output from a run of fast_protein_cluster and computes
the linkage matrix for use with scipy.cluster.hierarchy
"""

import argparse
import os
import sys
import numpy as np
from pdb import set_trace as tr

def validate_history(filename):
    """Examine the fast_protein_cluster output history file for validation.

       File contains lines of the format "int int int float". The last two lines
       must be of the form "int -1 -1 float".

       :param filename: the name of the history file

       :return: True only if file has valid format
    """
    if not os.path.exists(filename):
        return False
    lines = open(filename).readlines()

    items = [float(x) for x in lines[0].split()]
    # line must contain four numbers
    if(len(items)!=4):
        return False
    # the first three numbers must be integers
    for number in items[0:3]:
        if not number.is_integer():
            return False
    # Only the last two lines must contain "int -1 -1 float",
    # indicating clustering stopped when two clusters remained
    last_three_lines = lines[-3:]
    if [int(i) for i in last_three_lines[0].split()[1:3]] == [-1, -1]:
        return False
    for line in last_three_lines[1:]:
        if not [int(i) for i in line.split()[1:3]] == [-1, -1]:
            return False
    # We require that the ID's of the singletons are numbers running
    # from zero up to the number of singletons minus one
    nsingletons=len(lines)
    for line in lines:
        if not [int(i)<nsingletons for i in line.split()[1:3]] == [True, True]:
            return False
    # All tests passed
    return True

def generate_linkage(history_handle):
    """Create a linkage matrix for use with scipy.cluster.hierarchy

    The linkage matrix has shape (N-1,4) where N is number of singletons. Detailed
    explanation of the matrix in the scipy.cluster.hierarchy.linkage webpage. Also
    read the process by which the linkage matrix is created.

    :param history_handle: handle to the history file
    :type history_handle: file

    :return: linkage matrix
    :rtype: numpy.ndarray
    """
    # Figure out number of singletons
    lines = history_handle.readlines()
    nsingletons = len(lines)

    # Initialization of several variables
    remains = range(nsingletons) # every step we remove a singleton from the pool
    translation = dict(zip(remains, remains))  # translate fpc ID to scipy ID
    last_id = nsingletons  # scipy assigns cluster ID's above number of singletons
    linkage_matrix = list()
    cluster_size = [1,]*nsingletons
    distance=-1.0
    # process each line except the last two (see function validate_history)
    for line in lines[:-2]:
        items = [float(x) for x in line.split()]
        fpc_id1 = int(items[1])  # fast_protein_cluster ID
        fpc_id2 = int(items[2])
        distance = float(items[3])
        scp_id1 = translation[fpc_id1]  # scipy ID
        scp_id2 = translation[fpc_id2]
        new_size = cluster_size[scp_id1]+cluster_size[scp_id2]
        cluster_size.append(new_size)
        remains.remove(fpc_id1) # fast_protein_cluster discards ID of the first singleton
        linkage_matrix.append([scp_id1, scp_id2, distance, new_size])
        translation[fpc_id2] = last_id  # mapping to the ID of the newly created cluster
        last_id += 1
    # insert the joining of the last two clusters with a fictitious distance
    fpc_id1 = remains[0]
    fpc_id2 = remains[1]
    final_distance = 1.1*distance  # fictious distance when joining the two clusters
    linkage_matrix.append([translation[fpc_id1], translation[fpc_id2], final_distance, nsingletons])

    return np.array(linkage_matrix)  # shape is (nsingletons, 4)

def save_linkage(linkage_matrix, filename):
    """Save the linkage matrix to file

    Each row of the matrix is saved to a row in the output file

    :param linkage_matrix:
    :type linkage_matrix: numpy.ndarray
    :param filename: file name
    :type filename: basestring
    :return: success of the saving process
    :rtype: bool
    """
    buffer=''
    for i in range(len(linkage_matrix)):
        row=linkage_matrix[i]
        buffer += "{0:6d} {1:6d} {2:7.3f} {3:7d}\n".format(int(row[0]), int(row[1]), row[2], int(row[3]))
    try:
        open(filename, "w").write(buffer)
        return True
    except:
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create linkage matrix for scipy.hierarchy")
    parser.add_argument("history", type=str, help="agglomeration.history file from fast_protein_cluster")
    parser.add_argument("linkage", type=str, help="save linkage matrix to file linkage")
    args=parser.parse_args()
    if not validate_history(args.history):
        sys.stderr.write("Error: History file is not valid.\n")
        sys.exit(1)
    z = generate_linkage(open(args.history, "r"))
    if not save_linkage(z, args.linkage):
        sys.stderr.write("Error: Unable to save linkage matrix to file.\n")
        sys.exit(1)
    sys.exit(0)
