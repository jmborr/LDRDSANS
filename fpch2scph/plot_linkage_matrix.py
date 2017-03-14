# -*- coding: utf-8 -*-

"""
Matplotlib of the dendogram associated with the linkage matrix.

Thanks to Jorn's Blog
<https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/>
"""

# needed imports
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plots a dendogram from a scipy.cluster.hierarchy linkage matrix.")
    parser.add_argument("linkage", type=str, help="linkage matrix file, output from fpch2scph.py")
    parser.add_argument("--p", type=int, default=10, help="show only the last p merged clusters")
    args=parser.parse_args()
    Z=np.loadtxt(args.linkage)
    plt.title('Hierarchical Clustering Dendrogram (truncated)')
    plt.xlabel('sample index')
    plt.ylabel('RMSD (Angstroms)')
    dendrogram(
        Z,
        truncate_mode='lastp',  # show only the last p merged clusters
        p=args.p,  # show only the last p merged clusters
        show_leaf_counts=False,  # otherwise numbers in brackets are counts
        leaf_rotation=90.,
        leaf_font_size=12.,
        show_contracted=True,  # to get a distribution impression in truncated branches
    )
    plt.show()
    sys.exit(0)
