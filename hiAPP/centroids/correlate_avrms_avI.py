import iprofile
import math
import numpy as np
import scipy.cluster.hierarchy as sch
import navigate
import matplotlib.pyplot as plt
from mantid import simpleapi as msapi

# Load centroid names
cframe = tuple(str(int(l.strip()))
               for l in open('centroids.names').readlines())  # trajectory frame for each centroid
ncentroids = len(cframe)

# Populate rms and sim matrixes
ws = msapi.LoadNexus("centroids.sans.nxs")  # Load I(Q) for each centroid
rmsM = np.zeros(ncentroids ** 2).reshape(ncentroids, ncentroids)
simM = np.zeros(ncentroids ** 2).reshape(ncentroids, ncentroids)
counter = 0
for line in open('rms.matrix', 'r').readlines():
    id1, id2, rms = list(i for i in line.strip().split())
    id1 = int(id1)
    id2 = int(id2)
    if id1 < 0:
        break
    rmsM[id1][id2] = float(rms)
    rmsM[id2][id1] = rmsM[id1][id2]
    simM[id1][id2] = math.sqrt(iprofile.similarity(ws.dataY(id1), ws.dataY(id2)))
    simM[id2][id1] = simM[id1][id2]
    counter += 1
    if not counter%1000:
        print(counter)

# Calculate average rms and average sim between clusters
rmsl = list()
siml = list()
max_cluster_level = 5
Z = np.loadtxt('./centroids.linkage_matrix')
root_node = sch.to_tree(Z)
cad = navigate.clusters_afloat_depth(root_node, depth=max_cluster_level)
for level in range(1, len(cad)):
    clusters = cad[level]  # list of clusters at a particular level from the root node
    nclusters = len(clusters)
    for i in range(nclusters-1):
        for j in range(i+1, nclusters):
            rms = 0.0
            sim = 0.0
            ileafs = [leaf.id for leaf in navigate.leafs_under_node(clusters[i])]
            jleafs = [leaf.id for leaf in navigate.leafs_under_node(clusters[j])]
            for il in ileafs:
                for jl in jleafs:
                    rms += rmsM[il][jl]
                    sim += simM[il][jl]
            rmsl.append(rms/(len(ileafs)*len(jleafs)))
            siml.append(sim/(len(ileafs)*len(jleafs)))
np.savetxt('correlate_avrms_avI.dat',np.array([rmsl,siml]).transpose())
plt.scatter(rmsl, siml)
plt.show()


