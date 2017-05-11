import iprofile
import matplotlib.pyplot as plt
import os.path
import math
import numpy as np

"""
Read file centroids.agglomeration.history to fetch RMSD and centoid ID pairs
at each step of the clustering. Then calculate I(Q) for each centroid and
similarity beween the pairs. Finally, compare RMSD and similarity 
"""

# Location of sassena files containing I(Q)
sId = '/SNS/CAMM/users/jbq/development/LDRDSANS/scrkinase/sassena'

cframe = tuple(str(int(l.strip()))
               for l in open('centroids.names').readlines())  # trajectory frame for each centroid
rmsl = list()
siml = list()
for line in open('centroids.agglomeration.history', 'r').readlines():
    step, id1, id2, rms = list(i for i in line.strip().split())
    id1 = int(id1)
    id2 = int(id2)
    if id1 < 0:
        break
    rmsl.append(float(rms))
    p1 = iprofile.profile('{}/frame{}.h5.gz'.format(sId, cframe[id1])).i
    p2 = iprofile.profile('{}/frame{}.h5.gz'.format(sId, cframe[id2])).i
    siml.append(math.sqrt(iprofile.similarity(p1, p2)))
np.savetxt('centroid_rms_vs_i.dat',np.array([rmsl,siml]).transpose())
plt.scatter(rmsl, siml)
plt.show()

