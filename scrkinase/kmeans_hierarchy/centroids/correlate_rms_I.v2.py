import iprofile
import matplotlib.pyplot as plt
import os.path
import math
import numpy as np

"""
Reads file rms.matrix to fetch centroid RMSD and Similarity in I(Q) for all centroid pairs.
Finally, compare RMSD and Similarity
"""
# Location of sassena files containing I(Q)
sId = '/SNS/CAMM/users/jbq/development/LDRDSANS/scrkinase/sassena'

cframe = tuple(str(int(l.strip()))
               for l in open('centroids.names').readlines())  # trajectory frame for each centroid
rmsl = list()
siml = list()
counter = 0

for line in open('rms.matrix', 'r').readlines():
    id1, id2, rms = list(i for i in line.strip().split())
    id1 = int(id1)
    id2 = int(id2)
    if id1 < 0:
        break
    rmsl.append(float(rms))
    q, p1 = iprofile.profile('{}/frame{}.h5.gz'.format(sId, cframe[id1]))
    p2 = iprofile.profile('{}/frame{}.h5.gz'.format(sId, cframe[id2])).i
    siml.append(math.sqrt(iprofile.similarity(p1, p2)))
    counter += 1
    print(counter)
np.savetxt('all_rms_vs_i.dat',np.array([rmsl,siml]).transpose())
plt.scatter(rmsl, siml)
plt.show()

