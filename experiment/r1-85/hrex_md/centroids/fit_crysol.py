"""For the topmost clusters in the hierarchical tree, fit their
crysol-derived profiles to experiment.
Also load radius of gyration for the centroids and propagate up
the tree.
"""

from __future__ import print_function

import os
import numpy as np
from idpflex import cnextend, properties, bayes
import matplotlib.pyplot as plt
from pdb import set_trace as tr

# Load the clustering
tree = cnextend.Tree()
Z = np.loadtxt('centroids.linkage_matrix')
tree.from_linkage_matrix(Z)

# Load the SAXS profiles for each centroid and propagate up the tree nodes
saxs_props = list()
profdir = 'crysol_int'
for index in open('centroids.names', 'r').readlines():
    # SAXS property
    file_name = os.path.join(profdir,'{}00.int'.format(index.strip()))
    ps = properties.SaxsProperty(name='saxs')
    ps.from_crysol_int(file_name)
    ps.y *= 0.5e-05  # Make computed and experimental intensities similar
    saxs_props.append(ps)
properties.propagator_size_weighted_sum(saxs_props, tree)

# Load the Rg profiles for each centroid and propagate up the tree nodes
rg_props = list()
for rg in np.loadtxt('crysol_fit/rg_chi2.dat', usecols=1):
    ps = properties.ScalarProperty(name='Rg')
    ps.setScalar(rg)
    rg_props.append(ps)
properties.propagator_size_weighted_sum(rg_props, tree)

# Load experimental profile to a SaxsProperty and fit the topmost clusters
exp_prof = '../SAXS_SH4_UD.v2.dat'
exp_prop = properties.SaxsProperty()
exp_prop.from_ascii(exp_prof)
fits_out = bayes.do_fit_at_depth_tree(tree, exp_prop, 'saxs', max_depth=10)

# Output scaling parameters and chi2
for depth, fit_out in enumerate(fits_out):
    chi_squared = fit_out.OutputChi2overDoF
    scalings = list()
    p_table = fit_out.OutputParameters
    for irow in range(p_table.rowCount()):
        row = p_table.row(irow)
        if '.Scaling' in row['Name']:
            scalings.append(row['Value'])
    scalings = np.array(scalings) / sum(scalings)
    print('depth = ', depth, 'chi2 ={:5.2f}'.format(chi_squared),
          'scalings =',
          ' '.join(['{:5.3f}'.format(s) for s in list(scalings)]))
chis_squared = [fit_output.OutputChi2overDoF for fit_output in fits_out]
print(' '.join('{:5.2f}'.format(c) for c in chis_squared))
buffer = '# depth  Chi2\n'
for i, chi2 in enumerate(chis_squared):
    buffer += '{:3d} {:5.2f}\n'.format(i, chi2)
open('chi2_vs_depth.dat', 'w').write(buffer)
plt.plot(range(len(fits_out)), chis_squared)
plt.show()

"""
depth =  0 chi2 =19.38 scalings = 1.000
depth =  1 chi2 = 9.05 scalings = 0.999 0.001
depth =  2 chi2 =10.99 scalings = 0.000 0.787 0.213
depth =  3 chi2 = 5.32 scalings = 0.613 0.036 -0.000 0.351
depth =  4 chi2 = 5.52 scalings = 0.514 0.145 0.018 0.323 0.000
depth =  5 chi2 = 6.96 scalings = 0.424 -0.000 0.043 0.190 0.065 0.278
depth =  6 chi2 = 7.26 scalings = 0.147 0.354 0.000 0.036 0.154 0.053 0.256
depth =  7 chi2 = 6.94 scalings = 0.153 0.120 0.237 0.314 0.000 0.029 0.109 0.038
depth =  8 chi2 = 6.56 scalings = 0.182 0.140 0.100 0.189 0.242 0.019 0.000 0.023 0.105
depth =  9 chi2 = 6.64 scalings = 0.106 0.172 0.118 0.084 0.090 0.172 0.224 0.013 0.000 0.021
"""
# Chi2 values at increasing depths, starting with zero
# 19.38  9.05 10.99  5.32  5.52  6.96  7.26  6.94  6.56  6.64

# Output radius of gyration and cluster size for the clusters at the level
# that minimizes chi2 from the previous fit.
depth_optimal = 3  # minimizes chi2
nodes = tree.clusters_at_depth(depth_optimal)
print(' '.join('{:5.2f}'.format(node['Rg'].y) for node in nodes ))  #

print(' '.join('{:3d}'.format(node.count) for node in nodes ))  # 22 270 330 378
