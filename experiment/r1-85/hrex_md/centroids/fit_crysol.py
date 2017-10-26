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
for rg in [float(r) for r in open('crysol_fit/shell_Rg.dat').readlines()]:
    ps = properties.ScalarProperty(name='Rg')
    ps.setScalar(rg)
    rg_props.append(ps)
properties.propagator_size_weighted_sum(rg_props, tree)

# Fit to experimental profile
exp_prof = '../SAXS_SH4_UD.v2.dat'
exp_prop = properties.SaxsProperty()
exp_prop.from_ascii(exp_prof)
fits_out = bayes.do_fit_at_depth_tree(tree, exp_prop, 'saxs', max_depth=10)

# Output scaling parameters
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
plt.plot(range(len(fits_out)), chis_squared)
#plt.show()

"""
depth =  0 chi2 =29.21 scalings = 1.000
depth =  1 chi2 =17.63 scalings = 1.000 0.000
depth =  2 chi2 =19.99 scalings = 0.000 0.775 0.225
depth =  3 chi2 =11.15 scalings = 0.644 0.027 0.000 0.329
depth =  4 chi2 =11.69 scalings = 0.535 0.123 0.023 0.319 0.000
depth =  5 chi2 =14.25 scalings = 0.422 0.000 0.054 0.175 0.075 0.274
depth =  6 chi2 =14.68 scalings = 0.152 0.351 0.000 0.045 0.142 0.062 0.249
depth =  7 chi2 =14.40 scalings = 0.147 0.126 0.230 0.308 0.000 0.038 0.104 0.048
depth =  8 chi2 =13.41 scalings = 0.189 0.129 0.104 0.189 0.244 0.023 0.000 0.029 0.092
depth =  9 chi2 =13.92 scalings = 0.096 0.173 0.114 0.082 0.094 0.168 0.223 0.022 0.000 0.028

"""
#  Chi2 values: 29.21 17.63 19.99 11.15 11.69 14.25 14.68 14.40 13.41 13.92
depth_optimal = 3  # minimizes chi2
nodes = tree.clusters_at_depth(depth_optimal)
print(' '.join('{:5.2f}'.format(node['Rg'].y) for node in nodes ))  # 26.66 22.66 20.70 28.53
tr()
print(' '.join('{:3d}'.format(node.count) for node in nodes ))  # 