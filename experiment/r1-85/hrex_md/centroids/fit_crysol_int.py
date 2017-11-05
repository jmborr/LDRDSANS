"""For each centroid, do a fit of its crysol profile to experiment
"""

from __future__ import print_function

import os
import numpy as np
from scipy.cluster.hierarchy import linkage
from idpflex import cnextend, properties, bayes
import matplotlib.pyplot as plt
from pdb import set_trace as tr

# Load experimental profile to a SaxsProperty
root_dir = '/home/jbq/repositories/LDRDSANS/experiment/r1-85/hrex_md'
exp_prof = os.path.join(root_dir, 'SAXS_SH4_UD.v2.dat')
exp_prop = properties.SaxsProperty()
exp_prop.from_ascii(exp_prof)

# Load the SAXS profiles for each centroid and fit to experiment with Mantid
chis_squared = list()
scalings = list()
backgrounds = list()
centroid_indexes = [i.strip() for i in open('centroids.names', 'r').readlines()]
for index in centroid_indexes:
    print('index =', index)
    # Load the centroid's crysol profile to a SaxsProperty
    file_name = os.path.join('crysol_int', '{}00.int'.format(index))
    ps = properties.SaxsProperty(name='saxs')
    ps.from_crysol_int(file_name)
    ps.y *= 1.0e-05  # Make computed and experimental intensities similar
    node = cnextend.ClusterNodeX(0)
    node.add_property(ps)
    # Create the fit model with the centroid's profile and fit
    f = bayes.model_at_node(node, 'saxs')
    fit_output = bayes.do_fit(f, exp_prop)
    # Retrieve chi2, scaling and background parameters
    chis_squared.append(fit_output.OutputChi2overDoF)
    p_table = fit_output.OutputParameters
    for irow in range(p_table.rowCount()):
        row = p_table.row(irow)
        if '.Scaling' in row['Name']:
            scalings.append(row['Value'])
        elif '.A0' in row['Name']:
            backgrounds.append(row['Value'])
    #SaveNexus(InputWorkspace=fit_output,
    #          Filename='/tmp/junk.nxs')

# Save chi2, scaling and background to a file
i_centroid_indexes = [int(i) for i in centroid_indexes]
cxs = np.array([i_centroid_indexes, scalings, backgrounds, chis_squared])
np.savetxt('fit_crysol_int.dat', cxs.transpose(),
           fmt='%05d           %7.3f         %7.3f            %5.2f',
           header='frame-number  scaling(1E-5)  background(1E-5)  chi_squared')
# Scatter plot the the chi2 values
plt.scatter(i_centroid_indexes, chis_squared)
plt.show()


