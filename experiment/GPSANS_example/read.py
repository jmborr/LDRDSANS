from __future__ import print_function, absolute_import

import numpy as np

from sas.sascalc.dataloader.readers.ascii_reader import Reader
from sas.sascalc.data_util.qsmearing import smear_selection
from sasmodels.models import tabulated, lorentz
from sasmodels.resolution import Pinhole1D

def fermi(q_values, nu=None, kt=None):
    """Silly simulated I(Q) profile using a Fermi-Dirac distribution"""
    if kt is None:
        kt = np.mean(q_values)
    if nu is None:
        nu = np.std(q_values)
    return 1.0/(1 + np.exp((q_values - nu) / kt))

# Read in experimental reduced file
r = Reader()
exp_file = 'HiResSANS_exp9_scan0013_0001_Iq.txt'
exp_data = r.read(exp_file)[0]
q_exp = exp_data.x
nu_0 = 2 * np.mean(q_exp)
kt_0 = nu_0 / 4.0

# Generate an array of Q values for the simulated profile, and
# parameters defining the shape of the simulated intensity lines
q_sim = np.arange(0.0, 2 * np.max(q_exp), 0.5 * np.max(q_exp) / len(q_exp))
tabulated.initialize(q=q_sim)
x = 1 - 0.5 * np.random.sample(10)  # generate 10 simulated profiles
nu = nu_0 * x
kt = kt_0 * x

# Iterate generating one simulated profile at a time, and
# convolving this profile with the experimental resolution
resolution = Pinhole1D(exp_data.x, exp_data.dx)
convolved_profiles = list()
for (a, b) in zip(nu, kt):
    intensities = fermi(q_sim, nu=a, kt=b)
    tabulated.initialize(I=intensities)
    i_sim = tabulated.Iq(resolution.q_calc)
    convolved_profiles.append(resolution.apply(i_sim))


for profile in convolved_profiles:
    buf = '# X  Y  dY  dX\n'
    for i, y in enumerate(profile):
        buf += '{} {}\n'.format(q_exp[i], y)
    open('/tmp/sim_convolved_{}.dat'.format(i), 'w').write(buf)




