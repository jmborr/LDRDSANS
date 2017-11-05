from __future__ import print_function

import re
from pdb import set_trace as tr

summary = open('crysol_summary.txt').read()

# Extract the hits
patterns = dict(name=r'Model\:\s+(\d+)\.pdb',
                rg=r'Rg:\s+(\d+\.\d+)\s',
                chi2=r'Chi\^2:\s+(\d+\.\d+)\n')
to_type = dict(name=str, rg=float, chi2=float)
hits = {key: [to_type[key](h) for h in re.findall(pattern, summary)]
        for key, pattern in patterns.items()}

# Validation
lengths = [len(hits[key]) for key in hits]
if lengths[1:] != lengths[:-1]:
    raise ValueError('Not all hits have the same length')

# Save to file
buffer = '#  centroid  Rg    Chi2\n'
for i in range(len(hits['name'])):
    vals = {key: hits[key][i] for key in hits}
    buffer += '    {name:6s}  {rg:5.2f}  {chi2:5.2f}\n'.format(**vals)
open('rg_chi2.dat','w').write(buffer)
