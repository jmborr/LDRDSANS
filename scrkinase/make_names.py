nframes=4630
outfilename="HA.names"
from pdb import set_trace as tr
buf=""
for i in range(nframes):
    #tr()
    buf += "f{0:05d}\n".format(i)
open(outfilename, 'w').write(buf)
