import h5py
import gzip
import shutil
import tempfile
import numpy as np
import os
from collections import namedtuple
import matplotlib.pyplot as plt

qi = namedtuple('qi', 'q i')

def profile(gzipfile, log10=False):
    """
    Return momentum transfer and intensities 
    :param gzipfile: gzipped sassena file
    :param log10: return log10 of the intensities
    :return: namedtuple qi
    """
    col1 = np.s_[:,0]  # first column
    ft, ftname = tempfile.mkstemp(dir='/tmp')
    ft2 = open(ftname, 'wb')
    gt = gzip.open(gzipfile)
    shutil.copyfileobj(gt, ft2)
    gt.close()
    ft2.close()
    with h5py.File(ftname) as f:
        i = f['fq'][col1]
        q = f['qvectors'][col1]
        reorder = np.argsort(q)
        q = q[reorder]
        i = i[reorder]
        if log10:
            i = np.log10(i)
        f.close()
    os.close(ft)
    os.remove(ftname)
    return qi(q, i)

def similarity(profile1, profile2, weights = None):
    """
    A measure of similarity between two I(Q) profiles
        1/len(profile1) * Sum ((profile1-profile2)/(profile1+profile2))**2)
    :param profile1: first profile
    :param profile2: second profile
    :return: similarity measure
    :except: profiles of different length   
    """
    if weights is None:
        weights = np.ones(len(profile1))
    if len(profile1) != len(profile2):
        raise IndexError("profiles have different length")
    return np.sum(weights*((profile1-profile2)/(profile1+profile2))**2)/len(profile1)

if __name__ == '__main__':
    """
    For debugging, calculate average intensities
    """
    nframe = 4628
    p = list()
    for iframe in range(1, 1+nframe):
        p.append(profile('frame{}.h5.gz'.format(iframe)).i)
        if iframe%100 == 0:
            print(iframe)
    p = np.array(p)
    y = np.mean(p, axis=0)
    e = np.std(p, axis=0)
    x = profile('frame1.h5.gz').q
    plt.errorbar(x, y, yerr=e)
    plt.show()




             

