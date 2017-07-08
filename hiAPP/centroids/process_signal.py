import os
cwd = os.getcwd()
LoadSassena(Filename=os.path.join(cwd, 'signal.h5'), OutputWorkspace='sans')
icq=Transpose(InputWorkspace='sans_fqt.Re')
# Recall that LoadSassena simmetrize time to negative values. We have to remove these
# as they have no meaning
nh = icq.getNumberHistograms()
icq=ExtractSpectra(icq,StartWorkspaceIndex=int((nh-1)/2), EndWorkspaceIndex=nh-1)
SaveNexus(icq, os.path.join(cwd, 'centroids.sans.nxs'))
SaveAscii(icq, os.path.join(cwd, 'centroids.sans.dat'),
          ColumnHeader=False, WriteSpectrumID=False,
          Separator='Space')

