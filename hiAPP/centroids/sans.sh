#!/bin/bash

#####################################################################
#  SASSENA run produces SANS scattering for each centroid in fqt.Re
#  Problem: sassena in my machine does not work. I had to use the one in camm2.sns.gov
#####################################################################
module load sassena
mpirun -np 8 sassena --config scatter_sans.xml &> scatter_sans.log  # produces signal.h5
python process_signal.h5  # saves SANS profiles from signal.h5 in centroids.sans.h5 and centroids.sans.dat



