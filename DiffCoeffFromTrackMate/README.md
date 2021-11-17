Description:
This script will extract the track data and metadata from the output .xml 
file from the TrackMate plug-in of Image J. It will then use that extracted
data to determine the 2D diffusion coefficient for each track, by fitting
into a linear model. The best number of points to include in the fit is
determined following the algorithm in Michalet, Phys Rev. E., 2010, 82: 041914
(https://doi.org/10.1103/PhysRevE.82.041914). The outputs are a histogram
showing all the measured diffusion coefficients, and a histogram showing
only the diffusion coefficients of those above the immobile cutoff.
All of the data is saved as a .MAT file.

Code written by Bob Rawle, Williams College, 2021
Published online in conjunction with the manuscript by Lam et al. (2021).
