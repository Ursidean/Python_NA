"""
This is the base file for running the Neighbourhood Analyser in Python.
The neighbourhood analyser requires a Python interpreter to be running
a version of Python 2. Python 3 will NOT be able to read in the data.
There are several modules that must be installed:
- Numpy
- gdal
The anaconda file is recommended for this:
https://www.continuum.io/downloads
"""

# Modules
from read_map import read_map
import numpy as np
from considered_distances import considered_distances
from EF import EF
import math
import csv


# Example application
# Specify the data paths.
map_path_1 = (
    "C:\\Users\\charl\\OneDrive\\Documents\\Python_NA\\Example_data\\lu1989.asc"
)
map_path_2 = (
    "C:\\Users\\charl\\OneDrive\\Documents\\Python_NA\\Example_data\\lu2000.asc"
)
mask_path = (
    "C:\\Users\\charl\\OneDrive\\Documents\\Python_NA\\Example_data\\region.asc"
)

# Read in the maps.
omap = read_map(map_path_1)
amap = read_map(map_path_2)
mask = read_map(mask_path)

# Analyse the maps.
mapshape = np.shape(omap)
rows = mapshape[0]
cols = mapshape[1]

# Specify the land-use class names.
luc_names = [
    "Other agriculture", "Pastures", "Arable land", "Greenhouses",
    "Housing low density", "Housing high density", "Industry", "Services",
    "Socio cultural uses", "Forest", "Extensive grasslands", "Nature",
    "Recreation areas", "Airport", "Fresh water", "Marine water"
]
# Determine the number of land-use classes.
luc = len(luc_names)
# Set the number of passive land-use classes.
pas = 3
# Set the number of feature land-use classes.
fea = 3
# Calculate the number of active land-use classes.
act = luc - (fea + pas)

# Specify the maximum neighbourhood distance considered.
dmax = 9
# Partition the distances into discrete rings.
x = considered_distances(dmax)
cd = x[0]
cdl = x[1]
# Generate a list of maximum possible neighbourhood sizes.
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
# Partition to size of distances considered
N = []
for c in range(0, dmax):
    N.append(N_all[c])

# Calculate the enrichment factor values.
data_ef = EF(luc, dmax, cdl, cd, N, omap, amap, mask, rows, cols)

# Log scale the enrichment factor values.
log_ef = np.zeros(shape=(dmax, luc, luc))
for p in range(0, luc):
    for q in range(0, luc):
        for c in range(0, dmax):
            if data_ef[c, p, q] == 0:
                log_ef[c, p, q] = -9999
            else:
                log_ef[c, p, q] = math.log(data_ef[c, p, q], 10)

# Write the output.
output_file_path = (
    "C:\\Users\\charl\\OneDrive\\Documents\\Python_NA\\Example_output\\NA_results.csv"
)
store = [0]*4
with open(output_file_path, "wb") as csv_file:
    writer = csv.writer(csv_file)
    values = ["from", "to", "d", "EF"]
    writer.writerow(values)
    for i in range(0, luc):
        for j in range(0, act):
            for c in range(0, dmax):
                store[0] = luc_names[i]
                store[1] = luc_names[j + pas]
                store[2] = c
                store[3] = log_ef[c, j + pas, i]
                writer.writerow(store)

# Analyser completed running.
