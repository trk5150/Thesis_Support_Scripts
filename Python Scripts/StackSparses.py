# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 15:06:14 2024

@author: tik105
"""

import numpy as np
from scipy import sparse
from scipy.io import mmread
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack, save_npz

# Directory where the matrices are stored
directory = "C://Users//tik105//Desktop//Gut and UC SC SOM analysis//"

# Initialize a list to hold all the matrices
matrices = []

# Loop through the matrix files (assuming they are named matrix.mtx, matrix2.mtx, ..., matrix18.mtx)
for i in range(1, 19):
    # Construct the file path
    file_path = f"{directory}matrix{i}.mtx" if i > 1 else f"{directory}matrix.mtx"
    
    # Read the matrix and append to the list
    matrices.append(mmread(file_path))

# Stack all the matrices horizontally
xy = sparse.hstack(matrices)

print(xy.shape)
# (10, 20)

output_file_mtx = "C://Users//tik105//Desktop//Gut and UC SC SOM analysis//stacked_matrix.mtx"
mmwrite(output_file_mtx, xy)

print(type(xy))
# <class 'scipy.sparse.coo.coo_matrix'>