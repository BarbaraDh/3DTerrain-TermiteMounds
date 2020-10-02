# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:58:51 2020

@author: bgdhont
"""

#%% DIT IS HEM


import os
import numpy as np
from sklearn.neighbors import KDTree
import laspy
import copy
import time

def keep_minima(dataset, k):
    """ 
    Based on the X Y values, the nearest neighbours are identified for each point. Only the point with the lowest Z value is kept.
    """
    datasetxy = dataset[:,:2] 
    tree= KDTree(datasetxy)
    ind = tree.query(datasetxy, k = k, return_distance = False) #works! gives int 
    f = dataset[ind]
    min_vals = np.amin(f[:,:,2], axis = 1)
    idmins = np.array([np.where(f[i,:,2] == min_vals[i])[0][0] for i in range(len(min_vals))]) #duurt lang
    minima = np.array([f[i, idmins[i]] for i in range(len(idmins))])
    unique_rows = np.unique(minima, axis=0)
    return unique_rows

def save_as_las(inFile,dataset, input_file, neighbours):
    hdr = copy.copy(inFile.header)
    directory, name  = os.path.split(input_file)
    path2 = os.path.join(directory, 'low_points_knn' + str(neighbours) + '_' + name)
    outfile = laspy.file.File(path2, mode="w",vlrs = inFile.header.vlrs,  header=hdr)
    outfile.x = dataset.transpose()[0]
    outfile.y = dataset.transpose()[1]
    outfile.z = dataset.transpose()[2]
    outfile.close()

def calc_lower_points(input_file, neighbours):
    """
    main function.
    1) las data is read and converted to numpy array
    2) local minima are calculated
    3) result is converted to .las file
    
    """
    t0 = time.process_time()
    print('Searching for minima...')
    inFile = laspy.file.File(input_file, mode = "r") #1
    dataset = np.vstack([inFile.x, inFile.y, inFile.z]).transpose() #1
    minima = keep_minima(dataset, neighbours) #2
    save_as_las(inFile, minima, input_file, neighbours) #3
    inFile.close()
    t1 = time.process_time()
    len_minima = len(minima)
    len_dataset = len(dataset)
    perc = np.round(len_minima/len_dataset*100,2)
    print(str(len_minima) + ' of ' + str(len_dataset) + ' points kept (' + str(perc) + ' %)')
    print("Time spent: " + str(t1 - t0) + ' s' )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = 'Only keep low points')
    parser.add_argument("input", help = 'input las file')
    parser.add_argument("--neighbours", help = 'NN. Default: 30')
    args = parser.parse_args()
    input_file = args.input
    if args.neighbours:
        neighbours = int(args.neighbours)
    else:
        neighbours = 30
    calc_lower_points(input_file, neighbours)
    