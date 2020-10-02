# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 10:06:46 2020
Filters termite mounds from las data.

PACKAGES NEED TO BE INSTALLED:
   os, numpy, pandas, sklearn, matplotlib, laspy, copy, time

STEPS:
1) Load data into CloudCompare; use 'Cloth Simulation Filter' (settings: 
filter Scene: Steep slope, cloth resolution: 1.0). Save as .las
2) python: keep lowest points. input_file_CSF: 
filtered point cloud in .las format (result of the cloth simulation filter
(CloudCompare), filter the cloud)
3) CALCULATE NORMALS IN CC (radius = .75); convert to DIP DIRECTION/DIP DEGREE. Use quadric! Save it as GF_[name].txt                                                      
4) python: filter based on dip direction/dip degree
5) CALCULATE NORMALS IN CC (radius = .30) of the 'csf_to_calculate_normals.txt' file; convert to DIP DIRECTION/DIP DEG.
Save it as GF_csf_to_calculate_normals.txt
6) python: extra filter based on dip direction/dip degree

@author: bgdhont
"""
import keep_bottom_points as kbp
mport mounddetection2
import os

#%% STEP 2
input_file_CSF = 'C:/Users/bgdhont/Documents/Data/TEMP/Litchfield_F2_UTM52S_without_sl15_16_thinned_CSF_p1.las'
neighbours = 30 #amount of nearest neighbours used to calculate 
directory, name  = os.path.split(input_file_CSF) # 2
kbp.calc_lower_points(input_file_CSF, neighbours)
file_low_p = os.path.join(directory, 'low_points_knn' + str(neighbours) + '_' + name)
directory, name  = os.path.split(file_low_p)

#%% STEP 4
input_file_dips_GF = os.path.join(directory, 'GF_' + name.split('.las')[0] + '.txt')
output_dir = directory
mounddetection1.filter1(input_file_dips_GF, input_file_CSF, output_dir) 
file_to_calc_normals = output_dir + '/csf_to_calculate_normals.txt'
#%% STEP 6
directory, name  = os.path.split(file_to_calc_normals)
input_file_dips_GF2 = os.path.join(directory, 'GF_' + name)
mounddetection2.filter2(input_file_dips_GF2, output_dir)










