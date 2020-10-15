# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 10:06:46 2020
Filters termite mounds from las data.

PACKAGES NEED TO BE INSTALLED:
   os, numpy, pandas, sklearn, matplotlib, laspy, copy, time

STEPS:
1) Load data into CloudCompare; use 'Cloth Simulation Filter' 
    SETTINGS: filter Scene: Steep slope, cloth resolution: 1.0, max. iterations: 500, classification threshold: 1).
    OUTPUTS: Save as [name].las
    EXPLANATION: ground filtering algorithm, should comprise the complete termite mounds. Under story can be included.
    https://www.cloudcompare.org/doc/wiki/index.php?title=CSF_(plugin)
    
2) python: keep lowest points.
    INPUTS: input_file_CSF: filtered point cloud in .las format (result of the cloth simulation filter (CloudCompare)
            neighbours: default 30
    OUTPUTS:low_points_knn[neighbours]_[name].las
    EXPLANATION: For each point, neighbours are taken into consideration and only the lowest point remains. Duplicates are removed.
        
3) Calculate normals in CloudCompare
    INPUTS: low_points_knn[neighbours]_[name].las
    SETTINGS: local surface model: Quadric; radius: default 0.75 m; Orientation: Use mininum Spanning tree - knn = 6
    convert to DIP DIRECTION/DIP DEGREE. 
    OUTPUTS: Save as GF_[name].txt  
    EXPLANATION: https://www.cloudcompare.org/doc/wiki/index.php?title=Normals%5CConvert_to_Dip_and_Dip_direction_SFs
                                                   
4) python: filter based on dip direction/dip degree
    INPUTS:  GF_low_points_knn[neighbours]_[name].txt  
            th_dip_degree: default (7,87) (degrees)
            th_NN: default 50 
            th_clust: default 30
            th_distance: default .75 m 
            th_RMSE: default 1
            extra_info: default False
    OUTPUTS: csf_to_calculate_normals.txt
            extra_info: see mounddetection1 - function optimize_center 
    EXPLANATION: 1. points are filtered: all points between th_dip_degree remain
                2. only points that have > th_NN points within 2 meter remain
                3. Data is clustered using agglomerative clustering (distance thershold = th_clust)
                4. Cluster is checked if it has a conical shape: center is optimized based on dip direction. 
                If optimized center point is > th_distance from mean (x,y), or RMSE of the optimalisation
                > th_RMSE, cluster is discarded.
            
5)  Calculate normals in CloudCompare
    INPUTS: csf_to_calculate_normals.txt
    SETTINGS: local surface model: Quadric; radius: default 0.30 m; Orientation: Use mininum Spanning tree - knn = 6
    convert to DIP DIRECTION/DIP DEGREE.    
    OUTPUTS: Save as GF_csf_to_calculate_normals.txt
    EXPLANATION: https://www.cloudcompare.org/doc/wiki/index.php?title=Normals%5CConvert_to_Dip_and_Dip_direction_SFs
    
6) python: extra filter based on dip degree. 
    INPUTS: GF_csf_to_calculate_normals.txt
            th1: default 77 (%)
            th2: default 77 (degrees)
            n: default 0.15 (m)
    OUTPUTS: termite_mounds.txt: labelled termite mounds
             termite_mounds.png: map of the termite mounds
    EXPLANATION: all cluster of which > th1 % of points (above n m ground surface) have a 
            dip degree > th2° are discarded.
                                                                       
                                                                                     
@author: bgdhont                                                          
"""
import keep_bottom_points as kbp
import mounddetection1
import mounddetection2
import os

#%% STEP 1: NOT POSSIBLE VIA PYTHON. NEEDS TO BE DONE MANUALLY IN CLOUDCOMPARE
#%% STEP 2: keep lowest points
#parameters:
neighbours = 30 

input_file_CSF = r'C:\Users\bgdhont\Desktop\temp\F1_HR_CSF_p2.las'
#amount of nearest neighbours used to calculate low points 
directory, name  = os.path.split(input_file_CSF) # 2
kbp.calc_lower_points(input_file_CSF, neighbours)
file_low_p = os.path.join(directory, 'low_points_knn' + str(neighbours) + '_' + name)
directory, name  = os.path.split(file_low_p)


#%%  STEP 3:calculate normals (MANUAL OR VIA COMMAND LINE)
rad = .75

PC = file_low_p
directory, name  = os.path.split(file_low_p)
PC_output = os.path.join(directory, 'GF_' + name.split('.las')[0] + '.txt')

command = "Cloudcompare -SILENT -C_EXPORT_FMT ASC -ADD_HEADER -SEP COMMA -AUTO_SAVE OFF -O -GLOBAL_SHIFT AUTO " + PC + " -OCTREE_NORMALS " + str(rad) +  " -NORMALS_TO_DIP -SAVE_CLOUDS FILE " + PC_output
os.system("start \"\" cmd /K \"cd C:\\Program Files\\CloudCompare & " + command + "\"" )
#open  


#%% STEP 4: Label possible locations
#parameters:
th_dip_degree = (7,87)
th_NN = 50 
th_clust = 30
th_distance =  .75 
th_RMSE = 1
extra_info =  False

input_file_dips_GF = os.path.join(directory, 'GF_' + name.split('.las')[0] + '.txt')
output_dir = directory
mounddetection1.filter1(input_file_dips_GF, input_file_CSF, output_dir,
                        th_dip_degree, th_NN,th_clust, th_distance, th_RMSE, extra_info)
file_to_calc_normals = output_dir + '\\csf_to_calculate_normals.txt'
#%% STEP 5: calculate normals (MANUAL OR VIA COMMAND LINE)
rad = .3

PC = file_to_calc_normals
directory, name  = os.path.split(file_to_calc_normals)
PC_output = os.path.join(directory, 'GF_' + name)


command = "Cloudcompare -SILENT -C_EXPORT_FMT ASC -ADD_HEADER -SEP COMMA -AUTO_SAVE OFF -O -GLOBAL_SHIFT AUTO " + PC + " -OCTREE_NORMALS " + str(rad) +  " -NORMALS_TO_DIP -SAVE_CLOUDS FILE " + PC_output
os.system("start \"\" cmd /K \"cd C:\\Program Files\\CloudCompare & " + command + "\"" )
 

#%% STEP 6
th1 = 77
th2 = 77
n = 0.15

directory, name  = os.path.split(file_to_calc_normals)
input_file_dips_GF2 = os.path.join(directory, 'GF_' + name)
mounddetection2.filter2(input_file_dips_GF2, output_dir, th1,th2, n)










