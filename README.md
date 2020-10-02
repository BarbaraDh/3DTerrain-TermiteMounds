# 3DTerrain-TermiteMounds

main script: filter_termites_CC.py

Filters termite mounds from UAV las data.

PACKAGES NEED TO BE INSTALLED:
   os, numpy, pandas, sklearn, matplotlib, laspy, copy, time

STEPS:
1) Load data into CloudCompare; use 'Cloth Simulation Filter' (settings: 
filter Scene: Steep slope, cloth resolution: 1.0). Save as .las
2) python: keep lowest points. input_file_CSF: filtered point cloud in .las format (result of the cloth simulation filter (CloudCompare), filter the cloud)
3) CALCULATE NORMALS IN CC (radius = .75); convert to DIP DIRECTION/DIP DEGREE. Use quadric! Save it as GF_[name].txt                                                      
4) python: filter based on dip direction/dip degree
5) CALCULATE NORMALS IN CC (radius = .30) of the 'csf_to_calculate_normals.txt' file; convert to DIP DIRECTION/DIP DEG. Save it as GF_csf_to_calculate_normals.txt
6) python: extra filter based on dip dip degree


Extra information

step 2: file 'keep_bottom_points'. Is possible to run from the command line. Based on the X Y values, the nearest neighbours (default: 30) are identified for each point. Only the point with the lowest Z value is kept.

step 4: file 'mounddetection1'. 
    1) data is filtered based on dip degree (7-86 degrees)
    2) data is filtered: points with less than 50 points within a 2 m radius are discarded
    3) data is clustered using agglomerative clustering
    4) extra points are added to each cluster, from the 'low_points' filtered data
    5) Based on the dip direction, certain clusters are discarded. Creates a file 'resultsoptimizingcenter.csv' which lists the RMSE before and optimizing the termite mound center, and the distance between the initial center and optimized. Creates a folder 'plots' in which the optimized center is visualised. Creates a file 'temporaryresultsmounds.txt'. 
    6) extra points are added to each cluster, coming from the 'CSF' data

step 6: file 'mounddetection2'. Clusters are filtered based on dip degree. All clusters that have > 77 % amount of points (above 15 cm ground surface) that are steeper than 77◦ are discarded.
