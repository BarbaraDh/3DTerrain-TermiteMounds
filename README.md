# 3DTerrain-TermiteMounds

Geometric approach to filters termite mounds from Unmanned Aerial Vehicle - LiDAR data.

main script: filter_termites_CC.py

PACKAGES NEED TO BE INSTALLED:
   os, numpy, pandas, sklearn, matplotlib, laspy, copy, time

STEPS:
1) Load point cloud into CloudCompare; use 'Cloth Simulation Filter' plugin (settings: filter Scene: Steep slope, cloth resolution: 1.0). Save ground result as .las file.

2) Open 'filter_termites_CC.py' Go to the section 'STEP 2'. Change the 'input_file_CSF' to the filtered point cloud in .las format (result of the cloth simulation filter (CloudCompare)). Run section 2.
3) CALCULATE NORMALS in CloudCompare (quadric, radius = .75); convert to DIP DIRECTION/DIP DEGREE. Save the point cloud as GF_[name].txt                                                      
4) Run section 'STEP 4' in python script 'filter_termites_CC.py': filter based on dip direction/dip degree

5) CALCULATE NORMALS in Cloudcompare (quadric, radius = .30) of the 'csf_to_calculate_normals.txt' file; convert to DIP DIRECTION/DIP DEG. Save it as GF_csf_to_calculate_normals.txt

6) Run section 'STEP 6' in python script 'filter_termites_CC.py': extra filter based on dip dip degree


Extra information

step 2: file 'keep_bottom_points'. Based on the X Y values, the nearest neighbours (default: 30) are identified for each point. Only the point with the lowest Z value is kept.

step 4: file 'mounddetection1'. 
    1) data is filtered based on dip degree (7-86 degrees)
    2) data is filtered: points with less than 50 points within a 2 m radius are discarded
    3) data is clustered using agglomerative clustering
    4) extra points are added to each cluster, from the 'low_points' filtered data
    5) Based on the dip direction, certain clusters are discarded. Creates a file 'resultsoptimizingcenter.csv' which lists the RMSE before and optimizing the termite mound center, and the distance between the initial center and optimized. Creates a folder 'plots' in which the optimized center is visualised. Creates a file 'temporaryresultsmounds.txt'. 
    6) extra points are added to each cluster, coming from the 'CSF' data

step 6: file 'mounddetection2'. Clusters are filtered based on dip degree. All clusters that have > 77 % amount of points (above 15 cm ground surface) that are steeper than 77◦ are discarded.
