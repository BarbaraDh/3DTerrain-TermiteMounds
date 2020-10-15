# 3DTerrain-TermiteMounds

Geometric approach to filters termite mounds from Unmanned Aerial Vehicle - LiDAR data.

main script: filter_termites_CC.py

PACKAGES NEED TO BE INSTALLED:
   os, numpy, pandas, sklearn, matplotlib, laspy, copy, time

STEPS:

1) Load data into CloudCompare; use 'Cloth Simulation Filter' 

    SETTINGS:
    
            filter Scene: Steep slope
            
            cloth resolution: 1.0
            
            max. iterations: 500
            
            classification threshold: 1.
            
    
    OUTPUTS: Save as [name].las
    
    EXPLANATION: ground filtering algorithm, should comprise the complete termite mounds. Under story can be included.
    
    https://www.cloudcompare.org/doc/wiki/index.php?title=CSF_(plugin)
    
    
2) python: keep lowest points.

    INPUT: input_file_CSF: filtered point cloud in .las format (result of the cloth simulation filter (CloudCompare)
    
    PARAMETERS:
    
            neighbours: default 30
            
    OUTPUT:low_points_knn[neighbours]_[name].las
    
    EXPLANATION: For each point, neighbours are taken into consideration and only the lowest point remains. Duplicates are removed.
        
        
3) Calculate normals in CloudCompare, can be done via python

    INPUT: low_points_knn[neighbours]_[name].las
    
    SETTINGS: local surface model: Quadric; radius: default 0.75 m; Orientation: Use mininum Spanning tree - knn = 6
    
    convert to DIP DIRECTION/DIP DEGREE. 
    
    OUTPUT: Save as GF_[name].txt  
    
    EXPLANATION: https://www.cloudcompare.org/doc/wiki/index.php?title=Normals%5CConvert_to_Dip_and_Dip_direction_SFs
    
                                                   
4) python: filter based on dip direction/dip degree

    INPUT:  GF_low_points_knn[neighbours]_[name].txt 
    
    PARAMETERS:
    
            th_dip_degree: default (7,87) (degrees)
            
            th_NN: default 50 
            
            th_clust: default 30
            
            th_distance: default .75 m 
            
            th_RMSE: default 1
            
            extra_info: default False
            
    OUTPUT: csf_to_calculate_normals.txt
    
            
    EXPLANATION: 
    
                1. points are filtered: all points between th_dip_degree remain
    
                2. only points that have > th_NN points within 2 meter remain
                
                3. Data is clustered using agglomerative clustering (distance thershold = th_clust)
                
                4. Cluster is checked if it has a conical shape: center is optimized based on dip direction. 
                
                If optimized center point is > th_distance from mean (x,y), or RMSE of the optimalisation > th_RMSE, cluster is discarded.
            
5)  Calculate normals in CloudCompare, can be done via python

    INPUT: csf_to_calculate_normals.txt
    
    SETTINGS: local surface model: Quadric; radius: default 0.30 m; Orientation: Use mininum Spanning tree - knn = 6
    
    convert to DIP DIRECTION/DIP DEGREE.    
    
    OUTPUT: Save as GF_csf_to_calculate_normals.txt
    
    EXPLANATION: https://www.cloudcompare.org/doc/wiki/index.php?title=Normals%5CConvert_to_Dip_and_Dip_direction_SFs
    
    
6) python: extra filter based on dip degree. 

    INPUT: GF_csf_to_calculate_normals.txt
    
    PARAMETERS:
    
            th1: default 77 (%)
            
            th2: default 77 (degrees)
            
            n: default 0.15 (m)
            
    OUTPUT: 
    
             termite_mounds.txt: labelled termite mounds
    
             termite_mounds.png: map of the termite mounds
             
    EXPLANATION: all cluster of which > th1 % of points (above n m ground surface) have a dip degree > th2° are discarded.
            
                                                                       
Extra information

step 2: file 'keep_bottom_points'. Based on the X Y values, the nearest neighbours (default: 30) are identified for each point. Only the point with the lowest Z value is kept.

step 4: file 'mounddetection1'. 

    1) data is filtered based on dip degree (default 7-86 degrees)
    
    2) data is filtered: points with less than (default) 50 points within a 2 m radius are discarded
    
    3) data is clustered using agglomerative clustering
    
    4) extra points are added to each cluster, from the 'low_points' filtered data
    
    5) Based on the dip direction, certain clusters are discarded. Creates a file 'resultsoptimizingcenter.csv' which lists the RMSE before and optimizing the termite mound center, and the distance between the initial center and optimized. Creates a folder 'plots' in which the optimized center is visualised. Creates a file 'temporaryresultsmounds.txt'. 
    
    6) extra points are added to each cluster, coming from the 'CSF' data

step 6: file 'mounddetection2'. Clusters are filtered based on dip degree. All clusters that have > 77 % amount of points (above 15 cm ground surface) that are steeper than 77◦ are discarded (default).

