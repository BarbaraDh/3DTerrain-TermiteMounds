# 3DTerrain-TermiteMounds

filter_termites_CC.py

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
