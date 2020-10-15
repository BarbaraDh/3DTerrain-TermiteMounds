# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 11:32:27 2020
@author: bgdhont
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.neighbors import KDTree
import time
import laspy


def clustering(filtered, th_clust):
    """
    Clusters data using Agglomerative Clustering. Distance Threshold: 30 default

    """
    from sklearn.cluster import AgglomerativeClustering as AC
    agC = AC(n_clusters=None, distance_threshold=th_clust, memory=None)
    X = np.array(filtered.iloc[:, :3])
    agC.fit(X)
    labels = np.array([agC.labels_]).T
    #res = np.concatenate((X, np.array([agC.labels_]).T), axis = 1)
    amountclusters = len(set(agC.labels_))
    print('Amount of clusters: ' + str(amountclusters))
    return labels, amountclusters


def filter_density(hills, th_NN):
    """
    If less than th_NN (default: 50) points within a radius of 2m, the point is discarded.

    """
    idx_filtered = []
    i, rad, am = 0, 2, 50
    arr = np.array(hills.iloc[:, 0:3])
    tree = KDTree(arr)
    while i < len(arr):
        amount = tree.query_radius(arr[[i]], r=rad, count_only=True)
        if amount > am:
            idx_filtered.append(i)
        i = i + 1
    filtered = hills.iloc[idx_filtered, :].copy()
    return filtered

def add_all_points(filtered, data):
    """
    Adds extra points to the cluster that might have been discarded in the previous step.
    Those points exist mostly of the flat top of the termite mounds.
    Parameters
    ----------
    filtered : clustered data that was kept in previous step
    data : low points file
    """
    k = list(set(filtered['label'])) 
    #what if we add points to the clusters?
    for i in k:
        cluster = filtered.loc[filtered['label'] == i, :].copy()
        if len(cluster) > 15:    
            x0, y0 = [np.mean(cluster.iloc[:, 0]), np.mean(cluster.iloc[:, 1])] #f.max(axis=1)[0][2]
            minx=np.min(cluster.iloc[:,0])
            maxx=np.max(cluster.iloc[:,0])
            miny=np.min(cluster.iloc[:,1])
            maxy= np.max(cluster.iloc[:,1])
            if ((maxy-miny)/(maxx-minx) <3/2) & ((maxx-minx)/(maxy-miny) <3/2):
                rad = np.mean([(maxx - minx),(maxy-miny)])/2
                deel = data.loc[(data.iloc[:,0] <maxx) & (data.iloc[:,0] > minx) &
                          (data.iloc[:,1] >miny) & (data.iloc[:,1] < maxy) ,:].copy() 
                deel['label'] = i
                arr = np.array(deel.iloc[:, 0:2])
                tree = KDTree(arr)
                ind  = tree.query_radius([[x0,y0]], r=rad)
                to_append = deel.iloc[ind[0],:]
                filtered = pd.concat([filtered, to_append])
        return filtered
    
def overview_clusters(filtered, output_dir):
    """
    Creates a map of the locations of the clusters.

    """
    labels = filtered.loc[:,'label']
    x_means = []
    y_means = []
    z_means = []
    for i in np.unique(labels):
        x_t = np.mean(filtered.loc[filtered['label'] == i, '//X']) 
        y_t = np.mean(filtered.loc[filtered['label'] == i, 'Y']) 
        z_t = np.mean(filtered.loc[filtered['label'] == i, 'Z']) 
        x_means.append(x_t)
        y_means.append(y_t)
        z_means.append(z_t)
    
    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot()
    ax.scatter(x_means, y_means,c ='k', marker = '^');
    for i, txt in enumerate(np.unique(labels)):
        ax.annotate(txt, (x_means[i], y_means[i]), fontsize = 8)
    #ax.scatter(x, y,c ='b', marker = '.');
    ax.set_xlabel( 'X (m)')
    ax.set_ylabel( 'Y (m)')
    plt.grid()
    plt.savefig( output_dir  +  '/all_clusters.png', dpi = 800)
    
def optimize_center(filtered, output_dir, th_distance, th_RMSE, extra_info):
    """
    filters mounds based on dip degree and dip direction.
    
    Creates a file 'resultsoptimizingcenter.csv' which lists the RMSE before
    and optimizing the termite mound center, 
    and the distance between the initial center and optimized.
    
    Creates a folder 'plots' in which the optimized center is visualised.
    
    Creates a file 'temporaryresultsmounds.txt'. 
    """
    from scipy import optimize
    overview = pd.DataFrame({})
    k = list(set(filtered['label']))  # range(amountclusters)
    #k = [1]
    def angle_between(x0,y0): #p2 blijft vast 
        ang1 = np.arctan2(0,1) #y is 1, x is 0
        ang2 = np.arctan2(p2[:,0]-x0, p2[:,1]-y0)
        return  (ang2 - ang1)  % (2 * np.pi) #np.rad2deg kan je nog extra gebruiken ook
    
    def minimize(pars):
        dif_angle = np.array(angle_between(*pars) - angle)
        ids1 =  (dif_angle > 2*np.pi - np.pi/4) 
        ids2 =  (dif_angle <  (-2*np.pi +  np.pi/4))
        if len(set(ids1)) > 1:
            dif_angle[ids1] = dif_angle[ids1] - 2*np.pi
        if len(set(ids2)) > 1:
            dif_angle[ids2] = dif_angle[ids2] + 2*np.pi
        return list(dif_angle) 
    
    def calculateDistance(x1,y1,x2,y2):  
         dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)  
         return dist  
     
    for i in k:         
        cluster = filtered.loc[filtered['label'] == i, :].copy()
        if len(cluster) > 15:    
            x0, y0, z0 = [np.mean(cluster.iloc[:, 0]), np.mean(cluster.iloc[:, 1]),np.mean(cluster.iloc[:, 2])] #f.max(axis=1)[0][2]
            x0_init = np.mean(cluster.iloc[:, 0])
            y0_init = np.mean(cluster.iloc[:, 1])
            x = np.array(cluster.loc[:, '//X'])
            y = np.array(cluster.loc[:, 'Y'])
            z = np.array(cluster.loc[:, 'Z'])
            p2= np.array(cluster.loc[:, ['//X','Y']])
            cluster.loc[:,'Dip direction (rad)'] =  cluster.loc[:,'Dip direction (degrees)']/180*np.pi
            angle = np.array(cluster.loc[:,'Dip direction (rad)'])
            initial_pars = [x0, y0]
            pars = initial_pars
            pars2 = initial_pars
            nr_try = 0
            x_c = [0, -.2, -.2, .2, .2]
            y_c = [0, -.2, .2, -.2, .2]
            RMSEs = []
            while (nr_try < 5): #(x0 == initial_pars[0]) & (y0 == initial_pars[1]) & 
                initial_pars = [x0_init, y0_init]
                pars2[0] = initial_pars[0] + x_c[nr_try]
                pars2[1] = initial_pars[1] + y_c[nr_try]
                pars, ier = optimize.leastsq(minimize, pars2)
                RMSEs.append( np.mean((np.array(minimize(pars)))**2))
                nr_try = nr_try + 1
            if abs(min(RMSEs) - RMSEs[0]) <.005:
                idr = 0
            else:
                idr = RMSEs.index(min(RMSEs))
            initial_pars = [x0_init, y0_init]    
            pars2[0] = initial_pars[0] + x_c[idr]
            pars2[1] = initial_pars[1] + y_c[idr]
            pars, ier = optimize.leastsq(minimize, pars2)
            x0, y0 = pars
            cluster.loc[:,'based_on_mean'] =  angle_between(x0,y0)
            cols = cluster.loc[:,'Dip direction (rad)']
            levels = np.arange(0, .1, max(cols))
            if extra_info == True:
                output_dir_plots = output_dir + '/plots/plot' + str(i) + '.png' 
                if not os.path.exists(output_dir + '/plots'):
                    os.makedirs(output_dir + '/plots')
                fig = plt.figure(figsize=(5,4), num = str(i))
                ax = fig.add_subplot()
                plt.title(str(i))
                plt.tricontour(x, y, z, 10, linewidths=0.5, colors='grey')
                cmap = mpl.cm.viridis
                pl = ax.scatter(x, y,s=11, c = cols, vmin = min(cols), vmax = max(cols), cmap = cmap, edgecolor = 'k', linewidth=.4);
                pl3 = ax.scatter(x0, y0 ,s=20, c = 'red', edgecolor = 'k', linewidth=2, marker = 'x', label = 'optimized center');
                pl4 = ax.scatter(initial_pars[0], initial_pars[1] ,s=20, c = 'darkred', edgecolor = 'k', linewidth=2, marker = 'x', label = 'initial center');
                ax.set_xlabel('x (m)')
                ax.set_ylabel('y (m)')
                cbar = fig.colorbar(pl, shrink=1, aspect=15, ticks=[min(cols), np.pi/2, np.pi, 3/2*np.pi, max(cols)])
                cbar.ax.set_yticklabels(['0', '90', '180', '270', '360'])
                cbar.set_label('Azimuth angle (°)', fontsize = 8)
                ax.axis('equal')
                plt.legend()
                plt.savefig(output_dir_plots, dpi = 500)
            initial_pars = [x0_init, y0_init] 
            RMSE_before = np.mean((np.array(minimize(initial_pars)))**2)
            RMSE_opt =  np.mean((np.array(minimize(pars)))**2)
            distance = calculateDistance(initial_pars[0], initial_pars[1], x0,y0)
           
            overview.loc[i, 'RMSE_before'] = RMSE_before
            overview.loc[i, 'RMSE_opt'] = RMSE_opt
            overview.loc[i,'Distance'] = distance
            overview.loc[i,'id'] = i
            filtered.loc[filtered['label'] == i, 'RMSE'] = RMSE_opt
            filtered.loc[filtered['label'] == i, 'mound'] = (distance < th_distance ) & (RMSE_opt < th_RMSE)   
            filtered.loc[filtered['label'] == i, 'x0'] = x0
            filtered.loc[filtered['label'] == i, 'y0'] = y0
            
        else:
            overview.loc[i, 'RMSE_before'] = 10
            overview.loc[i, 'RMSE_opt'] = 10
            overview.loc[i,'Distance'] = 10
            overview.loc[i,'id'] = i
            filtered.loc[filtered['label'] == i, 'RMSE'] = 10
            filtered.loc[filtered['label'] == i, 'mound'] = False
            
    filtered.loc[:, 'mound'] = filtered.loc[:, 'mound'].astype('int')
    if extra_info == True:
        filtered.to_csv(output_dir + '/temporaryresultsmounds.txt', index=False)
        overview.to_csv(output_dir + '/resultsoptimizingcenter.csv')
    return filtered
    
def add_csf_points(filtered, input_file_CSF, output_dir):
    """
    adds points from the results of the 'Cloth Simulation Filter', and saves the file
    """
    #input_CSF = r'C:/Users/bgdhont/Documents/Data/RIEGL/LR_4PARTS/Litchfield_F2_UTM52S_without_sl15_16_thinned_CSF_p1.las'
    inFile = laspy.file.File(input_file_CSF, mode = "r")
    csf = pd.DataFrame(np.array((inFile.x,inFile.y,inFile.z)).transpose(), columns = ['//X', 'Y','Z']) #np.vstack([inFile.x, inFile.y, inFile.z]).transpose()
    dat = filtered.loc[filtered['mound']== True, :].copy()
    dat['csf'] = 0
    csf['csf'] = 1
    k = list(set(dat['label']))  
    
    for i in k:
        cluster = dat.loc[dat['label'] == i, :].copy()
        if len(cluster) > 15:    
            x0, y0 = [np.mean(cluster.iloc[:, 0]), np.mean(cluster.iloc[:, 1])] #f.max(axis=1)[0][2]
            minx=np.min(cluster.iloc[:,0])
            maxx=np.max(cluster.iloc[:,0])
            miny=np.min(cluster.iloc[:,1])
            maxy= np.max(cluster.iloc[:,1])
            rad = np.mean([(maxx - minx),(maxy-miny)])/2
            deel = csf.loc[(csf.iloc[:,0] <maxx) & (csf.iloc[:,0] > minx) &
                          (csf.iloc[:,1] >miny) & (csf.iloc[:,1] < maxy) ,:].copy() 
            deel['label'] = i
            arr = np.array(deel.iloc[:, 0:2])
            tree = KDTree(arr)
            ind  = tree.query_radius([[x0,y0]], r=rad)
            to_append = deel.iloc[ind[0],:]
            dat = pd.concat([dat, to_append])
    dat.to_csv(output_dir + '/csf_to_calculate_normals.txt', index=False, columns = ['//X', 'Y','Z', 'label','csf'])
    
    
def filter1(input_file_GF, input_file_CSF, output_dir, th_dip_degree, th_NN, th_clust, th_distance, th_RMSE, extra_info):
    """ 
    main function.
    1) data is filtered based on dip degree (7-86 degrees)
    2) data is filtered: points with less than 50 points within a 2 m radius are discarded
    3) data is clustered using agglomerative clustering
    4) extra points are added to each cluster, from the 'low_points' filtered data
    5) Based on the dip direction, certain clusters are discarded.
    6) extra points are added to each cluster, coming from the 'CSF' data
    """ 
    t0 = time.process_time()
    print('start')
    data = pd.read_csv(input_file_GF, delimiter=',')
    th1 = th_dip_degree[0]
    th2 = th_dip_degree[1]
    hills = data.loc[(data.loc[:, 'Dip (degrees)'] > th1) & (data.loc[:, 'Dip (degrees)'] < th2), :].copy() #1
    filtered = filter_density(hills,th_NN) #2
    filtered['label'], amountclusters = clustering(filtered,th_clust) #3
    
    filtered = add_all_points(filtered, data) #4
    overview_clusters(filtered, output_dir)
    filtered = optimize_center(filtered,output_dir, th_distance, th_RMSE, extra_info) #5
    add_csf_points(filtered, input_file_CSF,output_dir) #6
    t1 = time.process_time()
    print("Time spent: " + str(t1 - t0) + ' s' )
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = 'After CSF/lowpoints/Normals/dipdirectiondegrees(rad.75) ')
    parser.add_argument("input_GF", help = 'lowpoints with dip direction/degrees from CloudCompare (radius: .75)')
    parser.add_argument("input_CSF", help = 'CSF file from cloudcompare')
    parser.add_argument("output_dir", help= 'directory for the results')
    parser.add_argument("--th_dip_degree", help = 'default: (7,86)')
    parser.add_argument("--th_NN", help = 'Default: 50')
    parser.add_argument("--th_clust", help = 'Default: 30')
    parser.add_argument("--th_distance", help = 'Default: .75')
    parser.add_argument("--th_RMSE", help = 'Default: 1')
    parser.add_argument("--extra_info", help = 'NN. Default: False')
    args = parser.parse_args()
    input_file_GF = args.input_GF
    input_file_CSF = args.input_CSF
    output_dir = args.output_dir
    
    if args.th_dip_degree:
        th_dip_degree = int(args.th_dip_degree)
    else:
        neighbours = (7,87)
    if args.th_NN:
        th_NN = int(args.th_NN)
    else:
        th_NN = 50
    if args.th_clust:
        th_clust = int(args.th_clust)
    else:
        th_clust = 50   
    if args.th_distance:
        th_distance = int(args.th_distance)
    else:
        th_distance = .75
    if args.th_RMSE:
        th_RMSE = int(args.th_RMSE)
    else:
        th_RMSE = 1
    if args.extra_info:
        extra_info = bool(args.extra_info)
    else:
        extra_info = False
    filter1(input_file_GF, input_file_CSF, output_dir, th_dip_degree, 
            th_NN,th_clust, th_distance, th_RMSE, extra_info)
        
    