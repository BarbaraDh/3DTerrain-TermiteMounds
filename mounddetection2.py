# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 12:09:47 2020
@author: bgdhont
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time


def overview_clusters(filtered, output_dir):
    """
    creates a map with the locations of the termite mounds 
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
    ax.set_xlabel( 'X (m)')
    ax.set_ylabel( 'Y (m)')
    plt.grid()
    plt.savefig( output_dir  +  '/termite_mounds.png', dpi = 800)
    
def filter2(input_file, output_dir,th1,th2,n):
    """
    main function.
    Clusters are filtered based on dip degree.
    all clusters that have > th1 % amount of points (above n m ground surface) that are steeper than th2° are discarded.
    
    """
    t0 = time.process_time()
    print('start')
    data = pd.read_csv(input_file, delimiter=',')
    k = list(set(data['label'])) #'Scalar field' 
    for i in k:
        cluster = data.loc[data['label'] == i, :].copy()
        cluster['steep'] = (cluster.loc[:, 'Dip (degrees)'] > th1)
        uppercluster = cluster.loc[cluster['Z'] > (cluster['Z'].min() + n),:]
        percentage = uppercluster.steep.value_counts()[0]/ uppercluster.steep.count()
        data.loc[data['label'] == i, 'percentage'] = percentage
        mounds = data.loc[(data['percentage'] >th2/100), :]
    mounds.to_csv(output_dir + '/termite_mounds.txt', index=False)
    overview_clusters(mounds, output_dir)
    t1 = time.process_time()
    print("Time spent: " + str(t1 - t0) + ' s' )
    
    
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description = 'last step')
    parser.add_argument("input_CSF_GF", help = 'lowpoints with dip direction/degrees from CloudCompare (radius: .3)')
    parser.add_argument("output_dir", help= 'directory where results will be saved')
    parser.add_argument("--th1", help = 'th1. Default: 77 (%)')
    parser.add_argument("--th2", help = 'th2. Default: 77 (degrees)')
    parser.add_argument("--n", help = 'n. Default: 0.15 (m)')
    args = parser.parse_args()
    input_file = args.input_CSF_GF
    output_dir = args.output_dir
    if args.th1:
        th1 = int(args.th1)
    else:
        th1 = 77
    if args.th2:
        th2 = int(args.th2)
    else:
        th2 = 77
    if args.n:
        n = int(args.n)
    else: n = 0.15
    filter2(input_file, output_dir, th1,th2, n)
        
    
    