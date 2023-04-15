# coding = utf-8


'''
An example of plot the distance result in the paper
By angdi  and bing

'''


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.ndimage
import copy
from matplotlib.image import NonUniformImage
from scipy.interpolate import griddata 
import time


def vis_actin_surround_actin(actin_vect_map_lst, file_idx, actin_num, point_num, ):
    
    points_x = []
    points_y = []

    for idx, actin_vects in enumerate(actin_vect_map_lst):
        if len(actin_vects) >0:
            for vect in actin_vects:
                points_x.append(vect[0])
                points_y.append(vect[1])
                
    bins_ = 50
    fig = plt.figure(figsize=(8,6), dpi=1200)
    hist, xedges, yedges = np.histogram2d(points_x, points_y, bins=bins_, range =[[-1000,1000],[-1000,1000]] )

    values = hist.reshape(-1, 1)
    points = []
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            points.append([i,j])

    grid_x, grid_y = np.mgrid[1:hist.shape[0]:1000j, 1:hist.shape[1]:1000j]

    grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear',fill_value=0)

    grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.3) )] = np.nan
    grid_z0[np.where(grid_z0 < (np.max(grid_z0) * 0.45) )] = 0
    

    plt.imshow(grid_z0,interpolation='bilinear', extent=(-100,100,-100,100),cmap= 'RdYlBu_r')
    plt.title(f'{file_idx} Actin surrounding probability density ditribution testClim \n climmax = {np.max(grid_z0)}, actin num {actin_num}, points num {point_num}' )
    

    plt.colorbar()

