# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 11:22:16 2021

@author: Ejdrup
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 16:04:26 2021

@author: Ejdrup
"""

import pandas as pd
import numpy as np
import os, glob
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.spatial import cKDTree
import matplotlib.path as path
import seaborn as sns
# import storm_analysis.sa_utilities.hdf5_to_image as h5_image
from sklearn.neighbors import KernelDensity
from scipy import stats
from matplotlib.path import Path
from shapely.geometry import Polygon
import cv2 as cv2
from itertools import compress
from scipy.optimize import linprog
from scipy.spatial import Delaunay

def get_neighbours_list(voronoi):
    neighbours_list = [[] for point_index in  range(len(voronoi.points))]
    counter = 0
    for count, neighbour_pair in enumerate(voronoi.ridge_points):
        if count%int((len(voronoi.ridge_points)-1)/5) == 0:
            #print('Looping every voronoi region - progress: {}%'.format(20*counter))
            counter += 1
        neighbours_list[neighbour_pair[0]].append(neighbour_pair[1])
        neighbours_list[neighbour_pair[1]].append(neighbour_pair[0])
    return neighbours_list

def REp(ch1_points, ch2_points):

    # Compute the voronoi tessellation for each channel
    print("Computing 2D voronoi tessellation.")
    ch1_vor = Voronoi(ch1_points)
    ch2_vor = Voronoi(ch2_points)

    
    unfiltered_region = [ch1_vor.regions[i] for i in ch1_vor.point_region]
    bool_index = np.ones((len(unfiltered_region)), dtype=bool)
    for i in range(0,len(bool_index)):
        if -1 in unfiltered_region[i]:
            bool_index[i] = False
    sorted_region = list(compress(unfiltered_region, bool_index))
    sorted_vertices = [ch1_vor.vertices[i] for i in sorted_region]
    sorted_region_path = [path.Path(i) for i in sorted_vertices]
    sorted_region_area = [Polygon(i).area for i in sorted_vertices]
    sorted_points_ch1 = ch1_vor.points[bool_index]
    
    ckd_tree_ch2 = cKDTree(ch2_points) # Construct a knn tree
    # Figure out how to circumvent k = max
    if len(ch2_vor.points) < 1000:
            dist_ch2, idx_ch2 = ckd_tree_ch2.query(sorted_points_ch1, k=len(ch2_vor.points), n_jobs = -1) #query closest neighbor
            n_points_ch2_in_region = [np.count_nonzero(i.contains_points(ch2_vor.points[idx_ch2[j,:],:])) for j, i in enumerate(sorted_region_path)]
    else:
        dist_ch2, idx_ch2 = ckd_tree_ch2.query(sorted_points_ch1, k=1000, n_jobs = -1) #query closest neighbor
        n_points_ch2_in_region = [np.count_nonzero(i.contains_points(ch2_vor.points[idx_ch2[j,:],:])) for j, i in enumerate(sorted_region_path)]
        
    # Construct dictionary of points in adjacent regions for each point
    #print("Constructing dictionary of adjacent regions.")
    
    # Compute first order mean distance
    neighbours_list = get_neighbours_list(ch1_vor)
    neighbour_indices = list(compress(neighbours_list, bool_index))
    neighbour_points = [ch1_vor.points[i] for i in neighbour_indices]
    
    # Compute mean ditance to first order neighbors
    neighbour_distances = [np.sqrt(np.sum((np.array(j)-sorted_points_ch1[i])**2,axis=1)) for i, j in enumerate(neighbour_points)]
    first_order_mean_distance = [np.nanmean(i) for i in neighbour_distances]
    
    return np.array(n_points_ch2_in_region), np.array(sorted_region_area), np.array(first_order_mean_distance), bool_index

def in_hull(hull, p):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    https://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
    """
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return sum(hull.find_simplex(p)>=0)
            

def REp3D(ch1_points, ch2_points):

    # Compute the voronoi tessellation for each channel
    print("Computing 3D voronoi tessellation.")
    ch1_vor = Voronoi(ch1_points)
    ch2_vor = Voronoi(ch2_points)
    
    # # Construct dictionary of points in adjacent regions for each point
    # print("Constructing dictionary of adjacent regions.")
    # neighbours_dict = get_neighbours_dict(ch1_vor)
    # # ckd_tree_DAT = cKDTree(ch1_points) # Construct a knn tree of DAT molecules
    # # dist_DAT, idx_DAT = ckd_tree_DAT.query(ch1_points, k=4, n_jobs = -1) #query all the STX molecules for closest neighbor
    # # mean_neighbor = np.mean(dist_DAT[:,1:],axis=1)
    
    
    unfiltered_region = [ch1_vor.regions[i] for i in ch1_vor.point_region]
    
    bool_index = np.ones((len(unfiltered_region)), dtype=bool)
    for i in range(0,len(bool_index)):
        if -1 in unfiltered_region[i]:
            bool_index[i] = False
    sorted_region = list(compress(unfiltered_region, bool_index))
    
    sorted_vertices = [ch1_vor.vertices[i] for i in sorted_region]
    sorted_region_hull = [ConvexHull(i) for i in sorted_vertices]
    sorted_region_volume = [i.volume for i in sorted_region_hull]
    sorted_points_ch1 = ch1_vor.points[bool_index]
    
    # Construct and query a knn tree
    print("Constructing a knn tree")
    ckd_tree_ch2 = cKDTree(ch2_points)
    dist_ch2, idx_ch2 = ckd_tree_ch2.query(ch1_vor.points[bool_index], k=100, n_jobs = -1)
    
    # Compare the old and new hull for differences
    print("Comparing hulls")
    n_points_ch2_in_region = [in_hull(i.points,ch2_vor.points[idx_ch2[j,:],:]) for j, i in enumerate(sorted_region_hull)]
    
    # Compute first order mean distance
    print("Identifying neighboring points")
    neighbours_list = get_neighbours_list(ch1_vor)
    neighbour_indices = list(compress(neighbours_list, bool_index))
    neighbour_points = [ch1_vor.points[i] for i in neighbour_indices]
    
    # Compute mean distance to first order neighbors
    print("Computing mean distance to neighbors")
    neighbour_distances = [np.sqrt(np.sum((np.array(j)-sorted_points_ch1[i])**2,axis=1)) for i, j in enumerate(neighbour_points)]
    first_order_mean_distance = [np.nanmean(i) for i in neighbour_distances]
    
    return np.array(n_points_ch2_in_region), np.array(sorted_region_volume), np.array(first_order_mean_distance), bool_index

def bin_REp(n_points, areas_ch1, first_ord_dist, max_dist, step_size, total_volume = "None"):
    # Try to bin
    if total_volume is "None":
        threshold = np.percentile(areas_ch1,99.5)
        total_volume = np.sum(areas_ch1[areas_ch1 < threshold])
    sorted_density_ch1 = np.sort(1/areas_ch1)
    sorted_density_ch1_idx = np.argsort(1/areas_ch1)
    first_ord_dist_sorted = first_ord_dist[sorted_density_ch1_idx]
    maxbin = int(max_dist/step_size)
    dr = step_size
    density = sum(n_points)/total_volume
    pts_ratio = np.zeros((maxbin+1,1))
    no_regions = np.zeros((maxbin+1,1))
    no_loc_per_region = np.zeros((maxbin+1,1))
    area_per_bins = np.zeros((maxbin+1,1))
    for i in range(0,len(sorted_density_ch1)):
        if np.isfinite(first_ord_dist_sorted[i]):
            bin_no = int(np.log10(first_ord_dist_sorted[i])/dr)
            if (bin_no <= maxbin):
                obs_dens = n_points[sorted_density_ch1_idx[i]]/(1/sorted_density_ch1[i]*density)
                pts_ratio[bin_no] += obs_dens
                no_regions[bin_no] += 1
                no_loc_per_region[bin_no] += n_points[sorted_density_ch1_idx[i]]
                area_per_bins[bin_no] += areas_ch1[sorted_density_ch1_idx[i]]
                
    return np.divide(pts_ratio,no_regions), no_regions, no_loc_per_region, area_per_bins

def bin_REp_area(n_points, areas_ch1, first_ord_dist, max_dist, step_size, size_threshold = "None", total_volume = "None"):
    # Try to bin
    if total_volume is "None":
        threshold = np.percentile(areas_ch1,99.5)
        total_volume = np.sum(areas_ch1[areas_ch1 < threshold])
    sorted_density_ch1 = np.sort(1/areas_ch1)
    sorted_density_ch1_idx = np.argsort(1/areas_ch1)
    first_ord_dist_sorted = first_ord_dist[sorted_density_ch1_idx]
    area_sorted = areas_ch1[sorted_density_ch1_idx]
    maxbin = int(max_dist/step_size)
    dr = step_size
    density = sum(n_points)/total_volume
    pts_ratio = np.zeros((maxbin+1,1))
    no_regions = np.zeros((maxbin+1,1))
    no_loc_per_region = np.zeros((maxbin+1,1))
    area_per_bins = np.zeros((maxbin+1,1))
    for i in range(0,len(sorted_density_ch1)):
        if np.isfinite(first_ord_dist_sorted[i]):
            bin_no = int(np.log10(area_sorted[i])/dr)
            if (bin_no <= maxbin):
                obs_dens = n_points[sorted_density_ch1_idx[i]]/(1/sorted_density_ch1[i]*density)
                pts_ratio[bin_no] += obs_dens
                no_regions[bin_no] += 1
                no_loc_per_region[bin_no] += n_points[sorted_density_ch1_idx[i]]
                area_per_bins[bin_no] += areas_ch1[sorted_density_ch1_idx[i]]
                
    return np.divide(pts_ratio,no_regions), no_regions, no_loc_per_region, area_per_bins