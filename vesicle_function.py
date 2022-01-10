# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:21:27 2021

@author: Ejdrup
"""
import numpy as np

def generate_spherical(npoints, radius = 50, ndim = 3):
    vec = np.random.randn(ndim, npoints)
    # vec = np.linalg.norm(vec, axis=0)
    vec[0,:] = vec[0,:]*radius + np.random.randn(npoints)*10
    vec[1,:] = vec[1,:]*radius + np.random.randn(npoints)*10
    vec[2,:] = vec[2,:]*radius + np.random.randn(npoints)*25
    return vec

def generate_circle(npoints, radius = 50, zshift = -75):
    coords = np.zeros((3, npoints))
    r = radius * np.sqrt(np.random.random(npoints))
    theta = np.random.random(npoints) * 2 * 3.14
    coords[0,:] = r * np.cos(theta) + np.random.randn(npoints)*10
    coords[1,:] = r * np.sin(theta) + np.random.randn(npoints)*10
    coords[2,:] = np.random.randn(npoints)*25 + zshift
    return coords
   
def generate_vesicle_content(npoints, radius, ndim = 3):
    npoints = int(npoints * (1/0.524)) #Volume of sphere/cube
    coords = np.zeros((3, npoints))
    coords[0,:] = np.random.random(npoints)-0.5
    coords[1,:] = np.random.random(npoints)-0.5
    coords[2,:] = np.random.random(npoints)-0.5
    distance = np.sqrt((coords[0,:])**2 + (coords[1,:])**2 + (coords[2,:])**2)
    coords = coords[:,distance < 0.5]
    coords[0,:] = coords[0,:] + np.random.randn(coords.shape[1])*10
    coords[1,:] = coords[1,:] + np.random.randn(coords.shape[1])*10
    coords[2,:] = coords[2,:] + np.random.randn(coords.shape[1])*25
    return coords 

def vesicle_simulation(pseudo):
    # Simulate position on non-docked vesicles
    no_non_docked = 20
    non_docked = np.zeros((3,no_non_docked))
    non_docked[0,:] = np.random.random(no_non_docked) * 800 - 400
    non_docked[1,:] = np.random.random(no_non_docked) * 800 - 400
    non_docked[2,:] = np.random.random(no_non_docked) * 200
    
    # Simulate position of docked vesicle
    no_docked = 20
    docked = np.zeros((3,no_docked))
    docked[0,:] = np.random.random(no_docked) * 800 - 400
    docked[1,:] = np.random.random(no_docked) * 800 - 400
    
    # Gather both populations
    all_vesicles = np.hstack((non_docked,docked))
    
    # Visualize the locations of the vesicle centers
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d') 
    # ax.scatter(docked[0,:], docked[1,:], docked[2,:], s=10, c='r')
    # ax.scatter(non_docked[0,:], non_docked[1,:], non_docked[2,:], s=10, c='b')
    
    # Simulate individual localizations on vesicle pool
    vesicle_coords = np.zeros((3,1))
    for i in range(all_vesicles.shape[1]):
        new_vesicle = generate_spherical(100, radius = 30) + all_vesicles[:,i][:,None]
        vesicle_coords = np.hstack((vesicle_coords,new_vesicle))
    vesicle_coords = vesicle_coords[:,1:]
    
    # Simulate individual localizations on SNAREs
    
    SNARE_coords = np.zeros((3,1))    
    for i in range(docked.shape[1]):
        new_SNARE = generate_circle(50, radius = 30) + docked[:,i][:,None]
        SNARE_coords = np.hstack((SNARE_coords,new_SNARE))
    SNARE_coords = SNARE_coords[:,1:]
    
    # Simulate individual localizations of content
    
    content_coords = np.zeros((3,1))    
    for i in range(all_vesicles.shape[1]):
        new_content = generate_vesicle_content(50, radius = 20) + all_vesicles[:,i][:,None]
        content_coords = np.hstack((content_coords,new_content))
    content_coords = content_coords[:,1:]
    
    # Add noise to both channels
    noise = 300
    noise_SNARE = np.zeros((3,noise))
    noise_SNARE[0,:] = np.random.random(noise) * 800 - 400
    noise_SNARE[1,:] = np.random.random(noise) * 800 - 400
    noise_SNARE[2,:] = np.random.random(noise) * 400 - 200
    SNARE_coords = np.hstack((SNARE_coords,noise_SNARE))
    
    noise_vesicle = np.zeros((3,noise))
    noise_vesicle[0,:] = np.random.random(noise) * 800 - 400
    noise_vesicle[1,:] = np.random.random(noise) * 800 - 400
    noise_vesicle[2,:] = np.random.random(noise) * 400 - 200
    vesicle_coords = np.hstack((vesicle_coords,noise_vesicle))
    
    noise_content = np.zeros((3,noise))
    noise_content[0,:] = np.random.random(noise) * 800 - 400
    noise_content[1,:] = np.random.random(noise) * 800 - 400
    noise_content[2,:] = np.random.random(noise) * 400 - 200
    content_coords = np.hstack((content_coords,noise_content))
    
    return vesicle_coords, content_coords, SNARE_coords