import matplotlib.pyplot as plt
import pandas as pd
import Water_level as water_level
import numpy as np

import IS2_file as is2File


def p_is2(m,is2,laser,images_fp):

    #readjusting photon depths to the sea level
    f = is2.get_sea_level_function(laser)
    # print(f)
    lats = m.Latitude.drop_duplicates()
    sea = f(lats)
    mean_sea = np.mean(sea)
    sea = sea - np.mean(sea)

    #thresholds to classify photon as land or reef or water
    threshold = -0.3
    threshold_land = 2
    #plotting the sea, land and reef in different colours based on the thresholds
    sea_level = m.loc[(m.Photon_depth > threshold) & (m.Photon_depth < threshold_land) ]
    reef = m.loc[m.Photon_depth <= threshold]
    land = m.loc[m.Photon_depth > threshold_land]
    fig,ax = plt.subplots()
    plt.scatter(sea_level.Latitude, sea_level.Photon_depth, s = 0.1)
    plt.scatter(reef.Latitude, reef.Photon_depth, s= 0.1, color = 'green')
    plt.scatter(land.Latitude, land.Photon_depth, s =0.1, color = 'brown')
    plt.xlabel('Latitude')
    plt.ylabel('Height (m)')
    date = is2.get_date()
    track = is2.get_track()
    plt.title('Track: ' + str(track) + ' - Date: ' + str(date.date()) + ' - Laser: ' + laser)

    #plotting depth predictions on top of the photon data
    # print(sea)
    # print(len(sea))
    # print(lats)
    # print(len(lats))
    # print(mean_sea)
    plt.plot(lats, sea, linestyle='-', marker='o',color='blue', linewidth=3, markersize=2)

    depth_prof = plt.plot(m.Latitude, m.Predicted_depth, linestyle='-', marker='o',color='orange', linewidth=1, markersize=2,alpha = 0.4)
    plt.savefig(images_fp + '/photon_preds' + str(date.date())+laser + '.png', dpi=300)
    return ax
    
def plot_is2_depths(photons, is2, laser, f, bounds, titleprefix):
    #photons is ATL03 geolocated photon dataframe
    #is2 is ATL03 file object
    #laser is beam id ['gt1r']
    #f is sea level function
    #bounds = [-25, 5]
    #titleprefix is optional prefix to plot title
    
    #initialize plotting
    plt.figure(figsize=(16,6))
    plt.xlabel('Latitude', fontsize=16)
    plt.ylabel('Height (m)', fontsize=16)
    date = is2.get_date()
    is2_track = is2.get_track()
    titletext = 'Track: ' + str(is2_track) + ' - Date: ' + str(date.date()) + ' - Laser: ' + laser
    if titleprefix == 'False':
        plt.title(titletext, fontsize=16)
    else: 
        plt.title(titleprefix + ' - ' + titletext, fontsize=16)
    if bounds != 'False':
        plt.ylim(bounds)

    #plot photons
    plt.scatter(photons.Latitude, photons.Height, s = 0.1, color = 'black')
    
    #plotting sea level
    if f != 'False':
        maxlat = max(photons.lats)
        minlat = min(photons.lats)
        sea_lats = np.arange(minlat,maxlat,0.01).tolist()
        sea_level = f(sea_lats)
        plt.scatter(sea_lats, sea_level, s= 0.1, color = 'yellow')
        plt.plot(sea_lats, sea_level, linestyle='-', marker='o', color='yellow', linewidth=3, markersize=2)

    return
    
def plot_is2_depths_bathy(photons, is2, laser, depths, f, bounds, titleprefix):

    #initialize plotting
    plt.figure(figsize=(16,6))
    plt.xlabel('Latitude', fontsize=16)
    plt.ylabel('Height (m)', fontsize=16)
    date = is2.get_date()
    is2_track = is2.get_track()
    titletext = 'Track: ' + str(is2_track) + ' - Date: ' + str(date.date()) + ' - Laser: ' + laser
    if titleprefix == 'False':
        plt.title(titletext, fontsize=16)
    else: 
        plt.title(titleprefix + ' - ' + titletext, fontsize=16)
    if bounds != 'False':
        plt.ylim(bounds)

    #plot photons
    plt.scatter(photons.Latitude, photons.Height, s = 0.1, color = 'black')
    
    #plot depths
    plt.xlim(min(depths.Latitude), max(depths.Latitude))
    plt.plot(depths.Latitude, depths.Height, linestyle='-', marker='o',color='orange', linewidth=1, markersize=2,alpha = 0.4)

    #plotting sea level
    if f != 'False':
        maxlat = max(ph_lats)
        minlat = min(ph_lats)
        sea_lats = np.arange(minlat,maxlat,0.01).tolist()
        sea_level = f(sea_lats)
        plt.scatter(sea_lats, sea_level, s= 0.1, color = 'yellow')
        plt.plot(sea_lats, sea_level, linestyle='-', marker='o', color='yellow', linewidth=3, markersize=2)

    return