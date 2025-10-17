#plotting packages
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.colors import BoundaryNorm
import matplotlib as mpl
import matplotlib.pylab as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib import patches
import matplotlib.pyplot as plt
from rasterio.plot import show
import matplotlib.ticker as ticker

import math
import numpy as np
import os
import pandas as pd


def plot_s2_band_raw(meta):
#Plot all four bands in one plot.
    
    plt.rc('font', family='serif')
    imgs = meta['imgs']
    mask_thresh = meta['mask_thresh']
    infred = imgs[3]
    msk = infred>mask_thresh
    fig = plt.figure(figsize=(17, 17))
    
    #mask_recipe = lambda a: False if a < mask_thresh else True
    #msk = mask_recipe(np.ndarray.flatten(infred))
    #print(msk)

    #define x axis
    step_x = 500 #define ticks for plotting
    nx = imgs[0].shape[2]
    x_positions = np.arange(0, nx, step_x)
    x_positions = x_positions + (nx - np.max(x_positions))/2
    x_labels = (x_positions - np.mean(x_positions))*meta['xdim']/1000
    utm_x0 = meta['ulx'] + meta['xdim']*nx/2
    x_axis_text = 'UTM X - {:.0f} (km)'.format(utm_x0/1000)

    #define y axis
    step_y = 500
    ny = imgs[0].shape[1]
    y_positions = np.arange(0, ny, step_y)
    y_positions = y_positions + (ny - np.max(y_positions))/2
    y_labels = (y_positions - np.mean(y_positions))*meta['ydim']/1000
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0/1000)

    #set overall figure size
    fig = plt.figure(figsize=(17, 17))

    #plot Blue band      
    blue = np.copy(imgs[0])
    blue[msk] = 0
    ax = fig.add_subplot(2, 2, 1)
    ax.set_title('Sentinel-2 Band 02 (Blue)', fontsize=16)   
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    lowbound = np.nanmean(imgs[0]) - 2*np.nanstd(imgs[0])
    highbound = np.nanmean(imgs[0]) + 2*np.nanstd(imgs[0])
    norm = mpl.colors.Normalize(vmin = lowbound, vmax = highbound)
    imgplot = show(blue, ax = ax, cmap='Blues', norm=norm)   
   
    #plot Green band
    green = np.copy(imgs[1])
    green[msk] = 0
    ax = fig.add_subplot(2, 2, 2)
    ax.set_title('Sentinel-2 Band 03 (Green)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(green, ax = ax, cmap='Greens', norm=norm)   

    #plot Red band   
    red = np.copy(imgs[2])
    red[msk] = 0
    ax = fig.add_subplot(2, 2, 3)
    ax.set_title('Sentinel-2 Band 04 (Red)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(red, ax = ax, cmap='Reds', norm=norm)  

    #plot IR band        
    ax = fig.add_subplot(2, 2, 4)
    ax.set_title('Sentinel-2 Band 08 (Infrared)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[3], ax = ax, cmap='Greys', norm=norm)  
    
    #plot B1 coastal aerosol band
    aerosol = np.copy(imgs[4])
    aerosol[msk] = 0
    fig = plt.figure(figsize=(17, 17))
    ax = fig.add_subplot(2, 2, 1)
    ax.set_title('Sentinel-2 Band 01 (Coastal Aerosol)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(aerosol, ax = ax, cmap='Purples', norm=norm)  


def plot_s2_band_histogram(meta):

    np.seterr(divide='ignore', invalid='ignore')
    imgs = meta['imgs']

    mask_thresh = meta['mask_thresh']
    infred = imgs[3]
    msk = infred>mask_thresh
    msk = np.ndarray.flatten(msk)
    
    #extract band data
    blue_flat = np.ndarray.flatten(imgs[0])
    blue_flat[msk] = 0
    green_flat = np.ndarray.flatten(imgs[1])
    green_flat[msk] = 0
    red_flat = np.ndarray.flatten(imgs[2])
    red_flat[msk] = 0
    ir_flat = np.ndarray.flatten(imgs[3])
    aerosol_flat = np.ndarray.flatten(imgs[4])
    aerosol_flat[msk] = 0

    #set up plotting
    fig = plt.figure(figsize=(17, 5))

	#Plot histogram
    ax = fig.add_subplot(1, 1, 1)
    plt.hist([blue_flat, green_flat, red_flat, ir_flat, aerosol_flat], bins=np.arange(0,6000,10), histtype='step', stacked=False, color = ['Blue', 'Green', 'Red', 'Black','Purple'])
    ax.set_title('Sentinel-2 Band Intensity Distribution', fontsize=16)
    ax.set_xlabel('Intensity in Sentinel L2A units')
    ax.set_ylabel('N')
    ax.set_yscale('log')
    plt.setp(ax, xlim = [0, 6000], ylim = [1, 1000000])


def plot_s2_band_diff_histogram(meta):

    np.seterr(divide='ignore', invalid='ignore')
    imgs = meta['imgs']
    
    mask_thresh = meta['mask_thresh']
    infred = imgs[3]
    msk = infred>mask_thresh
    msk = np.ndarray.flatten(msk)

    #extract band data
    blue_flat = np.ndarray.flatten(imgs[0])
    blue_flat[msk] = 0
    green_flat = np.ndarray.flatten(imgs[1])
    green_flat[msk] = 0
    red_flat = np.ndarray.flatten(imgs[2])
    red_flat[msk] = 0
    ir_flat = np.ndarray.flatten(imgs[3])
    aerosol_flat = np.ndarray.flatten(imgs[4])
    aerosol_flat[msk] = 0

    #set up plotting
    fig = plt.figure(figsize=(17, 5))

	#Plot histogram
    ax = fig.add_subplot(1, 1, 1)
    plt.hist([np.log(aerosol_flat) - np.log(blue_flat)], bins=np.arange(-6,10,.02), histtype='step', stacked=False, color = ['Purple'])
    plt.hist([np.log(blue_flat) - np.log(green_flat)], bins=np.arange(-6,10,.02), histtype='step', stacked=False, color = ['Blue'])
    plt.hist([np.log(green_flat) - np.log(red_flat)], bins=np.arange(-6,10,.02), histtype='step', stacked=False, color = ['Red'])
    ax.set_title('Sentinel-2 Log Band Difference Distribution', fontsize=16)
    ax.set_xlabel('Log(Band A) - Log(Band B)')
    ax.set_ylabel('N')
    ax.set_yscale('log')
    plt.setp(ax, xlim = [-4, 8], ylim = [1, 1000000])


def plot_s2_band_diff(meta, flag):

    np.seterr(divide='ignore', invalid='ignore')
    imgs = meta['imgs']
    
    mask_thresh = meta['mask_thresh']
    infred = imgs[3]
    msk = infred>mask_thresh

    #define x axis
    step_x = 1000
    nx = imgs[0].shape[2]
    x_positions = np.arange(0, nx, step_x)
    x_positions = x_positions + (nx - np.max(x_positions))/2
    x_labels = x_positions - np.mean(x_positions)
    utm_x0 = meta['ulx'] + meta['xdim']*nx/2
    x_axis_text = 'UTM X - {:.0f} (m)'.format(utm_x0)

    #define y axis
    step_y = 1000
    ny = imgs[0].shape[1]
    y_positions = np.arange(0, ny, step_y)
    y_positions = y_positions + (ny - np.max(y_positions))/2
    y_labels = y_positions - np.mean(y_positions)
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    #set up plotting
    fig = plt.figure(figsize=(17, 7.5))

 	#plot Blue-Green
    ax = fig.add_subplot(1, 2, 1)
    blue = imgs[0].astype(float)
    blue[msk] = 0
    blue = np.log(blue)
    green = imgs[1].astype(float)
    green[msk] = 0
    green = np.log(green)
    logdiff = blue - green
    lowbound = np.nanmean(logdiff) - 2*np.nanstd(logdiff)
    highbound = np.nanmean(logdiff) + 2*np.nanstd(logdiff)
    print('Blue/Green plotbounds:', lowbound, highbound)
    if flag == 'green_high':
        #only plot pixels where green band is higher than blue band
        logdiff = np.where((logdiff > 0) & (blue < 6.6), 1, logdiff)    
    ax.set_title('Sentinel-2 Log(Blue) - Log(Green)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = lowbound, vmax = highbound)
    imgplot = show(logdiff, ax = ax, cmap='Blues', norm=norm)

 	#plot Green-Red
    ax = fig.add_subplot(1, 2, 2)
    red = imgs[2].astype(float)
    red[msk]=0
    red = np.log(red)
    logdiff = green - red
    lowbound = np.nanmean(logdiff) - 2*np.nanstd(logdiff)
    highbound = np.nanmean(logdiff) + 2*np.nanstd(logdiff)
    print('Green/Red plotbounds:', lowbound, highbound)    
    ax.set_title('Sentinel-2 Log(Green) - Log(Red)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = lowbound, vmax = highbound)
    imgplot = show(logdiff, ax = ax, cmap='Reds', norm=norm)
    
    
    fig = plt.figure(figsize=(17, 7.5))
    
    #plot Aerosol-Blue 
    ax = fig.add_subplot(1, 2, 1)
    aerosol = imgs[4].astype(float)
    aerosol[msk]=0
    aerosol = np.log(aerosol)
    blue_a = np.copy(blue)
    ind = np.isinf(aerosol)
    blue_a[ind] = aerosol[ind]
    logdiff = aerosol - blue_a
    lowbound = np.nanmean(logdiff) - 2*np.nanstd(logdiff)
    highbound = np.nanmean(logdiff) + 2*np.nanstd(logdiff)
    print('Aerosol/Blue plotbounds:', lowbound, highbound)    
    ax.set_title('Sentinel-2 Log(Aerosol) - Log(Blue)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = lowbound, vmax = highbound)
    imgplot = show(logdiff, ax = ax, cmap='Blues', norm=norm)#cmap='PuBu'
    
    #plot Aerosol-Red 
    ax = fig.add_subplot(1, 2, 2)
    red_a = np.copy(red)
    red_a[ind] = aerosol[ind]
    logdiff = aerosol - red_a
    lowbound = np.nanmean(logdiff) - 2*np.nanstd(logdiff)
    highbound = np.nanmean(logdiff) + 2*np.nanstd(logdiff)
    print('Aerosol/Red plotbounds:', lowbound, highbound)    
    ax.set_title('Sentinel-2 Log(Aerosol) - Log(Red)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = lowbound, vmax = highbound)
    imgplot = show(logdiff, ax = ax, cmap='Reds', norm=norm)#cmap='PuRd'


def plot_s2_band_vs_band(meta):
#Plot bands against each other
    
    imgs = meta['imgs']
    
    mask_thresh = meta['mask_thresh']
    infred = imgs[3]
    msk = infred>mask_thresh
    msk = np.ndarray.flatten(msk)
    
    #extract band data
    blue_flat = np.ndarray.flatten(imgs[0])
    blue_flat[msk] = 0
    green_flat = np.ndarray.flatten(imgs[1])
    green_flat[msk] = 0
    red_flat = np.ndarray.flatten(imgs[2])
    red_flat[msk] = 0
    ir_flat = np.ndarray.flatten(imgs[3])
    aerosol_flat = np.ndarray.flatten(imgs[4])
    aerosol_flat[msk] = 0

    #set up plotting
    fig = plt.figure(figsize=(17, 17))
    xlimit = [-20, 0]
    plt.rc('font', size=12)
    
    #Log(Blue) vs. Log(Aerosol)
    ax = fig.add_subplot(2, 2, 1)
    plt.title('Blue vs Aerosol Comparison')
    ax.set_xlabel('Log(Aerosol)', fontsize=13)
    ax.set_ylabel('Log(Blue)', fontsize=13)
    plt.setp(ax, xlim = [3, 9.5], ylim = [3, 9.5])
    plt.scatter(np.log(aerosol_flat), np.log(blue_flat), color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')

    #Log(Green) vs. Log(Blue) 
    ax = fig.add_subplot(2, 2, 2)
    plt.title('Green vs Blue Comparison')
    ax.set_xlabel('Log(Blue)', fontsize=13)
    ax.set_ylabel('Log(Green)', fontsize=13)
    plt.setp(ax, xlim = [3, 9.5], ylim = [3, 9.5])
    plt.scatter(np.log(blue_flat), np.log(green_flat), color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    #Log(Red) vs. Log(Green)
    ax = fig.add_subplot(2, 2, 3)
    plt.title('Red vs Green Comparison')
    ax.set_xlabel('Log(Green)', fontsize=13)
    ax.set_ylabel('Log(Red)', fontsize=13)
    plt.setp(ax, xlim = [3, 9.5], ylim = [3, 9.5])
    plt.scatter(np.log(green_flat), np.log(red_flat), color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    fig = plt.figure(figsize=(17, 17))
    xlimit = [-20, 0]
    plt.rc('font', size=12)
    
    #Log(Band) Aerosol - Blue against each other
    ax = fig.add_subplot(2, 2, 1)
    plt.title('Blue vs Aerosol Comparison')
    ax.set_xlabel('Log(Aerosol)', fontsize=13)
    ax.set_ylabel('Log(Blue) - Log(Red)', fontsize=13)
    plt.setp(ax, xlim = [3, 9.5], ylim = [-2.0, 2.5])
    plt.scatter(np.log(aerosol_flat), np.log(aerosol_flat) - np.log(blue_flat), color = 'blue', s = 2)
    
    #Log(Blue) - Log(Green) vs. Log(Blue)
    ax = fig.add_subplot(2, 2, 2)
    plt.title('Green vs Blue Comparison')
    ax.set_xlabel('Log(Blue)', fontsize=13)
    ax.set_ylabel('Log(Blue) - Log(Green)', fontsize=13)
    plt.setp(ax, xlim = [3, 9.5], ylim = [-2.0, 2.5])
    plt.scatter(np.log(blue_flat), np.log(blue_flat) - np.log(green_flat), color = 'blue', s = 2)
    
    #Log(Green) - Log(Red) vs. Log(Green)
    ax = fig.add_subplot(2, 2, 3)
    plt.title('Red vs Green Comparison')
    ax.set_xlabel('Log(Green)', fontsize=13)
    ax.set_ylabel('Log(Green) - Log(Red)', fontsize=13)
    plt.setp(ax, xlim = [3, 9.5], ylim = [-2.0, 2.5])
    plt.scatter(np.log(green_flat), np.log(green_flat) - np.log(red_flat), color = 'blue', s = 2)


def plot_s2_band_ratio(meta):

    np.seterr(divide='ignore', invalid='ignore')
    imgs = meta['imgs']

    #define x axis
    step_x = 500
    nx = imgs[0].shape[2]
    x_positions = np.arange(0, nx, step_x)
    x_positions = x_positions + (nx - np.max(x_positions))/2
    x_labels = x_positions - np.mean(x_positions)
    utm_x0 = meta['ulx'] + meta['xdim']*nx/2
    x_axis_text = 'UTM X - {:.0f} (m)'.format(utm_x0)

    #define y axis
    step_y = 500
    ny = imgs[0].shape[1]
    y_positions = np.arange(0, ny, step_y)
    y_positions = y_positions + (ny - np.max(y_positions))/2
    y_labels = y_positions - np.mean(y_positions)
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    #set up plotting
    fig = plt.figure(figsize=(15, 7.5))

 	#plot Blue-Green
    ax = fig.add_subplot(1, 2, 1)
    blue = imgs[0].astype(float)
    blue = np.log(blue)
    green = imgs[1].astype(float)
    green = np.log(green)
    logdiff = blue/green
    
    ax.set_title('Log(Blue)/Log(Green)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = .5, vmax = 1.5)
    imgplot = show(logdiff, ax = ax, cmap='Blues', norm=norm)

 	#plot Green-Red
    ax = fig.add_subplot(1, 2, 2)
    red = imgs[2].astype(float)
    red = np.log(red)
    logdiff = green/red
    
    ax.set_title('Log(Green)/Log(Red)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = 1, vmax = 2)
    imgplot = show(logdiff, ax = ax, cmap='Reds', norm=norm)

