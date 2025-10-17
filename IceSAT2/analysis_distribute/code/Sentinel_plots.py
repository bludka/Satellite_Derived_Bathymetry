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

    #set overall figure size
    fig = plt.figure(figsize=(15, 15))

    #plot Blue band        
    ax = fig.add_subplot(2, 2, 1)
    ax.set_title('Sentinel-2 Band 02 (Blue)', fontsize=16)   
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[0], ax = ax, cmap='Blues')
   
    #plot Green band        
    ax = fig.add_subplot(2, 2, 2)
    ax.set_title('Sentinel-2 Band 03 (Green)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[1], ax = ax, cmap='Greens')

    #plot Red band        
    ax = fig.add_subplot(2, 2, 3)
    ax.set_title('Sentinel-2 Band 04 (Red)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[2], ax = ax, cmap='Reds')

    #plot IR band        
    ax = fig.add_subplot(2, 2, 4)
    ax.set_title('Sentinel-2 Band 08 (Infrared)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[3], ax = ax, cmap='Greys')


def plot_s2_band_histogram(meta):

    np.seterr(divide='ignore', invalid='ignore')
    imgs = meta['imgs']

    #extract band data
    blue_flat = np.ndarray.flatten(imgs[0])
    green_flat = np.ndarray.flatten(imgs[1])
    red_flat = np.ndarray.flatten(imgs[2])
    ir_flat = np.ndarray.flatten(imgs[3])

    #set up plotting
    fig = plt.figure(figsize=(15, 5))

	#Plot histogram
    ax = fig.add_subplot(1, 1, 1)
    plt.hist([blue_flat, green_flat, red_flat, ir_flat], bins=np.arange(0,6000,10), histtype='step', stacked=False, color = ['Blue', 'Green', 'Red', 'Black'])
    ax.set_title('Sentinel-2 Band Intensity Distribution', fontsize=16)
    ax.set_xlabel('Intensity in Sentinel L2A units')
    ax.set_ylabel('N')
    ax.set_yscale('log')
    plt.setp(ax, xlim = [0, 6000], ylim = [1, 500000])


def plot_s2_band_diff_histogram(meta):

    np.seterr(divide='ignore', invalid='ignore')
    imgs = meta['imgs']

    #extract band data
    blue_flat = np.ndarray.flatten(imgs[0])
    green_flat = np.ndarray.flatten(imgs[1])
    red_flat = np.ndarray.flatten(imgs[2])
    ir_flat = np.ndarray.flatten(imgs[3])

    #set up plotting
    fig = plt.figure(figsize=(15, 5))

	#Plot histogram
    ax = fig.add_subplot(1, 1, 1)
    plt.hist([np.log(blue_flat) - np.log(green_flat)], bins=np.arange(-6,10,.02), histtype='step', stacked=False, color = ['Blue'])
    plt.hist([np.log(green_flat) - np.log(red_flat)], bins=np.arange(-6,10,.02), histtype='step', stacked=False, color = ['Red'])
    ax.set_title('Sentinel-2 Log Band Difference Distribution', fontsize=16)
    ax.set_xlabel('Log(Band A) - Log(Band B)')
    ax.set_ylabel('N')
    ax.set_yscale('log')
    plt.setp(ax, xlim = [-4, 8], ylim = [1, 500000])


def plot_s2_band_diff(meta):

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
    logdiff = blue - green
    
    ax.set_title('Sentinel-2 Log(Blue) - Log(Green)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = -0.25, vmax = 1.25)
    imgplot = show(logdiff, ax = ax, cmap='Blues', norm=norm)

 	#plot Green-Red
    ax = fig.add_subplot(1, 2, 2)
    red = imgs[2].astype(float)
    red = np.log(red)
    logdiff = green - red
    
    ax.set_title('Sentinel-2 Log(Green) - Log(Red)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = 0, vmax = 7)
    imgplot = show(logdiff, ax = ax, cmap='Reds', norm=norm)


def plot_s2_band_vs_band(meta):
#Plot bands against each other
    
    imgs = meta['imgs']

    #extract band data
    blue_flat = np.ndarray.flatten(imgs[0])
    green_flat = np.ndarray.flatten(imgs[1])
    red_flat = np.ndarray.flatten(imgs[2])
    ir_flat = np.ndarray.flatten(imgs[3])

    #set up plotting
    fig = plt.figure(figsize=(15, 7.5))
    xlimit = [-20, 0]
    plt.rc('font', size=12)

    #Log(Band) Blue - Green against each other
    ax = fig.add_subplot(1, 2, 1)
    plt.title('Blue vs Green Comparison')
    ax.set_xlabel('Log(Blue)', fontsize=13)
    ax.set_ylabel('Log(Green)', fontsize=13)
    plt.setp(ax, xlim = [3, 8.5], ylim = [3, 8.5])
    plt.scatter(np.log(blue_flat), np.log(green_flat), color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    #Log(Band) Green - Red against each other
    ax = fig.add_subplot(1, 2, 2)
    plt.title('Green vs Red Comparison')
    ax.set_xlabel('Log(Green)', fontsize=13)
    ax.set_ylabel('Log(Red)', fontsize=13)
    plt.setp(ax, xlim = [1, 8], ylim = [1, 8])
    plt.scatter(np.log(green_flat), np.log(red_flat), color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')


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

