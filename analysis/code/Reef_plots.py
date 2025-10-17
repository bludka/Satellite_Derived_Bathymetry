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
import cmocean

def band_histogram(df):
#method to plot sentinel histogram of just the training 

    fig = plt.figure(figsize=(15, 5))
    ax = fig.add_subplot(1, 1, 1)
    plt.hist([df['b2'], df['b3'], df['b4'], df['b8'],df['b1']], bins=np.arange(0,3000,10), histtype='step', stacked=False, color = ['Blue', 'Green', 'Red', 'Black','Purple'])
    ax.set_title('Sentinel-2 Band Intensity Distribution: Training Pixels Only', fontsize=16)
    ax.set_xlabel('Intensity in Sentinel L2A units')
    ax.set_ylabel('Frequency')
    plt.setp(ax, xlim = [0, 6000], ylim = [0, 2000])


def plot_median_variance_graph(medians, variances, median_threshold, variance_threshold,outpath):
    fig, ax = plt.subplots()
    sns.scatterplot(x = medians, y = variances, legend = 'full',ax = ax)
    plt.xlabel('Median pixel value for log difference')
    plt.ylabel('Variance pixel value for log difference')

    ax.add_patch(
        patches.Rectangle(
            xy=(-median_threshold, -variance_threshold),  # point of origin.
            width=2*median_threshold,
            height=2*variance_threshold,
            linewidth=1,
            color='red',
            fill=False
        )
    )

    fn = os.path.join(outpath,'median_vs_variance.png')
    plt.savefig(fn)


def plot_reefs(reef_depths, training_data, sf, line, reef, vmin, vmax):
#method to plot the reef depth histogram and a scatter plot of the same

    plt.rc('font', family='sans-serif', size=13)

    reef_name = reef.get_reef_name()
    dt = sf.get_date().strftime("%Y%m%d%H%M%S")

    #plot histogram of depths
    fig, ax = plt.subplots(1, 3, figsize = (28,10))
    fig.patch.set_facecolor('white')
    depth_histogram(ax[0], reef_depths, reef)
    #plot line fit to training data
    line_of_best_fit(ax[1], training_data, sf, line)    
    #plot reef depths
    reef_scatter(fig, ax[2], reef_depths, training_data, sf, reef, vmin, vmax)
    print('Done Plotting')

    #save plot
    fn = '{reef_name}-{date}.png'.format(reef_name = reef_name, date = dt)
    out = os.path.join(reef.get_reef_depth_plots(), fn)
    plt.savefig(out)
    print('Done Saving')
    plt.close(fig)
    return

def depth_histogram(ax, reef_depths, reef, suffix = ''):
    reef_name = reef.get_reef_name()
    reef_depths['Height'].plot.hist(bins = np.arange(-30,20,1), ax = ax)
    ax.set_xlabel('Height (m)')
    ax.set_ylabel('Frequency')
    ax.set_title('{reef_name} {suffix} Depth Histogram'.format(reef_name = reef_name, suffix = suffix))

def line_of_best_fit(ax, training_data, sf, line):
    r = 3
    sns.scatterplot(x = training_data['diff_b2_b3'], y = training_data.Height, color = 'blue', ax = ax)
    sns.scatterplot(x = training_data['diff_b3_b4'], y = training_data.Height, color = 'green', ax = ax)
    sns.lineplot(x = [-r,r], y = [line(-r),line(r)], color = 'black', ax = ax)
    sns.lineplot(x = [0,0], y = [line(-r),line(r)], color = 'black', dashes=True, ax = ax)
    ax.set_xlabel('Log(Blue Band) - Log(Green Band)')
    #ax.set_xlabel('Log(Green Band) - Log(Red Band)')
    ax.set_ylabel('Depth')
    xt = (list(ax.get_xticks()))[:-1]

    for i,x in enumerate(xt):
        xt[i] = np.round(x,2)
    ax.set_title(sf.get_date().strftime("%Y/%m/%d %H:%M:%S") + '. Tide = ' + str(sf.get_tide()) + 'm')
    xlim = (-r, r)
    ylim = ( -25,0)
    plt.setp(ax, xlim=xlim, ylim=ylim)


def reef_scatter(fig, ax, reef_depths, training_data, safe_file, reef, vmin, vmax, suffix = ''):

    reef_name = reef.get_reef_name()
    meta = safe_file.get_meta()
    imgs = meta['imgs']

    # Color map
    cmap = cmocean.tools.crop(cmocean.cm.topo, vmin, vmax, 1)
    bounds = np.arange(vmin, vmax, 1)
    norm = BoundaryNorm(bounds, ncolors = 256)

    ##Custom Color Map
    #cmap = cm.colors.ListedColormap(['black','darkblue','mediumblue','blue','royalblue','cornflowerblue',
    #                                 'skyblue','darkgreen','forestgreen','limegreen','greenyellow','yellow',
    #                                 'orange','tomato','red'])
    #bounds = np.arange(-24,7,2)
    #norm = BoundaryNorm(bounds,cmap.N)
    ## extract all colors from the .jet map
    #cmaplist = [cmap(i) for i in range(cmap.N)]
    ## create the new map
    #cmap = mpl.colors.LinearSegmentedColormap.from_list(
    #    'Custom cmap', cmaplist, cmap.N)

    # create a second axes for the colorbar
    ax2 = fig.add_axes([0.905, 0.125, 0.025, 0.755])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
        spacing='proportional', ticks=bounds, boundaries=bounds, format='%.1f')

    #Define x axis
    utm_x0 = np.min(reef_depths.x)
    utm_y0 = np.min(reef_depths.y)
    x_axis_text = 'UTM X - {:.0f} (m)'.format(utm_x0)
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    #Scatter plot of the predicted depths
    pts = ax.scatter(x = reef_depths.x - utm_x0, y = reef_depths.y - utm_y0, c = reef_depths.Height, s = 1, cmap = cmap, norm = norm)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    ax.set_title('{reef_name} {suffix} Depth Predictions (m)'.format(reef_name = reef_name, suffix = suffix))   
    
    #scatter plot of the track lines from the ICESAT 2 data
    ax.scatter(x = training_data.x - utm_x0, y = training_data.y - utm_y0, s = 3, c = 'black', label = 'ICESAT-2 tracks')
    
    custom_lines = [Line2D([0], [0], color='black', lw=4)]
    ax.legend(custom_lines, ['ICESAT-2 tracks'])


def aggregate_plot(data,df,sf, f):
    df['Height'] = df['median']
    reef_name = sf.get_reef_name()
    fig ,ax = plt.subplots(1,2, figsize = (28,12))
    depth_histogram(ax[0], df, reef_name, f)
    reef_scatter(fig,ax[1],df,reef_name,data,f)
    fn = '{reef_name}-{f}.png'.format(reef_name = reef_name, f = f)
    out = os.path.join(sf.get_img_path(), fn)
    plt.savefig(out)
    plt.close(fig)


def depth_histogram_plot(df):
    fig = plt.figure(figsize=(16, 5))

    ax = fig.add_subplot(1, 2, 1)
    plt.hist([np.log(df['b2']), np.log(df['b3'])], bins=np.arange(5,9,0.05), histtype='step', stacked=False, color = ['Blue', 'Green'])
    ax.set_xlabel('Log Raw Band Intensity')
    ax.set_ylabel('Frequency')

    ax = fig.add_subplot(1, 2, 2)
    plt.hist([np.log(df['b2']) - np.log(df['b3'])], bins=np.arange(-0.5,1.5,0.025), histtype='step', color ='Black')
    ax.set_xlabel('Log(Blue) - Log(Green)')
    ax.set_ylabel('Frequency')


def plot_sentinel_cal_separate_profiles(data):

    #set limits
    xlimit = [-20, 0]

    #set up plotting
    fig = plt.figure(figsize=(17, 11))
    plt.rc('font', size=12)

    #Raw Log(Bands)
    ax = fig.add_subplot(2, 3, 1)
    plt.title('Log intensity Sentinel-2 bands')
    ax.set_xlabel('Depth', fontsize=13)
    ax.set_ylabel('Log Band Intensity', fontsize=13)
    plt.setp(ax, xlim = xlimit, ylim = [0, 8.5])
    plt.scatter(data.Height, data.log_b2, marker='.', color = 'blue', s = 2)  
    plt.scatter(data.Height, data.log_b3, marker='.', color = 'green', s = 2)  
    plt.scatter(data.Height, data.log_b4, marker='.', color = 'red', s = 2)  

    #Log(Band) Blue/Green against each other
    ax = fig.add_subplot(2, 3, 2)
    plt.title('Comparison of log intensity of Blue vs Green')
    ax.set_xlabel('Blue', fontsize=13)
    ax.set_ylabel('Green', fontsize=13)
    plt.setp(ax, xlim = [3, 8.5], ylim = [3, 8.5])
    plt.scatter(data.log_b2, data.log_b3, color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    #Log(Band) differences    
    ax = fig.add_subplot(2, 3, 3)
    plt.setp(ax, xlim = xlimit, ylim = [-.75, 2.0])
    plt.title('Log band difference Blue vs. Green')
    ax.set_xlabel('Depth', fontsize=13)
    plt.scatter(data.Height, data.diff_b2_b3, marker='.', color = 'aqua', s = 2)  
    plt.plot([-20,0],[0,0],linestyle='dashed', color = 'black')

    #Log(Band) Green/Red against each other
    ax = fig.add_subplot(2, 3, 5)
    plt.title('Comparison of log intensity of Green vs Red')
    ax.set_xlabel('Green', fontsize=13)
    ax.set_ylabel('Red', fontsize=13)
    plt.setp(ax, xlim = [1, 8], ylim = [1, 8])
    plt.scatter(data.log_b3, data.log_b4, color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    #Log(Band) differences    
    ax = fig.add_subplot(2, 3, 6)
    plt.setp(ax, xlim = xlimit, ylim = [0, 8])
    plt.title('Log band difference Green vs. Red')
    plt.scatter(data.Height, data.diff_b3_b4, marker='.', color = 'aqua', s = 2)  
    plt.plot([-20,0],[0,0],linestyle='dashed', color = 'black')


def plot_sentinel_cal(data, slope1, intercept1, slope2, intercept2, slope3, intercept3):

    #set limits
    xlimit = [-20, 0]

    #set up plotting
    fig = plt.figure(figsize=(20, 20))
    plt.rc('font', size=12)

    #Raw Log(Bands)
    ax = fig.add_subplot(3, 3, 1)
    plt.title('Log intensity Sentinel-2 bands')
    ax.set_xlabel('Depth', fontsize=13)
    ax.set_ylabel('Log Band Intensity', fontsize=13)
    plt.setp(ax, xlim = xlimit, ylim = [3, 9])
    plt.scatter(data.Height, data.log_b2, marker='.', color = 'blue', s = 2)  
    plt.scatter(data.Height, data.log_b3, marker='.', color = 'green', s = 2)  
    plt.scatter(data.Height, data.log_b4, marker='.', color = 'red', s = 2)  

    #Log(Band) Blue/Green against each other
    ax = fig.add_subplot(3, 3, 2)
    plt.title('Comparison of log intensity of Blue vs Green')
    ax.set_xlabel('Blue', fontsize=13)
    ax.set_ylabel('Green', fontsize=13)
    plt.setp(ax, xlim = [3, 9], ylim = [3, 9])
    plt.scatter(data.log_b2, data.log_b3, color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    #Log(Band) differences Blue vs. Green   
    ax = fig.add_subplot(3, 3, 3)
    plt.setp(ax, xlim = xlimit, ylim = [-1.0, 2.0])
    plt.title('Log band difference Blue vs. Green')
    ax.set_xlabel('Depth', fontsize=13)
    plt.scatter(data.Height, data.diff_b2_b3, marker='.', color = 'aqua', s = 2)  
    plt.plot(np.array(xlimit), np.array(xlimit)*slope1 + intercept1, linestyle='solid', color = 'black')  
    plt.plot([-20,0],[0,0],linestyle='dashed', color = 'black')

    #Log(Band) Green/Red against each other
    ax = fig.add_subplot(3, 3, 5)
    plt.title('Comparison of log intensity of Green vs Red')
    ax.set_xlabel('Green', fontsize=13)
    ax.set_ylabel('Red', fontsize=13)
    plt.setp(ax, xlim = [3, 9], ylim = [3, 9])
    plt.scatter(data.log_b3, data.log_b4, color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    #Log(Band) differences Green vs. Red
    ax = fig.add_subplot(3, 3, 6)
    plt.setp(ax, xlim = xlimit, ylim = [0.0, 3.0])
    plt.title('Log band difference Green vs. Red')
    plt.scatter(data.Height, data.diff_b3_b4, marker='.', color = 'aqua', s = 2)  
    plt.plot(np.array(xlimit), np.array(xlimit)*slope2 + intercept2, linestyle='solid', color = 'black')  
    plt.plot([-20,0],[0,0],linestyle='dashed', color = 'black')
    
    #Log(Band) Aerosol/Blue against each other
    ax = fig.add_subplot(3, 3, 8)
    plt.title('Comparison of log intensity of Aerosol vs Blue')
    ax.set_xlabel('Aerosol', fontsize=13)
    ax.set_ylabel('Blue', fontsize=13)
    plt.setp(ax, xlim = [3, 9], ylim = [3, 9])
    plt.scatter(data.log_b1, data.log_b2, color = 'blue', s = 2)
    plt.plot([1,20],[1,20],linestyle='dashed', color = 'black')
    
    #Log(Band) differences Aerosol vs. Blue
    ax = fig.add_subplot(3, 3, 9)
    #plt.setp(ax, xlim = xlimit, ylim = [0.0, 3.0])
    plt.title('Log band difference Aerosol vs. Blue')
    plt.scatter(data.Height, data.diff_b1_b2, marker='.', color = 'aqua', s = 2)  
    plt.plot(np.array(xlimit), np.array(xlimit)*slope3 + intercept3, linestyle='solid', color = 'black')  
    plt.plot([-20,0],[0,0],linestyle='dashed', color = 'black')


def plot_sentinel_icesat(imgs, meta, df):

    plt.rc('font', family='serif')

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

    #set overall figure size
    fig = plt.figure(figsize=(8, 8))

    #plot Blue band        
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Sentinel-2 Band 02 (Blue) + Training Pixels', fontsize=16)   
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[0], ax = ax, cmap='Blues')
   
    #overplot Blue>Green photon locations on Blue band
    ax.scatter(x = df.image_coord_x, y = df.image_coord_y, c = 'blue', marker = 0, s = 10)


def plot_sentinel_icesat_split(imgs, meta, df_use, df_toss):

    plt.rc('font', family='serif')

    #define x axis
    step_x = 200
    nx = imgs[0].shape[2]
    x_positions = np.arange(0, nx, step_x)
    x_positions = x_positions + (nx - np.max(x_positions))/2
    x_labels = x_positions - np.mean(x_positions)
    utm_x0 = meta['ulx'] + meta['xdim']*nx/2
    x_axis_text = 'UTM X - {:.0f} (m)'.format(utm_x0)

    #define y axis
    step_y = 200
    ny = imgs[0].shape[1]
    y_positions = np.arange(0, ny, step_y)
    y_positions = y_positions + (ny - np.max(y_positions))/2
    y_labels = y_positions - np.mean(y_positions)
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    #set overall figure size
    fig = plt.figure(figsize=(10, 10))

    #plot Blue band        
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Band 02 (Blue)', fontsize=16)   
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[0], ax = ax, cmap='Blues')
   
    #overplot Blue>Green photon locations on Blue band
    photon_y = [Item[0] for Item in df_use.image_coords]
    photon_x = [Item[1] for Item in df_use.image_coords]
    ax.scatter(x = photon_x, y= photon_y, c = 'blue', marker = 0, s = 10)

    #overplot Green>Blue locations on Green band
    photon_y = [Item[0] for Item in df_toss.image_coords]
    photon_x = [Item[1] for Item in df_toss.image_coords]
    ax.scatter(x = photon_x, y= photon_y, c = 'green', marker = 1, s = 10)


def corr_plot(datum,reef_name,outpath):
    r = 3
    num_blocks = int(np.ceil(np.sqrt(len(datum))))
    fig, ax = plt.subplots(num_blocks,num_blocks, figsize = (20,24))
    xlim = (-r, r)
    ylim = ( -25,0)
    plt.setp(ax, xlim=xlim, ylim=ylim)

    axlist = []
    for axl in ax:
        for axl2 in axl:
            axlist.append(axl2)

    day_keys = list(datum.keys())
    for i,dict_item in enumerate(datum.items()):
        d = dict_item[1][2]
        line = dict_item[1][0]
        meta = dict_item[1][1]

        sns.scatterplot(x = d['diff_b2_b3'], y = d.Height, color = 'blue', ax = axlist[i])
        sns.lineplot(x = [-r,r], y = [line(-r),line(r)], color = 'black', ax = axlist[i])
        axlist[i].set_xlabel('Log(Blue Band) - Log(Green Band)')
        axlist[i].set_ylabel('Depth')
        axlist[i].set_title(str(meta['dt'].date()))
        xt = (list(axlist[i].get_xticks()))[:-1]
        for i,x in enumerate(xt):
            xt[i] = np.round(x,2)

    fn = os.path.join(outpath,'corr_plot.png')
    plt.savefig(fn)

