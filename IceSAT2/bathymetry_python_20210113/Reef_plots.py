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

#method to plot the reef depth histogram and a scatter plot of the same
def plot_reefs(fp,data,sf,line):
    df = pd.read_csv(fp)
    reef_name = sf.get_reef_name()
    dt = sf.get_date().strftime("%Y%m%d%H%M%S")
    #plot histogram of depths
    fig, ax = plt.subplots(1,3,figsize = (28,12))
    depth_histogram(ax[0], df, reef_name)
    line_of_best_fit(ax[1],data,sf,line)
    reef_scatter(fig,ax[2],df,reef_name,data)
    fn = '{reef_name}-{date}.png'.format(reef_name = reef_name, date = dt)
    out = os.path.join(sf.get_img_path(), fn)
    plt.savefig(out)
    plt.close(fig)
    return

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

def depth_histogram(ax,df,reef_name, suffix = ''):
    df['Height'].plot.hist(bins = np.arange(-30,20,1), ax = ax)
    ax.set_xlabel('Height (m)')
    ax.set_ylabel('Frequency')
    ax.set_title('{reef_name} {suffix} Depth Histogram'.format(reef_name = reef_name, suffix = suffix))

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


def band_histogram(df):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    plt.hist([df['b2'], df['b3'], df['b4'], df['b8']], bins=np.arange(0,3000,10), histtype='step', stacked=False, color = ['Blue', 'Green', 'Red', 'Black'])
    ax.set_xlabel('Raw Band Intensity')
    ax.set_ylabel('Frequency')
    plt.setp(ax, xlim = [0, 3000], ylim = [0, 1000])

def band_histogram_full(imgs, meta):

    np.seterr(divide='ignore', invalid='ignore')

    blue_flat = np.ndarray.flatten(imgs[0])
    green_flat = np.ndarray.flatten(imgs[1])
    red_flat = np.ndarray.flatten(imgs[2])
    ir_flat = np.ndarray.flatten(imgs[3])

    fig = plt.figure(figsize=(14, 4))
    ax = fig.add_subplot(1, 2, 1)
    plt.hist([blue_flat, green_flat, red_flat, ir_flat], bins=np.arange(0,6000,10), histtype='step', stacked=False, color = ['Blue', 'Green', 'Red', 'Black'])
    ax.set_xlabel('Raw Band Intensity')
    ax.set_ylabel('Frequency')
    ax.set_yscale('log')
    plt.setp(ax, xlim = [0, 6000], ylim = [1, 50000])

    ax = fig.add_subplot(1, 2, 2)
    plt.hist([np.log(blue_flat) - np.log(green_flat)], bins=np.arange(0,6,.02), histtype='step', stacked=False, color = ['Green'])
    plt.hist([np.log(green_flat) - np.log(red_flat)], bins=np.arange(0,6,.02), histtype='step', stacked=False, color = ['Red'])
    ax.set_yscale('log')
    plt.setp(ax, xlim = [0, 4], ylim = [1, 500000])

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
    y_positions = np.arange(0, nx, step_x)
    y_positions = y_positions + (nx - np.max(y_positions))/2
    y_labels = y_positions - np.mean(y_positions)
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(1, 2, 1)
    blue = imgs[0].astype(float)
    blue = np.log(blue)
    green = imgs[1].astype(float)
    green = np.log(green)
    logdiff = blue - green
    
    ax.set_title('Log(Blue) - Log(Green)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = 0, vmax = 2)
    imgplot = show(logdiff, ax = ax, cmap='Blues', norm=norm)

    ax = fig.add_subplot(1, 2, 2)
    red = imgs[2].astype(float)
    red = np.log(red)
    logdiff = green - red
    
    ax.set_title('Log(Green) - Log(Red)', fontsize=16)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    norm = mpl.colors.Normalize(vmin = 0, vmax = 2)
    imgplot = show(logdiff, ax = ax, cmap='Greens', norm=norm)

def reef_scatter(fig,ax,df,reef_name,data, suffix = None):
    #getting just depths between +- 45m
    df = df.loc[(df.Height <= 10) & (df.Height >= -25)]
    #creating a color scale at 5m intervals
    cmap = cm.colors.ListedColormap(['black','navy','mediumblue' ,'blue','royalblue', 'dodgerblue',
                                     'skyblue','limegreen',  'lime' , 'yellow'
                                      ,'orange','tomato','red','firebrick' ,'maroon'])
    bounds = np.arange(-25,11,2.5)
    norm = BoundaryNorm(bounds,cmap.N)
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)

    # create a second axes for the colorbar
    ax2 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
        spacing='proportional', ticks=bounds, boundaries=bounds, format='%.1f')

    #scatter plot of the predicted depths
    pts = ax.scatter(x = df.x, y = df.y, c = df.Height, s= 1, cmap = cmap, norm = norm)
    #scatter plot of the track lines from the ICESAT 2 data
    ax.scatter(x = data.x, y= data.y, s = 3, c = 'black', label = 'ICESAT-2 tracks')
    custom_lines = [Line2D([0], [0], color='black', lw=4)]
    ax.legend(custom_lines, ['ICESAT-2 tracks'])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('{reef_name} {suffix} Depth Predictions (m)'.format(reef_name = reef_name, suffix = suffix))

def plot_sentinel_bands(imgs, meta, df_use, df_toss):

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
    y_positions = np.arange(0, nx, step_x)
    y_positions = y_positions + (nx - np.max(y_positions))/2
    y_labels = y_positions - np.mean(y_positions)
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    #set overall figure size
    fig = plt.figure(figsize=(17, 17))

    #plot Blue band        
    ax = fig.add_subplot(2, 2, 1)
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
    ax.scatter(x = photon_x, y= photon_y, c = 'blue', marker = '.', s = 5)

    #plot Green band        
    ax = fig.add_subplot(2, 2, 2)
    ax.set_title('Band 03 (Green)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[1], ax = ax, cmap='Greens')

    #overplot Green>Blue locations on Green band
    photon_y = [Item[0] for Item in df_toss.image_coords]
    photon_x = [Item[1] for Item in df_toss.image_coords]
    ax.scatter(x = photon_x, y= photon_y, c = 'green', marker = '.', s = 5)

    #plot Red band        
    ax = fig.add_subplot(2, 2, 3)
    ax.set_title('Band 04 (Red)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[2], ax = ax, cmap='Reds')

    #plot IR band        
    ax = fig.add_subplot(2, 2, 4)
    ax.set_title('Band 08 (Infrared)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[3], ax = ax, cmap='Greys')

def plot_sentinel_bands_closeup(imgs, meta, df_use, df_toss):

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
    y_positions = np.arange(0, nx, step_x)
    y_positions = y_positions + (nx - np.max(y_positions))/2
    y_labels = y_positions - np.mean(y_positions)
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    #set overall figure size
    fig = plt.figure(figsize=(17, 17))

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

def plot_sentinel_cal(data, line_slope, line_intercept):

    #set limits
    xlimit = [-18, 0]

    #make line for fitting
    x = np.array(xlimit)
    y = line_slope*x + line_intercept
    
    #set up plotting
    fig = plt.figure(figsize=(17, 10))
    plt.rc('font', size=12)

    #Raw Log(Bands)
    ax = fig.add_subplot(2, 3, 1)
    plt.title('Log intensity Sentinel-2 bands')
    plt.setp(ax, xlim = xlimit, ylim = [0, 8])
    plt.scatter(data.Height, np.log(data.b2), marker='.', color = 'blue', s = 2)  
    plt.scatter(data.Height, np.log(data.b3), marker='.', color = 'green', s = 2)  
    plt.scatter(data.Height, np.log(data.b4), marker='.', color = 'red', s = 2)  

    #Log(Band) Blue/Green against each other
    ax = fig.add_subplot(2, 3, 2)
    plt.title('Comparison of log intensity of Blue vs Green')
    ax.set_xlabel('Blue', fontsize=13)
    ax.set_ylabel('Green', fontsize=13)
    plt.setp(ax, xlim = [5, 8], ylim = [5, 8])
    plt.scatter(np.log(data.b2), np.log(data.b3), color = 'blue', s = 2)
    scale = 1.00
    plt.plot([1,20],[1*scale,20*scale],linestyle='dashed', color = 'black')
    
    #Log(Band) differences    
    ax = fig.add_subplot(2, 3, 3)
    plt.setp(ax, xlim = xlimit, ylim = [-.5, 1.5])
    plt.title('Log band difference Blue vs. Green')
    plt.scatter(data.Height, np.log(data.b2) - np.log(data.b3), marker='.', color = 'aqua', s = 2)  
    #plt2 = plt.scatter(data.Height, np.log(data.b3) - np.log(data.b4), marker='.', color = 'purple')
    plt.plot(x, y, linestyle='dashed', color = 'black')  
    plt.plot([-20,0],[0,0],linestyle='dashed', color = 'black')
    #plt.scatter(np.log(data.b2), np.log(data.b4), color = 'red')
    #plt.scatter(np.log(data.b3), np.log(data.b4), color = 'purple')

    #Log(Band) Green/Red against each other
    ax = fig.add_subplot(2, 3, 5)
    plt.title('Comparison of log intensity of Green vs Red')
    ax.set_xlabel('Green', fontsize=13)
    ax.set_ylabel('Red', fontsize=13)
    plt.setp(ax, xlim = [1, 8], ylim = [1, 8])
    plt.scatter(np.log(data.b3), np.log(data.b4), color = 'blue', s = 2)
    scale = 1.00
    plt.plot([1,20],[1*scale,20*scale],linestyle='dashed', color = 'black')
    
    #Log(Band) differences    
    ax = fig.add_subplot(2, 3, 6)
    plt.setp(ax, xlim = xlimit, ylim = [0, 8])
    plt.title('Log band difference Green vs. Red')
    plt.scatter(data.Height, np.log(data.b3) - np.log(data.b4), marker='.', color = 'aqua', s = 2)  
    #plt2 = plt.scatter(data.Height, np.log(data.b3) - np.log(data.b4), marker='.', color = 'purple')
    plt.plot(x, y, linestyle='dashed', color = 'black')  
    plt.plot([-20,0],[0,0],linestyle='dashed', color = 'black')
    #plt.scatter(np.log(data.b2), np.log(data.b4), color = 'red')
    #plt.scatter(np.log(data.b3), np.log(data.b4), color = 'purple')

def plot_sentinel_bands_simple(imgs, meta):

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
    y_positions = np.arange(0, nx, step_x)
    y_positions = y_positions + (nx - np.max(y_positions))/2
    y_labels = y_positions - np.mean(y_positions)
    utm_y0 = meta['uly'] + meta['ydim']*ny/2
    y_axis_text = 'UTM Y - {:.0f} (m)'.format(utm_y0)

    #set overall figure size
    fig = plt.figure(figsize=(17, 17))

    #plot Blue band        
    ax = fig.add_subplot(2, 2, 1)
    ax.set_title('Band 02 (Blue)', fontsize=16)   
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[0], ax = ax, cmap='Blues')
   
    #plot Green band        
    ax = fig.add_subplot(2, 2, 2)
    ax.set_title('Band 03 (Green)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[1], ax = ax, cmap='Greens')

    #plot Red band        
    ax = fig.add_subplot(2, 2, 3)
    ax.set_title('Band 04 (Red)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[2], ax = ax, cmap='Reds')

    #plot IR band        
    ax = fig.add_subplot(2, 2, 4)
    ax.set_title('Band 08 (Infrared)', fontsize=16)
    ax.set_ylabel(y_axis_text, fontsize=13)
    ax.set_xlabel(x_axis_text, fontsize=13)
    plt.xticks(x_positions, x_labels, fontsize=12)
    plt.yticks(y_positions, y_labels, fontsize=12)
    ax.set_xticklabels(['{:,}'.format(int(x)) for x in x_labels])    
    ax.set_yticklabels(['{:,}'.format(int(y)) for y in y_labels])    
    imgplot = show(imgs[3], ax = ax, cmap='Greys')


def line_of_best_fit(ax,data,sf,line):
    r = 3
    sns.scatterplot(x = data['diff_b2_b3'], y = data.Height, color = 'blue', ax = ax)
    sns.lineplot(x = [-r,r], y = [line(-r),line(r)], color = 'black', ax = ax)
    ax.set_xlabel('Log(Blue Band) - Log(Green Band)')
    ax.set_ylabel('Depth')
    xt = (list(ax.get_xticks()))[:-1]

    for i,x in enumerate(xt):
        xt[i] = np.round(x,2)
    ax.set_title(sf.get_date().strftime("%Y/%m/%d %H:%M:%S") + ' -> tide - ' + str(sf.get_tide()) + 'm')
    xlim = (-r, r)
    ylim = ( -25,0)
    plt.setp(ax, xlim=xlim, ylim=ylim)

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
