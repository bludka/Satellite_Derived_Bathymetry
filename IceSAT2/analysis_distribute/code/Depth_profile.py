import os
import pandas as pd
import numpy as np
import math
import geopandas as gpd
import pyproj as proj

from sklearn.cluster import DBSCAN

import Water_level as water_level
import Coral_Reef as coral_reef
import ICESat2_plots as is2_plot
import ICESat2_file as is2_file

#Changes to Karan's file:
#	remove .src from imports, since we are starting notebook in src directory

def create_photon_df(photon_data):
    """
    Creates dataframe from photon data
    Params - 1. photon_data (list) - photon data extracted from h5 file
    Return - DataFrame - containing photon data
    """
    df = pd.DataFrame(photon_data).T
    #sets the column names
    df.columns = ['Height','Latitude','Longitude','Conf_land','Conf_ocean','Conf_inlandwater']
    return df


def individual_confidence(df):
    """
    Gets the confidence of photons for land and ocean
    Params - 1. df (DataFrame) - photon dataframe
    Return - DataFrame - added columns for photon confidence
    """
    # print(df)
    # 0 - land, 1 - ocean, 2 - sea ice, 3 - land ice, 4 - inland water
    df['Conf_land'] = df.apply(lambda x: x.Confidence[0], axis = 1)
    df['Conf_ocean'] = df.apply(lambda x: x.Confidence[1], axis = 1)
    df['Conf_inlandwater'] = df.apply(lambda x: x.Confidence[3], axis = 1)
    return df


def convert_h5_to_csv(is2_file,laser,out_fp):
    """
    Converts photon data from h5 file to a csv file
    Corrects the photon depths making them relative to sea level, adjusted for 
    speed of light in water and tide on the day the satellite orbits over reef.

    Params - 1. is2_file (IS2_file) - object representing icesat file
             2. laser (str) - laser of satellite we want photons for
             3. out_fp (str) - path to store csv file
    Return - DataFrame - photon data over reef
    """
    #creates dataframe with photon data
    photon_data = is2_file.get_photon_data(laser)
    df_laser = create_photon_df(photon_data)

    #Clip dataframe to reef bounding shape
    reef_polygon = reef.return_reef_polygon()
    gdf = gpd.GeoDataFrame(df_laser, geometry=gpd.points_from_xy(df_laser.Longitude, df_laser.Latitude))
    poly_gdf = gpd.GeoDataFrame([1], geometry=[reef_polygon])
    gdf_clipped = gpd.clip(gdf, poly_gdf)
    df = pd.DataFrame(gdf_clipped)
  
    #if there are photons in clipped dataframe
    if len(df) != 0:
        #adjust for sea level, speed of light in water, tide	
        df,f = water_level.normalise_sea_level(df)  #f contains coefficients of polynomial representing sea level
        if len(df) == 0:
            print('No photons')
            return df
        is2_file.set_sea_level_function(f,laser)
        df = water_level.adjust_for_speed_of_light_in_water(df)
        #writes dataframe to a file containing just the photon data that is required
        df.to_csv(out_fp)
    return df


def apply_DBSCAN(df, eps_val, min_samp, water_thresh_low, water_thresh_high, sf):

    empty_df = pd.DataFrame()
    meta = sf.get_meta()
    
    #threshold below sea level for which we consider reefs
    df = df.loc[(df.Height > water_thresh_low) & (df.Height < water_thresh_high)].copy()
    print(len(df))
    if len(df) == 0:
        return empty_df, empty_df

    #getting just high confidence land photons
    #df = df.loc[(df.Conf_land == 4)]
    #if len(df) == 0:
    #    return empty_df, empty_df

    # setup projections
    crs_wgs84 = proj.CRS('EPSG:4326') # assuming you're using WGS84 geographic
    crs_local = proj.CRS(meta['crs'])  # use CRS for sentinel 2 image
    transformer = proj.Transformer.from_crs(crs_wgs84, crs_local)
    df['x'], df['y'] = transformer.transform(df.Latitude.values, df.Longitude.values)

    #initialize dbscan object with algorithm parameters (see https://scikit-learn.org/stable/modules/clustering.html#dbscan)
    dbscan = DBSCAN(eps = eps_val, min_samples = min_samp)

    x = df[['Height','x']]
    if len(x) == 0:
        return empty_df, empty_df
 
    model = dbscan.fit(x)
    labels = model.labels_
    df['labels'] = labels
    clustered_photons = df.loc[df.labels >= 0]
    unclassified_photons = df.loc[df.labels < 0]

    return clustered_photons, unclassified_photons


def depth_profile_adaptive(df, out_path,is2,laser):
    """
    Cleans noisy photons from ICESAT output using an adaptive rolling window
    Params - 1. df (DataFrame) - contains all icesat photons
             2. out_path - path to save cleaned photon data
    Return - DataFrame - contains cleaned icesat photons
    """
    #threshold below sea level for which we consider reefs
    water_level_thresh = 0.5

    #generating arrays for the height, latitude
    h,l = [], []

    #iterating through the latitudes 0.0005 degrees at a time
    start = df.Latitude.min()
    end = df.Latitude.max()
    dx = 0.0005

    #getting just high confidence land photons
    df = df.loc[(df.Conf_land == 4)]
    if len(df) == 0:
        return
    df = df.astype({'Latitude': 'float', 'Longitude':'float'})
    #sorting photons by latitude
    df = df.sort_values('Latitude')

    #getting line of best fit for longitudes using ICESAT photon data
    lon_samples = df.sample(min(1000, len(df)))
    z1 = np.polyfit(lon_samples.Latitude, lon_samples.Longitude,1)
    lon_func = np.poly1d(z1)

    f = is2.get_sea_level_function(laser)
    lats = df.Latitude.drop_duplicates()
    sea = f(lats)
    mean_sea = np.mean(sea)

    while start <= end:
        #subsetting photons that fall within window of size dx
#         temp = df.loc[(df.Latitude >= start-((dx/2))) & (df.Latitude <start+((dx/2)))]
        temp = df.loc[(df.Latitude >= start) & (df.Latitude <start+dx)]
        #getting the midpoint latitude of the window
        mean_lat = (temp.Latitude.max() + temp.Latitude.min()) / 2
        #subsetting coral reef photons
        temp = temp.loc[(temp.Height < f(mean_lat) - mean_sea  - water_level_thresh)]
        if len(temp) == 0:
            start += dx
            continue
        #getting the IQR of photons
        uq = temp["Height"].quantile(0.75)
        lq = temp["Height"].quantile(0.25)
        temp = temp.loc[(temp.Height >= lq) & (temp.Height < uq) ]

        #if IQR contains more than 3 photons we proceed with the depth analysis
        if temp.shape[0] > 3:
            #getting depths that we will iterate throguh
            min_depth = math.ceil(temp.Height.min()-1)
            max_depth = min(0,math.ceil(temp.Height.max()))
            median_depth = pd.DataFrame()

            #iterating through intervals of 0.5m at a time
            for x in range(min_depth,max_depth):
                for y in range(2):
                    #subsetting photons within each 0.5m interval
                    depth_intervals = temp.loc[(temp.Height >= x + (y/2)) & (temp.Height < x+ ((y+1)/2))]
                    #if the interval contains one or more photons we will store the photon information for future calculations
                    if depth_intervals.shape[0] >=1:
                        median_depth = pd.concat([depth_intervals,median_depth])

            #if more than 2 photons are saved from the previous step we set the median to be the predicted height else nan
            if median_depth.shape[0] >= 2:
                h.append(median_depth.Height.median())
            else:
                h.append(np.nan)

            #append the mid point of the latitude to represent the latitude for the calculated height
            l.append((start + start+dx)/2)
            #move to the next window
            start += dx

        #if the IQR does not contain more than 3 photons we use an adaptive window
        else:
            #we have already completed the first iteration by checking if the window has more than 3 photons
            i = 2
            #boolean flag to check if we have met the requirements to make a prediction
            bool_check = False
            #saving the starting latitude
            ts = start
            #adaptive window check done at max 4 times
            while i <= 4:
                #subset data in the window
                temp = df.loc[(df.Latitude >= start-(i*(dx/2))) & (df.Latitude <start+(i*(dx/2)))]
                #get the midpoint of the latitudes in the window
                mean_lat = (temp.Latitude.max() + temp.Latitude.min()) / 2
                #get coral reef photons
                temp = temp.loc[(temp.Height < f(mean_lat) - mean_sea  - water_level_thresh)]

                #setting counter to move to the next adaptive window
                i+=1
                #check if there are more than 30 photons in the window
                if temp.shape[0] > 30:

                    #find depths through which we will iterate
                    min_depth = math.ceil(temp.Height.min()-1)
                    max_depth = min(0,math.ceil(temp.Height.max()))
                    median_depth = pd.DataFrame()
                    #iterate through depths of 0.5m at a time
                    for x in range(min_depth,max_depth):
                        for y in range(2):
                            depth_intervals = temp.loc[(temp.Height >= x + (y/2)) & (temp.Height < x+ ((y+1)/2))]
                            #if any depth interval has 2 or more photons, we will store that information for future use
                            if depth_intervals.shape[0] >=1:
                                median_depth = pd.concat([depth_intervals,median_depth])

                    #if we had more than 2 photons saved from the previous step, we will preduct the median of those photons as the height
                    if median_depth.shape[0] >= 2:
                        #set the boolean flag to true
                        bool_check = True
                        h.append(median_depth.Height.median())
                        #latitude for the given depth prediction is represented by the latitude midpoint
                        l.append((start + start+dx)/2)
                        i = 5

            #if we did not meet any of the criteria we will predict nan for the height and the midpoint of the latitude for lat
            if bool_check == False:
                h.append(np.nan)
                l.append((start + start+dx)/2)
            #move to the next window
            start = ts + dx

    #remove noise
    if h:
        h = remove_depth_noise(h)
        #creating dataframe with the depth and latitudes
        depth = pd.DataFrame([h,l,lon_func(l)]).T
        depth.columns = ['Height','Latitude','Longitude']


        #disregards files with less than 10 depth predictions
        if depth.dropna().shape[0] >= 15:
            depth.to_csv(out_path)
            return depth
    return pd.DataFrame()


def remove_depth_noise(depths):
    """
    Removes noise from cleaned data
    Params - 1. depths (DataFrame) - cleaned photon predictions with some noise
    Return - DataFrame - cleaned from noisy predictions
    """
    prev = depths[0]
    curr = depths[1]
    #only retains photons that have a prediction before and after
    for i in range(1, len(depths) -1):
        next = depths[i+1]
        if np.isnan(prev) and np.isnan(next):
            depths[i] = np.nan
        prev,curr = curr,next
    return depths


def process_h5(reef, is2_file):
    """
    Calculate depth predictions for a single ICESAT-2 file
    Params - 1. reef (Coral_Reef) - reef ICESAT-2 is orbiting over
             2. is2_file (IS2_file) - ICESAT-2 file we are processing
    """
    #gets directories to save outfiles to
    
    icesat_fp, proc_fp, images_fp,data_plots_path = reef.get_file_drectories()

    #get local variables
    reef_name = reef.get_reef_name()
    is2_file_tag = is2_file.get_file_tag()
    strong_beams = is2_file.get_strong_lasers()

    #looping through each strong laser
    for laser in strong_beams:

        #path for csv file containing raw photon data
        photon_fn = '{reef_name}_photons_{h5_fn}_{laser}.csv'.format(reef_name=reef_name, h5_fn=is2_file_tag, laser=laser)
        photons_path = os.path.join(reef.get_icesat_photons_path(), photon_fn)

        #loading raw photon data if it exists, else extracting it from h5 file
        if not os.path.exists(photons_path):
            df = convert_h5_to_csv(is2_file, laser, photons_path)
        else:
            df = pd.read_csv(photons_path)
            is2_file.metadata = is2_file.load_json()

        #if length of photons is zero, move onto next lase  	
        if len(photons) == 0:
            continue
        print('Number of ICESAT-2 Photons in {laser} is {reef_length}'.format(laser=laser, reef_length=str(len(df))))

        #Classify seafloor photons
        eps_val = 0.75
        min_samp = 6
        water_thresh_low = -50
        water_thresh_high = -0.4
        bathy, non_bathy = depth.apply_DBSCAN(df, eps_val, min_samp, water_thresh_low, water_thresh_high)
        if len(bathy) == 0:
            print('No bathymetric photons found for {}'.format(laser))
            continue
        print('Number bathy photons in {laser} is {reef_length}'.format(laser=laser, reef_length=str(len(bathy))))
            
        #saves bathy data
        depths_fn =  '{reef_name}_{h5_fn}_{laser}.csv'.format(reef_name=reef_name, h5_fn=is2_file_tag, laser=laser)
        bathymetry_output_path = os.path.join(reef.get_icesat_bathmetry_path(), depths_fn)
        bathy.to_csv(bathymetry_output_path)

        #Plot depths
        is2_plot.plot_is2_depths_bathy(df, is2_file, laser, bathy, reef.get_icesat_images_path(), 'False', [-25, 5], 'Corrected Bathymetric Profile')


def get_depths(reef):
    """
    Wrapper function that takes in the reef and outputs the depth profile of each reef
    Params - 1. reef (Coral_Reef) - reef over which the ICESAT satelitte is orbitting
    """
    h5_dir = reef.get_icesat_rawdata_path()
    #looping through each h5 file and generating cleaned photon data
    for h5_fn in os.listdir(h5_dir):
        if h5_fn.endswith('.h5'):
            print(h5_fn)
            is2 = is2_file.IS2_file(reef, h5_fn)
            process_h5(reef, is2)