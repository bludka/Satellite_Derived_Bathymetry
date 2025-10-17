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


def apply_DBSCAN(df, eps_val, min_samp, water_thresh_low, water_thresh_high, sf, smooth):
    """
    if smooth is set, it is in units of decimal degrees and is the distance over which photons are aggregated
    """
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
    bathy_photons = df.loc[df.labels >= 0]
    unclassified_photons = df.loc[df.labels < 0]
    
    if smooth != 0:
        
        dx = smooth
        h,lat = [], []

        #get line of best fit for longitudes using ICESAT photon data
        lon_samples = bathy_photons.sample(min(1000, len(bathy_photons)))
        z1 = np.polyfit(lon_samples.Latitude, lon_samples.Longitude,1)
        lon_func = np.poly1d(z1)
      
        #sorting photons by latitude
        bathy_photons = bathy_photons.astype({'Latitude': 'float', 'Longitude':'float'})
        bathy_photons = bathy_photons.sort_values('Latitude')

        #iterate through the latitudes dx degrees at a time
        start = bathy_photons.Latitude.min()
        end = bathy_photons.Latitude.max()

        while start <= end:
        
            #subset photons that fall within window of size dx
            temp = bathy_photons.loc[(bathy_photons.Latitude >= start) & (bathy_photons.Latitude <start+dx)]

            #if analysis cell contains more than 3 photons we proceed with the depth analysis
            if temp.shape[0] >= 3:
                h.append(temp["Height"].median())
                lat.append((start + start+dx)/2)
            #else:
            #    h.append(np.nan)
            #lat.append((start + start+dx)/2)
            start += dx
            
        bathy_smoothed = pd.DataFrame([h,lat,lon_func(lat)]).T
        bathy_smoothed.columns = ['Height','Latitude','Longitude']

        #Apply DBSCAN one more time to be sure we don't get outliers in smoothed bathymetry
        
        #re-initialize dbscan object [optional]
        dbscan = DBSCAN(eps = 3, min_samples = 2)
        
        bathy_smoothed['x'], bathy_smoothed['y'] = transformer.transform(bathy_smoothed.Latitude.values, bathy_smoothed.Longitude.values)
        x = bathy_smoothed[['Height','x']]
        model = dbscan.fit(x)
        labels = model.labels_
        bathy_smoothed['labels'] = labels
        depth = bathy_smoothed.loc[bathy_smoothed.labels >= 0]
        return depth
    else:
        return bathy_photons, unclassified_photons


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