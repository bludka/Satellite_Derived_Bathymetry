import geopandas as gpd
import glob
import lxml
import numpy as np
import os
import rasterio
import re
import requests

from bs4 import BeautifulSoup
from datetime import datetime
from pathlib import Path
from rasterio import mask
from shapely.geometry import box

#importing other python files
import Tide_API as tide

class Sentinel2_image():
    """
    Class representing a Sentinel-2 image
    """
    def __init__(self,safe_file_path,coords):
        """
        Initialise a Sentinel-2 Image
        Params - 1. safe_file_path (str) - path of safe file
                 2. coords (str) - coords of bounding box around coral reef
        """
        #saving input params as instance variables
        self.bbox_coords = coords
        self.safe_file_path = safe_file_path
        #extracting information about reef and paths
        self.reef_path = os.path.dirname(safe_file_path)
        self.safe_file = os.path.basename(safe_file_path)
        self.reef_name = os.path.basename(self.reef_path)
        
        #generating outpaths
        self.predictions_path = os.path.join(self.reef_path,'Output', 'Depth_Predictions')
        self.imgs_path = os.path.join(self.predictions_path, 'Imgs')
        self.depth_preds_path = os.path.join(self.predictions_path, 'Reef_depth_predictions')
        self.training_data_path = os.path.join(self.predictions_path,'Training_data')
        self.dirs = [self.predictions_path, self.imgs_path, self.depth_preds_path,\
                        self.training_data_path]
        self.create_directories()
        
        #generating path of metadata
        fn = 'MTD_TL.xml'
        self.meta_path = list(Path(self.reef_path).glob('**/' + self.safe_file + '/**/' + fn))[0]
        self.meta = self.get_metadata()
        
        #calculating tide on day of sentinel-2 image
        self.tide_level = tide.get_tide(coords,self.get_date())
        
    def create_directories(self):
        """
        Create output directories
        """
        #if directory does not exist, create it
        for dir in self.dirs:
            if not os.path.exists(dir):
                os.mkdir(dir)

    def get_coords(self):
        """
        Return - [min-x, min-y, max-x, max-y] - bounding box of coral reef
        """
        return self.bbox_coords

    def get_tide(self):
        """
        Return - int - tide on day of sentinel-2 image
        """
        return self.tide_level

    def get_safe_file(self):
        """
        Return - str - name of safefile
        """
        return self.safe_file

    def get_img_path(self):
        """
        Return - str - path to save images
        """
        return self.imgs_path

    def get_file_directories(self):
        """
        Return - list - paths to save outfiles
        """
        return [self.imgs_path, self.depth_preds_path, self.training_data_path]

    def get_reef_name(self):
        """
        Return - str - name of coral reef
        """
        return self.reef_name

    def get_metadata(self):
        """
        Return - dict - containing relevant data (date, crs, dimestions of img)
                        from sentinel-2 image
        """
        #load in metadata of sentinel-2 image
        metadata_file = open(self.meta_path, 'r')
        contents = metadata_file.read()
        #print(contents)

        soup = BeautifulSoup(contents, 'xml')
        meta = {}

        #getting the time of the image and creating a datetime object
        dt = (soup.find('SENSING_TIME').text.split('.')[0].replace('T',''))
        meta['dt'] = datetime.strptime(dt, '%Y-%m-%d%H:%M:%S')

        #getting the crs of the image
        geo_info = soup.find('n1:Geometric_Info')
        meta['crs'] = geo_info.find('HORIZONTAL_CS_CODE').text.lower()
        
        #getting the step of the image in the x and y dircetions
        geo_pos = geo_info.find('Geoposition' , {'resolution':"10"})
        meta['xdim'] = int(geo_pos.find('XDIM').text)
        meta['ydim'] = int(geo_pos.find('YDIM').text)

        metadata_file.close()
        return meta


    def get_date(self):
        """
        Return - datetime - date of the sentinel-2 image
        """
        if 'dt' in self.meta:
            return self.meta['dt']
        else:
            return None

    def get_crs(self):
        """
        Return - dict - crs of sentinel-2 image
        """
        return self.meta['crs']

    def read_gjson(self):
        """
        Loading in geojson of coral reef in the crs of sentinel-2 image
        Return - geodataframe - geometry object of reef
        """
        #load in geojson boundary of coral reef
        fp = os.path.join(self.reef_path, os.path.basename(self.reef_path) +'.geojson')
        df = gpd.read_file(fp)
        #set the current crs to WGS84 lat/lon
        df = df.set_crs(epsg=4326)
        #print(df['geometry'].bounds)
        #transform coords to that of S2 image
        df = df.to_crs(self.meta['crs'])
        #print(df['geometry'].bounds)
        return df

    def get_meta(self):
        """
        Get metadata of Sentinel-2 image
        Return - dict - containing metadata
        """
        return self.meta

    def load_sentinel(self):
        """
        Load in sentinel images
        Return - list where each index holds a band's image
        """
        #get boundary of reef
        geom = self.read_gjson()['geometry']
        #print(geom)
        #get bounding box coordinates of coral reef in [min-x, min-y, max-x, max-y]
        bb = geom.bounds
        #print(bb)
        #get upper left coordinates of bounding box
        self.meta['ulx'] = bb.minx[0]
        self.meta['uly'] = bb.maxy[0]
 
        #select the bands that we want
        bands = ['B02','B03','B04','B08']
        imgs = []
        #loops through the bands
        for b in bands:
            #getting paths for each band image
            img_dir = os.path.dirname(self.meta_path)
            img_path = list(Path(img_dir).glob('**/' + 'IMG_DATA' + '/**/*'+b+'_10m.jp2'))[0]
            #reads in image and crops out the reef
            band = rasterio.open(img_path, driver = 'JP2OpenJPEG')
            #print(band.crs)
            #print(band.bounds)
            out_image, out_transform = mask.mask(band, geom, crop=True)
            imgs.append(out_image)
        self.meta['imgs'] = imgs
        return imgs

    def load_sentinel_full(self):
        """
        Load in sentinel images
        Return - list where each index holds a band's image
        """
        #get boundary of reef
        geom = self.read_gjson()['geometry']
        #print(geom)
        #get bounding box coordinates of coral reef in [min-x, min-y, max-x, max-y]
        bb = geom.bounds
        #print(bb)
        #get upper left coordinates of bounding box
        self.meta['ulx'] = bb.minx[0]
        self.meta['uly'] = bb.maxy[0]
 
        #print(bb.minx[0], bb.miny[0], bb.maxx[0], bb.maxy[0])
        bb_polygon = box(bb.minx[0], bb.miny[0], bb.maxx[0], bb.maxy[0])
        bb_gpd = gpd.GeoDataFrame({"id":1,"geometry":[bb_polygon]})
        #print(bb_gpd.geometry)

        #select the bands that we want
        bands = ['B02','B03','B04','B08']
        imgs = []
        #loops through the bands
        for b in bands:
            #getting paths for each band image
            img_dir = os.path.dirname(self.meta_path)
            img_path = list(Path(img_dir).glob('**/' + 'IMG_DATA' + '/**/*'+b+'_10m.jp2'))[0]
            #reads in image and crops out the reef
            band = rasterio.open(img_path, driver = 'JP2OpenJPEG')
            #print(band.crs)
            #print(band.bounds)
            out_image, out_transform = mask.mask(band, bb_gpd.geometry, crop=True)
            imgs.append(out_image)
        self.meta['imgs'] = imgs
        return imgs
