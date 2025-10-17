import geopandas as gpd
import glob
import lxml
import numpy as np
import numpy.matlib
import os
import rasterio
from rasterio.transform import Affine
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
    def __init__(self, reef, safe_file_path, coords):
        """
        Initialise a Sentinel-2 Image
        Params - 1. safe_file_path (str) - path of safe file
                 2. coords (str) - coords of bounding box around coral reef
        """
        #saving input params as instance variables
        self.reef_path = reef.get_reef_dir()
        self.reef_name = reef.get_reef_name()
        self.safe_file_path = safe_file_path
        self.output_path = reef.get_output_dir()
        self.bbox_coords = coords  #coords in WGS84
        
        #extracting information about reef and paths
        self.s2_rawdata_path = os.path.dirname(safe_file_path)
        self.safe_filename = os.path.basename(safe_file_path)
        
        #generating path of metadata
        fn = 'MTD_TL.xml'
        self.meta_path = list(Path(self.s2_rawdata_path).glob('**/' + self.safe_filename + '/**/' + fn))[0]
        self.meta = self.get_metadata()
        
        #load bands
        self.load_sentinel()
        
        #calculating tide on day of sentinel-2 image
        self.tide_level = tide.get_tide(coords,self.get_date(), reef)
        

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

    def get_safe_filename(self):
        """
        Return - str - name of safefile
        """
        return self.safe_filename

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
        meta['crs'] = geo_info.find('HORIZONTAL_CS_CODE').text.lower()	#UTM Zone
        
        #getting the step of the image in the x and y dircetions
        geo_pos = geo_info.find('Geoposition' , {'resolution':"10"})
        meta['xdim'] = int(geo_pos.find('XDIM').text)
        meta['ydim'] = int(geo_pos.find('YDIM').text)
        
        #also grab step of image for 60m resolution for B1 band. Other info the same as the other bands?
        geo_pos = geo_info.find('Geoposition' , {'resolution':"60"})
        meta['xdim60'] = int(geo_pos.find('XDIM').text)
        meta['ydim60'] = int(geo_pos.find('YDIM').text)

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
        fp = os.path.join(self.reef_path, 'metadata', os.path.basename(self.reef_path) +'.geojson')
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
        
        #this was helpful:
        #https://rasterio.readthedocs.io/en/latest/quickstart.html
        
        #get boundary of reef
        geom = self.read_gjson()['geometry']
        #get bounding box coordinates of coral reef in [min-x, min-y, max-x, max-y]
        bb = geom.bounds

        #select the bands that we want
        bands = ['B02','B03','B04','B08','B01'] #note B01 band is 60m resolution, while other bands are 10m res
        imgs = []
        imgs_info = []
        #loops through the bands
        for b in bands:
            #getting paths for each band image
            img_dir = os.path.dirname(self.meta_path)
            if b == 'B01':
                img_path = list(Path(img_dir).glob('**/' + 'IMG_DATA' + '/**/*'+b+'_60m.jp2'))[0]
                #reads in image and crops out the reef
                band = rasterio.open(img_path, driver = 'JP2OpenJPEG')
                out_image, out_transform60 = mask.mask(band, geom, crop=True, all_touched=True)
                imgs.append(out_image)
            else:
                img_path = list(Path(img_dir).glob('**/' + 'IMG_DATA' + '/**/*'+b+'_10m.jp2'))[0]
                #reads in image and crops out the reef
                band = rasterio.open(img_path, driver = 'JP2OpenJPEG')
                out_image, out_transform = mask.mask(band, geom, crop=True, all_touched=True)
                imgs.append(out_image)
                driver = band.driver
                crs = band.crs
        
        #get upper left coordinates of 10m resolution image grid 
        ulx, uly = out_transform*(0,0)
        self.meta['ulx'] = ulx
        self.meta['uly'] = uly
        
        #resample 60m resolution B01 band to 10m res 
        resrat = self.meta['xdim60']/self.meta['xdim'] #resrat = 60/10 = 6
        img4sqz = np.squeeze(imgs[4])
        img4_10m = np.repeat(np.repeat(img4sqz,resrat,axis=0),resrat,axis=1)
        
        #define affine tranform for resampled image
        ulx60, uly60 = out_transform60*(0,0)
        transform = Affine.translation(ulx60, uly60) * Affine.scale(10, 10)
                
        #write new image file for resampled B01 band 
        #where to put this?
        test_path = os.path.join(self.output_path,'test.jp2')
        new_dataset = rasterio.open(test_path, 'w', driver= driver, height=img4_10m.shape[0], width=img4_10m.shape[1], count=1, dtype=img4_10m.dtype, crs=crs, transform=transform)
        new_dataset.write(img4_10m, 1)
        new_dataset.close()
        
        #crop again for a higher resolution crop
        band = rasterio.open(test_path, driver = driver)
        out_image, out_transform60to10 = mask.mask(band, geom, crop=True, all_touched=True)
        
        #replace 60m resolution image with image resampled/recropped to 10m
        imgs[4] = out_image
        
        self.meta['imgs'] = imgs
        return imgs

    def load_sentinel_full(self):
        # old and obsolete?
        """
        Load in sentinel images
        Return - list where each index holds a band's image
        """
        #get boundary of reef
        geom = self.read_gjson()['geometry']
        #get bounding box coordinates of coral reef in [min-x, min-y, max-x, max-y]
        bb = geom.bounds
        #get upper left coordinates of bounding box
        self.meta['ulx'] = bb.minx[0]
        self.meta['uly'] = bb.maxy[0]
 
        #print(bb.minx[0], bb.miny[0], bb.maxx[0], bb.maxy[0])
        bb_polygon = box(bb.minx[0], bb.miny[0], bb.maxx[0], bb.maxy[0])
        bb_gpd = gpd.GeoDataFrame({"id":1,"geometry":[bb_polygon]})

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
            out_image, out_transform = mask.mask(band, bb_gpd.geometry, crop=True)
            imgs.append(out_image)
        self.meta['imgs'] = imgs
        return imgs
