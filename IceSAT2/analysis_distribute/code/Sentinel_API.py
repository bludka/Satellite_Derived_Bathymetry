from sentinelsat import SentinelAPI, read_geojson, geojson_to_wkt
from datetime import date
import os
import geopandas as gpd
import zipfile


def get_sentinel_images(reef, start_date, end_date, num_images, user,password, cloud_cover_percentage):
    """
    Method to download Sentinel-2 images using Sentinel API
    Params - 1. reef (str) - Coral reef object
             2. start_date (str) - starting date of sentinel images
             3. end_date (str) - end date of sentinel images
             4. num_images (int) - number of sentinel-2 images to download
             5. user (str) - username on scihub.copernicus.eu
             6. password (str) - password on scihub.copernicus.eu
    """

    #load geojson of reef
    reef_gjson_fp = os.path.join(reef.get_metadata_dir(), reef.get_reef_name()+'.geojson')
    reef_footprint = geojson_to_wkt(read_geojson(reef_gjson_fp))

    #query sentinel sat api
    api = SentinelAPI(user, password, 'https://scihub.copernicus.eu/dhus')
    products = api.query(reef_footprint,date = (start_date, end_date), platformname = 'Sentinel-2',\
                         area_relation = 'Intersects',processinglevel = 'Level-2A',\
                         cloudcoverpercentage = cloud_cover_percentage, order_by = 'cloudcoverpercentage')  
    print('Number of scenes: {}'.format(len(products)))

    #creating folder for saving sentinel images
    sentinel_path = os.path.join(reef_path, 'SAFE_files')
    if not os.path.exists(sentinel_path):
        os.makedirs(sentinel_path)

    #downloading num_images
    for i,x in enumerate(products.items()):
        k,v = x[0],x[1]
        if i < num_images:
            print('Downloading image {} of {}'.format(i, len(products)))
            safe_folder = os.path.join(sentinel_path, v['title'] + '.SAFE')
            if not os.path.exists(safe_folder):
                print('Downloading: {}'.format(safe_folder))
                api.download(k, directory_path = sentinel_path)
            else:
                print('{} exists!'.format(safe_folder))
    
    #unzipping files
    for file in os.listdir(sentinel_path):
        if file.endswith('.zip'):
            file_path = os.path.join(sentinel_path, file)
            out_path = os.path.join(sentinel_path, file.split('.')[0])

            if os.path.exists(file_path) and not os.path.exists(out_path):
                with zipfile.ZipFile(file_path,"r") as zip_ref:
                    zip_ref.extractall(sentinel_path)
                os.remove(file_path)
