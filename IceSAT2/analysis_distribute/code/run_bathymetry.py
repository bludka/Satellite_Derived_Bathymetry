#import ATL03_API as is2
#import Depth_profile as dp
#import Tide_API as ta
#import Pixel_transformation as pt
#import Sentinel_API as sa
#import Coral_Reef as coral_reef

import json

def run_pipeline(working_directory):

    params_f = open(working_directory + '/metadata/project-params.json')
    params = json.load(params_f)
    reef_name = params['reef_name']
    start_date = params['sentinel_start_date']
    end_date = params['sentinel_end_date']
    num_imgs = params['num_sentinel_imgs']
    redownload_is2 = params['redownload_is2']
    download_sentinel = params['download_sentinel']
    earthdata_login = params['earthdata_login']
    earthdata_password = params['earthdata_password']
    sentinel_username = params['sentinel_username']
    sentinel_password = params['sentinel_password']
    params_f.close()
    
    print(reef_name)  
    
    #reef = coral_reef.Coral_Reef(reef_name, working_directory)
    #data_dir = reef.get_data_dir()
    #if redownload_is2:
    #     dp.get_depths(reef)
    #if download_sentinel:
    #    sa.get_sentinel_images(reef, start_date, end_date, num_imgs, sentinel_username, sentinel_password)
    #pt.all_safe_files(reef)
    
