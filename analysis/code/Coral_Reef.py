import os
import geopandas as gpd

#CLASS CORAL REEF
class Coral_Reef():
    """
    Class to represent a coral reef
    """
    def __init__(self, reef_name, working_directory):
        """
        Initialise a coral reef object
        Params 1. reef_name (str) - name of the coral reef
               2. working_directory (str) - project directory for this reef
        """
        #stores input variables as class variables
        self.reef_name = reef_name
        self.working_directory = working_directory
        self.reef_dir = os.path.join(self.working_directory, self.reef_name)

        #data/metadata directories
        self.data_dir = os.path.join(self.reef_dir, 'data')
        self.metadata_dir = os.path.join(self.reef_dir, 'metadata')
        self.output_dir = os.path.join(self.reef_dir, 'output')

        #data subdirectories
        #former data_cleaning_path is now just data directory
        self.icesat_rawdata_path = os.path.join(self.data_dir, 'is2_rawdata')
        self.icesat_photons_path = os.path.join(self.data_dir, 'is2_photons')
        self.icesat_bathymetry_path = os.path.join(self.data_dir, 'is2_bathymetry')
        self.icesat_validation_path = os.path.join(self.data_dir, 'is2_validation')
        self.icesat_images_path = os.path.join(self.output_dir, 'is2_plots')
        self.reef_depth_plots = os.path.join(self.output_dir, 'reef_depth_plots')
        self.sentinel_rawdata_path = os.path.join(self.data_dir, 's2_rawdata')
        self.sentinel_depths_path = os.path.join(self.data_dir, 's2_depth_estimates')
        self.training_data_path = os.path.join(self.data_dir, 'training_data')

        self.dirs = [self.data_dir, self.output_dir, \
                       self.icesat_rawdata_path, self.icesat_photons_path, \
                       self.icesat_bathymetry_path, self.icesat_images_path, \
                       self.reef_depth_plots, self.sentinel_rawdata_path, \
                       self.sentinel_depths_path, self.training_data_path, self.icesat_validation_path]
        self.create_directories()

        #stores coordinates of bounding box
        self.bbox_coords = self.get_bounding_box()
        self.reef_polygon = self.get_reef_polygon()


    def get_reef_name(self):
        """
        Return - str - name of coral reef 
        """
        return self.reef_name
        
        
    ## METHODS FOR FILE STRUCTURE
    def get_reef_dir(self):
        """
        Return - str - path of the main reef project directory
        """
        return self.reef_dir

    def get_data_dir(self):
        """
        Return - str - path of the reef data directory
        """
        return self.data_dir

    def get_metadata_dir(self):
        """
        Return - str - path of the reef metadata directory
        """
        return self.metadata_dir

    def get_output_dir(self):
        """
        Return str - path of reef output directory
        """
        return self.output_dir

    def get_icesat_rawdata_path(self):
        """
        Return - str - path for raw ICESAT data
        """
        return self.icesat_rawdata_path

    def get_icesat_photons_path(self):
        """
        Return - str - path for processed ICESAT photons
        """
        return self.icesat_photons_path

    def get_icesat_bathymetry_path(self):
        """
        Return - str - path for processed ICESAT output
        """
        return self.icesat_bathymetry_path

    def get_icesat_images_path(self):
        """
        Return - str - path for processed ICESAT output
        """
        return self.icesat_images_path
        
    def get_reef_depth_plots(self):
        """
        Return - str - path for processed ICESAT output
        """
        return self.reef_depth_plots
        
    def get_sentinel_rawdata_path(self):
        """
        Return - str - path for processed ICESAT output
        """
        return self.sentinel_rawdata_path

    def get_sentinel_depths_path(self):
        """
        Return - str - path for processed ICESAT output
        """
        return self.sentinel_depths_path
        
    def get_training_data_path(self):
        """
        Return - str - path for processed ICESAT output
        """
        return self.training_data_path
        
    def create_directories(self):
        """
        Creates directories for all outfiles
        """
        #checks if directory exists, if not it is created
        for dir in self.dirs:
            if not os.path.exists(dir):
                os.mkdir(dir)
                
    def get_icesat_validation_path(self):
        """
        Return - str - path for processed ICESAT output
        """
        return self.icesat_validation_path


    ## METHODS FOR REEF BOUNDS
    def return_reef_polygon(self):
        """
        Return - str - polygon outlining reef shape
        """
        return self.reef_polygon

    def get_bounding_box(self):
        """
        Gets the coordinates of bounding box around coral reef
        Return - [min-x min-y max-x max-y]
        """
        #loads in geojson of reef into geopandasa
        geojson_fp = os.path.join(self.metadata_dir, self.reef_name + '.geojson')
        reef_gjson = gpd.read_file(geojson_fp)
        #returns coordinates of bounding box around coral reef
        reef_polygon = reef_gjson.geometry[0]
        coords = reef_polygon.bounds
        return coords

    def get_reef_polygon(self):
        """
        Gets the coordinates of bounding box around coral reef
        Return - [min-x min-y max-x max-y]
        """
        #loads in geojson of reef into geopandasa
        geojson_fp = os.path.join(self.metadata_dir, self.reef_name + '.geojson')
        reef_gjson = gpd.read_file(geojson_fp)
        reef_polygon = reef_gjson.geometry[0]
        return reef_polygon
