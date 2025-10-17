import pandas as pd
import numpy as np
import Depth_profile as depth


def get_water_level(df):
    """
    Calculate function to represent sea level
    Params - 1. df (Dataframe) - photon depth, lat, lon
    Return np.poly1d - function representing the sea level
    """

    water,lat = [],[]
    #gets just ocean photons
    df = df.loc[df.Conf_ocean == 4]
    if len(df) == 0:
        return None

    #extract photons associated with ocean surface
    water_thresh_high = 1.25
    water_thresh_low = -1.25
    surface = df.loc[(df.Height > df.Height.median() + water_thresh_low) & (df.Height < df.Height.median() + water_thresh_high)]
    #eps_val = 1
    #min_samp = 5
    #surface, non_surface = depth.apply_DBSCAN(df, eps_val, min_samp, -1.25, 1.25)

    #creating a df with just the latitude and height
    sea_level = pd.DataFrame([surface.Height, surface.Latitude]).T.dropna()
    sea_level.columns = ['water','latitude']

    #fitting line to remaining points
    z = np.polyfit(sea_level.latitude, sea_level.water,1)
    f = np.poly1d(z)

    #getting absolute error for each point
    sea_level['abs_diff'] = np.abs(sea_level.water - f(sea_level.latitude))
    #retaining only points with absolute error less than 2
    sea_level = sea_level.loc[sea_level.abs_diff < 2]
    #fitting a parabolic function to the remaining points
    z2 = np.polyfit(sea_level.latitude, sea_level.water,2)
    f2 = np.poly1d(z2)

    #getting absolute error for each point
    sea_level['abs_diff'] = np.abs(sea_level.water - f2(sea_level.latitude))
    #retaining only points with absolute error less than 2
    sea_level = sea_level.loc[sea_level.abs_diff < 1]
    #fitting a parabolic function to the remaining points
    z3 = np.polyfit(sea_level.latitude, sea_level.water,2)
    f3 = np.poly1d(z3)

    return f3
    

def adjust_for_speed_of_light_in_water(df):
    """
    Adjust photon depth to account for change in the speed of photons
    in air and water.
    Params - 1. df(Dataframe) - photon depth, lat, lon
    Return Dataframe - photon depth, lat, lon
    """
    speed_of_light_air = 299710		#assumes 532 nm wavelength, ir = 1.000273027 (https://emtoolbox.nist.gov/wavelength/ciddor.asp)
    speed_of_light_water = 218826	#assumes 532 nm wavelength, index of refraction in water of 1.37 (Mangano et al., 2013)
    coef = speed_of_light_water / speed_of_light_air
    df['Height'] = df['Height'] * coef
    return df


def normalise_sea_level(df):
    """
    Adjusting photons to reference local instantaneous sea level
    Params - 1. df (Dataframe) - photon depth, lat, lon
    Return - df: Dataframe - photon depth adjusted, lat, lon, sea_surface
             f: np.poly1d - function representing sea level
    """
    #calculating sea level
    f = get_water_level(df)
    if not f:
        return pd.DataFrame(), None

    #adjust photons to sea level
    sea = np.array(f(df.Latitude))
    df = df.assign(sea_surface=sea)
    df['Height']=df['Height']-df['sea_surface']
    
    #restrict depths to less than 10 meters above ocean [WHY?]
    df = df.loc[df.Height < 10]

    return df,f
