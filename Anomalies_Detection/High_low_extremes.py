# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:23:52 2022

# IDENTIFICATION OF HIGH AND LOW EXTREMES FROM SEVERAL VARIABLES 
# INPUTS: catchment ID and THR_dry and THR_wet

@author: amo232
"""
# Set libraries
import os 
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.stats import percentileofscore
import xarray as xr
from shapely.geometry import mapping
from operator import itemgetter
import re
import logging
import logging.handlers
import time
import sys

# Import Modules
from Streamflow import *
from Precipitation import *
from SoilMoisture import *
from SurfaceWater import *
from Groundwater import *
from NDVI import *


def set_logger(arg_start, basin_name, verbose=True):
    """
    Set-up the logging system, exit if this fails
    """
    # assign logger file name and output directory
    datelog = time.ctime()
    datelog = datelog.replace(':', '_')
    reference = f'Basin_loop_number-name_{arg_start}-{basin_name}'


    logfilename = ('logger' + os.sep + reference + '_logfile_' + '.log')

    # create output directory if not exists
    if not os.path.exists('logger'):
        os.makedirs('logger')

    # create logger and set threshold level, report error if fails
    try:
        logger = logging.getLogger(reference)
        logger.setLevel(logging.DEBUG)
    except IOError:
        sys.exit(f'IOERROR: Failed to initialize logger with: {logfilename}')

    # set formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s -'
                                  '%(levelname)s - %(message)s')

    # assign logging handler to report to .log file
    ch = logging.handlers.RotatingFileHandler(logfilename,
                                              maxBytes=10*1024*1024,
                                              backupCount=5)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # assign logging handler to report to terminal
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # start up log message
    logger.info(f'File logging to {logfilename}')

    return logger, ch

#%%
# Set sub-function
## Function of the percentile computed for the whole time series (complete) 
## -> percentile excluding ]0, 1] -> : min values is never attributed to 0 as it is assumed 
## that another min exist which is not contained in the analysed time series)
def Percentile_comp(df, variable):
    arr = df[variable]
    arr_sorted = sorted(arr)
    percentiles = arr.apply(lambda x: percentileofscore(arr_sorted, x))
    percentiles.rename("%s_percentiles" % variable, inplace=True)
    # Add percentile values to your dataframe according to index
    df_m = pd.merge(df, percentiles, left_index=True, right_index=True)
    return(df_m) # return the whole dataframe with the additional column of the percentile values

## Function of the percentile computed grouping the time series per month; "excluding percentile"
def Percentile_month(df, variable):
    percentiles_month = []

    for x, y in (df[variable].groupby(df.date.dt.month)):
        a = y.sort_values()
        # print(y.apply(lambda y: percentileofscore(a, y)))
        percentiles_month.append(y.apply(lambda y: percentileofscore(a, y)))

    df_percentiles_month = pd.merge(df[["date", variable]], pd.concat(
        percentiles_month, axis=0), left_index=True, right_index=True)
    df["%s_percentiles_month" % variable] = df_percentiles_month["%s_y" % variable] 
    return (df["%s_percentiles_month" % variable]) # return the column of the percentile values, ordered according to the index of the provided dataframe 

#%%   
# Set inputs 
## Set the datasets' path
path = '/projects/0/FWC2/Perfect_Storm/'

      
## -> Extract only the selected GSIM station (according to the quality of the basin delineation. at least 20 years of GOOD timeseries and last year greater or equal than 1980)
# Import csv file with basin metadata already filtered according to basin quality and lenght of the GOOD time series 
GSIM_metadata_path = os.path.join(path,'GSIM/GSIM_analysis_out/GSIM_stations_High_Medium_Catch_&_great_240_timeSteps_with_AU.csv')
GSIM_metadata =  pd.read_csv(GSIM_metadata_path)
GSIM_metadata ['last.time.'] = pd.to_datetime(GSIM_metadata ['last.time.'].str.strip(), format='%m/%d/%Y')
GSIM_metadata = GSIM_metadata [GSIM_metadata['last.time.'].dt.year >= 1980] ## Otherwise it would never combine with the other variables 
GSIM_metadata.reset_index(drop=True, inplace = True)

GSIM_ID = GSIM_metadata["gsim.no"].values

# ## Filter GSIM_ID

arg_start = int(sys.argv[1])
GSIM_ID_filt = [GSIM_ID [arg_start]]
logger, ch = set_logger(arg_start, GSIM_ID_filt[0])

## Define start and end year of the analysis 
Start_year = 0 # 1973: time period used in the paper  #1950
End_year = 0 # 2005: time period used in the paper #2020

## Define threshold for dry and wet extremes
thr_dry = 15
thr_wet = 85 

#%%
# Import dataset 

# ------ Surface water -------

path_HydroLakes = os.path.join(path,'Hydro_GR/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp')
ne_HydroLakes = gpd.read_file(path_HydroLakes)
ne_HydroLakes_v2 = ne_HydroLakes[(ne_HydroLakes['Lake_type'] == 1) & (ne_HydroLakes['Lake_area'] >= 1) & (ne_HydroLakes['Lake_area'] < 2)]   

path_HydroGR = os.path.join(path,'Hydro_GR/Hydro_GR_shp/HR_GRanD_WP.shp')
ne_HydroGR = gpd.read_file(path_HydroGR)

# --------- NDVI -------------
# IMPORT CROP AND URBAN MASKS and SAMPLE AND RE-PROJECT 
## We will mask urban areas and crop areas in order to better identify the drougth signal -> crops are masked because we will have a low NDVI 
## after the harvesting so given the fact that we are not considering only the crop growing season we preferred to cut out the crops
logger.info("Import masks and sample for the NDVI analysis")
# Import crops and sample them according to netcdf file (from 1 km to 4 km)
cropland =  os.path.join(path,'NDVI/LandUse/asap_mask_crop_v02.tif')
# Read mask rasters with xarray
rs_cropland = xr.open_rasterio (cropland) #import rioxarray as rio -> rio.open_rasterio(cropland)
# Open one monthly NDVI netcdf file to re-sample the mask accordingly
ds = xr.open_dataset(os.path.join(path,"NDVI/NDVI_month/NDVI_1981-08.nc"), engine="netcdf4")
# Re-sample
rs_cropland_int = rs_cropland.interp_like(ds, method='linear')

# Urban Area
urbanarea =  os.path.join(path,'NDVI/LandUse/URBAN_AREAS/exposure_base_2010_urbkm2.nc')
# Read mask rasters with xarray
rs_urbanarea = xr.open_rasterio (urbanarea)
# Re-sample
rs_urbanarea_int = rs_urbanarea.interp_like(ds, method='linear')
logger.info("Used masks are: %s and %s"%(cropland,urbanarea))

#%%
# Start Loop
for station_GSIM in GSIM_ID_filt: 
    logger.info(station_GSIM)
    
    # Run Modules
    # STREAMFLOW ------------------------------------------------------------------
    logger.info ("----------------Streamflow analysis----------------")
    dt_s, dt_anomaly_s = Streamflow(path, station_GSIM, Start_year, End_year, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger)
    
    ## Export input values and respective percentiles  to excel   
    path_file_inputs = os.path.join(path, 'Output/Streamflow/Streamflow_%s.xlsx'%station_GSIM) 
    dt_s.to_excel(path_file_inputs)
    
    ## Export dt_anomaly as excel
    path_file_outputs = os.path.join(path, 'Output/DF_anomalies/DF_anomalies_%s.xlsx'%station_GSIM) 
    dt_anomaly_s.to_excel(path_file_outputs)
    
    # PRECIPITATION---------------------------------------------------------------
    logger.info ("----------------Precipitation analysis----------------")
    dt_p, dt_anomaly_p = Precipitation(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger) #mm/month 
    
    ## Export Precipitation inputs (extracted according to the analysed catchment)
    path_file_inputs = os.path.join(path, 'Output/Precipitation/Precipitation_%s.xlsx'%station_GSIM)
    dt_p.to_excel(path_file_inputs)
    
    ## Export outputs 
    ## Merge with the main dataset
    path_file_outputs = os.path.join(path,'Output/DF_anomalies/DF_anomalies_%s.xlsx'%station_GSIM) 
    main_dt = pd.read_excel(path_file_outputs)
    ## Add a code line below in order to be sure about the format of date as datetime 
    main_dt['date'] = main_dt.date.astype(str).str[0:10].astype('datetime64[ns]')
    ## Merge dataset
    new_main_dt = pd.merge(main_dt, dt_anomaly_p, how= "left", on="date")
    ## Export
    new_main_dt.to_excel(path_file_outputs)
    
    # SOIL MOISTURE---------------------------------------------------------------
    logger.info ("----------------Soil Moisture analysis----------------")
    dt_sm, dt_anomaly_sm = SoilMoisture(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger) #mm/day 
    
    if len(dt_sm)>0: # Check if there are SM values 
        ## Export Precipitation inputs (extracted according to the analysed catchment)
        path_file_inputs = os.path.join(path, 'Output/SoilMoisture/SoilMoisture_%s.xlsx'%station_GSIM)
        dt_sm.to_excel(path_file_inputs)
        
        ## Export outputs 
        ## Merge with the main dataset
        path_file_outputs = os.path.join(path,'Output/DF_anomalies/DF_anomalies_%s.xlsx'%station_GSIM) 
        main_dt = pd.read_excel(path_file_outputs)
        ## Add a code line below in order to be sure about the format of date as datetime 
        main_dt['date'] = main_dt.date.astype(str).str[0:10].astype('datetime64[ns]')
        ## Merge dataset
        new_main_dt = pd.merge(main_dt, dt_anomaly_sm, how= "left", on="date")
        ## Export
        new_main_dt.to_excel(path_file_outputs)
    
    # SURFACE WATER---------------------------------------------------------------
    logger.info ("----------------Surface Water analysis----------------")
    
    dt_sw, dt_anomaly_sw = SurfaceWater(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, ne_HydroGR, ne_HydroLakes_v2, logger) # [m2]
    
    if len(dt_sw)>0: # Check if there are reservoirs 
        ## Export Precipitation inputs (extracted according to the analysed catchment)
        path_file_inputs = os.path.join(path, 'Output/SurfaceWater/SurfaceWater_%s.xlsx'%station_GSIM)
        dt_sw.to_excel(path_file_inputs)
        
        ## Export outputs 
        ## Merge with the main dataset
        path_file_outputs = os.path.join(path,'Output/DF_anomalies/DF_anomalies_%s.xlsx'%station_GSIM) 
        main_dt = pd.read_excel(path_file_outputs)
        ## Add a code line below in order to be sure about the format of date as datetime 
        main_dt['date'] = main_dt.date.astype(str).str[0:10].astype('datetime64[ns]')
        ## Merge dataset
        new_main_dt = pd.merge(main_dt, dt_anomaly_sw, how= "left", on="date")
        ## Export
        new_main_dt.to_excel(path_file_outputs)
    
    
    # GROUNDWATER---------------------------------------------------------------
    logger.info ("----------------Groundwater analysis----------------")
    dt_gw, dt_anomaly_gw = Groundwater(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger) #index [-]

    if len(dt_gw)>0: # Check if there are GW values    
        ## Export Precipitation inputs (extracted according to the analysed catchment)
        path_file_inputs = os.path.join(path, 'Output/Groundwater/Groundwater_%s.xlsx'%station_GSIM)
        dt_gw.to_excel(path_file_inputs)
        
        ## Export outputs 
        ## Merge with the main dataset
        path_file_outputs = os.path.join(path,'Output/DF_anomalies/DF_anomalies_%s.xlsx'%station_GSIM) 
        main_dt = pd.read_excel(path_file_outputs)
        ## Add a code line below in order to be sure about the format of date as datetime 
        main_dt['date'] = main_dt.date.astype(str).str[0:10].astype('datetime64[ns]')
        ## Merge dataset
        new_main_dt = pd.merge(main_dt, dt_anomaly_gw, how= "left", on="date")
        ## Export
        new_main_dt.to_excel(path_file_outputs)
    
    
    # NDVI---------------------------------------------------------------
    logger.info ("----------------NDVI analysis----------------")
    dt_ndvi, dt_anomaly_ndvi = NDVI(path, station_GSIM, thr_dry, thr_wet, Percentile_month, rs_urbanarea_int, rs_cropland_int, logger) #index [-]
    
    # Check if there are cells in the basin that are not urban or crops: 
    if len(dt_ndvi)>0: # Check if there are reservoirs
        ## Export NDVI inputs (extracted according to the analysed catchment)
        path_file_inputs = os.path.join(path, 'Output/NDVI/NDVI_%s.xlsx'%station_GSIM)
        dt_ndvi.to_excel(path_file_inputs)
        
        ## Export outputs 
        ## Merge with the main dataset
        path_file_outputs = os.path.join(path,'Output/DF_anomalies/DF_anomalies_%s.xlsx'%station_GSIM) 
        main_dt = pd.read_excel(path_file_outputs)
        ## Add a code line below in order to be sure about the format of date as datetime 
        main_dt['date'] = main_dt.date.astype(str).str[0:10].astype('datetime64[ns]')
        ## Merge dataset
        new_main_dt = pd.merge(main_dt, dt_anomaly_ndvi, how= "left", on="date")
        ## Export
        new_main_dt.to_excel(path_file_outputs)
    
    







