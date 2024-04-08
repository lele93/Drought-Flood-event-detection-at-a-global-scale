# -*- coding: utf-8 -*-
"""
From NDVI_analysis_loop_SMN_4km -> SMN : No noise (smoothed) Normalized Difference Vegetation Index

@author: amo232

NDVI loop analysis - 4 km NDVI monthly - extracted from: SMN https://www.star.nesdis.noaa.gov/smcd/emb/vci/VH/vh_browse.php
URBAN AREA: https://ghsl.jrc.ec.europa.eu/download.php?ds=bu
            https://www.nature.com/articles/s41467-020-15788-7 https://dataverse.harvard.edu/dataverse/geospatial_human_dimensions_data

"""
import xarray as xr 
import numpy as np
import os
import glob
import geopandas
from shapely.geometry import mapping
import pandas as pd
#%% 

def NDVI (path, station_GSIM, thr_dry, thr_wet, Percentile_month, rs_urbanarea_int, rs_cropland_int, logger):
    
    #%%    
    # Import Netcdf file names
    nc_files = glob.glob(os.path.join(path,"NDVI/NDVI_month/*.nc"))
    
    # Import basin shapefile 
    filename_bsn = os.path.join(path,'GSIM/GSIM_metadata/GSIM_catchments/%s.shp'%station_GSIM.lower())
    shp = geopandas.read_file(filename_bsn, crs="epsg:4326")
    
    ## Create an empty dataframe 
    ds_ndvi_urb_crop_out = pd.DataFrame(columns=['date', 'NDVI_hist', 'NDVI_MEDIAN', 'NDVI_MAX', 'NDVI_MIN'])
    
    # Loop over each netcdf file, mask it according to the mask layers, extract NDVI values and store them in a dataframe togteher with the datetime information
    logger.info("loop over NDVI files to clip and mask")
    i = 0
    for nc_file in nc_files:  
        logger.info(nc_file)
        # Open NDVI netcdf file
        ds_nc = xr.open_dataset(nc_file, engine="netcdf4")
        date = ds_nc['time'].values[0]
        
        # Clip netcdf file according to basin shapefile
        ds_nc = ds_nc.rio.write_crs(4326)
        ds_nc_clip = ds_nc.rio.clip(shp.geometry.apply(mapping), shp.crs)
        # MASK
        ds_nc_clip_urb = ds_nc_clip.where(rs_urbanarea_int<0.3) # mask out urban area -> if the cell has a urban coverage greater than 30 % it is masked out
        ds_nc_clip_urb_cr = ds_nc_clip_urb.where(rs_cropland_int<50) # rs_cropland_int==0) # mask out crop ->crop layer varies betwen 0 and 200 and represent 0-100% physical value
        # ds_nc_clip_urb_cr.NDVI_mean.plot()
        
        # Convert to dataframe 
        dts_ndvi = ds_nc_clip_urb_cr.to_dataframe()
        
        if ~dts_ndvi.NDVI_mean.isnull().values.all(): # if at least one value is no NAN
            logger.info ('NDVI values have been extracted') 
        
            # Extract one NDVI value rapresentative of the whole basin 
            count, division = pd.np.histogram(dts_ndvi["NDVI_mean"], bins=12, range=(dts_ndvi.NDVI_mean.min(), dts_ndvi.NDVI_mean.max()))
            # Identify the index of the highest count
            idx_max = np.where(count == np.max(count))
            # Identify value
            NDVI_hist = ((division[idx_max] + division[idx_max[0][0] +1])/2)[0]
            NDVI_median = dts_ndvi["NDVI_mean"].median()
            NDVI_max = dts_ndvi["NDVI_mean"].max()
            NDVI_min = dts_ndvi["NDVI_mean"].min()
            
        else:
            logger.info ('the basin is all crop area and/or urban area - hence no NDVI values have been extracted')
            NDVI_hist = np.nan
            NDVI_median = np.nan
            NDVI_max = np.nan
            NDVI_min = np.nan
        
        # Save in the dataframe
        ds_ndvi_urb_crop_out.loc [i] = [date, NDVI_hist, NDVI_median, NDVI_max, NDVI_min]
        i = i+1
      
    # Reading/Interpreting the NDVI values: Values greater than .1 generally denote increasing degrees
    # in the greenness and intensity of vegetation. Values between 0 and .1 are commonly 
    # characteristic of rocks and bare soil, and values less than 0 sometimes indicate water/snow.  
    
    # Check if dataframe is all made of nan or not
    if (~ds_ndvi_urb_crop_out.drop(columns='date').isnull().values.all()): #if no all values are NAN proceed 
        # Compute anmolay according to variable threshold - both for extreme dry and extreme "wet"
        ds_ndvi_urb_crop_out['date'] = pd.to_datetime(ds_ndvi_urb_crop_out['date']) #, format='%Y%m').dt.to_period('M')
        ds_ndvi_urb_crop_out['NDVI_MEDIAN_percentiles_month'] = Percentile_month(ds_ndvi_urb_crop_out, 'NDVI_MEDIAN')
        ds_ndvi_urb_crop_out["Low_NDVI_int"] = np.where(ds_ndvi_urb_crop_out['NDVI_MEDIAN_percentiles_month']<=thr_dry,ds_ndvi_urb_crop_out['NDVI_MEDIAN_percentiles_month'] - thr_dry, np.nan)
        ds_ndvi_urb_crop_out["High_NDVI_int"] = np.where(ds_ndvi_urb_crop_out['NDVI_MEDIAN_percentiles_month']>=thr_wet,ds_ndvi_urb_crop_out['NDVI_MEDIAN_percentiles_month'] - thr_wet, np.nan)
        
        # # PLOT
        # fig, ax = plt.subplots(figsize=(40,10))
        # ax.plot(ds_ndvi_urb_crop_out.date, ds_ndvi_urb_crop_out["NDVI_MEDIAN"], color='black', ls="-", label="NDVI median [-] monthly")
        # ax.plot(ds_ndvi_urb_crop_out.date, ds_ndvi_urb_crop_out['MEDIAN_percentile_month'], color='grey', ls=":", label="NDVI median percentile [-] monthly")
        # ax.fill_between(ds_ndvi_urb_crop_out.date, ds_ndvi_urb_crop_out["MEDIAN_percentile_month"], thr_dry, where=(ds_ndvi_urb_crop_out["MEDIAN_percentile_month"]<= thr_dry), color='red',alpha=0.5, interpolate=True)  
        # ax.fill_between(ds_ndvi_urb_crop_out.date, ds_ndvi_urb_crop_out["MEDIAN_percentile_month"], thr_wet, where=(ds_ndvi_urb_crop_out["MEDIAN_percentile_month"]>= thr_wet), color='red',alpha=0.5, interpolate=True) 
        
        # ax.legend(bbox_to_anchor=(1, 1.05), prop={'size': 20})
        # ax.xaxis.label.set_size(20)
        # plt.xticks(fontsize=30 )
        # plt.yticks(fontsize=30 )
        # ax.set_ylabel("NDVI anomaly [-]", fontsize=22)
        # plt.title('NDVI',  fontsize=25)
        # ax.tick_params(axis='both', labelsize=25)
        # plt.show()
        
        ## Add NDVI anomaly values to the time period defined by the streamflow 
        main_streamflow =  pd.read_excel(os.path.join(path, 'Output/Streamflow/Streamflow_%s.xlsx'%station_GSIM))
        ## Add a code line below in order to be sure about the format of date as datetime 
        main_streamflow['date'] = main_streamflow.date.astype(str).str[0:10].astype('datetime64[ns]')
        ## Merge dataset
        dataset = pd.merge(main_streamflow["date"], ds_ndvi_urb_crop_out.fillna(-9999), how="left", left_on =main_streamflow.date.dt.to_period('M'), right_on=ds_ndvi_urb_crop_out.date.dt.to_period('M'))
        # Remember: nan values mean that we do not have input values, -9999 means that there is no anomaly (this is why we transform nan value that represents no anomalies in -9999)
        ## Clean dataset
        dataset.drop(columns = ["date_y", "key_0"], inplace= True)
        dataset.rename(columns={"date_x" : "date"}, inplace = True)
        dt = dataset[["date", 'NDVI_hist', 'NDVI_MEDIAN', 'NDVI_MAX', 'NDVI_MIN', "NDVI_MEDIAN_percentiles_month"]]
        dt_anomaly = dataset [["date", "Low_NDVI_int", "High_NDVI_int"]]
        
    else:
        logger.info('No NDVI value has been extracted')
        dt = []
        dt_anomaly = []
    
    return(dt, dt_anomaly)

#%%
# # Run function
# dt_ndvi, dt_anomaly_ndvi = NDVI(path, station_GSIM, thr_dry, thr_wet, Percentile_month) # [-]

# ## Export NDVI inputs (extracted according to the analysed catchment)
# path_file_inputs = os.path.join(path, 'Output/NDVI/NDVI_%s.xlsx'%station_GSIM)
# dt_ndvi.to_excel(path_file_inputs)

# ## Export outputs 
# ## Merge with the main dataset
# # --- (check with the results of the previous script)
# # path_file_outputs = os.path.join(path,'Output/D&F_analysis_%s.xlsx'%station_GSIM)
# # main_dt = pd.read_excel(path_file_outputs)
# # main_dt.rename(columns ={"key_0": "date"}, inplace = True)
# # new_main_dt = pd.merge(main_dt, dt_anomaly_ndvi, how= "left", left_on="date", right_on = dt_anomaly_ndvi.date.dt.to_period("M").astype(str))
# # --
# path_file_outputs = os.path.join(path,'Output/D&F_analysis/D&F_analysis_%s.xlsx'%station_GSIM) 
# main_dt = pd.read_excel(path_file_outputs)
# new_main_dt = pd.merge(main_dt, dt_anomaly_ndvi, how= "left", on="date")
# ## Export
# new_main_dt.to_excel(path_file_outputs)
