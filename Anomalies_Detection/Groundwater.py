# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 09:49:53 2022

@author: amo232

From: groundwater_GRACE_anomaly.py
GRACE analysis _Total water storage anomaly - Standardized value 
Dataset: https://edo.jrc.ec.europa.eu/gdo/php/index.php?id=2112
https://data.jrc.ec.europa.eu/dataset/0fd62e28-241f-472c-8966-98744920e181#dataaccess

50 km

"""
import xarray as xr
import pandas as pd
import geopandas
import os 
from os import listdir
from os.path import isfile, join
from shapely.geometry import mapping
import numpy as np 
from operator import itemgetter
from scipy.stats import percentileofscore


#%%
# Define main function 
def Groundwater (path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger):
    ## Import related basin
    ## Define the path of the shapefiles
    filename = os.path.join(path,'GSIM/GSIM_metadata/GSIM_catchments/%s.shp'%station_GSIM.lower())
    basin_Shape = geopandas.read_file(filename, crs="epsg:4326")
    bounds = basin_Shape.geometry.apply(lambda x: x.bounds).tolist()
    # Compute area of the catchment 
    area = basin_Shape['geometry'].to_crs({'proj':'cea'}).area/10**6 #km2
    minx, miny, maxx, maxy = min(bounds, key=itemgetter(0))[0], min(bounds, key=itemgetter(1))[1], max(bounds, key=itemgetter(2))[2], max(bounds, key=itemgetter(3))[3] #https://stackoverflow.com/questions/13145368/find-the-maximum-value-in-a-list-of-tuples-in-python
    
    # Extract list of TWS_GRACE netcdf file in folder
    path_folder = os.path.join(path, 'GRACE/GFZ_GRAVIS_files') #/twsan_m_wld_20220101_20220401_m.nc
    onlyfiles = [f for f in listdir(path_folder) if isfile(join(path_folder, f))]
    gw = pd.DataFrame(columns=['date', 'GW_MEDIAN', 'GW_MEAN', 'GW_MAX', 'GW_MIN'])
    
    median = pd.DataFrame()
    mean = pd.DataFrame()
    maximum = pd.DataFrame()
    minimum = pd.DataFrame()
    
    for namefile in onlyfiles:
        logger.info (namefile)

        path_tmp = os.path.join(path_folder, namefile)
        ds=xr.open_dataset(path_tmp)
        # Change longitude coordinates from 0 - 360 degrees to -180 -- +180
        ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
        ds = ds.sortby(ds.lon)
        # Mask netcdf file according to shapefile basin 
        ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        ds.rio.write_crs("epsg:4326", inplace=True)
    
        
        ## Check if catchment size is larger or smaller than cell size 
        if area.values[0] <=  2500:
            # Catchment is smaller than the cell size
            # Use lat and lon to identify related cell in soil moisture raster
            # Change longitude coordinates from 0 - 360 degrees to -180 -- +180
            # ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            # ds = ds.sortby(ds.lon)
            
            tolerance = 0.5
            clipped = ds.sel(lat=slice(int(maxy)+tolerance, int(miny)-tolerance), lon =slice(int(minx)-tolerance, int(maxx)+tolerance)) # lat goes from the largest to the smallest while lon goes from the smallest to the largest 
        else:
            try: 
                clipped = ds.rio.clip(basin_Shape.geometry.apply(mapping), basin_Shape.crs, drop=True) 
            
            # except:
            except Exception as e:
    
                print('1. basin area > 2500, however:--- %s -----,  we then clip according to sel function providing max and min of lat and lon'%e)
                
                clipped = ds.sel(lat=slice((maxy), (miny)), lon =slice((minx), (maxx))) # lat goes from the largest to the smallest while lon goes from the smallest to the largest
                if (len(clipped.lat)==0) or (len(clipped.lon)==0):
                    print('2. no pixel centroids have been captured with min/max lat/lon -> exception:--- %s -----,  we then use all_touched'%e) #https://gis.stackexchange.com/questions/390777/no-data-found-in-bounds-using-rioxarray-to-clip-geotiff-to-shapefile
                    ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
                    ds.rio.write_crs("epsg:4326", inplace=True)
                    clipped = ds.rio.clip(basin_Shape.geometry.apply(mapping), basin_Shape.crs, all_touched=True, drop=True)
                    

        df = clipped.to_dataframe()
        df = df.reset_index()
        
        median = pd.concat([median, df['tws'].groupby(df.time.dt.to_period('M')).median()])
        mean = pd.concat([mean, df['tws'].groupby(df.time.dt.to_period('M')).mean()])
        maximum = pd.concat([maximum, df['tws'].groupby(df.time.dt.to_period('M')).max()])
        minimum = pd.concat([minimum, df['tws'].groupby(df.time.dt.to_period('M')).min()])

    if (len(median)>0 & ((median.notna()).any())).all(): # check if there are values for at least one time step and if the values are not all nan 
        # Save variables in one dataframe 
        gw["date"] = pd.to_datetime(median.index.astype(str))
        gw["GW_MEDIAN"] = median.values
        gw["GW_MEAN"] = mean.values
        gw["GW_MIN"] = minimum.values
        gw["GW_MAX"] = maximum.values
        
        # Delete nan values
        gw = gw[~gw["GW_MEDIAN"].isnull()]
           
        # Define anomaly and identify GW drought conditions
        gw["Median_percentile_months"] = Percentile_month(gw, "GW_MEDIAN")
        gw["Median_percentile"] = Percentile_comp(gw, "GW_MEDIAN")["GW_MEDIAN_percentiles"]
        gw["Low_gw"] = np.where(gw["Median_percentile_months"] <= thr_dry, 1, np.nan)
        gw["Low_gw_int"] = np.where(gw["Low_gw"] == 1,gw["Median_percentile_months"]- thr_dry, np.nan)
        gw ["High_gw"] = np.where(gw["Median_percentile"]>=thr_wet, 1, np.nan)
        gw ["High_gw_int"] = np.where(gw ["High_gw"] == 1, gw["Median_percentile"] - thr_wet, np.nan) 
      
        # # Filter GW data
        # # Filter the dataframe according to date start and date end 
        # Start_year= 1980
        # End_year = 2016
        # gw_filt= gw[(gw.date.dt.year >= Start_year) & (gw.date.dt.year <= End_year)]
        
        # # Plot
        # fig, ax = plt.subplots(figsize=(40,10))
        # ax.plot(gw_filt.date, gw_filt["GW_MEDIAN"], color='black', ls="-", label="GTotal water storage anomaly [-] monthly")
        # ax2 = ax.twinx()
        # ax2.plot(gw_filt.date, gw_filt["Median_percentile_months"], color='red', ls="-", label="Percentile month")
        # ax2.plot(gw_filt.date, gw_filt["Median_percentile"], color='blue', ls="-", label="Percentile")
        # ax2.fill_between(gw_filt.date,gw_filt["Median_percentile_months"], thr_dry, where=(gw_filt["Median_percentile_months"]<= thr_dry), color='red',alpha=0.5, interpolate=True)  
        # ax2.fill_between(gw_filt.date,gw_filt["Median_percentile"], thr_wet, where=(gw_filt["Median_percentile"]>= thr_wet), color='blue',alpha=0.5, interpolate=True)  
        
        # ax.legend(bbox_to_anchor=(1, 1.05), prop={'size': 20})
        # ax.xaxis.label.set_size(20)
        # plt.xticks(fontsize=30 )
        # plt.yticks(fontsize=30 )
        # ax.set_ylabel("Total water storage anomaly [-]", fontsize=22)
        # ax.tick_params(axis='both', labelsize=25)
        # ax.xaxis.set_major_locator(mdates.YearLocator()) 
        # plt.show()
        
        ## Add groundwater values to the time period defined by the streamflow 
        main_streamflow =  pd.read_excel(os.path.join(path, 'Output/Streamflow/Streamflow_%s.xlsx'%station_GSIM))
        ## Add a code line below in order to be sure about the format of date as datetime 
        main_streamflow['date'] = main_streamflow.date.astype(str).str[0:10].astype('datetime64[ns]')
        ## Merge dataset
        dataset = pd.merge(main_streamflow["date"], gw.fillna(-9999), how="left", left_on =main_streamflow.date.dt.to_period('M'), right_on=gw.date.dt.to_period('M'))
        # Remember: nan values mean that we do not have input values, -9999 means that there is no anomaly 
        ## Clean dataset
        dataset.drop(columns = ["date_y", "key_0"], inplace= True)
        dataset.rename(columns={"date_x" : "date"}, inplace = True)
        dt = dataset[["date", "GW_MIN", "GW_MAX", "GW_MEAN", "GW_MEDIAN", "Median_percentile_months", "Median_percentile"]]
        dt_anomaly = dataset [["date", "Low_gw_int", "High_gw_int"]]
        
    else:
        logger.info('no GW values have been detected')
        dt = []
        dt_anomaly = []

    return (dt, dt_anomaly)
#%%
# # Run function
# dt_gw, dt_anomaly_gw = Groundwater(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp) #index [-]

# ## Export Precipitation inputs (extracted according to the analysed catchment)
# path_file_inputs = os.path.join(path, 'Output/Groundwater/Groundwater_%s.xlsx'%station_GSIM)
# dt_gw.to_excel(path_file_inputs)

# ## Export outputs 
# ## Merge with the main dataset
# # --- (check with the results of the previous script)
# # path_file_outputs = os.path.join(path,'Output/D&F_analysis_%s_Jul_17.xlsx'%station_GSIM)
# # main_dt = pd.read_excel(path_file_outputs)
# # main_dt.rename(columns ={"key_0": "date"}, inplace = True)
# # new_main_dt = pd.merge(main_dt, dt_anomaly_gw, how= "left", left_on="date", right_on = dt_anomaly_gw.date.dt.to_period("M").astype(str))
# # --
# path_file_outputs = os.path.join(path,'Output/D&F_analysis/D&F_analysis_%s.xlsx'%station_GSIM) 
# main_dt = pd.read_excel(path_file_outputs)
# new_main_dt = pd.merge(main_dt, dt_anomaly_gw, how= "left", on="date")
# ## Export
# new_main_dt.to_excel(path_file_outputs)
