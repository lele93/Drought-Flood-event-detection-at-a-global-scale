# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:29:52 2022

@author: amo232

From: SoilMoisture_analysis.py

soilmoisture analyisis with GLEAM v3.5a monthly timescale (26 km x 26 km)

Variable and fixed threshold analysis (respectively for drought and flood)
"""
import xarray as xr
# import rioxarray
import pandas as pd
import geopandas
import os 
# from os import listdir
# from os.path import isfile, join
from shapely.geometry import mapping
import numpy as np 
# import matplotlib.pyplot as plt
from operator import itemgetter
from scipy.stats import percentileofscore



#%%%
def SoilMoisture (path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger):
    # Identify extreme dry and extreme wet events through the threshold analysis 
    
    path_sm = os.path.join(path, 'GLEAM_soilmoisture/SMsurf_1980-2020_GLEAM_v3.5a_MO.nc')
    path_sm_root = os.path.join(path, 'GLEAM_soilmoisture/SMroot_1980-2020_GLEAM_v3.5a_MO.nc')
    
    var = ["surf", "root"]
    path_x = [path_sm, path_sm_root]
    
    ## Import related basin
    ## Define the path of the shapefiles
    filename = os.path.join(path,'GSIM/GSIM_metadata/GSIM_catchments/%s.shp'%station_GSIM.lower())
    basin_Shape = geopandas.read_file(filename, crs="epsg:4326")
    bounds = basin_Shape.geometry.apply(lambda x: x.bounds).tolist()
    # Compute area of the catchment 
    area = basin_Shape['geometry'].to_crs({'proj':'cea'}).area/10**6 #km2
    ## Extract max and min of both lat and lon
    minx, miny, maxx, maxy = min(bounds, key=itemgetter(0))[0], min(bounds, key=itemgetter(1))[1], max(bounds, key=itemgetter(2))[2], max(bounds, key=itemgetter(3))[3] #https://stackoverflow.com/questions/13145368/find-the-maximum-value-in-a-list-of-tuples-in-python
    
    ## Define empty dataframe
    SoilMoisture = pd.DataFrame()
    
    # Extract Soil Moisture values for the catchment under analysis, respectively for soil moisture surf and soil moisture root
    for i in range(0, len(path_x)):
        logger.info("%s"%var[i])
        ds = xr.open_dataset(path_x[i])
        ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        ds.rio.write_crs("epsg:4326", inplace=True)
        if area.values[0] <=  625:
            # catchment is smaller than the cell size
            # Use lat and lon to identify related cell in soil moisture raster
            tolerance = 0.2
            clipped = ds.sel(lat=slice(int(maxy)+tolerance,int(miny)-tolerance), lon =slice(int(minx)-tolerance, int(maxx)+tolerance))
        else:
            try:
                clipped = ds.rio.clip(basin_Shape.geometry.apply(mapping), basin_Shape.crs, drop=True)

            except Exception as e:
    
                print('1. basin area > 625, however:--- %s -----,  we then clip according to sel function providing max and min of lat and lon'%e)
                
                clipped = ds.sel(lat=slice((maxy), (miny)), lon =slice((minx), (maxx))) # lat goes from the largest to the smallest while lon goes from the smallest to the largest
                if (len(clipped.lat)==0) or (len(clipped.lon)==0):
                    print('2. no pixel centroids have been captured with min/max lat/lon -> exception:--- %s -----,  we then use all_touched'%e) #https://gis.stackexchange.com/questions/390777/no-data-found-in-bounds-using-rioxarray-to-clip-geotiff-to-shapefile
                    ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
                    ds.rio.write_crs("epsg:4326", inplace=True)
                    clipped = ds.rio.clip(basin_Shape.geometry.apply(mapping), basin_Shape.crs, all_touched=True, drop=True)
        
        df = clipped.to_dataframe()
        df = df.reset_index()
    
        # Extract min, max, mean and median for the respective soil moisture variable
        SM_var = df.groupby(df.time).agg({'SM%s'%var[i]:['median', 'mean', 'max', 'min']})
        SM_var = SM_var.reset_index()
        
        # Extreme dry events according to variable threshold analysis
        ## Create 2D dataframe 
        SM_var_x = SM_var ["SM%s"%var[i]] 
        SM_var_x["date"] = SM_var["time"]
        # Delete nan values -----------------------------------------------------------
        SM_var_x = SM_var_x[~SM_var_x["median"].isnull()]
        
        # Check if you have values
        if len(SM_var_x) > 0: 
            print(len(SM_var_x))
            
            ## Convert SM variable in percentile 
            SM_var_x["median_percentiles_month"] = Percentile_month(SM_var_x,"median")
    
            
            SM_var_x["Low_sm"] = np.where(SM_var_x ["median_percentiles_month"]<= thr_dry, 1, np.nan)
            SM_var_x["Low_sm_int"] = np.where(SM_var_x["Low_sm"] == 1, SM_var_x ["median_percentiles_month"]- thr_dry, np.nan)
            
            # Extreme wet
            SM_var_x["median_percentiles"] = Percentile_comp(SM_var_x,"median") ["median_percentiles"]
            SM_var_x["High_sm"] = np.where(SM_var_x["median_percentiles"]>= thr_wet, 1, np.nan)
            SM_var_x["High_sm_int"] = np.where(SM_var_x["High_sm"] == 1, SM_var_x["median_percentiles"] - thr_wet, np.nan)
       
        
            # # Filter the dataframe according to date start and date end 
            # Start_year= 1980
            # End_year = 2016
            # SM_var_filt= SM_var_x[(SM_var_x.date.dt.year >= Start_year) & (SM_var_x.date.dt.year <= End_year)]
            
            # # Plot
            # fig, ax = plt.subplots(figsize=(40,10))
            # ax.plot(SM_var_filt.date, SM_var_filt ["median"], color='black', ls="-", label="SoilMOisture %s median [mm/day] avg month"%var[i])
            # ax2 = ax.twinx()
            # ax2.plot(SM_var_filt.date, SM_var_filt["median_percentiles_month"], color='red', ls=":", label="Percentile months")
            # ax2.plot(SM_var_filt.date, SM_var_filt["median_percentiles"], color='blue', ls=":", label="Percentile")
            
            # ax.fill_between(SM_var_filt.date.values, SM_var_filt ["median_percentiles_month"], thr_dry, where=(SM_var_filt ["median_percentiles_month"]<= thr_dry), color='red',alpha=0.5, interpolate=True)  
            # ax.fill_between(SM_var_filt.date.values, SM_var_filt ["median_percentiles"], thr_wet, where=(SM_var_filt ["median_percentiles"]>= thr_wet), color='blue',alpha=0.5, interpolate=True)  
            
            # ax.legend(bbox_to_anchor=(1, 1.05), prop={'size': 20})
            # ax.xaxis.label.set_size(20)
            # plt.xticks(fontsize=30 )
            # plt.yticks(fontsize=30 )
            # ax.set_ylabel("Soilmoisture [mm/day]", fontsize=22)
            # plt.title('Extreme dry and wet soil moisture',  fontsize=25)
            # ax.tick_params(axis='both', labelsize=25)
            # # adjust the ticks for the secondary y axes
            # plt.show()
            
            # Save inputs
            SoilMoisture ["median_%s"%var[i]] = SM_var_x["median"]
            SoilMoisture ["mean_%s"%var[i]] = SM_var_x["mean"]
            SoilMoisture ["max_%s"%var[i]] = SM_var_x["max"]
            SoilMoisture ["min_%s"%var[i]] = SM_var_x["min"]
            SoilMoisture ["median_%s_percentile_month"%var[i]] = SM_var_x["median_percentiles_month"]
            SoilMoisture ["median_%s_percentile"%var[i]] = SM_var_x["median_percentiles"]
            
            # Save analysis
            SoilMoisture ["Low_sm_%s"%var[i]] = SM_var_x["Low_sm"]
            SoilMoisture ["High_sm_%s"%var[i]] = SM_var_x["High_sm"]
            SoilMoisture ["Low_sm_%s_int"%var[i]] = SM_var_x["Low_sm_int"]
            SoilMoisture ["High_sm_%s_int"%var[i]] = SM_var_x["High_sm_int"]
        
        
    #%%
    if ((SoilMoisture.columns == "median_surf").any() or (SoilMoisture.columns == "median_root").any()).any():

        ## Add precipitation values to the time period defined by the streamflow 
        main_streamflow =  pd.read_excel(os.path.join(path, 'Output/Streamflow/Streamflow_%s.xlsx'%station_GSIM))
        ## Add a code line below in order to be sure about the format of date as datetime 
        main_streamflow['date'] = main_streamflow.date.astype(str).str[0:10].astype('datetime64[ns]')
        ## Merge dataset
        SoilMoisture["date"] = SM_var_x["date"]
        dataset = pd.merge(main_streamflow["date"], SoilMoisture.fillna(-9999), how="left", left_on =main_streamflow.date.dt.to_period('M'), right_on=SoilMoisture.date.dt.to_period('M'))
        # Remember: nan values mean that we do not have input values, -9999 means that there is no anomaly 
        ## Clean dataset
        dataset.drop(columns = ["date_y", "key_0"], inplace= True)
        dataset.rename(columns={"date_x" : "date"}, inplace = True)
        dt = dataset[['date', 'median_surf', 'mean_surf', 'max_surf', 'min_surf', 'median_root', 'mean_root', 'max_root', 'min_root', "median_surf_percentile_month", "median_root_percentile_month", "median_surf_percentile", "median_root_percentile" ]]
        dt_anomaly = dataset [['date', "Low_sm_surf_int", "Low_sm_root_int", "High_sm_surf_int", "High_sm_root_int"]]
    else:
        print('no SM values have been detected. If (%d > 0) & (%d > 0), means that pixels have been detected but the values are nan for the whole time series' %(len(clipped.lat), len(clipped.lon)))
        dt = []
        dt_anomaly = []   
    
    return (dt, dt_anomaly)
#%%

# # Run function
# dt_sm, dt_anomaly_sm = SoilMoisture(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp) #mm/day 

# ## Export Precipitation inputs (extracted according to the analysed catchment)
# path_file_inputs = os.path.join(path, 'Output/SoilMoisture/SoilMoisture_%s.xlsx'%station_GSIM)
# dt_sm.to_excel(path_file_inputs)

# ## Export outputs 
# ## Merge with the main dataset
# #--- (check with the results of the previous script)
# # path_file_outputs = os.path.join(path,'Output/D&F_analysis_%s_Jul_17.xlsx'%station_GSIM)
# # main_dt = pd.read_excel(path_file_outputs)
# # main_dt.rename(columns ={"key_0": "date"}, inplace = True)
# # new_main_dt = pd.merge(main_dt, dt_anomaly_sm, how= "left", left_on="date", right_on = dt_anomaly_sm.date.dt.to_period("M").astype(str))
# # --
# path_file_outputs = os.path.join(path,'Output/D&F_analysis/D&F_analysis_%s.xlsx'%station_GSIM) 
# main_dt = pd.read_excel(path_file_outputs)
# new_main_dt = pd.merge(main_dt, dt_anomaly_sm, how= "left", on="date")
# ## Export
# new_main_dt.to_excel(path_file_outputs)



