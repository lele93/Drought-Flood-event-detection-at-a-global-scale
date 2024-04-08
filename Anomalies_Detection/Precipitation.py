# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:40:16 2022 ---> FOR SNELLIUS

@author: amo232

Respective script: Precipitation_MSWEP_analysis2

Analysis of monthly MSWEP precipitation time series per catchment, applying variable and fixed threshold to the precipitation percentile values
Resolution:  11km X 11km

"""
import xarray as xr
# import rioxarray
import pandas as pd
import geopandas as gpd
import os 
from os import listdir
from os.path import isfile, join
from shapely.geometry import mapping
import numpy as np 
# import matplotlib.pyplot as plt
from operator import itemgetter
from scipy.stats import percentileofscore
# #%%


#%%
def Precipitation (path,station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger):
    ## Import related basin
    ## Define the path of the shapefiles (GSIM basin)
    filename = os.path.join(path,'GSIM/GSIM_metadata/GSIM_catchments/%s.shp'%station_GSIM.lower())
    basin_Shape = gpd.read_file(filename, crs="epsg:4326")
    bounds = basin_Shape.geometry.apply(lambda x: x.bounds).tolist()
    ## Compute area of the catchment 
    area = basin_Shape['geometry'].to_crs({'proj':'cea'}).area/10**6 #km2
    ## Identify max and min for lat and lon
    minx, miny, maxx, maxy = min(bounds, key=itemgetter(0))[0], min(bounds, key=itemgetter(1))[1], max(bounds, key=itemgetter(2))[2], max(bounds, key=itemgetter(3))[3] #https://stackoverflow.com/questions/13145368/find-the-maximum-value-in-a-list-of-tuples-in-python
    ## Import netcdf precipitation file
    path_folder = os.path.join(path, "MSWEP_precipitation/Monthly/") #mm/month 
    
    
    # Extract list of MSWEP netcdf file in folder
    onlyfiles = [f for f in listdir(path_folder) if isfile(join(path_folder, f))]
    precip = pd.DataFrame(columns=['date', 'precip_SUM', 'precip_MEDIAN', 'precip_MEAN', 'precip_MAX', 'precip_MIN'])
    
    ## Loop over the netcdf files (each netcdf file correspond to a certain month and year) to etract one value for the analysed catchment per time period
    i = 0
    for namefile in onlyfiles:
        logger.info (namefile)
        path_tmp = os.path.join(path_folder, namefile)
        ds=xr.open_dataset(path_tmp)
            
        # Clip netcdf file according to shapefile basin 
        ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        ds.rio.write_crs("epsg:4326", inplace=True)
        ## Check catchment size: if below MSWEP resolution than just take the value of the corresponding pixel
        if area.values[0] <=  121:
            # catchment is smaller than the cell size
            # Use lat and lon to identify related cell in precipitation raster
            tolerance = 0.5
            clipped = ds.sel(lat=slice(int(maxy)+tolerance,int(miny)-tolerance), lon =slice(int(minx)-tolerance, int(maxx)+tolerance))
            ## REMEMBER: if after clipping df is empty could be because the order of the max and min lat or lon is wrong. Check ds to see if lat and lon are stored in ascending or descending order!!!     
        else:
            try:
                clipped = ds.rio.clip(basin_Shape.geometry.apply(mapping), basin_Shape.crs, drop=True)

            except Exception as e:
                print('1. basin area > 121, however:--- %s -----,  we then clip according to sel function providing max and min of lat and lon'%e)
                clipped = ds.sel(lat=slice((maxy), (miny)), lon =slice((minx), (maxx))) # lat goes from the largest to the smallest while lon goes from the smallest to the largest
                
                if (len(clipped.lat)==0) or (len(clipped.lon)==0):
                    print('2. no pixel centroids have been captured with min/max lat/lon -> exception:--- %s -----,  we then use all_touched'%e) #https://gis.stackexchange.com/questions/390777/no-data-found-in-bounds-using-rioxarray-to-clip-geotiff-to-shapefile
                    ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
                    ds.rio.write_crs("epsg:4326", inplace=True)
                    clipped = ds.rio.clip(basin_Shape.geometry.apply(mapping), basin_Shape.crs, all_touched=True, drop=True)
        
        df = clipped.to_dataframe()
        df = df.reset_index()
    
        ## Extract min, max and mean precipitation of the catchment under analysis 
        df_sum= df['precipitation'].sum()
        df_median = df['precipitation'].median()
        df_mean = df['precipitation'].mean()
        df_max =  df['precipitation'].max()
        df_min =  df['precipitation'].min()
        df_time = df['time'].drop_duplicates(keep = 'first') [0]
        ## Save variables
        precip.loc [i] = [df_time, df_sum, df_median, df_mean, df_max, df_min]
        i = i+1
        
    # Delete NaN values
    precip_na = precip[~precip["precip_MEDIAN"].isnull()]
    # # Sort values according to ascending date
    # precip = precip_na.sort_values(by="date", ascending=False) 
    # # reset index and delete previous one 
    # precip = precip.reset_index(drop=True)
    
    
    #%%
    ## Extreme dry events according to variable threshold analysis  -- anomaly as negative value
    ## Tranform precipitation values in percentile 
    precip["precip_SUM_percentile_months"] = Percentile_month(precip,"precip_SUM")
    precip["Low_prec"] = np.where(precip["precip_SUM_percentile_months"]<= thr_dry, 1, np.nan)
    precip ["Low_prec_int"] = np.where(precip["Low_prec"] == 1, precip["precip_SUM_percentile_months"] - thr_dry, np.nan )
    
    # Extreme wet according to fixed threshold -- anomaly as positive value
    precip["precip_SUM_percentile"] = Percentile_comp(precip,"precip_SUM") ["precip_SUM_percentiles"]
    precip["High_prec"] = np.where(precip["precip_SUM_percentile"]>= thr_wet, 1, np.nan)
    precip ["High_prec_int"] = np.where(precip["High_prec"] == 1, precip["precip_SUM_percentile"] - thr_wet, np.nan )
    
    ## PLOT
    # # Filter the dataframe according to date start and date end 
    # Start_year= 1980
    # End_year = 2016
    # precip_filt= precip[(precip.date.dt.year >= Start_year) & (precip.date.dt.year <= End_year)]
    
    # # Plot
    # fig, ax = plt.subplots(figsize=(40,10))
    # ax.plot(precip_filt.date, precip_filt["precip_SUM"], color='black', ls="-", label="Precipitation area [mm/month]")
    # ax2 = ax.twinx()
    # ax2.plot(precip_filt.date, precip_filt["precip_SUM_percentile_months"], color='red', ls=":", label="Precipitation_SUM percentile month")
    # ax2.plot(precip_filt.date, precip_filt["precip_SUM_percentile"], color='blue', ls=":", label="Precipitation_SUM percentile")
    
    # ax2.fill_between(precip_filt.date.values, precip_filt["precip_SUM_percentile_months"], thr_dry, where=(precip_filt["precip_SUM_percentile_months"]<= thr_dry), color='red',alpha=0.5, interpolate=True)  
    # ax2.fill_between(precip_filt.date.values, precip_filt["precip_SUM_percentile"], thr_wet, where=(precip_filt["precip_SUM_percentile"]>= thr_wet), color='blue',alpha=0.5, interpolate=True)  
    
    # ax.legend(bbox_to_anchor=(1, 1.05), prop={'size': 20})
    # ax.xaxis.label.set_size(20)
    # plt.xticks(fontsize=30 )
    # plt.yticks(fontsize=30 )
    # ax.set_ylabel("Precipitation [mm/month]", fontsize=22)
    # plt.title('Extreme dry and wet weather events',  fontsize=25)
    # ax.tick_params(axis='both', labelsize=25)
    # # adjust the ticks for the secondary y axes
    # plt.show()
    
    #%% 
    ## Add precipitation values to the time period defined by the streamflow 
    main_streamflow =  pd.read_excel(os.path.join(path, 'Output/Streamflow/Streamflow_%s.xlsx'%station_GSIM))
    ## Add a code line below in order to be sure about the format of date as datetime 
    main_streamflow['date'] = main_streamflow.date.astype(str).str[0:10].astype('datetime64[ns]')
    ## Merge dataset
    dataset = pd.merge(main_streamflow["date"], precip.fillna(-9999), how="left", left_on =main_streamflow.date.dt.to_period('M'), right_on=precip.date.dt.to_period('M'))
    # Remember: nan values mean that we do not have input values, -9999 means that there is no anomaly 
    ## Clean dataset
    dataset.drop(columns = ["date_y", "key_0"], inplace= True)
    dataset.rename(columns={"date_x" : "date"}, inplace = True)
    dt = dataset[['date', 'precip_SUM', 'precip_MEDIAN', 'precip_MEAN', 'precip_MAX', 'precip_MIN', "precip_SUM_percentile_months", "precip_SUM_percentile" ]]
    dt_anomaly = dataset [['date', "Low_prec_int", "High_prec_int"]]
    
    return (dt, dt_anomaly)
#%%

# # Run function
# dt_p, dt_anomaly_p = Precipitation(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp) #mm/month 

# ## Export Precipitation inputs (extracted according to the analysed catchment)
# path_file_inputs = os.path.join(path, 'Output/Precipitation/Precipitation_%s.xlsx'%station_GSIM)
# dt_p.to_excel(path_file_inputs)

# ## Export outputs 
# ## Merge with the main dataset
# path_file_outputs = os.path.join(path,'Output/D&F_analysis/D&F_analysis_%s.xlsx'%station_GSIM) 
# main_dt = pd.read_excel(path_file_outputs)
# new_main_dt = pd.merge(main_dt, dt_anomaly_p, how= "left", on="date")
# ## Export
# new_main_dt.to_excel(path_file_outputs)
