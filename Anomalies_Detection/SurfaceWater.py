# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 18:14:59 2022

@author: amo232

From: Surface_Water_explorer.py
Surface water explorer used to detect drought conditions; from 1984 to 2019 (included) 
"""

import hydroengine as he
# import ogr
import numpy as np
# from openpyxl import load_workbook
import pandas as pd
# import matplotlib.pyplot as plt
import geopandas as gpd
import os 
from operator import itemgetter
from scipy.stats import percentileofscore

# #%%

#%%
def SurfaceWater (path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp, ne_HydroGR, ne_HydroLakes_v2, logger): 
    # Identify ID of the reservoirs contained in the basin 
    ## Import basin GSIM
    filename = os.path.join(path,'GSIM/GSIM_metadata/GSIM_catchments/%s.shp'%station_GSIM.lower())
    basin_Shape = gpd.read_file(filename, crs="epsg:4326")
    bounds = basin_Shape.geometry.apply(lambda x: x.bounds).tolist()
    # Compute area of the basin 
    area = basin_Shape['geometry'].to_crs({'proj':'cea'}).area/10**6 #km2
    minx, miny, maxx, maxy = min(bounds, key=itemgetter(0))[0], min(bounds, key=itemgetter(1))[1], max(bounds, key=itemgetter(2))[2], max(bounds, key=itemgetter(3))[3] #https://stackoverflow.com/questions/13145368/find-the-maximum-value-in-a-list-of-tuples-in-python
    
    # Open the layer with the outline of the lakes/reservoirs
    # path_HydroGR =os.path.join(path,'Hydro_GR/Hydro_GR_shp/HR_GRanD_WP.shp')
    # ne_HydroGR = gpd.read_file(path_HydroGR)
    
    # Are both layers (GSIM basin and lakes shapefiles) in the same CRS?
    # if (ne_HydroGR.crs == basin_Shape.crs):
        # print("Both layers are in the same crs!", ne_HydroGR.crs, basin_Shape.crs)
    
    # Simplify the geometry of the clip extent for faster processing -> Use this with caution as it modifies your data.
    basin_boundary_sim = basin_Shape.simplify(.2, preserve_topology=True)
    
    # Clip data -> lakes layer clipped respect to the basin 
    ne_HydroGR_clip = gpd.clip(ne_HydroGR, basin_boundary_sim)
    
    # Ignore missing/empty geometries
    ne_HydroGR_clip = ne_HydroGR_clip[~ne_HydroGR_clip.is_empty]
    
    # Repeat the same process for HydroLakes
    # Clip data -> lakes layer clipped respect to the basin 
    ne_HydroLK_clip = gpd.clip(ne_HydroLakes_v2, basin_boundary_sim)
    
    # Ignore missing/empty geometries
    ne_HydroLK_clip = ne_HydroLK_clip[~ne_HydroLK_clip.is_empty]
    
    logger.info("In the basin, there are: %d reservoirs" %(len(ne_HydroGR_clip)))
    logger.info("In the basin, there are: %d lakes with area around 1 km2" %(len(ne_HydroLK_clip)))
    
    # Filter the reservoirs that are used for water supply 
    ne_HydroGR_clip_ws = ne_HydroGR_clip[ne_HydroGR_clip['G_MAIN_USE']=="Water supply"]
    
    # Filter the reservoirs that are used for irrigation
    ne_HydroGR_clip_irr = ne_HydroGR_clip[ne_HydroGR_clip['G_MAIN_USE']=="Irrigation"]
    
    logger.info("In the basin, there are: %d reservoirs for water supply and %d reservoirs used for irrigation" %(len(ne_HydroGR_clip_ws), len(ne_HydroGR_clip_irr)))
    
    # Extract the Hydro_ID for the reservoir with the largest volume
    Hylak_id = []
    typology = []
    
    if len(ne_HydroGR_clip_ws) > 0:
        Hylak_id.append(ne_HydroGR_clip_ws.loc[ne_HydroGR_clip_ws['Vol_total'].idxmax()] ['Hylak_id'])
        typology.append("ws")
        
    if len(ne_HydroGR_clip_irr) > 0:
        Hylak_id.append(ne_HydroGR_clip_irr.loc[ne_HydroGR_clip_irr['Vol_total'].idxmax()] ['Hylak_id'])
        typology.append("irr")
        
    # Extract the  Hylak_id of the 3 largest lakes in terms of volume
    if len(ne_HydroLK_clip) > 0:
        extract = ne_HydroLK_clip.loc[ne_HydroLK_clip.sort_values(by = 'Vol_total') [:3].index] ['Hylak_id'].values
        for id_L in extract: 
            Hylak_id.append(id_L)
            typology.append("Lake")

    #%%
    # Compute anomaly of surface water extent
    k = 0
    
    if len(Hylak_id)>0:
       reservoirs = pd.DataFrame()
       lakes = pd.DataFrame()
       
       for i in Hylak_id:
           print(i)
           # Extract surface water time series 
           time_series=he.get_lake_time_series(int(i),'water_area')
           area_series=time_series['water_area'] #[m2]
           df_area_series = pd.DataFrame (area_series, columns = ['area_series'])
   
           # Compute the percentile only if you have at least 20 years of no nan data -> 240 months 
           if len(df_area_series[df_area_series['area_series'] != 0]) > 240: 
               if(( typology[k] == 'ws') or (typology[k] == 'irr')):
                   
                   # Save in dataframe
                   sw_x = pd.DataFrame(time_series)
                   sw_x['date'] = pd.to_datetime(sw_x['time'], unit='ms')
                   sw_x.drop(columns="time", inplace=True)
                   reservoirs ["date"] = sw_x['date'] 
           
                   # Clean the dataset (0 values are NaN) --> CHECK in your previous work with Deltares
                   sw_x = sw_x[sw_x["water_area"] != 0]
               
                   # Estimation of high surface and low surface 
                   # Low surface        
                   sw_x["Percentile_months_%d"%i] = Percentile_month(sw_x,"water_area")
                   sw_x["Low_sw_extent_%d"%i] = np.where(sw_x["Percentile_months_%d"%i]<= thr_dry, 1, np.nan)
                   sw_x["Low_sw_extent_%d_%s_int"%(i,typology[k])] = np.where(sw_x["Low_sw_extent_%d"%i]==1, sw_x["Percentile_months_%d"%i]- thr_dry, np.nan)
                   # I have only put low and not, low and equal because otherwise it will also consider 
                   # low the months in which the reservoir is always empty ##CHECK 
                   
                   # High surface (or max capacity)    
                   sw_x["Percentile_%d"%i] = Percentile_comp(sw_x,"water_area") ["water_area_percentiles"]
                   sw_x["High_sw_extent_%d"%i] = np.where(sw_x["Percentile_%d"%i]>= thr_wet, 1, np.nan)
                   sw_x["High_sw_extent_%d_%s_int"%(i,typology[k])] = np.where(sw_x["High_sw_extent_%d"%i] ==1, sw_x["Percentile_%d"%i] - thr_wet, np.nan)
                   
                   # Merge the dataset to the main one
                   reservoirs = pd.merge(reservoirs, sw_x.fillna(-9999), how="left", on="date")
                   # Remember: nan values mean that we do not have input values, -9999 means that there is no anomaly 
                   reservoirs.rename(columns={"water_area": "water_area_%d"%i}, inplace=True)
                   
               else: 
                   # Save in Lake dataframe 
                   lk_x = pd.DataFrame(time_series)
                   lk_x['date'] = pd.to_datetime(lk_x['time'], unit='ms')
                   lk_x.drop(columns="time", inplace=True)
                   lakes ["date"] = lk_x['date'] 
           
                   # Clean the dataset (0 values are NaN) --> CHECK in your previous work with Deltares
                   lk_x = lk_x[lk_x["water_area"] != 0]
               
                   # Estimation of high surface and low surface 
                   # Low surface        
                   lk_x["Percentile_months_%d"%i] = Percentile_month(lk_x,"water_area")
                   lk_x["Low_sw_extent_%d"%i] = np.where(lk_x["Percentile_months_%d"%i]<= thr_dry, 1, np.nan)
                   lk_x["Low_sw_extent_%d_%s_int"%(i,typology[k])] = np.where(lk_x["Low_sw_extent_%d"%i]==1, lk_x["Percentile_months_%d"%i]- thr_dry, np.nan)
                   # low the months in which the reservoir is always empty ##CHECK 
                   
                   # High surface (or max capacity)      
                   lk_x["Percentile_%d"%i] = Percentile_comp(lk_x,"water_area") ["water_area_percentiles"]
                   lk_x["High_sw_extent_%d"%i] = np.where(lk_x["Percentile_%d"%i]>= thr_wet, 1, np.nan)
                   lk_x["High_sw_extent_%d_%s_int"%(i,typology[k])] = np.where(lk_x["High_sw_extent_%d"%i] ==1, lk_x["Percentile_%d"%i] - thr_wet, np.nan)
                   
                   # Merge the dataset to the main one
                   lakes = pd.merge(lakes, lk_x.fillna(-9999), how="left", on="date")
                   # Remember: nan values mean that we do not have input values, -9999 means that there is no anomaly 
                   lakes.rename(columns={"water_area": "water_area_%d"%i}, inplace=True)
               
           else: 
               logger.info('LENGTH time series: %d'%(len(df_area_series[df_area_series['area_series'] != 0])))
           
           k = k+1
           
       # Compute two average values respectively for the resevoir and the lake surface water extent anomalies
       # Identify columns that have the followin patterns 
          
       columns_pattern = ["Low_sw_extent", "High_sw_extent"] 
       reservoirs_mean = pd.DataFrame()
       lakes_mean = pd.DataFrame()
   
       for col in columns_pattern: 
           reservoirs1 = reservoirs.loc[:, reservoirs.columns.str.startswith(col)]
           reservoirs2 = reservoirs1.loc[:, reservoirs1.columns.str.endswith("_int")]
           
           # In the below lines we consider the mean between anomalies and in case one of them does not have anomaly (-9999) then write no anomaly (-9999)
           # reservoirs_mean ['%s_reservoirs_int'%col] = reservoirs2[reservoirs2 != -9999].mean(axis = 1).values
           # reservoirs_mean ['%s_reservoirs_int'%col]= np.where((reservoirs2 == -9999).any(axis=1), -9999, reservoirs_mean['%s_reservoirs_int'%col])
       
           # In the below lines we consider the mean between anomalies and in case one of them does not have anomaly (-9999) then we compute half of the anomaly 
           reservoirs_mean ['%s_reservoirs_int'%col] = reservoirs2.replace(-9999, 0).mean(axis = 1).values
           reservoirs_mean ['%s_reservoirs_int'%col]= np.where((reservoirs2 == -9999).all(axis=1), -9999, reservoirs_mean['%s_reservoirs_int'%col])
           reservoirs_mean ['%s_reservoirs_int'%col]= np.where(((reservoirs2 == -9999).any(axis=1)) & ((reservoirs2.isna()).any(axis=1)), -9999, reservoirs_mean['%s_reservoirs_int'%col])
       
           lakes1 = lakes.loc[:, lakes.columns.str.startswith(col)]
           lakes2 = lakes1.loc[:, lakes1.columns.str.endswith("_int")]
       
       
           # lakes_mean ['%s_lakes_int'%col] = lakes2[lakes2 != -9999].mean(axis = 1).values
           # lakes_mean ['%s_lakes_int'%col] = np.where((lakes2 == -9999).any(axis=1), -9999, lakes_mean['%s_lakes_int'%col])
       
           # In the below lines we consider the mean between anomalies and in case one of them does not have anomaly (-9999) then we compute half of the anomaly 
           lakes_mean ['%s_lakes_int'%col] = lakes2.replace(-9999, 0).mean(axis = 1).values
           lakes_mean ['%s_lakes_int'%col] = np.where((lakes2 == -9999).all(axis=1), -9999, lakes_mean ['%s_lakes_int'%col])
           lakes_mean ['%s_lakes_int'%col] = np.where(((lakes2 == -9999).any(axis=1)) & ((lakes2.isna()).any(axis=1)), -9999, lakes_mean ['%s_lakes_int'%col])
       
       
       ## Add SurfaceWater values to the time period defined by the streamflow 
       main_streamflow =  pd.read_excel(os.path.join(path, 'Output/Streamflow/Streamflow_%s.xlsx'%station_GSIM))
       ## Add a code line below in order to be sure about the format of date as datetime 
       main_streamflow['date'] = main_streamflow.date.astype(str).str[0:10].astype('datetime64[ns]')
       ## Merge dataset
       # Check if the dataframe are not empty:
       if (len(reservoirs) > 0) & (len(lakes) > 0): 
           
           ## add time to the dataframes
           reservoirs_mean['date'] = reservoirs['date']
           lakes_mean['date'] = lakes['date']
           ## Add SurfaceWater values to the time period defined by the streamflow 
           dataset_0 = pd.merge(main_streamflow["date"], reservoirs_mean, how="left", left_on =main_streamflow.date.dt.to_period('M'), right_on=reservoirs_mean.date.dt.to_period('M'), copy=False)
           dataset_0.drop(columns = ["date_y", "key_0"], inplace= True)
           dataset = pd.merge(dataset_0, lakes_mean, how="left", left_on =dataset_0.date_x.dt.to_period('M'), right_on=lakes_mean.date.dt.to_period('M'))
           
           dataset.drop(columns = ["key_0", 'date'], inplace= True)
           dataset.rename(columns={"date_x" : "date"}, inplace = True)
           dt_anomaly = dataset.filter(like='int')
           dt_anomaly["date"] = dataset ["date"]
           
           dt = pd.merge(reservoirs, lakes, how='left', on= 'date' )
           
           
           
       elif (len(reservoirs) > 0): 
           reservoirs_mean['date'] = reservoirs['date']
           dataset = pd.merge(main_streamflow["date"], reservoirs_mean, how="left", left_on =main_streamflow.date.dt.to_period('M'), right_on=reservoirs_mean.date.dt.to_period('M'), copy=False)
           
           dataset.drop(columns = ["date_y", "key_0"], inplace= True)
           dataset.rename(columns={"date_x" : "date"}, inplace = True)
           dt_anomaly = dataset.filter(like='int')
           dt_anomaly["date"] = dataset ["date"]
           
           dt = reservoirs
           
       elif (len(lakes) > 0): 
           lakes_mean['date'] = lakes['date']
           dataset = pd.merge(main_streamflow["date"], lakes_mean, how="left", left_on =main_streamflow.date.dt.to_period('M'), right_on=lakes_mean.date.dt.to_period('M'), copy=False)
       
           dataset.drop(columns = ["date_y", "key_0"], inplace= True)
           dataset.rename(columns={"date_x" : "date"}, inplace = True)
           dt_anomaly = dataset.filter(like='int')
           dt_anomaly["date"] = dataset ["date"]    
       
           dt = lakes
           
       else:
           logger.info('The lakes and reservoirs detected in the basins do not fit the requirements for long timeseries')
           dt = []
           dt_anomaly = []

    else: 
        logger.info('There are no reservoirs in the basin')
        dt = []
        dt_anomaly = []
      
    return (dt, dt_anomaly)

#%%
# # Export the results 
# # Run function
# dt_sw, dt_anomaly_sw = SurfaceWater(path, station_GSIM, thr_dry, thr_wet, Percentile_month, Percentile_comp) # [m2]

# ## Export Precipitation inputs (extracted according to the analysed catchment)
# path_file_inputs = os.path.join(path, 'Output/SurfaceWater/SurfaceWater_%s.xlsx'%station_GSIM)
# dt_sw.to_excel(path_file_inputs)

# ## Export outputs 
# ## Merge with the main dataset
# # --- (check with the results of the previous script)
# # path_file_outputs = os.path.join(path,'Output/D&F_analysis_%s.xlsx'%station_GSIM)
# # main_dt = pd.read_excel(path_file_outputs)
# # main_dt.rename(columns ={"key_0": "date"}, inplace = True)
# # new_main_dt = pd.merge(main_dt, dt_anomaly_sw, how= "left", left_on="date", right_on = dt_anomaly_sw.date.dt.to_period("M").astype(str))
# # --
# path_file_outputs = os.path.join(path,'Output/D&F_analysis/D&F_analysis_%s.xlsx'%station_GSIM) 
# main_dt = pd.read_excel(path_file_outputs)
# new_main_dt = pd.merge(main_dt, dt_anomaly_sw, how= "left", on="date")
# ## Export
# new_main_dt.to_excel(path_file_outputs)
