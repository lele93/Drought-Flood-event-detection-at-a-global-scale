# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 13:37:03 2022

# Analyse anomalies of variables and identify Flood, Drought to Flood, Drought and Flood and Drought events per each station and 
report them in one excel file 

FOR EACH OF THIS EVENT, EXTRACT:
    - STATION
    - MAX/MIN anomalies (dependent on high/low variable)
    - Events type (Flood, Drought to Flood, Drought and Flood or Drought events)
    - TIME: Dstart, Dend, Fstart, Fend
    - Duration drought, flood

@author: amo232
"""

# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from os.path import isfile, join

# # Set path
# path = r"C:/Users/amo232/surfdrive/Documents/Paper2/Dataset"

# # import the dataset
# station_GSIM = 'AU_0002536'  # "GB_0000057_MEAN_low_flow" #"AU_0002536_Jul_17" #'ZA_0000257' #"AU_0000292_MEAN_low_flow" #"BR_0000004" #'ZA_0000257' #"AU_0000292" # "GB_0000057" #"AU_0002536" #'ZA_0000257'


def Detect_D_F_events (path, station_GSIM, logger):
    
        # dataset = pd.read_excel(os.path.join(path, 'Output/D&F_analysis_OLD/D&F_analysis_%s.xlsx' % station_GSIM))
        dataset = pd.read_excel(os.path.join(path, 'Output/DF_anomalies/DF_anomalies_%s.xlsx' % station_GSIM))
        dataset = dataset.drop(dataset.filter(regex='Unnamed').columns, axis=1)
        dataset['time'] = pd.to_datetime(dataset.date)
        dataset["ID"] = range(0, len(dataset))
        
        # Check which column variables exist
        columns_names_orig = ['Low_prec_int', 'Low_flow_int', 'Low_sm_surf_int', 'Low_sm_root_int', 'Low_sw_extent_reservoirs_int', 'Low_sw_extent_lakes_int', 'Low_gw_int', 'Low_NDVI_int']
        columns_names_dataset = dataset.columns
        ## Identify the variables that are not present in the dataset (because they would be all nan)
        columns_nanValues = list(set(columns_names_orig).difference(columns_names_dataset))
        # ADD these variables as columns with NaN values in order to not have error in the following part of the code
        dataset [columns_nanValues] = np.nan 

        
        # COLUMN INDICATING THE MISSING VARIABLES
        # Create a column indicatig the missing variables
        dataset['missing_variables'] = np.nan
        dataset['missing_variables'] = dataset['missing_variables'].astype(object)
        dataset["no_missing"] = np.nan
        dataset_drop = dataset.drop(columns=['missing_variables', 'dry/wet_month', 'no_missing'])
        
        for index, row in dataset.iterrows():
            missing_variables = dataset_drop.columns[np.unique(np.where(dataset_drop[index:index+1].T.isna())[0])]
            list_missing_values = [x for x in missing_variables.values.tolist() if x.startswith('Low')]
            dataset["missing_variables"][index] = list_missing_values
            dataset["no_missing"][index] = len(list_missing_values)
        
        # REPLACE -9999 with NAN
        dataset_na = pd.merge(dataset.drop(columns='missing_variables').replace(-9999, np.nan), dataset[['missing_variables', 'date']], how='left', on='date')
        
        
        ## DEFINE DROUGHT / FLOOD -------------------------------------------------------------------------------------------------------------------------------
        analysis = dataset_na.copy()
        analysis.reset_index(inplace=True)
        

        # Define consecutive drougth conditions accounting for multiple variables
        compound_drought_conditions = (analysis['Low_flow_int'] < 0)*1 + (analysis['Low_sm_surf_int'] < 0) * 1 + (analysis['Low_sm_root_int'] < 0) * 1 + (analysis['Low_gw_int'] < 0) * 1  # + (analysis['Low_NDVI_int'] <0) *1
        # Substitute 0 with NAN so that in the following part of the script if there are low conditions + low precipitation then add 1 otherwise is NaN (in this way we do not consider the events that have only low precipitation or only low NDVI)
        compound_drought_conditions.replace(0, np.nan, inplace=True)
        analysis["compound_drought_conditions"] = (compound_drought_conditions > 0) * compound_drought_conditions + (analysis['Low_prec_int'] < 0) * 1 + (analysis['Low_NDVI_int'] < 0) * 1  + (analysis['Low_sw_extent_lakes_int'] < 0) * 1 + (analysis['Low_sw_extent_reservoirs_int'] < 0) * 1 
        analysis["0_1"] = np.where(analysis["compound_drought_conditions"] > 0, 1, np.nan)
        # Define consecutive number
        analysis["0_1_cons"] = np.nan  # DROUGHT CONDITION PERIOD (NO. OF MONTHS)
        
        for t in range(1, len(analysis['0_1'])):
            if ((analysis['0_1'][t] > 0)):
                analysis['0_1_cons'][t] = analysis['0_1_cons'].fillna(0)[t-1] + 1
                analysis['0_1_cons'][t-1] = np.nan
            else:
                analysis['0_1_cons'][t] = np.nan
        
        # Compound high flow conditions
        analysis["High_flow_cons"] = np.nan
        for t in range(1, len(analysis['High_flow_cons'])):
            if analysis['High_flow_int'][t] > 0:
                analysis['High_flow_cons'][t] = analysis['High_flow_cons'].fillna(0)[t-1] + 1
                analysis['High_flow_cons'][t-1] = np.nan
            else:
                analysis['High_flow_cons'][t] = np.nan
        
        ## INITIALIZE DATASET -------------------------------------------------------------------------------
        # Initiate list 
        Station = []
        Missing_variables = []
        Event_type = []
        Dstart = []
        Dend = []
        Fstart = []
        Fend = []
        
        dataset_anomalies = pd.DataFrame(columns = ['Low_prec_int', 'Low_flow_int', 'Low_sm_surf_int', 'Low_sm_root_int', 'Low_sw_extent_reservoirs_int', 'Low_sw_extent_lakes_int', 'Low_gw_int',
                                                   'Low_NDVI_int', 'High_prec_int', 'High_flow_int', 'High_sm_surf_int', 'High_sm_root_int', 'High_sw_int',
                                                   'High_gw_int', 'High_NDVI_int', 'D_duration', 'F_duration', 'no_missing_variables', 'Index_anom'])
        
        # Shift High flow up of TIME_LAG row -> for the identification of cascading Drought TO Flood
        time_lag = 1 # CHANGE THIS VALUE TO HAVE DIFFERENT TIME LAG # if euqal to 1 means that you are looking at event where high flow end occurred just after the drought period 
        analysis['High_flow_shift_%d'%(time_lag-1)] = analysis["High_flow_cons"].shift(- time_lag) 
        
        z = 0
        
        for index, row in analysis.iterrows():  # example: [237:238]
            # print(index)
            
            ## DROUGHT CONDITIONS COMPOUNDING WITH HIGH FLOW--------------------------------------------------------------------------------------------------
            # Identify period of high flow with low conditions (considering low conditions characteristics from the moment low condition start until high flow)
            # Create a new column with typology of the dry and wet period when extreme dry and wet occurr in the same month
            
            # Check if drought condition end is positive
            if (row['0_1_cons'] > 0):
                # Identify start of the drought condition
                period_dr = row['0_1_cons']
                start_indx_dr = int(index - period_dr + 1)
                # Check if between the end and start of the drought condition there has been high flow
                High_flow_indx = analysis[start_indx_dr: index + 1][(analysis[start_indx_dr: index+1]["High_flow_cons"]).notna()].index.values
                # Filter columns in order to have only low conditions
                filter_col = [col for col in analysis if col.startswith('Low_')]
                low_ds = analysis[filter_col]
                filter_col_h = [col for col in analysis if col.startswith('High_')]
                high_ds = analysis[filter_col_h]
                # so that high_flow is not reported twice in the new column characteristic
                high_ds.drop(columns=["High_flow_cons", 'High_flow_shift_%d'%(time_lag-1)], inplace=True)
        
                # If there is high flow extract the drought conditions that occurred between the start of the drought and the start of the high flow and then extract the high condition that correspond to the high flow period
                if len(High_flow_indx) > 0:
                    for i in High_flow_indx:
                        # print(i)
                        # Identify the maximum anomalies for the low conditions
                        # Find all the low conditions that are not nan in period between the start of the drought condition and the the end of the drought
                        # you need to traspose the dataframe as you are not using row as variable
                        col_names_low = low_ds.columns[np.unique(np.where(low_ds[start_indx_dr: i + 1].T.notna())[0])] # this is to know only from the start of the drought to the end of the flood, if you want from the start of the drought until the end of the drought than substitute i with index
                        an_low = analysis[col_names_low.values][start_indx_dr: index + 1].min()
                        # Do the same for high values
                        col_names_high = high_ds.columns[np.unique(np.where(high_ds[int(i - analysis['High_flow_cons'][i]) : i+1].T.notna())[0])] # if you instead want to know only from the start of the drought to the start of the flood than substitute index with i 
                        an_high = analysis[col_names_high.values][int(i - analysis['High_flow_cons'][i]) : i+1].max()
                        # # Check if all values in the selected columns are not nan using the dataset dataframe to see if there is any gap in the data
                        # nan_values = (dataset [column_names.values][start_indx_dr: i + 1].isna()).any()
                        # no_nan_values = nan_values[nan_values == "True"].count()
                            
                        # Check if there is any missing value in the period of the identified event 
                        columns_nan_variables_1 = dataset.loc[:, dataset.columns.str.startswith("High")][start_indx_dr : i+1].isna().all() # i: until the end of the flood 
                        
                        # Add -9999 to the variable if we have the data but there is no anomaly 
                        col_anom_var_low_DandF = (dataset.loc[:, dataset.columns.str.startswith("Low")][start_indx_dr: i + 1]).isin([-9999]).all() # i: until the end of the flood
                        col_anom_var_high_DandF = (dataset.loc[:, dataset.columns.str.startswith("High")][int(i - analysis['High_flow_cons'][i]) : i+1]).isin([-9999]).all()
                        col_anom_DandF = np.concatenate([col_anom_var_low_DandF[col_anom_var_low_DandF == True].index.values, col_anom_var_high_DandF[col_anom_var_high_DandF == True].index.values])
                        
                        # Add information dataset 
                        dataset_anomalies.loc[z] = [np.nan] * 19 # 19 is the number of columns in dataset_anomalies
                        
                        if len(col_anom_DandF) >0:
                            for col_a in col_anom_DandF:
                                dataset_anomalies.loc[z] [col_a] =  - 9999
                    
                        for col in col_names_low:
                            dataset_anomalies.loc[z] [col] = an_low [col]
                            
                        for col in col_names_high:
                            dataset_anomalies.loc[z] [col] = an_high [col]
                            
                        Station.append(station_GSIM)
                        Event_type.append('D&F')
                        Dstart.append(analysis.date [start_indx_dr])
                        Dend.append(row ['date'])
                        dataset_anomalies.loc[z] ['D_duration'] =  row['0_1_cons'] 
                        dataset_anomalies.loc[z] ['F_duration'] =  analysis ['High_flow_cons'] [i]
                        Fstart.append(analysis.date [i - analysis['High_flow_cons'][i]])
                        Fend.append(analysis.date [i])
                        Missing_variables.append([i.replace('High_', '') for i in columns_nan_variables_1 [columns_nan_variables_1 == True].index.values])
                        dataset_anomalies.loc[z] ['no_missing_variables'] = columns_nan_variables_1[columns_nan_variables_1 == True].count()
                        dataset_anomalies.loc[z] ['Index_anom'] = i
        
                        z = z+1 
                       
          
            ## DROUGHT CONDITIONS CASCADING WITH HIGH FLOW ---------------------------------------------------------------------------------------------------------
            
            # Look for cascading drought conditions and high flow and write down the respective characteristics (with different time lags)   
            # Start with the high flow 
            
            # Check if high flow shifted according time lag is positive
            if (row['High_flow_shift_%d'%(time_lag-1)]>0):
                # Check if end of drought concide with start of high flow shifted according time lag
                idx_check = index - int(row['High_flow_shift_%d'%(time_lag-1)])+1 # start of the high flow + the time lag in order to see if there is a drought ending in that period
                if (analysis.iloc[idx_check] ['0_1_cons'] > 0) and (idx_check>=0):
                    end_indx_dr = idx_check 
                    start_indx_dr = end_indx_dr - int(analysis.iloc[idx_check] ['0_1_cons'])+1 # since when the period of the drought is equal to 1 it means is in that month 
                    end_indx_fl = index + time_lag
                    start_indx_fl = end_indx_fl - int(row['High_flow_shift_%d'%(time_lag-1)])+1
                    
                   #  Filter columns in order to have only high and low conditions
                    filter_col = [col for col in analysis if col.startswith('Low_')]
                    low_ds = analysis[filter_col]
                    filter_col_h = [col for col in analysis if col.startswith('High_')]
                    high_ds = analysis[filter_col_h]
                    high_ds.drop(columns=["High_flow_cons", 'High_flow_shift_%d'%(time_lag-1)], inplace=True)# so that high_flow is not reported twice in the new column characteristic
        
                    # Identify min and max anomalies (according respectively to low and high conditions)
                    column_names_low = low_ds.columns[np.unique(np.where(low_ds [start_indx_dr : end_indx_dr+1].T.notna())[0])] 
                    anomalies_low = analysis[column_names_low.values][start_indx_dr : end_indx_dr+1].min()
                    # Do the same for high values
                    column_names_high = high_ds.columns[np.unique(np.where(high_ds[(start_indx_fl):(end_indx_fl+1)].T.notna())[0])] 
                    anomalies_high = analysis[column_names_high.values][(start_indx_fl):(end_indx_fl+1)].max()
                    
                    # Check if there is any missing value in the period of the identified event 
                    columns_nan_variables_2 = dataset.loc[:, dataset.columns.str.startswith("High")][start_indx_dr : (end_indx_fl+1)].isna().all() # if a variabe is missing in High then it is missing also in Low so you can check just in one of the two 
                    
                    # Add -9999 to the variable if we have the data but there is no anomaly 
                    col_anom_var_low = (dataset.loc[:, dataset.columns.str.startswith("Low")][start_indx_dr : (end_indx_dr+1)]).isin([-9999]).all()
                    col_anom_var_high = (dataset.loc[:, dataset.columns.str.startswith("High")][start_indx_fl : (end_indx_fl+1)]).isin([-9999]).all()
                    col_anom = np.concatenate([col_anom_var_low[col_anom_var_low == True].index.values, col_anom_var_high[col_anom_var_high == True].index.values])
                    
                    # Add information dataset 
                    dataset_anomalies.loc[z] = [np.nan] * 19
                    
                    if len(col_anom) >0:
                        for col_a in col_anom:
                            dataset_anomalies.loc[z] [col_a] =  - 9999
        
                    for col in column_names_low:
                        dataset_anomalies.loc[z] [col] = anomalies_low [col]
                        
                    for col in column_names_high:
                        dataset_anomalies.loc[z] [col] = anomalies_high [col]
                        
                    Station.append(station_GSIM)
                    Event_type.append('DtoF')
                    Dstart.append(analysis.date [start_indx_dr])
                    Dend.append(analysis.date [end_indx_dr])
                    dataset_anomalies.loc[z] ['D_duration'] =  analysis ['0_1_cons'] [end_indx_dr] 
                    dataset_anomalies.loc[z] ['F_duration'] =  row['High_flow_shift_%d'%(time_lag-1)] 
                    Fstart.append(analysis.date [start_indx_fl])
                    Fend.append(analysis.date [end_indx_fl])
                    Missing_variables.append([i.replace('High_', '') for i in columns_nan_variables_2 [columns_nan_variables_2 == True].index.values])
                    dataset_anomalies.loc[z] ['no_missing_variables'] = columns_nan_variables_2[columns_nan_variables_2 == True].count()
                    dataset_anomalies.loc[z] ['Index_anom'] = end_indx_fl
                    
                    z = z+1
        
        # analysis.drop(analysis.filter(regex='shift').columns, axis=1, inplace=True)     # you need to delete this column before setting another time lag
        
        
            ## FLOOD --------------------------------------------------------------------------------------------------------------------------------------
            # Identify flood events that do not compound or are not preceeded by drought 
            
            # Start with the flood
            # Check if high flow  cons is positive
            if (row['High_flow_cons']>0):
                # Extract start and end of the flood
                strt_flood_indx = int(index - row['High_flow_cons'] + 1)
                end_flood_indx = index
                
                # Check if there is NOT a drought condition within or just before the flood (-1) PS (+1) has been added just to be able to inclue the flood end index
                if analysis['compound_drought_conditions'] [strt_flood_indx - 1 : end_flood_indx +1].isna().all(): 
                    
                    filter_col_h = [col for col in analysis if col.startswith('High_')]
                    high_ds = analysis[filter_col_h]
                    high_ds.drop(columns=["High_flow_cons", 'High_flow_shift_%d'%(time_lag-1)], inplace=True)# so that high_flow is not reported twice in the new column characteristic
        
                    # Identify max anomaly 
                    column_names_high = high_ds.columns[np.unique(np.where(high_ds[(strt_flood_indx):(end_flood_indx+1)].T.notna())[0])] 
                    anomalies_high = analysis[column_names_high.values][(strt_flood_indx):(end_flood_indx+1)].max()
                    
                    # Check if there is any missing value in the period of the identified event 
                    columns_nan_variables_3 = dataset.loc[:, dataset.columns.str.startswith("High")][strt_flood_indx : (end_flood_indx+1)].isna().all()
                    
                    # Add -9999 to the variable if we have the data but there is no anomaly 
                    col_anom_var_high = (dataset.loc[:, dataset.columns.str.startswith("High")][strt_flood_indx : (end_flood_indx+1)]).isin([-9999]).all() # Because if there is an anomaly it will be reported and if it is nan than we will leave it nan
                    col_anom_var_low = (dataset.loc[:, dataset.columns.str.startswith("Low")][strt_flood_indx : (end_flood_indx+1)]).notna().all() # Here we won't report low anomalies so we need to check if there are no nan in the dataset 
                    col_anom = np.concatenate([col_anom_var_low[col_anom_var_low == True].index.values, col_anom_var_high[col_anom_var_high == True].index.values])
                   
                    # Add information dataset 
                    dataset_anomalies.loc[z] = [np.nan] * 19
                    
                    if len(col_anom) >0:
                        for col_a in col_anom:
                            dataset_anomalies.loc[z] [col_a] =  - 9999
                    
                    # Add to dataset 
                    for col in column_names_high:
                        dataset_anomalies.loc[z] [col] = anomalies_high [col]
        
                    Station.append(station_GSIM)
                    Event_type.append('F')
                    Dstart.append(np.nan)
                    Dend.append(np.nan)
                    dataset_anomalies.loc[z] ['D_duration'] =  np.nan 
                    dataset_anomalies.loc[z] ['F_duration'] =  row['High_flow_cons'] 
                    Fstart.append(analysis.date [strt_flood_indx])
                    Fend.append(analysis.date [end_flood_indx])
                    
                    Missing_variables.append([i.replace('High_', '') for i in columns_nan_variables_3 [columns_nan_variables_3 == True].index.values])
                    dataset_anomalies.loc[z] ['no_missing_variables'] = columns_nan_variables_3[columns_nan_variables_3 == True].count()
                    
                    dataset_anomalies.loc[z] ['Index_anom'] = end_flood_indx
                    
                    z = z+1
        
        
            ## DROUGHT --------------------------------------------------------------------------------------------------------------------------------------
            # Identify drought periods that do not compound or are followed by floods
        
            # Start with the drought 
            # Check if drought condition end is positive
            
            if (row['0_1_cons'] > 0):
                # Extract the start and end of the drought 
                strt_dr_indx = int(index - row['0_1_cons'] + 1)
                end_dr_indx = index
        
                # Check if no flood compound or follow the drought compound (+1 because to consider following flood +1 to account for that index in the slicing process)
                if analysis['High_flow_int'] [strt_dr_indx : end_dr_indx +2].isna().all(): 
                    
                    filter_col_l = [col for col in analysis if col.startswith('Low_')]
                    low_ds = analysis[filter_col_l]
        
                    # Identify min anomaly 
                    column_names_low = low_ds.columns[np.unique(np.where(low_ds[(strt_dr_indx):(end_dr_indx+1)].T.notna())[0])] 
                    anomalies_low = analysis[column_names_low.values][(strt_dr_indx):(end_dr_indx+1)].min()
                    
                    # Check if there is any missing value in the period of the identified event 
                    columns_nan_variables = dataset.loc[:, dataset.columns.str.startswith("High")][strt_dr_indx : (end_dr_indx+1)].isna().all()
                    
                    # Add -9999 to the variable if we have the data but there is no anomaly 
                    col_anom_var_high = (dataset.loc[:, dataset.columns.str.startswith("High")][strt_dr_indx : (end_dr_indx+1)]).notna().all()
                    col_anom_var_low = (dataset.loc[:, dataset.columns.str.startswith("Low")][strt_dr_indx : (end_dr_indx+1)]).isin([-9999]).all()
                    col_anom = np.concatenate([col_anom_var_low[col_anom_var_low == True].index.values, col_anom_var_high[col_anom_var_high == True].index.values])
                   
                    # Add information dataset 
                    dataset_anomalies.loc[z] = [np.nan] * 19
                    
                    if len(col_anom) >0:
                        for col_a in col_anom:
                            dataset_anomalies.loc[z] [col_a] =  - 9999
                    
                    # Add to dataset 
                    for col in column_names_low:
                        dataset_anomalies.loc[z] [col] = anomalies_low [col]
        
                    Station.append(station_GSIM)
                    Event_type.append('D')
                    Dstart.append(analysis.date [strt_dr_indx])
                    Dend.append(analysis.date [end_dr_indx])
                    dataset_anomalies.loc[z] ['D_duration'] =  row['0_1_cons']
                    dataset_anomalies.loc[z] ['F_duration'] =  np.nan  
                    Fstart.append(np.nan)
                    Fend.append(np.nan)
                    
                    Missing_variables.append([i.replace('High_', '') for i in columns_nan_variables [columns_nan_variables == True].index.values])
                    dataset_anomalies.loc[z] ['no_missing_variables'] = columns_nan_variables[columns_nan_variables == True].count()
                    
                    dataset_anomalies.loc[z] ['Index_anom'] = end_dr_indx
                    
                    z = z+1
                    
                
                
         
                
        # Inster lists in dataframe
        dataset_anomalies.insert (0, "Station", Station) 
        dataset_anomalies.insert (18, "Dstart", Dstart) 
        dataset_anomalies.insert (19, "Dend", Dend) 
        dataset_anomalies.insert (20, "Fstart", Fstart) 
        dataset_anomalies.insert (21, "Fend", Fend) 
        dataset_anomalies.insert (23, "Missing_variables", Missing_variables) 
        dataset_anomalies.insert (24, "Event_type", Event_type) 
        dataset_anomalies.insert (25, 'no_months_from1980', len(analysis[analysis.date.dt.year>=1980]))
        
        
        # PRINT some information
        ## How many Drought & Flood events have been found in the catchment and over which time frame 
        print ('The no. of D_&_F found is: %d over %d months' %(len(dataset_anomalies[dataset_anomalies['Event_type'] == 'D&F']), len(analysis[analysis['dry/wet_month'].notna()])))

        ## How many Drought to Flood events have been found in the catchment and over which time frame 
        print ('The no. of D_to_F found is: %d  over %d months' %(len(dataset_anomalies[dataset_anomalies['Event_type'] == 'DtoF']), len(analysis[analysis['dry/wet_month'].notna()])))
        
        ## How many Flood events have been found in the catchment and over which time frame 
        print ('The no. of individual F found is: %d  over %d months' %(len(dataset_anomalies[dataset_anomalies['Event_type'] == 'F']), len(analysis[analysis['dry/wet_month'].notna()])))
        
        return(dataset_anomalies)
    
