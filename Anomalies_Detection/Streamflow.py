# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 10:01:53 2022
(Respective script: High_low_flow_GSIM_Marjolein_modified.py)


Application of a generic method for high and low flow identification using streamflow GSIM values. The moethod is taken from (Huijgevoort et al., 2012) :
                                    www.hydrol-earth-syst-sci.net/16/2437/2012 
This general method is a combination of the variable threshold method and the consecutive dry period method
both applied to monthly streamflow values. The methods allow the analysis of ephemeral rivers (usually present in arid climates) 
as well as the analysis of rivers in cold areas where runoff is equal to zero when temperature is below zero). We applied a few changes to the oroginal code:
    
    - For checking which methods to use (combined or TLM), here I use the percentile grouped by month and then check if it is above or below 5
      (in the original code,  checking of 0 values was done considering the whole time series) 
    - The "perc_positive_runoff_month" has been computed by accounting only for the positive streamflow, while in the original script all discharge values have been considered. 
    - Fdry will be multiplied to percentile values corresponding to 0 discharge, while Fwet will be multiplied to percentile values corresponding to discharge value greater than zero
      In the original script instead, Fdry was multiplied to the percentile values correspondng to the combined serie values equal to 0. 
      Fwet, instead, was multiplied to values corresponding to combined serie values grater than 0

    
    
    In general:
        - EVERYTHING IS FIRST CONVERTED TO PERCENTILE AND THEN HIGH AND LOW FLOW HAVE BEEN IDENTIFIED (LOW FLOW ACCORDING TO THE COMBINED PERCENTILE IN CASE OF EPHEMERAL RIVER)
        This means we are using the excluding percentile ]0,1]
 
---       
Dataset used: GSIM streamflow data 
---
        
@author: amo232
"""
# Set libraries
import os.path
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore

# #%%
### INPUTS 

#%% Create the streamflow function

def Streamflow (path, station_GSIM,Start_year, End_year, thr_dry, thr_wet, Percentile_month, Percentile_comp, logger):
    
    #%%
    # Import GSIM dataset' path
    GSIM_path = os.path.join(path, 'GSIM/GSIM_indices/TIMESERIES/monthly')
    
    # Import GSIM dataset
    GSIM = pd.read_csv(os.path.join(GSIM_path, '%s.mon' %
                       station_GSIM), skiprows=21)
    GSIM.columns.values
    # Clean dataset and convert object to float
    if (GSIM['\t"MEAN"'].astype(str).str[0:1] == '\t').any():
        GSIM['\t"MEAN"'] = pd.to_numeric(GSIM['\t"MEAN"'].map(
            lambda x: x.lstrip('\t')), errors='coerce').astype('float')
        GSIM['\t"MAX"'] = pd.to_numeric(GSIM['\t"MAX"'].map(
            lambda x: x.lstrip('\t')), errors='coerce').astype('float')
        GSIM['\t"MIN"'] = pd.to_numeric(GSIM['\t"MIN"'].map(
            lambda x: x.lstrip('\t')), errors='coerce').astype('float')
    # Convert date from object to datetime
    GSIM['date'] = pd.to_datetime(GSIM['date'], format='%Y-%m-%d')
    
    # Delete rows with NaN values
    GSIM_cl = GSIM.dropna(subset=['\t"MEAN"'])
    GSIM_cl.reset_index(drop=True, inplace=True)
    
    # %%
    # Filter the data for the time period that you are interested
    ## if we assigned 0 value to the start and end year then assign the respecitive initial and final year  
    if (Start_year == 0) & (End_year == 0):
        Start_year = GSIM_cl.date.dt.year.values[0]
        End_year = GSIM_cl.date.dt.year.values[-1]
    
    GSIM_filt = GSIM_cl[(GSIM_cl.date.dt.year >= Start_year)
                        & (GSIM_cl.date.dt.year <= End_year)]
    ## Reset index
    GSIM_filt.reset_index(drop=True, inplace=True)
    
    #%% 
    # Determine each month if it is more wet or dry
    
    ## Compute the mean monthly discharge regime: the mean discharge per month throughout the whole timeseries
    GSIM_regime = GSIM_filt[['\t"MEAN"', "date"]].groupby(
        GSIM_filt.date.dt.month).aggregate(({'date': 'first', '\t"MEAN"': 'mean'}))
    
    ## Add the results to the whole timeseries
    mean_year = pd.merge(GSIM_filt["date"], GSIM_regime, how="left",
                         left_on=GSIM_filt.date.dt.month, right_on=GSIM_regime.index)
    GSIM_filt['mean_years'] = mean_year['\t"MEAN"'] # we can directly put them equal as they have the same index
    
    ## Compute the difference between min/ max and the mean - the one with the lower difference will tell if dry or wet month
    dry = GSIM_filt['\t"MEAN"'] - GSIM_filt['\t"MIN"']
    wet = GSIM_filt['\t"MAX"'] - GSIM_filt['\t"MEAN"']
    GSIM_filt["dry/wet_month"] = np.where((dry <= wet) | (GSIM_filt['\t"MAX"']
                                        < GSIM_filt['mean_years']), "dry", "wet")  
    # We also added the condition that the maximum values has to be lower than the average monthly discharge for that month: 
        # this addition is becouse MEAN might be close to the min but still the min be a value above the average for that year 
    
    # # %% -------------------------
    # # Define functions
    
    # # Function of the percentile computed for the whole time series (complete) 
    # # -> percentile excluding ]0, 1] -> : min values is never attributed to 0 as it is assumed 
    # # that another min exist which is not contained in the analysed time series)
    # def Percentile_comp(df, variable):
    #     arr = df[variable]
    #     arr_sorted = sorted(arr)
    #     percentiles = arr.apply(lambda x: percentileofscore(arr_sorted, x))
    #     percentiles.rename("%s_percentiles" % variable, inplace=True)
    #     # Add percentile values to your dataframe according to index
    #     df_m = pd.merge(df, percentiles, left_index=True, right_index=True)
    #     return(df_m) # return the whole dataframe with the additional column of the percentile values
    
    # # Function of the percentile computed grouping the time series per month; "excluding percentile"
    # def Percentile_month(df, variable):
    #     percentiles_month = []
    
    #     for x, y in (df[variable].groupby(df.date.dt.month)):
    #         a = y.sort_values()
    #         # print(y.apply(lambda y: percentileofscore(a, y)))
    #         percentiles_month.append(y.apply(lambda y: percentileofscore(a, y)))
    
    #     df_percentiles_month = pd.merge(df[["date", variable]], pd.concat(
    #         percentiles_month, axis=0), left_index=True, right_index=True)
    #     df["%s_percentiles_month" % variable] = df_percentiles_month["%s_y" % variable] 
    #     return (df["%s_percentiles_month" % variable]) # return the column of the percentile values, ordered according to the index of the provided dataframe 
    
    #%%
    # -- STEP 1 (of the combined TLM and CDPM method)- Transform values to percentiles 
    
    ## Define discharge variable 
    var = '\t"MEAN"' 
    # var_descrip = "mean discharge" 
    
    #%% -------------------------------------------------------------------------------------------------
    # Transform the time series in percentile (computed grouping per month) and then check which months has percentile values equal or greater than 5 
    GSIM_filt["%s_percentiles_month" % var] = Percentile_month(GSIM_filt, var)
    GSIM_filt["month"] = GSIM_filt.date.dt.month
    
    # Extract percentiles corresponding to zero values 
    Percentiles_0 = GSIM_filt.loc[GSIM_filt[var] == 0,
                                ("%s_percentiles_month" %var, "month")]
    # If values are greater than 5 percentile extract the month
    CDPM_month = np.unique(
        Percentiles_0.loc[Percentiles_0["%s_percentiles_month" %var] >= 5, "month"].values)
    
    #%%
    # -- STEP 2 - check if less than 5 percent of the timeseries contains a value of zero 
    
    # If CDPM_month is empty than apply the TLM method to the whole timeseries
    if len(CDPM_month) == 0:
        logger.info("Less than 5 percent of the time series contains a value of zero,\nhence the TLM method is applied to the whole timeseries")
        check = 1
        # Low flow events
        # Use the percentile of the MEAN streamflow grouped per month to assess if above/below threshold, according to the TLM analysis
        GSIM_filt['Low_flow'] = np.where(GSIM_filt["%s_percentiles_month" % var] <= thr_dry, 1, np.nan)
        # Compute severity of low flow
        GSIM_filt['Low_flow_int'] = np.where(GSIM_filt['Low_flow'] == 1, GSIM_filt["%s_percentiles_month" %var]-thr_dry, np.nan)
    
    else: 
        logger.info ('More than 5 percent of the time series contains zero values, \nhence both TLM and CDPM methods need to be applied')
        logger.info (Percentiles_0.drop_duplicates(keep='first'))
        check = 0
        
        # STEP 4. All positive values (x>0) are transformed in their corresponding percentile statistics (grouping per month)
        ## Extract positive data values from time series
        GSIM_tlm = GSIM_filt[(GSIM_filt[var] != 0)]
        ## Resert index
        GSIM_tlm.reset_index(drop=True, inplace=True)
        ## Compute percentile of the positive streamflow grouped per month
        GSIM_tlm ["perc_positive_runoff_month"] = Percentile_month(GSIM_tlm, var)
        ## Add the column of the percentile with positive flow in the main dataframe
        GSIM_merged = pd.merge(GSIM_filt[["date", "dry/wet_month", var,'\t"MAX"']], GSIM_tlm[[
                                   "date", "perc_positive_runoff_month"]], how='left', on="date") 
    
        
        # STEP 6a: Compute Fdry and Fwet -- grouped per month
        ## Identify the fraction of positive flow values per month throughout the whole timeseries
        Fwet = (GSIM_merged[GSIM_merged[var] > 0][var].groupby(GSIM_merged.date.dt.month).count(
        )/GSIM_merged[var].groupby(GSIM_merged.date.dt.month).count())
        Fwet.rename("Fwet_month", inplace=True)
        Fdry = 1 - Fwet
        Fdry.rename("Fzero_month", inplace=True)
        ## Create dataframe Fraction (containing Fdry and Fwet as columns)
        Fractions = Fwet.to_frame()
        Fractions["Fzero_month"] = Fdry
        ## Add the column to the main dataframe, to order the values according to the index
        GSIM_merged_f = pd.merge(GSIM_merged, Fractions, how='left',
                               left_on=GSIM_filt.date.dt.month, right_on=Fwet.index)
      
        # Write Fzero in relation to discharge values equal to zero and Fwet when values are above 0
        GSIM_merged["Fraction"] = np.where(GSIM_merged[var] > 0, GSIM_merged_f["Fwet_month"], GSIM_merged_f["Fzero_month"] )
        
        # STEP 3 and 5. Periods of positive runoff that experience a drought are combined with the zero runoff observations to obtain a new serie
        ## The following lines do the same computations of Marjolein's code but faster
        ## Write 1 when the value is 0
        GSIM_merged["dry_0"]= np.where(GSIM_merged[var] == 0, 1, np.nan)
        ## Compute the cumulative value of consecutive zero flow
        GSIM_merged["drymonths"] = GSIM_merged["dry_0"] * (GSIM_merged["dry_0"].groupby((GSIM_merged["dry_0"] != GSIM_merged["dry_0"].shift()).cumsum()).cumcount() + 1)    
        ##  Write 1 when the TLM analysis on the positive flow has identified a drouht otherwise write the values of dry_0 column
        GSIM_merged["dry_0_TLM"]= np.where((GSIM_merged[var] > 0 ) & (GSIM_merged["perc_positive_runoff_month"] <= thr_dry), 1, GSIM_merged["dry_0"])
        ## Compute the cumulative value of consecutive TLM drought and zero flow
        GSIM_merged["combined"] = GSIM_merged["dry_0_TLM"] * (GSIM_merged["dry_0_TLM"].groupby((GSIM_merged["dry_0_TLM"] != GSIM_merged["dry_0_TLM"].shift()).cumsum()).cumcount() + 1)    
        ## Delete column: dry_0 and dry_0_TLM
        GSIM_merged.drop(columns=["dry_0","dry_0_TLM"], inplace=True)
                            
        # Calculate percentile of 0 flows (no grouped per month) using both the dry months and combined time series (in the combined time series we also consider when a drought is happening according to the TLM and then is followed by 0 flow values)
        percentile_drymonths = Percentile_comp(GSIM_merged[(GSIM_merged["drymonths"] > 0)], "drymonths") # we won't use it in the following part of the code
        percentile_combined =  Percentile_comp(GSIM_merged[(GSIM_merged["combined"] > 0)], "combined")
        
               
        # STEP 6: Rescale percentile of zero values (combined percentile) and percentile of positive values (percentile_positive runoff) and Merge percentiles from months with flow to months without flow to get an overall flow percentile timeseries     
        ## Merge the combined percentile to the main dataframe 
        GSIM_merged_fin = pd.merge(GSIM_merged,percentile_combined [["date",'combined_percentiles']], how= 'left', left_on = GSIM_merged.date, right_on = percentile_combined.date )
        
        ####### CHECK---------------
        ## Rescale: if zero flow than multiply the percentile for the Fzero fruction otherwise multiply for the Fwet fruction and add the respective Fdry
        GSIM_merged_fin ["percentile_fin"] = np.where (GSIM_merged_fin[var] == 0, GSIM_merged_fin ['Fraction'] * (100 - GSIM_merged_fin ['combined_percentiles']), (1- GSIM_merged_fin ['Fraction'])*100 + GSIM_merged_fin ['Fraction'] * GSIM_merged_fin ["perc_positive_runoff_month"])
        # GSIM_merged_fin ["percentile_fin"] = np.where (GSIM_merged_fin["combined"] == 0, GSIM_merged_fin ['Fraction'] * (100 - GSIM_merged_fin ['combined_percentiles']), (1- GSIM_merged_fin ['Fraction'])*100 + GSIM_merged_fin ['Fraction'] * GSIM_merged_fin ["perc_positive_runoff_month"])
        ####### ---------------
        
        # Identify low and high flow and severity
        # Low flow events
        GSIM_merged_fin['Low_flow'] = np.where(GSIM_merged_fin["percentile_fin"] <= thr_dry, 1, np.nan)
        GSIM_merged_fin['Low_flow_int'] = np.where(GSIM_merged_fin['Low_flow'] == 1, GSIM_merged_fin["percentile_fin"]-thr_dry, np.nan)    
       
        
    # High flow events 
    # Percentile of the MAX streamflow over the whole timeseries and TLM analysis
    var_p ='\t"MAX"'
    perc_max = Percentile_comp(GSIM_filt, var_p)
    # df_perc_max = pd.merge(GSIM_filt[["date", var_p]], perc_max, on="date")
    GSIM_filt["%s_percentiles" %var_p] = perc_max["%s_percentiles" % var_p]
    
    GSIM_filt['High_flow'] = np.where( GSIM_filt["%s_percentiles" %var_p] >= (thr_wet), 1, np.nan)
    GSIM_filt['High_flow_int'] = np.where(GSIM_filt['High_flow'] == 1, (GSIM_filt["%s_percentiles" %var_p])-(thr_wet), np.nan)
    
    #%%
    # Report values in the dataframe with complete time series to indicate when there are NAN values (in order to not confond no-anomalies with NAN values)
    
    
    # Merge
    if check == 1: # PERENNIAL
        dt = pd.merge(GSIM['date'], GSIM_filt[["date",'\t"MEAN"', '\t"MAX"', '\t"MIN"', "dry/wet_month", "%s_percentiles_month" %var,"%s_percentiles" %var_p]], how= "left", on= "date")
        dt ["river_type"] = "perennial"
        ## Replace NaN values (representing no anomalies) in GISM_filt to -9999 to not confuse with the Nan values that we will have when there is no flow value
        GSIM_filt_na = GSIM_filt.fillna(-9999)
        dt_anomaly = pd.merge(dt[['date', 'dry/wet_month']], GSIM_filt_na[["date", 'Low_flow_int', 'High_flow_int']], how= "left", on= "date")
    else:
        GSIM_filt["%s_percentiles_month" %var] = GSIM_merged_fin ["percentile_fin"] # we could put them equal because they have the same index 
        dt = pd.merge(GSIM['date'], GSIM_filt[["date",'\t"MEAN"', '\t"MAX"', '\t"MIN"', "dry/wet_month", "%s_percentiles_month" %var,"%s_percentiles" %var_p]], how= "left", on= "date")
        dt ["river_type"] = "ephemeral"
        GSIM_filt['Low_flow_int'] = GSIM_merged_fin['Low_flow_int']
        GSIM_filt_na = GSIM_filt.fillna(-9999)
        dt_anomaly = pd.merge(dt[['date', 'dry/wet_month']], GSIM_filt_na[["date", 'Low_flow_int', 'High_flow_int']], how= "left", on= "date")
    
    
    return (dt, dt_anomaly)

#%%    
# dt, dt_anomaly = Streamflow(path, station_GSIM,Start_year, End_year, thr_dry, thr_wet, Percentile_month, Percentile_comp)
# # Export input values and respective percentiles  to excel   
# path_file_inputs = os.path.join(path, 'Output/Streamflow/Streamflow_%s.xlsx'%station_GSIM) 
# dt.to_excel(path_file_inputs)

# # Export dt_anomaly as excel
# path_file_outputs = os.path.join(path, 'Output/D&F_analysis/D&F_analysis_%s.xlsx'%station_GSIM) 
# dt_anomaly.to_excel(path_file_outputs)