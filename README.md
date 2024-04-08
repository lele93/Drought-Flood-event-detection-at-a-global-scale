# Compound and Consecutive Drought-Flood event detection at a global scale 

Welcome to our repository! The scripts provided here are designed to detect drought and flood events, as well as their consecutive and compound occurrences, through the analysis of various hydro-meteorological and biophysical variables. These scripts are applied to a dataset comprising 8255 catchments.
Data Sources:

The codes utilize time-series data from various sources:

    Precipitation (P): Data sourced from MSWEP (Beck et al., 2019).
    Streamflow (STR) : Data sourced from GSIM (Do et al. 2018; Gudmundsson et al. 2018). 
    Soil Moisture (SM): Soil moisture readings near the surface (0 – 5 cm depth) and in the root zone (0 – 250 cm depth) are obtained from GLEAM (Martens et al., 2017).
    Total Water Storage (TWS): Data sourced from GRACE (Boergens et al., 2019).
    Normalized Difference Vegetation Index (NDVI): Data sourced from NOAA (Vermote et al., 2014).
    Surface Water Extent: Time-series data retrieved from Donchyts et al. (2016).

Description:

    Anomalies Detection:
        Within the "Anomalies_Detection" folder, high and low anomalies have been computed for different hydrometeorological time series using both variable and fixed threshold level methods (refer to Matanó, A.; Berghuijs, W.; Mazzoleni, M.; de Ruiter, M.C.; Ward, P.J., Van Loon, A.F.: Compound and consecutive drought-flood events at a global scale, for more detail on the methodology)



    Detect_D_F_events.py Script:
        This script, Detect_D_F_events.py, is the core of our event detection process. It identifies single drought and flood events, as well as their cascading and compounding occurrences. Notably, it incorporates drought propagation into the detection of drought periods for more accurate results.

Usage:

To utilize these scripts effectively, follow these steps:

    1. Ensure you have the necessary data sources available or replace them with your desired datasets.
    2. Navigate to the "Anomalies_Detection" folder.
    3. Run the script High_low_extremes.py. This script automatically calls other necessary functions contained within the folder's other scripts.
    4. Once the anomalies detection process is complete, proceed to run Detect_D_F_events.py to detect drought and flood events.

Feel free to explore and modify these scripts to adapt them to your research or analysis requirements.

