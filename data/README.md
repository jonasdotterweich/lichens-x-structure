Data Availability Statement

The raw data files for this project are not hosted in this GitHub repository, as it is sensitive biological field work data whis is not made public at this stage.

To run the analysis scripts in the /code folder, you must manually place the required datasets into this directory following the structure below.

Structural requirements for the /data folder:

csv file containing the raw data species data
This file should optimally contain teh following columns:
- "Plot_ID": Unique identifier for each sample plot 
- "Species": The species observed in the plot
- "Count": The count of individuals of the species in the plot
- "Location": The geographical location of the plot
- "Date": The date when the data was collected 


csv file containing the raw data structural data
This file should optimally contain teh following columns:
- "Plot_ID": Unique identifier for each sample plot 
- "Location": The geographical location of the plot
- "Date": The date when the data was collected #
- various structural variables recorded, each in a single column, such as e.g.:
  - "Canopy_Height": The height of the canopy in meters
  - "Basal_Area": The basal area of trees in square meters per hectare
  - "Leaf_Area_Index": The leaf area index for the plot

