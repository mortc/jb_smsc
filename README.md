# jb_smsc

### SNOTEL 1 April anomalies

Compute 1 April 2024 SWE anomalies relative to 1991 - 2020 from SNOTEL data.

This script calculates the 1 April SWE climatology for 1991-2020 for SNOTEL sites. It then calculates the 1 April SWE anomaly for 2024 relative to 1991-2020 climatology.

For each sites, instances of SWE == 0 are identified and written to the main dataframe as are the number of years during 1991-2020 with SWE observations.

It includes some code for an ugly map. The map only displays sites with observations in all 30yrs (and 2024).

The SNOTEL data used are directly from NRCS (i.e. not the BCQC dataset). These data will soon be included in v2 of https://zenodo.org/records/10287093 , which should make it easier to reference.

### .csv output [1April2024_SNOTEL_anom.csv]
time: 1 April 2024
station_id: SNOTEL sation ID appended to 'SNOTEL_'
station_name: SNOTEL station name
snw: 1 April 2024 snow water equivalent, units = mm or kg/m-2
snd: 1 April 2024 snow depth, units = m
lat: station latitude, degrees North, from SNOTEL metadata
lon: station longitude, degrees West, from SNOTEL metadata
Elevation_NTK: station elevation (m) taken from the SNOTEL metadata
clm: 30yr (1991-2020) climatological SWE. This value is calculated for all available sites regardless of the number of years with data.
nyrs: number of years with valid SWE observations
anom: 1 April 2024 SWE anomaly relative to the 30yr climatology. This value is calculated for all available sites regardless of the number of years with data.
snw0: years when swe = 0 on 1 April for the corresponding site. If NaN then all years with observations (all nyrs) have swe >0.
