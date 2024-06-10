# jb_smsc

### SNOTEL 1 April anomalies

Compute 1 April 2024 SWE anomalies relative to 1991 - 2020 from SNOTEL data.

This script calculates the 1 April SWE climatology for 1991-2020 for SNOTEL sites. It then calculates the 1 April SWE anomaly for 2024 relative to 1991-2020 climatology.

For each sites, instances of SWE == 0 are identified and written to the main dataframe as are the number of years during 1991-2020 with SWE observations.

It includes some code for an ugly map. The map only displays sites with observations in all 30yrs (and 2024).

The SNOTEL data used are directly from NRCS (i.e. not the BCQC dataset). These data will soon be included in v2 of https://zenodo.org/records/10287093 , which should make it easier to reference.
