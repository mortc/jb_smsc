#jbsmsc 1 April snotel

import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gp
from shapely.geometry import Point

def to_gdf(df,lat,lon,epsg):
    geoms = [Point(*xy) for xy in df[[lon, lat]].values]
    gdf = gp.GeoDataFrame(data=df, geometry=geoms, crs=epsg) 
    return gdf

# data and dirs
# currently using an internal version of the snotel data. Soon these data will be available as v2 of the dataset https://zenodo.org/records/10287093 .
sitestore = '/fs/site6/eccc/crd/ccrp/users/com002/'
f_snt = '/home/vvi001/store6/Eval_Data/SNOTEL_data/data_all_US_1979_2024/netcdf/snotel_US_1979_2024_QC.nc'

# get the snotel data and do some basic filtering on day and month
snt = xr.open_dataset(f_snt)
df = snt[['time','station_id','station_name','snw','snd','lat','lon','station_id','Elevation_NTK']].to_dataframe().reset_index()
df_snt = df.loc[(df.snd>=0)|(df.snw>=0)].copy()
apr = df_snt.loc[(df_snt.time.dt.month == 4)&(df_snt.time.dt.day==1)]

# mean and number of years with obs 1991-2020
apr9120 = apr.loc[(apr.time.dt.year>=1991)&(apr.time.dt.year<=2020)].groupby(['station_id','lat','lon','station_name','Elevation_NTK'])['snw'].mean().reset_index().rename(columns = {'snw':'clm'})
apr9120c = apr.loc[(apr.time.dt.year>=1991)&(apr.time.dt.year<=2020)].groupby(['station_id','lat','lon','station_name','Elevation_NTK'])['snw'].count().reset_index().rename(columns = {'snw':'nyrs'})
apr_info = pd.merge(apr9120[['station_id','clm']],apr9120c[['station_id','nyrs']],on = 'station_id')

# grab 1 April 2024 for the sites with obs in all 30yrs
apr1 = pd.merge(apr.loc[(apr.time.dt.year==2024)],apr_info,on='station_id')

# compute anom as 2024 minus 30yr climatology
apr1['anom'] = apr1.snw - apr1.clm

# check zeroes
# identify instances in the 30yr climatology where swe == 0 on 1 april.
df0 = pd.DataFrame()
for sta in apr.station_id.unique():
    df_sub = apr.loc[(apr.station_id==sta)&(apr.time.dt.year>=1991)&(apr.time.dt.year<=2020)]
    sub0 = df_sub.loc[df_sub.snw==0]
    if len(sub0)>0:
        #print(f'{sta}: {sub0.time.unique()}')
        df0 = pd.concat([df0,sub0])

# write the years with 0swe to the corresponding site in the dataframe
def get_unique(df):
    return df.time.dt.year.unique()

df_zeroes = df0.groupby(['station_id']).apply(get_unique).reset_index().rename(columns = {0:'snw0'})
df_apr = pd.merge(apr1,df_zeroes,on = 'station_id',how = 'left')

# save
df_apr.to_csv(f'{sitestore}jb_smsc/1Apil2024_SNOTEL_anom.csv',index=False)

# map anoms for sites with obs in all yrs (ugly map)
from matplotlib import pyplot as plt 
import matplotlib as mpl

# to geodataframe for plotting
apr1_gdf = to_gdf(df_apr.loc[df_apr.nyrs == 30],'lat','lon','EPSG:4326')

usa = gp.read_file(f'{sitestore}shapefiles/US/cb_2018_us_state_20m/cb_2018_us_state_20m.shp').to_crs('EPSG:4326')
usa = usa.loc[~usa.STUSPS.isin(['HI','PR'])]
conus = usa.loc[~usa.STUSPS.isin(['AK','HI','PR'])]

fig, ax = plt.subplots(nrows = 1, ncols = 2)
ax1,ax2 = ax
usa.loc[usa.STUSPS=='AK'].plot(ax=ax1,linewidth=0.5, edgecolor='grey', color='whitesmoke', alpha = 0.8)
usa.loc[usa.STUSPS!='AK'].plot(ax=ax2,linewidth=0.5, edgecolor='grey', color='whitesmoke', alpha = 0.8)
gp.sjoin(apr1_gdf,usa.loc[usa.STUSPS=='AK'],how = 'inner').plot(ax = ax1,c = 'grey',markersize = 2)
gp.sjoin(apr1_gdf,usa.loc[usa.STUSPS!='AK'],how = 'inner').plot(ax = ax2,c = 'grey',markersize = 2)
gp.sjoin(apr1_gdf,usa.loc[usa.STUSPS=='AK'],how = 'inner').plot(ax = ax1,column = 'anom',cmap = 'bwr',markersize = 0.75,vmin = -100,vmax = 100)
gp.sjoin(apr1_gdf,usa.loc[usa.STUSPS!='AK'],how = 'inner').plot(ax = ax2,column = 'anom',cmap = 'bwr',markersize = 0.75,vmin = -100,vmax = 100)
ax1.set_xlim([-180,-130]),ax2.set_xlim([-125,-100])
cbar_ax = fig.add_axes([0.14, 0.12, 0.3, 0.04])
norm = mpl.colors.Normalize(vmin=-100, vmax=100)
sm = plt.cm.ScalarMappable(cmap='bwr', norm=norm)
sm.set_array([])
plt.colorbar(sm, cax = cbar_ax,ticks=np.linspace(-100, 100, 5).round(0),
    orientation = 'horizontal',extend = 'both',label = 'SWE anom. (mm)')#, label = 'FMA obs'
title_ax = fig.add_axes([0.12, 0.75, 0.5, 0.5])
title_ax.text(0,0, 'SNOTEL 1 April 2024 SWE\nanomaly relative to\n1991 - 2020')
title_ax.axis('off')
fig.savefig(f'{sitestore}jb_smsc/test.png')
plt.close()
