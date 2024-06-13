import pandas as pd
import xarray as xr
import numpy as np
import geopandas as gp
from shapely.geometry import Point

def to_gdf(df,lat,lon,epsg):
    geoms = [Point(*xy) for xy in df[[lon, lat]].values]
    gdf = gp.GeoDataFrame(data=df, geometry=geoms, crs=epsg) 
    return gdf

sitestore = '/fs/site6/eccc/crd/ccrp/users/com002/'
f = '/fs/site5/eccc/crd/ccrp/data/CanSWE/CanSWE_2024/Regional_Data/netcdf_merged/ONR_OPG-AB-BC.nc'

# get 2024 data (not yet in CanSWE)
ds = xr.open_dataset(f)
df1 = ds.copy().to_dataframe().reset_index()
df1a = df1.loc[(df1.snw>=0)(df1.time.dt.month==4)&(df1.time.dt.day==1)] #|(df1.snd>=0)
df2024 = df1a[['station_id','time','lat','lon','source','station_name','snw','snd']]

# canswe for 30yr climatology
f2 = '/fs/site6/eccc/crd/ccrp/data/CanSWE/CanSWE_2023/Updated_Data_Pub/CanSWE-CanEEN_1928-2023_v6.nc'
can = xr.open_dataset(f2)
df2 = can.copy().to_dataframe().reset_index()
df2 = df2.loc[(df2.snw>=0)|(df2.snd>=0)]
df3 = df2.loc[(df2.type_mes == 2)&(df2.time.dt.year>=1991)&(df2.time.dt.month==4)&(df2.time.dt.day==1)]
apr = df3[['station_id','time','lat','lon','elevation','source','station_name','snw','snd']]

# mean and number of years with obs 1991-2020
apr9120 = apr.loc[(apr.time.dt.year>=1991)&(apr.time.dt.year<=2020)].groupby(['station_id','lat','lon','station_name','elevation'])['snw'].mean().reset_index().rename(columns = {'snw':'clm'})
apr9120c = apr.loc[(apr.time.dt.year>=1991)&(apr.time.dt.year<=2020)].groupby(['station_id','lat','lon','station_name','elevation'])['snw'].count().reset_index().rename(columns = {'snw':'nyrs'})
apr_info = pd.merge(apr9120[['station_id','elevation','clm']],apr9120c[['station_id','nyrs']],on = 'station_id')

# grab 1 April 2024 for the sites with obs in all 30yrs
apr1 = pd.merge(df2024.loc[(df2024.time.dt.year==2024)],apr_info,on='station_id')

# compute anom as 2024 minus 30yr climatology
apr1['anom'] = apr1.snw - apr1.clm 
apr1['anom_pct'] = np.int16(np.round((apr1.anom/apr1.clm)*100,0))

# check zeroes
# identify instances in the 30yr climatology where swe == 0 on 1 april.
df0 = pd.DataFrame()
for sta in apr.station_id.unique():
    df_sub = apr.loc[(apr.station_id==sta)&(apr.time.dt.year>=1991)&(apr.time.dt.year<=2020)]
    sub0 = df_sub.loc[df_sub.snw==0]
    if len(sub0)>0:
        print(f'{sta}: {sub0.time.unique()}')
        df0 = pd.concat([df0,sub0])

# write the years with 0swe to the corresponding site in the dataframe
def get_unique(df):
    return df.time.dt.year.unique()

stuff = df0.groupby(['station_id']).apply(get_unique).reset_index().rename(columns = {0:'snw0'})
df_apr = pd.merge(apr1,stuff,on = 'station_id',how = 'left')

# save
df_apr.drop(['snd'],axis=1).to_csv(f'{sitestore}jb_smsc/1Apil2024_BCAB_anom.csv',index=False)

# map it!
from matplotlib import pyplot as plt 
import matplotlib as mpl

# to geodataframe for plotting
apr1_gdf = to_gdf(df_apr,'lat','lon','EPSG:4326')

world = gp.read_file(gp.datasets.get_path('naturalearth_lowres'))
nam = world.loc[world.iso_a3.isin(['USA','CAN'])].copy()

# grab snotel anoms
snt = to_gdf(pd.read_csv(f'{sitestore}jb_smsc/1Apil2024_SNOTEL_anom.csv'),'lat','lon','EPSG:4326')

n = 25 # minimum number of years with valid SWE obs
# swe anom limits for plotting
smin = -100 
smax = 100

fig, ax = plt.subplots(nrows = 1, ncols = 1)
ax1 = ax
nam.plot(ax=ax,linewidth=0.5, edgecolor='grey', color='whitesmoke', alpha = 0.8)
apr1_gdf.loc[apr1_gdf.nyrs>=n].plot(ax = ax1,c = 'grey',markersize = 2)
snt.loc[snt.nyrs>=n].plot(ax = ax1,c = 'grey',markersize = 2)
apr1_gdf.loc[apr1_gdf.nyrs>=n].plot(ax = ax1,column = 'anom_pct',cmap = 'bwr_r',markersize = 0.75,vmin = smin,vmax = smax)
snt.loc[snt.nyrs>=n].plot(ax = ax1,column = 'anom_pct',cmap = 'bwr_r',markersize = 0.75,vmin = smin,vmax = smax)
ax1.set_xlim([-176,-100])
ax1.set_ylim([27,74])
ax1.set_title(f'1 April 2024 SWE anomaly relative to 1991-2020\nSites with min {n}yrs during 1991-2020 shown')

cbar_ax = fig.add_axes([0.29, 0.27, 0.25, 0.04]) #0.53
norm = mpl.colors.Normalize(vmin=smin, vmax=smax)
sm = plt.cm.ScalarMappable(cmap='bwr_r', norm=norm)
sm.set_array([])
plt.colorbar(sm, cax = cbar_ax,ticks=np.linspace(smin, smax, 5).round(0),
    orientation = 'horizontal',extend = 'both',label = 'SWE anom. (%)')

fig.savefig(f'{sitestore}jb_smsc/1Apr2024anom_NA_{n}yrs.png')
plt.close()
