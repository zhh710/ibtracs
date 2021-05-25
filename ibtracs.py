# -*- coding: utf-8 -*-
"""
5/12/2021

画台风路径
数据来自ibtracs
"""
from netCDF4 import Dataset
import numpy as np
import pandas as pd

filename = 'IBTrACS.last3years.v04r00.nc'
ibtracs = Dataset(filename,'r')
# dimensions
nstorm = ibtracs.dimensions['storm'].size
ndate_time = ibtracs.dimensions['date_time'].size
# variables, 
iso_time = ibtracs.variables['iso_time']
name = ibtracs.variables['name']
track_type = ibtracs.variables['track_type']
lat = ibtracs.variables['usa_lat']
lon = ibtracs.variables['usa_lon']
wind = ibtracs.variables['usa_wind']
pres = ibtracs.variables['usa_pres']
status = ibtracs.variables['usa_status']
basin = ibtracs.variables['basin']

## convert iso_time to pandas.Timestamp
def bytes_to_timestamp(bytes0):
    '''
    convert [b'2', b'0', b'1', b'8', b'-', b'0', b'1', b'-', b'0',
                   b'3', b' ', b'1', b'8', b':', b'0', b'0', b':', b'0',
                   b'0']
    to pandas.Timestamp
    '''
    return  pd.to_datetime(''.join(bytes0[:].astype('str')),format='%Y-%m-%d %H:%M:%S')
#print(bytes_to_timestamp(iso_time[0,0,:]))

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.transforms import offset_copy

def make_basemap_of_Atlantic(figure=None):
    '''
    '''
    if figure is not None:
        fig = figure
    else:
        fig = plt.figure(figsize=[16,16])
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_title('Atlantic')
    # extent = [-140, -20, 0, 90]
    extent = [-100, -40, 0, 40]
    ax.set_extent(extent)
    #ax.add_feature(dataset)
    ax.stock_img()
    ax.coastlines()
    ax.gridlines(x_inline=False, draw_labels=True)
    #plt.savefig("Atlantic.png")
    return ax,fig

def str_name(name0):
    '''
    return the length of name
    '''
    name0 = name0.astype('str')
    for count in range(len(name0)):
        try:
            tmp = ''.join(name0[count])
        except:
            break
    return ''.join(name0[0:count])
def get_storm_id_baseon_basin(basions=["NA","SA"]):
    '''
    EP=East_Pacific 
    NA=North_Atlantic 
    NI=North_Indian 
    SA=South_Atlantic 
    SI=South_Indian 
    SP=South_Pacific 
    WP=Western_Pacific
    '''
    idx = []
    for i in range(nstorm):
        basin0 = ''.join(basin[i,0,:].astype('str'))
        if basin0 in basions:
            name0 = str_name(name[i])
            idx.append(i)
        else:
            continue
    return len(idx),idx

# get storm in Atlantic basin
region_name = 'atlantic'
_,idx_atlantic = get_storm_id_baseon_basin(basions=["NA","SA"])

dict_color = {
'DB':'deepskyblue',
'TD':'aqua',
'TS':'bisque',
'TY':'gold',
'ST':'orange',
'TC':'coral',
'HU':'orangered',
'HR':'orangered',

'SD':'grey',
'SS':'grey',
'EX':'grey',
'PT':'grey',
'IN':'grey',
'DS':'grey',
'SD':'grey',
'LO':'grey',
'WV':'grey',
'ET':'grey',
'MD':'grey',
'XX':'grey',
}
def status_color(status0):
    '''
    return color:
     DB - disturbance,  deepskyblue
     TD - tropical depression,  aqua

     TS - tropical storm,  bisque
     TY - typhoon,  gold
     ST - super typhoon,  orange
     TC - tropical cyclone,  coral
     HU, HR - hurricane, orangered

     SD - subtropical depression,  
     SS - subtropical storm,  
     EX - extratropical systems,  
     PT - post tropical,  

     IN - inland,  
     DS - dissipating,  

     LO - low,  
     WV - tropical wave,  
     ET - extrapolated,  
     MD - monsoon depression,  
     XX - unknown.
    '''
    clr =  dict_color[''.join(status0[:].astype('str'))]
    if clr is not None:
        return clr
    else:
        return 'black'
#
from matplotlib.path import Path
import matplotlib.patches as patches
def add_storm_legend(axes):
    '''
    '''
    ax1 = axes.inset_axes([0.7,0.001,0.3,0.1],frameon=False)
    ax1.yaxis.set_tick_params(size=0)
    ax1.xaxis.set_tick_params(size=0)

    ax1.set_xticklabels('')
    ax1.set_yticklabels('')

    # patches
    bin0 = 1./7.
    verts_DB = [(0.,0.),(0.,0.4),(bin0,0.4),(bin0,0),(0.,0.)]
    verts_TD = [(bin0,0.),(bin0,0.4),(2*bin0,0.4),(2*bin0,0),(bin0,0.)]
    verts_TS = [(2*bin0,0.),(2*bin0,0.4),(3*bin0,0.4),(3*bin0,0),(2*bin0,0.)]
    verts_TY = [(3*bin0,0.),(3*bin0,0.4),(4*bin0,0.4),(4*bin0,0),(3*bin0,0.)]
    verts_ST = [(4*bin0,0.),(4*bin0,0.4),(5*bin0,0.4),(5*bin0,0),(4*bin0,0.)]
    verts_TC = [(5*bin0,0.),(5*bin0,0.4),(6*bin0,0.4),(6*bin0,0),(5*bin0,0.)]
    verts_HR = [(6*bin0,0.),(6*bin0,0.4),(7*bin0,0.4),(7*bin0,0),(6*bin0,0.)]
    codes = [
    Path.MOVETO,
    Path.LINETO,
    Path.LINETO,
    Path.LINETO,
    Path.CLOSEPOLY,]
    verts = [verts_DB,verts_TD,verts_TS,verts_TY,verts_ST,verts_TC,verts_HR]
    colors=[dict_color['DB'],dict_color['TD'],dict_color['TS'],dict_color['TY'],dict_color['ST'],dict_color['TC'],dict_color['HR']]
    labels=['DB',"TD",'TS','TY','ST','TC','HR']

    paths = [Path(vert0,codes) for vert0 in verts]
    i=1
    for l0,c0,p0 in zip(labels,colors,paths):
        pt = patches.PathPatch(p0,facecolor=c0,lw=0)
        ax1.add_patch(pt)
        ax1.text(i*bin0-bin0/2.,0.2,
            l0,
            horizontalalignment='center',
            verticalalignment='center',
            color='w',
            fontsize=12,
            fontweight='bold',
        )
        i+=1

def text_transform(axes):
    '''
    '''
    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(axes)
    return   offset_copy(geodetic_transform, units='dots', x=-10)
# plot
def plot_multistorm_onmap():
    '''

    '''
    ax , _ = make_basemap_of_Atlantic()
    for i in range(len(idx_atlantic)):

        idx0 = idx_atlantic[i]
        # get time length of storm
        time0 = []
        for it in range(ndate_time):
            try:
                time0.append(bytes_to_timestamp(iso_time[idx0,it,:]))
            except:
                break
        nt0 = it
        #
        lon0 = lon[idx0,0:nt0]
        lat0 = lat[idx0,0:nt0]
        status0 = status[idx0,:,:]
        name0 = str_name(name[idx0,:])
        #
        print(name0)
        for it in range(nt0):
            try:
                ax.plot(lon0[it:it+2],lat0[it:it+2],
                    linewidth=3,
                    transform=ccrs.Geodetic(),
                    color=status_color(status0[it,:]))
            except:
                break
    add_storm_legend(ax)    
    plt.savefig(region_name+'.png')

def plot_one_storm():
    '''
    '''
    fig = plt.figure(figsize=[16,16])
    for i in range(len(idx_atlantic)):
        idx0 = idx_atlantic[i]
        # get time length of storm
        time0 = []
        for it in range(ndate_time):
            try:
                time0.append(bytes_to_timestamp(iso_time[idx0,it,:]))
            except:
                break
        nt0 = it
        #
        lon0 = lon[idx0,0:nt0]
        lat0 = lat[idx0,0:nt0]
        status0 = status[idx0,:,:]
        name0 = str_name(name[idx0,:])
        time0_fs = time0[0].strftime(format='%Y%m%d%H%M%S')
        png_name0 = name0 + "_" + time0_fs
        t0_str = time0[0].strftime(format='%Y-%m-%d %H:%M:%S')
        t1_str = time0[-1].strftime(format='%Y-%m-%d %H:%M:%S')
        name0 = '{} {} - {}'.format(name0,t0_str,t1_str)
        #
        if name0 != 'ZETA':
            pass
        ax , _ = make_basemap_of_Atlantic(figure=fig)
        print("{} have {} points begining at {}".format(name0,nt0,time0[0]))
        # line
        ax.plot(lon0,lat0,
            linewidth=1.0,
            color='k',
            transform=ccrs.Geodetic())
        #text
        text_trans =  text_transform(ax)
        for it in range(nt0):
            #Point
            try:
                ax.plot(lon0[it],lat0[it],
                    marker="o",
                    markersize=4.5,
                    transform=ccrs.Geodetic(),
                    markeredgecolor=status_color(status0[it,:]),
                    markerfacecolor=status_color(status0[it,:]))
            except:
                break
            #text
            if time0[it].strftime(format='%H') == '00':
                ax.text(lon0[it],lat0[it],
                    time0[it].strftime(format='%d'),
                    verticalalignment='center', 
                    horizontalalignment='right',
                    transform=text_trans,
                    )
                ax.plot(lon0[it],lat0[it],
                    marker="o",
                    markersize=8,
                    transform=ccrs.Geodetic(),
                    markeredgecolor=status_color(status0[it,:]),
                    markerfacecolor=status_color(status0[it,:]))

        print(it+1,"ploted")
        ax.set_title(name0,fontsize=20,fontweight='bold')
        add_storm_legend(ax)
        #
        plt.savefig(png_name0+'.png')
        plt.cla()
    plt.close(fig)

#plot_multistorm_onmap()
plot_one_storm()

# Close file
ibtracs.close()
