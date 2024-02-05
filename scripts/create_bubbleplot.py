# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:25:33 2016

@author:    Sam Lamont
            NOAA National Water Center
            Tuscaloosa, AL 35404
            samuel.lamont@noaa.gov
            (205) 347-1454
"""
 
#import matplotlib
#matplotlib.use('WX') # 'GTKAgg', 'TkAgg', 'WX' ???

import netCDF4
import timeit
import sys
import  pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
import psycopg2 as psy2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import gc
import folium
import glob
import os
from folium import IFrame
from folium import plugins
#import gzip
#import shutil
from bs4 import BeautifulSoup
import math
import base64
import urllib2

#import rasterio
from rasterio.warp import calculate_default_transform, reproject #, RESAMPLING, transform_bounds
from rasterio import crs, transform

#import subprocess


pd.set_option('use_inf_as_null', True)

plt.style.use('ggplot')

def connect_ini():
    return psy2.connect(dbname='ini', user='sam', password='kinder') 
    
def connect_oneone():
    return psy2.connect(dbname='oneone', user='sam', password='kinder')    
    
'''

HRAP HrapToLatLong(HRAP hrap)
{
  double raddeg = 57.29577951;
  double earthrad = 6371.2;
  double stdlon = 105.;
  double mesh_len = 4.7625;
  double tlat, rr, gi, ang, x, y;
  HRAP ll;

  tlat=60./raddeg;

  x = hrap.x - 401.;
  y = hrap.y - 1601.;
  rr = x*x + y*y;
  gi = ((earthrad * (1+sin(tlat)))/mesh_len);
  gi=gi*gi;
  ll.y = asin((gi-rr)/(gi+rr))*raddeg;
  ang = atan2(y,x)*raddeg;
  if (ang < 0) ang=ang+360.;
  ll.x = 270+stdlon-ang;
  if (ll.x < 0) ll.x=ll.x+360;
  if (ll.x > 360) ll.x=ll.x-360;
  return ll;
}
'''

'''
        "#04e9e7",  # 0.01 - 0.10 inches
        "#019ff4",  # 0.10 - 0.25 inches
        "#0300f4",  # 0.25 - 0.50 inches
        "#02fd02",  # 0.50 - 0.75 inches
        "#01c501",  # 0.75 - 1.00 inches
        "#008e00",  # 1.00 - 1.50 inches
        "#fdf802",  # 1.50 - 2.00 inches
        "#e5bc00",  # 2.00 - 2.50 inches
        "#fd9500",  # 2.50 - 3.00 inches
        "#fd0000",  # 3.00 - 4.00 inches
        "#d40000",  # 4.00 - 5.00 inches
        "#bc0000",  # 5.00 - 6.00 inches
        "#f800fd",  # 6.00 - 8.00 inches
        "#9854c6",  # 8.00 - 10.00 inches
        "#a122b7"   # 10.00+
'''

def get_index(max_val):
    if max_val > 0.01 and max_val <= 0.1:
        ind=0
    elif max_val > 0.1 and max_val <= 0.25:
        ind=1
    elif max_val > 0.25 and max_val <= 0.5:
        ind=2        
    elif max_val > 0.5 and max_val <= 0.75:
        ind=3
    elif max_val > 0.75 and max_val <= 1.0:
        ind=4
    elif max_val > 1.0 and max_val <= 1.5:
        ind=5
    elif max_val > 1.5 and max_val <= 2.0:
        ind=6
    elif max_val > 2.0 and max_val <= 2.5:
        ind=7
    elif max_val > 2.5 and max_val <= 3.0:
        ind=8
    elif max_val > 3.0 and max_val <= 4.0:
        ind=9
    elif max_val > 4.0 and max_val <= 5.0:
        ind=10     
    elif max_val > 5.0 and max_val <= 6.0:
        ind=11
    elif max_val > 6.0 and max_val <= 8.0:
        ind=12
    elif max_val > 8.0 and max_val <= 10.0:
        ind=13
    elif max_val > 10.0:
        ind=14        
    return ind
        
# #############################################################################
##  Format a folium bubble plot by adding info in a <div>...
# #############################################################################        
def format_folium_map(str_folium_map, start_time, end_time, str_where, str_prodid):

    # Get the soup...
    with open(str_folium_map, 'r') as i_map:        
        soup = BeautifulSoup(i_map, 'lxml')

    new_title = soup.new_tag('title')
    new_title.string = '{} Summary Plots'.format(str_prodid)
    new_div = soup.new_tag('div')
    new_div['style'] = 'color:#5E6E7B; font-weight:bold; text-align:center; background:#CDD2D4' # font-style:italic; 
    new_h = soup.new_tag('h')

    new_div.insert(1, new_h)
    
    new_h.string = '{} | {} -- {} | {}'.format(str_prodid, start_time[:-9], end_time[:-9], str_where)
        
    soup.body.insert(1, new_div)
    soup.head.insert(1, new_title)
        
    # Write the file again to save edits...
    with open(str_folium_map, 'w') as out_map:
        out_map.write(str(soup))
        
# #############################################################################
##  Create a folium map with pre-existing pop-up files...
# #############################################################################                
def create_folium_map_png(str_popup_dir, start_time, end_time, str_where, str_prodid, str_preciprastername, yr, mon, day):
    
    print('\n >> Building map...')

    engine = create_engine('postgresql://', creator=connect_ini, poolclass=NullPool)

    # Popup setup...
    resolution=80 # originally set at 75 (try lower for size??)
    width=6.5
    height=4.5    
    
#    str_map_title='GagesII-Ref above minor flood stage'
    str_map_title='PBias Error (Red=Under simulated, Blue=Over simulated)'
    
#    str_stat = 'tot_vol_err'
    
    # Get all desired lid's from the location_nwm table...
    str_sql = 'select distinct gage_id, lid, usgs_lat, usgs_lon from nwm_location {};'.format(str_where)
    df_loc = pd.read_sql_query(str_sql, engine)
    
    # Watch out for missing lat/lon values...
    df_loc.usgs_lat = df_loc.usgs_lat.str.strip()
    df_loc.usgs_lon = df_loc.usgs_lon.str.strip()
    
    df_loc.usgs_lat.replace('', np.nan, inplace=True)
    df_loc.usgs_lon.replace('', np.nan, inplace=True)
    
    df_loc.dropna(subset=['usgs_lat','usgs_lon'], inplace=True)
    
    df_loc['usgs_lat'] = df_loc['usgs_lat'].astype(np.float)
    df_loc['usgs_lon'] = df_loc['usgs_lon'].astype(np.float)
    
    if df_loc.empty:
        print('problem with location table')
        
    # Get mean lat/long for map setup??
    mean_lat = df_loc['usgs_lat'].mean()
    mean_long = df_loc['usgs_lon'].mean()
    
    # Folium map...
    map_f = folium.Map(location=[mean_lat, mean_long], zoom_start=5, tiles='Cartodb Positron', control_scale=True)
    
    # Add NHD layer...
    folium.WmsTileLayer(url='https://basemap.nationalmap.gov/arcgis/services/USGSHydroCached/MapServer/WmsServer',
                            layers='0',
                            name='USGS-NHD',
                            format='image/png').add_to(map_f)    
    
    # Generate IN string...
    str_lid_in=''
    for tpl in df_loc.itertuples():
        str_lid_in += '\'{}\', '.format(tpl.lid.strip())        
    # Remove last two characters...
    str_lid_in = str_lid_in[:-2]
    
    # << Stats>> for specified time and group by Mean...
    str_sql = 'select lid, pbias, tot_vol_err from nwm_fcst_stats_sam where product_id = \'{}\' and basistime >= \'{}\' and basistime <= \'{}\' and lid in ({});'.format(str_prodid, start_time, end_time, str_lid_in)
    df_stats = pd.read_sql_query(str_sql, engine)  
    df_stats['lid'] = df_stats['lid'].str.strip()
    df_stats.set_index('lid', inplace=True)   
    
    df_stats.dropna(inplace=True)
#    df_stats = df_stats[df_stats['pbias'].notnull()]
    
#    gpd_abs_stats = df_stats.abs().groupby(df_stats.index) # WHY AM I DOING THIS?? Mean absolute error
    gpd_stats = df_stats.groupby(df_stats.index)
    
    df_mean_abs_stats = gpd_stats.agg(np.mean) 

    # Join df_loc, df_mean_stats based on lid to include lat/long and other attributes...
    df_loc.set_index('lid', inplace=True)
    df_final_abs_mean = pd.merge(df_loc, df_mean_abs_stats, left_index=True, right_index=True)
    
    ## Classifying based on quantiles...
    df_final_abs_mean['stat_cat']= df_final_abs_mean['pbias']
    df_final_abs_mean['stat_cat']=pd.qcut(np.abs(df_final_abs_mean['stat_cat']), 6, labels=[3,6,10,15,20,30])
    
            
    # <<< Get list of filenames and paths in selected directory >>>
    lst_filenames = glob.glob(str_popup_dir + '*.png') 
    cntr = 0
    fg_markers = folium.FeatureGroup(name='Stations') # Add markers to this layer
    for this_path in lst_filenames:

#        ###################
#        if cntr == 10:            
#            break
#        ###################
                    
        # Get the lid from the plot filename...
        str_file = os.path.basename(this_path)        
        n_lid = os.path.splitext(str_file)[0]
        
        print('\t{} -- lid: {}'.format(cntr, n_lid))
        
#        if n_lid=='MACF1':
#            print('pause')
    
        try:
            gage_id = df_final_abs_mean.loc[n_lid]['gage_id']
            gage_id = str(gage_id).strip()
            lat = df_final_abs_mean.loc[n_lid]['usgs_lat']
            lon = df_final_abs_mean.loc[n_lid]['usgs_lon']
            rad = df_final_abs_mean.loc[n_lid]['stat_cat']
        except:
            print('Selecting from df_final_abs_mean on index: {} failed'.format(n_lid))
            continue
                             
        if math.isnan(rad):
            print('nan!')
            rad=15
            # Skip it for now...?
#            continue                             
                             
        encoded = base64.b64encode(open(this_path, 'rb').read())        
        html = '<img src="data:image/png;base64,{}">'.format
        iframe = IFrame(html(encoded), width=(width*resolution)+75, height=(height*resolution)+50)
        popup = folium.Popup(iframe, max_width=2650)
        
        if df_final_abs_mean.loc[n_lid]['tot_vol_err'] < 0:
#        if srs_diff.loc[n_lid] < 0: 
            fillcolor = 'red' # v1.1 has less error than v1.0
        else:
            fillcolor = 'blue'
               
        fg_markers.add_child(folium.RegularPolygonMarker(location=[lat, lon], popup=popup, number_of_sides=50, fill_opacity=0.6, fill_color=fillcolor, color='#7e7e7e', radius=rad))        
        
        cntr+=1

    # NOTE: Could use topojson to reduce file size in future; need to figure out how to install
    print('Adding layers...')  
    print('RFC boundary layer')                  
    fg_json = folium.FeatureGroup(name='RFC Boundaries') 
    fg_json.add_child(folium.GeoJson(data=open('/home/share/nwm/stats_plots/programs/bubbleplots/rfc_boundaries.geojson'),
                                     name='RFC Boundaries',
                                     style_function=lambda feature: {
                                    'fillColor': '#ffff00',
                                    'color': 'gray',
                                    'weight': 2,
                                    'dashArray': '5, 5',
                                    'fillOpacity': 0.05
                                }))              
    
    fg_precip = folium.FeatureGroup(name='Daily Precip Sum')
    
#    # ==================HRAP======================   
    print('HRAP - downloading, reprojecting, and adding to map...')
    str_url = 'https://water.weather.gov/precip/p_download_new/{}/{:02}/{:02}/nws_precip_conus_{}{:02}{:02}.nc'.format(yr,mon,day,yr,mon,day)
    filename = '/home/share/sams_scripts/bubbleplots/temp_hrap.nc'
    
    data = urllib2.urlopen(str_url)
    
    with open(filename, 'wb') as f:
        f.write(data.read())
        
    # ==== HRAP to Sphere ====   
    ds_test = netCDF4.Dataset(filename)            
    in_arr = ds_test.variables['amountofprecip'][:]   
    
    in_arr = np.flipud(in_arr) # NOTE: THIS IS KEY, rotate on x-axis before reprojecting!
    
    in_arr = in_arr/2540. # convert to inches
    
    max_val = in_arr.max()
 
    # HRAP original...
    src_transform = transform.from_bounds(-387, -1591, 664, -778, 1051, 813)
    src_crs = crs.CRS.from_string('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-105 +k=1 +x_0=0 +y_0=0 +a=6371200 +b=6371200 +to_meter=4762.5 +no_defs')
    
    # HRAP Sphere...
    dst_crs = crs.CRS.from_string('+proj=longlat +a=6371200 +b=6371200 +no_defs')   
    dst_transform, dst_width, dst_height = calculate_default_transform(src_crs, dst_crs, 1051, 813, -387, -1591, 664, -778)
    
    # Reproject...
    out_arr_sphere = np.zeros((dst_height, dst_width), dtype=np.float64)
    reproject(source=in_arr, destination=out_arr_sphere, src_crs=src_crs, dst_crs=dst_crs,
              src_transform=src_transform, dst_transform=dst_transform) #, resampling=RESAMPLING.bilinear)     
              
    # ==== WGS84 to Web Mercator ====
    # WGS84...
    src_transform = dst_transform
    src_crs = {'init': 'EPSG:4326'}  

    # Get src bounds...
    src_bnds = transform.array_bounds(out_arr_sphere.shape[0], out_arr_sphere.shape[1], src_transform)
    dst_crs = {'init': 'EPSG:3857'}
    dst_transform, dst_width, dst_height = calculate_default_transform(src_crs, dst_crs, out_arr_sphere.shape[1], out_arr_sphere.shape[0], src_bnds[0], src_bnds[1], src_bnds[2], src_bnds[3])
    
    # Reproject...
    out_arr_merc = np.zeros((dst_height, dst_width), dtype=np.float64)
    reproject(source=out_arr_sphere, destination=out_arr_merc, src_crs=src_crs, dst_crs=dst_crs,
              src_transform=src_transform, dst_transform=dst_transform) #, resampling=RESAMPLING.bilinear)   
              
#    max_val = out_arr_merc.max()
#    out_arr_merc = out_arr_merc/2540.
    out_arr_merc[out_arr_merc<0.1] = np.nan  

    nws_precip_colors = [
        "#04e9e7",  # 0.01 - 0.10 inches
        "#019ff4",  # 0.10 - 0.25 inches
        "#0300f4",  # 0.25 - 0.50 inches
        "#02fd02",  # 0.50 - 0.75 inches
        "#01c501",  # 0.75 - 1.00 inches
        "#008e00",  # 1.00 - 1.50 inches
        "#fdf802",  # 1.50 - 2.00 inches
        "#e5bc00",  # 2.00 - 2.50 inches
        "#fd9500",  # 2.50 - 3.00 inches
        "#fd0000",  # 3.00 - 4.00 inches
        "#d40000",  # 4.00 - 5.00 inches
        "#bc0000",  # 5.00 - 6.00 inches
        "#f800fd",  # 6.00 - 8.00 inches
        "#9854c6",  # 8.00 - 10.00 inches
        "#a122b7"   # 10.00+
    ]
    
    ind_val = get_index(max_val)
    
#    lst_indexmap=[0.1, 0.25, 0.5, 0.75, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    
    colormap = mpl.colors.ListedColormap(nws_precip_colors[:ind_val]) # slicing reverses the list  
    colormap.set_under('k', alpha=0)    
    
    fg_precip.add_child(plugins.ImageOverlay(
        image=out_arr_merc,
        bounds=[[20.76866172, -131.4470], [55.3702, -64.5174350]],
        origin='upper',
        colormap=colormap,
#        colormap=plt.cm.jet, # or cm.s3pcpn
        mercator_project=False,
        opacity=0.7,
    )) 
    
    # =====================================================

    map_f.add_child(fg_precip)
    map_f.add_child(fg_json)               
    map_f.add_child(fg_markers)
    
    map_f.add_child(folium.LayerControl())

    str_html_filename = '{}_{}_GT10cfs.html'.format(str_prodid, start_time[0:10])
    
    print('Saving map...')
    map_f.save(str_popup_dir + str_html_filename)    
    format_folium_map(str_popup_dir + str_html_filename, start_time, end_time, str_map_title, str_prodid)
    
#    # Zip...
#    with open(str_popup_dir + str_html_filename, 'rb') as f_in, gzip.open(str_popup_dir + str_html_filename + '.gz', 'wb') as f_out:
#        shutil.copyfileobj(f_in, f_out)
    
    # Delete plots...
    print('Deleting pop-ups and temporary files...')
    for this_path in lst_filenames:
        os.remove(this_path)
        
#    os.remove('/home/share/nwm/stats_plots/programs/bubbleplots/temp_hrap.asc')
    os.remove(filename) # remove the ahps netcdf file
        
    return
    
# #############################################################################
##  QUERY from forecast vs. observed plotting forecast progression
##   and saving as .png
# #############################################################################
def create_fcst_prog_png_plots_noprecip(str_prodid, str_utc_start, str_utc_stop, str_where, str_plot_dir):
    
    try:
        print('--> Creating subplots')
           
        # Create sqlalchemy connection engine... 
        engine = create_engine('postgresql://', creator=connect_ini, poolclass=NullPool, execution_options=dict(stream_results=True))
       
        plt.ioff()    
        lst_lid=[] 
           
        # Get all desired lid's from the location_nwm table...
        print('\tAccessing location table...')
        str_sql = 'select distinct gage_id, lid, nws_name, da from nwm_location {};'.format(str_where)
        df_loc = pd.read_sql_query(str_sql, engine)
    #    df_loc['da'] = df_loc['da'].astype(np.float)       
        if df_loc.empty:
            print('problem with location table')
        
        # Generate IN string...
        print('\tGenerating IN string...')
        lst_lids=[]
        str_lid_in=''
        for tpl in df_loc.itertuples():
            str_lid_in += '\'{}\', '.format(tpl.lid.strip()) 
            lst_lids.append(tpl.lid.strip())
        str_lid_in = str_lid_in[:-2]
        
        # Get all observations for selected stations within time window...
        # You need to add time to the end of the query depending on which forecast was selected...
        if str_prodid=='NWM_medium':
            deltaT = pd.to_timedelta('10 days')
        elif str_prodid=='NWM_short':
            deltaT = pd.to_timedelta('18 hours')
        else:
            str_obs_utc_stop = str_utc_stop   # NOTE: UPDATE THIS FOR LONG RANGE
                       
        tstamp = pd.to_datetime(str_utc_stop) + deltaT                       
        str_obs_utc_stop = tstamp.strftime('%Y-%m-%d %H:%m:%S')
            
        # Get all USGS OBS for selected stations within time window...
        print('\tGetting usgs streamflow observations...')       
        str_sql = 'select lid, value, obstime from discharge where obstime >= \'{}\' and obstime <= \'{}\' and extract (minute from obstime)=0 and lid in ({}) order by obstime asc;'.format(str_utc_start, str_obs_utc_stop, str_lid_in)
        df_usgsQ = pd.read_sql_query(str_sql, engine)
        if df_usgsQ.empty:
            print('df_usgsQ does not contain any USGS data for this time period')      
        df_usgsQ.set_index('obstime', inplace=True)
        df_usgsQ['value'] = df_usgsQ['value'].astype(np.float) #/35.3147 # convert to cms              
        df_usgsQ_gpd = df_usgsQ.groupby('lid')
            
        # Get all FCSTs...
    #    start_time = timeit.default_timer()
        print('\tGetting fcsts...')       
        str_sql = 'select lid, value, basistime, validtime from fcstdischarge where basistime >= \'{}\' and basistime <= \'{}\' and product_id=\'{}\' and lid in ({});'.format(str_utc_start, str_utc_stop, str_prodid, str_lid_in)
        gen_simflow = pd.read_sql_query(str_sql, engine, chunksize=25000)       
        lst_dfsim=[]
        for df_chunk in gen_simflow:
            try:
                lst_dfsim.append(df_chunk)  
            except: # Exception, e:
                continue            
        df_chunks = pd.concat(lst_dfsim)
        df_chunks.sort_values('validtime', inplace=True)
        
        df_chunks_gpd = df_chunks.groupby('lid') # returns an iterable          
#        # ORDER BY VALIDTIME       
       
        # Get AnA...
        print('\tGetting AnA...')    
        str_sql = 'select lid, value, basistime, validtime from fcstdischarge where basistime >= \'{}\' and basistime <= \'{}\' and product_id=\'NWM_AnA\' and lid in ({}) order by validtime asc;'.format(str_utc_start, str_utc_stop, str_lid_in)
        df_AnA = pd.read_sql_query(str_sql, engine)
        gp_AnA = df_AnA.groupby('lid', sort=False)        
       
        # For now...
#        lst_skip = ['NW440','BREF1','MLBA4','LNGS1','CEDG1','SHAV2','IRCT1','EDNT2','CHIN7','IVNC1']
        i_small = 0
        for tpl in df_loc.itertuples(index=True):
            
            lst_dfsim=[]
            
            indx=tpl.Index
            i_lid=tpl.lid
            
#            if i_lid in lst_skip:
#                continue # skip this one
            
#            if indx>5:
#                break  
    
            print('\t{} -- {}'.format(indx,i_lid))   
               
            lst_lid.append(i_lid)                
                          
            # Initialize plot figure...        
            fig, ax = plt.subplots(1)
            
            try:                    
                df_simflow = df_chunks_gpd.get_group(i_lid)                  
                df_ana = gp_AnA.get_group(i_lid)
            except:
                print('lid=\'{}\' Does not contain any NWM data for this time period'.format(i_lid))
                continue 
                
            # Append AnA...
            df_simflow = df_simflow.append(df_ana)                
            df_ana = df_ana.set_index('validtime') 
              
                
            # Group by basistime and plot all fcsts...
            grouped = df_simflow.groupby('basistime')
        
            for i_basistime, df_fcst in grouped:
                
                df_fcst = df_fcst.set_index('validtime')  
                df_fcst.sort_index(axis=0, inplace=True)
               
                 # Plotting...    
                ax.plot(df_fcst.index, df_fcst['value'], label=np.str(pd.to_datetime(i_basistime)),linewidth=1.5)
                ax.plot(df_ana.index, df_ana['value'], marker='o', markersize=5, fillstyle='none', color='0.5', linestyle='none') 
                
                
            # Get observed time...
            try:           
                df_usgs = df_usgsQ_gpd.get_group(i_lid)
                                
                if df_usgs.value.mean() < 10.0:
                    print('\tSkipping this small stream!')
                    i_small += 1
                    continue
                
                df_usgs.sort_index(axis=0, inplace=True)
                # Plot observed flow...
                ax.plot(df_usgs.index, df_usgs['value'], linewidth=3.0, color='black', linestyle=':', label=tpl.gage_id)      
            except:
                print('\tNo USGS flow found for lid: {}'.format(i_lid))
                   
            # Formatting and annotations...
            ax.set_title('{} | {} | Drainage Area: {} $mi^2$'.format(i_lid.strip(),tpl.nws_name.strip(), tpl.da.strip()), fontsize=10)
            
            fig.autofmt_xdate()
        
            ax.set_ylabel('Q (cfs)')     
            
            # Save plots as png...
            figname = '{}{}.png'.format(str_plot_dir, i_lid)
            fig.savefig(figname, dpi=90, bbox_inches='tight')        
            
            fig.clf()
            plt.close(fig)        
            plt.close('all')
            
        print('\t {} small streams skipped'.format(i_small))

    except Exception as e:
        print('\r\nError creating png plots single-version! | Exception: {} \n'.format(e))
        
    finally:
        gc.collect()     
        
# ##############################################################################
#   << Query builder for MAIN >>
# ##############################################################################
def build_lid_query_over_timespan(str_proid, str_start_time, str_stop_time):
    
    # (Limit analysis based on some criteria)
    engine_ini = create_engine('postgresql://', creator=connect_ini, poolclass=NullPool)
    
    print('\tBuilding the whereclause...')
    str_sql = 'select distinct lid from nwm_location where ref=\'Y\' and atv=\'Y\';'
    df_lids = pd.read_sql_query(str_sql, engine_ini) 
    df_lids['lid'] = df_lids['lid'].str.strip()       
    
    if df_lids.empty:
        print('\t<< No lids found! >>\n')
        sys.exit()
        
    df_lids.drop_duplicates('lid', inplace=True)
    
    # Generate IN string...
    print('\tGenerating IN string...\n')
    str_lid_in=''
    for tpl in df_lids.itertuples():
        str_lid_in += '\'{}\', '.format(tpl.lid.strip())           
    str_in = str_lid_in[:-2]    
     
    str_where = 'where lid in ({})'.format(str_in) 
        
    return str_where

## #############################################################################
##  <<< MAIN >>>
## #############################################################################    
def main():
    
    print('\n<<< Start >>>\r\n')
    start_time_0 = timeit.default_timer()  
    
    # ===========  TESTING ========================
#    yr = 2017  
#    mon = 3
#    day = 20
#    str_prod_id = 'NWM_short'
    # =============================================
    
    lst_invars = sys.argv[1:]
   
    str_prod_id=lst_invars[0]
    str_yr = lst_invars[1]
    str_mon = lst_invars[2]
    str_day = lst_invars[3]
    
    yr = int(str_yr)  
    mon = int(str_mon) 
    day = int(str_day) 
    
    # =============================================

    str_preciprastername = 'test' # not used
        
    # ---------------- << Paths >> ----------------
#    str_subplot_dir = '/home/share/nwm/stats_plots/programs/bubbleplots/'
    str_subplot_dir = '/home/share/sams_scripts/bubbleplots/'
      
    str_start_time = '{}-{:02}-{:02} 00:00:00'.format(yr,mon,day)
    str_stop_time = '{}-{:02}-{:02} 23:00:00'.format(yr,mon,day) 

    
    print('Forecast type:  {}'.format(str_prod_id))
    print('Start time:  {}'.format(str_start_time))
    print('Stop time:  {}'.format(str_stop_time))
    
    
    # << QUERY >> 
    # From function...
    str_whereclause = build_lid_query_over_timespan(str_prod_id, str_start_time,  str_stop_time)       
      
    # << POP-UP PLOTS >> 
    create_fcst_prog_png_plots_noprecip(str_prod_id, str_start_time, str_stop_time, str_whereclause, str_subplot_dir)
        
    # << FOLIUM MAP >>
    create_folium_map_png(str_subplot_dir, str_start_time, str_stop_time, str_whereclause, str_prod_id, str_preciprastername, yr, mon, day)   


    print('\n<<< Done! >>>\r\n')
    print('   Total elapsed time:  ' + str(timeit.default_timer() - start_time_0)) ## TIMER 
    
    gc.collect()

## #############################################################################
##  <<< END MAIN >>>
## #############################################################################
    
if __name__ == "__main__":
    main()       
