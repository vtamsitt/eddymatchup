# ### function to match up oceanographic data with Chelton eddy database
# Author: Veronica Tamsitt, August 2021
#
# The altimetric Mesoscale Eddy Trajectories Atlas (META2.0) was produced by SSALTO/DUACS and distributed by AVISO+ (https://aviso.altimetry.fr ) with support from CNES, in collaboration with Oregon State University with support from NASA.
#
# To do:
# * make sure eddy database download works
# * add other aviso eddy trajectory products
# * Add error message for times outside of eddy database

#import dependent modules
import numpy as np
import xarray as xr
import ftplib
import getpass
import datetime
import glob


def match(lons,lats,datetimes,database='META3.2_allsat',download_database=True,
          latmin=-90,latmax=-35,lonmin=0,lonmax=360,hourrange=24,radiusrange=1):
    
    """
    Inputs:
    lons (array of floats): required; longitude of data to match up with eddy database, -180 to 180 or 0 to 360
    lats (array of floats): required; latitude of data to match up with eddy database, -90 to 90
    datetimes (array of python datetime): required; dates and times of data to match up with eddy database
    database (string): optional; which eddy tracking database/algorithm to use, default is META3.2_allsat
    download_database (True/False): optional, specify whether to download AVISO eddy database file or use existing saved file. 
    latmax (float): optional; maximum latitude to consider for matchups
    hourrange (int): optional; maximum time difference (in hours) to consider for eddy matchup, default is +/- 24 hours
    radiusrange (int or float): optional; maximum distance from eddy centre to include in matchups, such that the max distance considered is radiusrange*(radius from eddy centre to maximum speed). Default value is 1. 
    """
    
    #define constant for calculating distance
    dr = 2*np.pi*(6371.136*1000)/360.0

    #download and load eddy database - requires aviso login (FOR NOW JUST USE PRE-DOWNLOADED FILE)
    if database == 'META3.2_allsat':
        if download_database: 
            filepath = '/value-added/eddy-trajectory/delayed-time/META3.2_DT_allsat/' 
            IP = 'ftp-access.aviso.altimetry.fr'

            #input username and password
            USERNAME = input('Enter your username: ')
            PASSWORD = getpass.getpass('Enter your password: ')

            ftp = ftplib.FTP(IP) 
            ftp.login(USERNAME,PASSWORD) 
            ftp.cwd(filepath)
            fnames = ftp.nlst('*long*')
            for f in fnames:
                print(f)
                ftp.retrbinary("RETR " + f ,open(f, 'wb').write)
            ftp.quit()
        
        
        filename_cyclonic = glob.glob('META3.2_DT_allsat_Cyclonic_long_*.nc')
        filename_anticyclonic = glob.glob('META3.2_DT_allsat_Anticyclonic_long_*.nc')
        ed_c = xr.load_dataset(filename_cyclonic[0])
        ed_ac = xr.load_dataset(filename_anticyclonic[0])

        #reduce eddy database to south of latmax, within time frame of datetimes
        datetimes_n = datetimes.astype('datetime64')
        #for n,d in enumerate(datetimes):
        #    datetimes_n[n] = np.datetime64(d)
        datemin = datetimes_n.min()
        datemax = datetimes_n.max()
     
        #check longitude go from 0 to 360
        lons_pos = lons.copy() 
        if np.any(lons_pos<0):
            lons_pos[lons_pos<0] = lons_pos[lons_pos<0]+360.
            

        eddy_c = ed_c.where((ed_c.time>datemin) & (ed_c.time<datemax),drop=True) 
        eddy_c = eddy_c.where((eddy_c.latitude>latmin) & (eddy_c.latitude<latmax),drop=True)
        eddy_c = eddy_c.where((eddy_c.longitude>lonmin) & (eddy_c.longitude<lonmax),drop=True)

        eddy_ac = ed_ac.where((ed_ac.time>datemin) & (ed_ac.time<datemax),drop=True) 
        eddy_ac = eddy_ac.where((eddy_ac.latitude>latmin) & (eddy_ac.latitude<latmax),drop=True)
        eddy_ac = eddy_ac.where((eddy_ac.longitude>lonmin) & (eddy_ac.longitude<lonmax),drop=True)

        
        #now match up data set with eddy locations
   
        #set time window to search for eddies  
        dt = np.timedelta64(hourrange,'h')
    
        #loop through all obs data times and compare locations to eddy database
        match_type = np.zeros((len(datetimes)))
        match_track = np.empty((len(datetimes)))
        match_track[:] = np.nan
        match_time = np.copy(match_track)
        match_lat = np.copy(match_track)
        match_lon = np.copy(match_track)
        match_amp = np.copy(match_track)
        match_speed = np.copy(match_track)
        match_radius = np.copy(match_track)
        match_age = np.copy(match_track)
        match_dist = np.copy(match_track)
        
        for t in range(len(datetimes)):
            
            #run for both cyclonic and anticyclonic, then combine
            
            #find all eddies within +/- dt of the obs time
            eddy_c_day = eddy_c.where((eddy_c.time>datetimes_n[t]-dt) & (eddy_c.time<datetimes_n[t]+dt),drop=True)
            eddy_ac_day = eddy_ac.where((eddy_ac.time>datetimes_n[t]-dt) & (eddy_ac.time<datetimes_n[t]+dt),drop=True)
   
            #compute distance from eddy centre lats/lons
            dy = (eddy_c_day.latitude.values-lats[t])*dr
            dx = (eddy_c_day.longitude.values-lons_pos[t])*dr*np.cos(np.deg2rad(lats[t]))
            dists = np.sqrt(dx*dx + dy*dy)
            eddy_c_day['distance'] = xr.DataArray(dists,dims=['obs'])
        
            #repeat for anticyclonic
            dy = (eddy_ac_day.latitude.values-lats[t])*dr
            dx = (eddy_ac_day.longitude.values-lons_pos[t])*dr*np.cos(np.deg2rad(lats[t]))
            dists = np.sqrt(dx*dx + dy*dy)
            eddy_ac_day['distance'] = xr.DataArray(dists,dims=['obs'])
        
            #find if distance less than radiusrange*eddy radius                 
            if np.any(eddy_c_day.distance < radiusrange*eddy_c_day.speed_radius) or \
                np.any(eddy_ac_day.distance < radiusrange*eddy_ac_day.speed_radius): 
                
                eddy_c_match = eddy_c_day.where(eddy_c_day.distance<(radiusrange*eddy_c_day.speed_radius),drop=True)
                eddy_ac_match = eddy_ac_day.where(eddy_ac_day.distance<(radiusrange*eddy_ac_day.speed_radius),drop=True)
                              
                #if 1 match save it 
                if len(eddy_c_match.obs)+len(eddy_ac_match.obs) == 1:
                    if len(eddy_c_match.obs):
                        match_type[t] = -1 #-1 is cyclonic
                        match_track[t] = eddy_c_match.track
                        match_time[t] = eddy_c_match.time
                        match_lat[t] = eddy_c_match.latitude
                        match_lon[t] = eddy_c_match.longitude
                        match_amp[t] = eddy_c_match.amplitude
                        match_speed[t] = eddy_c_match.speed_average
                        match_radius[t] = eddy_c_match.speed_radius
                        match_age[t] = eddy_c_match.observation_number
                        match_dist[t] = eddy_c_match.distance
                    elif len(eddy_ac_match.obs):
                        match_type[t] = 1 #-1 is cyclonic
                        match_track[t] = eddy_ac_match.track
                        match_time[t] = eddy_ac_match.time
                        match_lat[t] = eddy_ac_match.latitude
                        match_lon[t] = eddy_ac_match.longitude
                        match_amp[t] = eddy_ac_match.amplitude
                        match_speed[t] = eddy_ac_match.speed_average
                        match_radius[t] = eddy_ac_match.speed_radius
                        match_age[t] = eddy_ac_match.observation_number
                        match_dist[t] = eddy_ac_match.distance
                #else >1 match find closest match and save it
                else: 
                    if len(eddy_c_match.obs)>=1:
                        closest_ind_c = eddy_c_match.distance.argmin()
                    if len(eddy_ac_match.obs)>=1:
                        closest_ind_ac = eddy_ac_match.distance.argmin()
                    if len(closest_ind_c) and len(closest_ind_ac):
                        if eddy_c_match.distance[closest_ind_c]<eddy_ac_match.distance[closest_ind_ac]:
                            match_type[t] = -1
                            match_track[t] = eddy_c_match.track[closest_ind_c]
                            match_time[t] = eddy_c_match.time[closest_ind_c]
                            match_lat[t] = eddy_c_match.latitude[closest_ind_c]
                            match_lon[t] = eddy_c_match.longitude[closest_ind_c]
                            match_amp[t] = eddy_c_match.amplitude[closest_ind_c]
                            match_speed[t] = eddy_c_match.speed_average[closest_ind_c]
                            match_radius[t] = eddy_c_match.speed_radius[closest_ind_c]
                            match_age[t] = eddy_c_match.observation_number[closest_ind_c]
                            match_dist[t] = eddy_c_match.distance[closest_ind_c]
                        else:
                            match_type[t] = 1
                            match_track[t] = eddy_ac_match.track[closest_ind_ac]
                            match_time[t] = eddy_ac_match.time[closest_ind_ac]
                            match_lat[t] = eddy_ac_match.latitude[closest_ind_ac]
                            match_lon[t] = eddy_ac_match.longitude[closest_ind_ac]
                            match_amp[t] = eddy_ac_match.amplitude[closest_ind_ac]
                            match_speed[t] = eddy_ac_match.speed_average[closest_ind_ac]
                            match_radius[t] = eddy_ac_match.speed_radius[closest_ind_ac]
                            match_age[t] = eddy_ac_match.observation_number[closest_ind_ac]
                            match_dist[t] = eddy_ac_match.distance[closest_ind_ac]
                    elif len(closest_ind_c):
                        match_type[t] = -1
                        match_track[t] = eddy_c_match.track[closest_ind_c]
                        match_time[t] = eddy_c_match.time[closest_ind_c]
                        match_lat[t] = eddy_c_match.latitude[closest_ind_c]
                        match_lon[t] = eddy_c_match.longitude[closest_ind_c]
                        match_amp[t] = eddy_c_match.amplitude[closest_ind_c]
                        match_speed[t] = eddy_c_match.speed_average[closest_ind_c]
                        match_radius[t] = eddy_c_match.speed_radius[closest_ind_c]
                        match_age[t] = eddy_c_match.observation_number[closest_ind_c]
                        match_dist[t] = eddy_c_match.distance[closest_ind_c]
                    elif len(closest_ind_ac):
                        match_type[t] = 1
                        match_track[t] = eddy_ac_match.track[closest_ind_ac]
                        match_time[t] = eddy_ac_match.time[closest_ind_ac]
                        match_lat[t] = eddy_ac_match.latitude[closest_ind_ac]
                        match_lon[t] = eddy_ac_match.longitude[closest_ind_ac]
                        match_amp[t] = eddy_ac_match.amplitude[closest_ind_ac]
                        match_speed[t] = eddy_ac_match.speed_average[closest_ind_ac]
                        match_radius[t] = eddy_ac_match.speed_radius[closest_ind_ac]
                        match_age[t] = eddy_ac_match.observation_number[closest_ind_ac]
                        match_dist[t] = eddy_ac_match.distance[closest_ind_ac]                            
    
    if database == 'META2.0':
        
        if download_database:
            filepath = '/value-added/eddy-trajectory/META2.0/' 
            IP = 'ftp-access.aviso.altimetry.fr'

            #input username and password
            USERNAME = input('Enter your username: ')
            PASSWORD = getpass.getpass('Enter your password: ')

            ftp = ftplib.FTP(IP) 
            ftp.login(USERNAME,PASSWORD) 
            ftp.cwd(filepath)
            filename = ftp.retrlines('LIST *eddy_trajectory*') 
            ftp.retrbinary("RETR " + filename ,open(filename, 'wb').write)
            ftp.quit()  
        
        filename_load = glob.glob('./data/eddy_trajectory_dt_2.0*.nc')
        ed = xr.load_dataset(filename_load[0])

        #reduce eddy database to south of latmax, within time frame of datetimes
        datemin = datetimes.min()
        datemax = datetimes.max()

        eddy = ed.where((ed.time>datemin) & (ed.time<datemax),drop=True) 
        eddy = eddy.where((eddy.latitude>latmin) & (eddy.latitude<latmax),drop=True)
     
        #check longitude go from 0 to 360
        lons_pos = lons.copy() 
        if np.any(lons_pos<0):
            lons_pos[lons_pos<0] = lons_pos[lons_pos<0]+360.
    
        #now match up data set with eddy locations
   
        #set time window to search for eddies  
        dt = np.timedelta64(hourrange,'h')
    
        #loop through all obs data times and compare locations to eddy database
        match_type = np.zeros((len(datetimes)))
        match_track = np.empty((len(datetimes)))
        match_track[:] = np.nan
        match_time = np.copy(match_track)
        match_lat = np.copy(match_track)
        match_lon = np.copy(match_track)
        match_amp = np.copy(match_track)
        match_speed = np.copy(match_track)
        match_radius = np.copy(match_track)
        match_age = np.copy(match_track)
        match_dist = np.copy(match_track)
        for t in range(len(datetimes)):
            #find all eddies within +/- dt of the obs time
            eddy_day = eddy.where((eddy.time>datetimes.iloc[t]-dt) & (eddy.time<datetimes.iloc[t]+dt),drop=True)
   
            #compute distance from eddy centre lats/lons
            dy = (eddy_day.latitude.values-lats.iloc[t])*dr
            dx = (eddy_day.longitude.values-lons_pos.iloc[t])*dr*np.cos(np.deg2rad(lats.iloc[t]))
            dists = np.sqrt(dx*dx + dy*dy)
            eddy_day['distance'] = xr.DataArray(dists,dims=['obs'])
        
            #find if distance less than radiusrange*eddy radius                 
            if np.any(eddy_day.distance < radiusrange*eddy_day.speed_radius): 
                eddy_match = eddy_day.where(eddy_day.distance<(radiusrange*eddy_day.speed_radius),drop=True)
            
                #if 1 match save it 
                if len(eddy_match.obs) == 1:
                    match_type[t] = eddy_match.cyclonic_type
                    match_track[t] = eddy_match.track
                    match_time[t] = eddy_match.time
                    match_lat[t] = eddy_match.latitude
                    match_lon[t] = eddy_match.longitude
                    match_amp[t] = eddy_match.amplitude
                    match_speed[t] = eddy_match.speed_average
                    match_radius[t] = eddy_match.speed_radius
                    match_age[t] = eddy_match.observation_number
                    match_dist[t] = eddy_match.distance
                #else >1 match find closest match and save it
                else: 
                    closest_ind = eddy_match.distance.argmin()
                    match_type[t] = eddy_match.cyclonic_type[closest_ind]
                    match_track[t] = eddy_match.track[closest_ind]
                    match_time[t] = eddy_match.time[closest_ind]
                    match_lat[t] = eddy_match.latitude[closest_ind]
                    match_lon[t] = eddy_match.longitude[closest_ind]
                    match_amp[t] = eddy_match.amplitude[closest_ind]
                    match_speed[t] = eddy_match.speed_average[closest_ind]
                    match_radius[t] = eddy_match.speed_radius[closest_ind]
                    match_age[t] = eddy_match.observation_number[closest_ind]
                    match_dist[t] = eddy_match.distance[closest_ind]
    
    #index of times with obs matchups
    match_ind = np.flatnonzero(match_type)
    
    nmatches = len(match_ind)
    print('Total eddy matchups = ' + str(nmatches))
    unique_eddies,unique_ind = np.unique(match_track[match_ind],return_index=True)
    neddies = len(unique_eddies)
    print('Total unique eddy IDs = ' + str(neddies))

 
    #return matchups as xarray Dataset to be merged with original data (use data time as dimensions?)
    # define data with variable attributes
    data_vars = dict(
                    eddy_type = (["obs"],match_type),
                    eddy_ID = (["obs"],match_track),
                    eddy_lat = (["obs"],match_lat),
                    eddy_lon = (["obs"],match_lon),
                    eddy_time = (["obs"],match_time),
                    eddy_amplitude = (["obs"],match_amp),
                    eddy_vmax = (["obs"],match_speed),
                    eddy_rad_to_vmax = (["obs"],match_radius),
                    eddy_age = (["obs"],match_age),
                    eddy_dist_to_ctr = (["obs"],match_dist)
    )
                    
    # define coordinates
    coords = dict(
                 time = (["obs"],datetimes),
                 lat = (["obs"],lats),
                 lon = (["obs"],lons)
    )

    # define global attributes
    attrs = dict(creation_date = datetime.datetime.now())

    # create dataset
    matchups = xr.Dataset(data_vars=data_vars, 
                coords=coords, 
                attrs=attrs)
    
    return matchups
