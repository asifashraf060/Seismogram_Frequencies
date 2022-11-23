#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 13:16:30 2022

@author: asifashraf
"""

#asif, OCT 20
##IMPORT
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import urllib.request, json
import os
#import libcomcat
# Local imports
from libcomcat.dataframes import get_detail_data_frame
from libcomcat.search import search,get_event_by_id
import pandas as pd
#import obspy
from obspy.core import read
from obspy.core import UTCDateTime
from obspy.clients.fdsn.client import Client
#To write data in a csv file
import csv

plt.close('all')

#INPUT Parameters to get eq (Variable)
starttime    = datetime(2014, 8, 1, 0, 0) #YY/MM/DD
endtime      = datetime(2014, 8, 30, 0, 0)
minlatitude  = 37
maxlatitude  = 40
minlongitude = -123
maxlongitude = -121
out_dir      = '/Users/asifashraf/Documents/Output/Miniseed_EMDexp/'
out_name     = 'EQ_test'
#INPUT Parameters to get eq (Consistent)
minmagnitude = 1
producttype  = 'finite-fault'
#INPUT Parameters to filter eq
dist_max     = 300 #for recorded stations in km
intn_min     = 3 #minimum intensity
##Length of the miniseed file in second
len_Sec      = 200


##auto file creation
glob_dir     = out_dir+out_name
seis_dir     = glob_dir+'/seismogram/' #directory to save seismogram figures
download_dir = glob_dir+'/download/' #directory to save downloaded miniseed files
try:
    os.mkdir(out_dir)
except:
    print('Out directory already exists')
    pass
try:
    os.mkdir(glob_dir)#MAKING OUTPUT DIRECTORY
except:
    print('Global directory already exists')
    pass
try:
    os.mkdir(seis_dir)
except:
    print('Seismogram folder already exists')
    pass
try:
    os.mkdir(download_dir)
except:
    print('download directory already exists')
    pass

##Search for EQs
eq = search(starttime=starttime, endtime=endtime, minlatitude=minlatitude, maxlatitude=maxlatitude, 
                        minlongitude=minlongitude, maxlongitude=maxlongitude,
                           orderby='time', minmagnitude=minmagnitude,producttype=producttype)
print(len(eq), 'earthquake is found')
print('eq info:', eq)

if len(eq) == 0:
    raise Exception('No Earthquake can be found')
if len(eq)>1:
    i = input('which earthquake you want to investigate? enter their serial no (top = 1): ')
    earthquake = eq[i-1]
if len(eq) == 1:
    earthquake = eq[0]

plt.figure(1)
plt.subplot(131)
plt.scatter(earthquake.longitude, earthquake.latitude)
plt.xlabel('longitude')
plt.ylabel('latitude')
plt.title('location of eq')
plt.subplot(133)
plt.scatter(earthquake.magnitude, earthquake.depth)
plt.xlabel('magnitude')
plt.ylabel('depth (km)')
plt.title('magnitude vs depth')
plt.savefig(glob_dir+'/'+'eq_info.png',format = 'png' )


##Search for stations recorded these eqs
#now check if sites have shakemap with sites
Nsta_out = [] # no. of station
for kev in range(len(eq)): #looping through number of events
    id = eq[kev].id 
    print('%d of %d, %s' % (kev,len(eq),id))
    
    event = get_event_by_id(eventid=id, includedeleted=True) #don't understand get event by id
    try:
        sm = event.getProducts(product_name='shakemap')[0]
        #get url to shakempa
        url_shake_map = sm.getContentURL('stationlist.json')

        if url_shake_map != None:
            #read station list to json
            with urllib.request.urlopen(url_shake_map) as url:
                sites_data = json.load(url)

            #Loop over sites
            Nsta = 0
            for ksta in range(len(sites_data['features'])):
                if sites_data['features'][ksta]['properties']['station_type'] == 'seismic':
                    Nsta += 1
        else:
            Nsta = 0          
    except:   
        Nsta = 0  
    print(Nsta, 'stations been found')
    Nsta_out.append(Nsta)
    
    
## Filtering No. of stations
#based on having the shakemap
kev=len(eq)-1
id = eq[kev].id
event = get_event_by_id(eventid=id, includedeleted=True)
sm = event.getProducts(product_name='shakemap')[0]
url_station_list = sm.getContentURL('stationlist.json')
with urllib.request.urlopen(url_station_list) as url:
    sites_data = json.load(url)
    
## FILTER the number of stations
## get the distances from all the stations
dist  = []
for idist in range(Nsta):
    Y = (sites_data.get('features')[idist])
    Z = Y.get('properties')
    dist.append(Z.get('distance'))
##filter out stations based on desired distance
ind_dist = [i for i,x in enumerate(dist) if x<dist_max]
print('No. of station is', (len(ind_dist)), 'after distance filtering')

## get the intensity from all the stations
intn = []
for idist in range(Nsta):
    Y = (sites_data.get('features')[idist])
    Z = Y.get('properties')
    intn.append(Z.get('intensity'))
## select a desired intensity and filter out stations
ind_rmv_intn = [i for i,x in enumerate(intn) if x == 'null']
for irmv in range(len(ind_rmv_intn)):
    yy = ind_rmv_intn[irmv]
    intn[yy] = 0
ind_intn = [i for i,x in enumerate(intn) if x>intn_min]
print('No. of station is ', (len(ind_intn)), 'after intensity filtering')

ind_dist_set  = set(ind_dist)
intersection1 = list(ind_dist_set.intersection(ind_intn))
print('No. of station is ', (len(intersection1)), 'after distance + intensity filtering')

## gather information from all these stations
codes         = []
names         = []
net_types     = []
channel_types = []
distances     = []
for i in range(len(intersection1)):
    k         = intersection1[i]
    code      = sites_data['features'][k]['properties']['code']
    codes.append(code)
    name      = sites_data['features'][k]['properties']['name']
    names.append(name)
    net_type  = sites_data['features'][k]['properties']['network']
    net_types.append(net_type)
    chan_type = sites_data['features'][k]['properties']['channels'][1]['name']
    channel_types.append(chan_type)
    distance  = sites_data['features'][k]['properties']['distance']
    distances.append(distance)
##exclude based on bad networks
ind_excld  = [i for i, j in enumerate(net_types) if j == '--']
for ind_excld in sorted(ind_excld, reverse = True):
    del codes[ind_excld]
    del names[ind_excld]
    del net_types[ind_excld]
    del channel_types[ind_excld]
    del distances[ind_excld]
ind_excld  = [i for i, j in enumerate(net_types) if j == 'NSMP']
for ind_excld in sorted(ind_excld, reverse = True):
    del codes[ind_excld]
    del names[ind_excld]
    del net_types[ind_excld]
    del channel_types[ind_excld]
    del distances[ind_excld]
print('After exclusion the number of station is', len(codes))

## INPUT--ObsPy
## set up the client to download waveforms
client = Client('IRIS')
ET       = str(earthquake.time)
shotTime = ET[0:10] + 'T' + ET[11:]
#
sTime = UTCDateTime(shotTime)
eTime = sTime+len_Sec
print("start: ", sTime)
print("end: ", eTime)
#check if the stations have data
codes_string       = ','.join(codes)
net_types_string   = ','.join(net_types)
stations_with_data = client.get_stations(starttime = sTime, endtime = eTime, network = net_types_string, station = codes_string, matchtimeseries = True)
#collect codes of those stations
st_codes_with_data = []
st_lats_with_data  = []
st_lons_with_data  = []
for i in range(len(stations_with_data)):
           A = stations_with_data[i]
           for j in range(len(A)):
                       st_code_with_data = A[j].code
                       st_codes_with_data.append(st_code_with_data)
                       st_lat_with_data  = A[j].latitude
                       st_lats_with_data.append(st_lat_with_data)
                       st_lon_with_data  = A[j].longitude
                       st_lons_with_data.append(st_lon_with_data)
codes_ind_with_data = []
for i in range(len(st_codes_with_data)):
               A = st_codes_with_data[i]     
               B = codes.index(A)
               codes_ind_with_data.append(B)

#manual_ind              = [1,2,3,4,5,6,7,8,9,10,11,33,41,42,43,44,45]                
ind                     = codes_ind_with_data
#ind                     = [codes_ind_with_data[i] for i in manual_ind]
#st_lats                 = [st_lats_with_data[i] for i in manual_ind]
st_lats                 = st_lats_with_data
#st_lons                 = [st_lons_with_data[i] for i in manual_ind]
st_lons                 = st_lons_with_data
codes_with_data         = [codes[i] for i in ind]
names_with_data         = [names[i] for i in ind]
net_types_with_data     = [net_types[i] for i in ind]
channel_types_with_data = [channel_types[i] for i in ind]
distances_with_data     = [distances[i] for i in ind]
print('only',len(codes_ind_with_data), 'stations have data from',len(codes), 'stations')

##correct for Instrument response
eq_miniseed   = []
ind_with_data = [] #collect indices for sites only with data
ids           = []
print('downloading earthquake from IRIS...')
for i in range(len(ind)):
    try:
        #network
        netcode = net_types_with_data[i]
        #station
        sta     = codes_with_data[i]
        #channel
        if len(channel_types_with_data[i]) == 6:
                channels = channel_types_with_data[i][3:6]
        else:
                channels = channel_types_with_data[i]
        st = client.get_waveforms(network = netcode, station = sta, location = '*', channel = channels, starttime = sTime, endtime = eTime,attach_response = True)
        eq_miniseed.append(st[0])
        ind_with_data.append(ind[i])
        ids.append(i)
    except:
        pass
print('eqs are downloaded')
print('Finally total no. of stations with data is', len(eq_miniseed))
print('Proceeding with data correction and applying EMD')

st_lons_final   = [st_lons_with_data[i] for i in ids]
st_lats_final   = [st_lats_with_data[i] for i in ids]
distances_final = [distances_with_data[i] for i in ids]
codes_final     = [codes_with_data[i] for i in ids ]
#plotting stations info
plt.figure(1234)
plt.plot(st_lons_final, st_lats_final, '.g')
plt.plot(earthquake.longitude, earthquake.latitude, '*r')
plt.title('stations relative to earthquake')
plt.xlabel('longitude')
plt.ylabel('latitude')
plt.savefig(glob_dir+'/'+'st_eq.png',format = 'png' )
plt.close(1234)
plt.figure(1235)
plt.plot(distances_final, '+')
plt.xlabel('No. of eqs')
plt.ylabel('distance(km)')
plt.title('Distances from eq')
plt.savefig(glob_dir+'/'+'st_dist.png',format = 'png' )
plt.close(1235)
#correct for instrument sensitivity
print('starting instrument sensitivity connection...')
for i in range(len(eq_miniseed)):
            st                  = eq_miniseed[i]
            sensitivity         = st.stats.response.instrument_sensitivity.value
            st_snst_applied     = st.data/sensitivity
            eq_miniseed[i].data = st_snst_applied
print('starting baseline correction...')
#baseline correction
for i in range(len(eq_miniseed)):
            A = eq_miniseed[i].data
            correct_factor = np.mean(A[1:10])*(-1)
            B = A+correct_factor
            eq_miniseed[i].data = B
#plotting seismograms
for i in range(len(eq_miniseed)):
            fig_no = i+99
            plt.figure(fig_no)
            plt.plot(eq_miniseed[i].data)
            title = codes_with_data[i]+' dist:' + str(distances_final[i])
            plt.title(title)
            plt.ylabel('acc')
            plt.xlabel('count')
            plt.savefig( seis_dir + 'spec'+str(fig_no)+'.png',format = 'png' )
            plt.close(fig = fig_no)
print('saving the downloaded miniseeds...')
#save the downloaded miniseeds into the folder
for i in range(len(eq_miniseed)):
            data_name = (str(i)+'_'+codes_final[i] + '_' + str(distances_final[i]) + 'km')
            eq_miniseed[i].write((download_dir + data_name + '.mseed'), format = 'MSEED')
            
print('Writing station info into a CSV file...')
codes_write = np.array(codes_final).T
lons_write  = np.array(st_lons_final).T
lats_write  = np.array(st_lats_final).T
dists_write = np.array(distances_final).T
combo  = np.stack([codes_write, lons_write, lats_write, dists_write], axis = 1)
header = ['Station Codes','Stations Lons','Station Lats','Station distances']
CV = zip(codes_final)
with open((glob_dir+'/info.csv'), 'w') as f:
    writer = csv.writer(f, delimiter = ',')
    writer.writerows(combo)
file = pd.read_csv((glob_dir+'/info.csv'))
file.to_csv((glob_dir+'/info_with_header(DontUse).csv'), header = header, index = False)
pd.read_csv('info_with_header.csv')

combo2  = np.stack([float(eq[0].longitude), float(eq[0].latitude),\
                    float(eq[0].depth), float(eq[0].magnitude)], axis = 0)
header2 = ['lon', 'lat', 'depth', 'magnitude']
with open((glob_dir+'/EQinfo.csv'), 'w') as f:
    writer = csv.writer(f, delimiter = ',')
    writer.writerows(map(lambda x: [x],combo2))
file = pd.read_csv(glob_dir+'/EQinfo.csv')

print('Data import complete!!!')







