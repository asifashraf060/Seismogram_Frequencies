#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 17:18:10 2022

@author: asifashraf
"""
import glob
import os
#import obspy
import obspy
from obspy.core import read
from obspy.core import UTCDateTime
from obspy.clients.fdsn.client import Client
#import EMD
from PyEMD import EMD
from hht_tools import spectrogram
import numpy as np
from scipy import signal
from scipy.signal import chirp
import scipy.fft
#import for plotting
import matplotlib as mtplt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#To write and read data in a csv file
import csv

plt.close('all')
#INPUT
inp_dir  = '/Users/asifashraf/Documents/Output/Miniseed_EMDexp/EQ_2021AntelopeValley/' #directory of the experiment (see 'get_eq_USGS_IRIS.py' script)

#Auto path calculation
download_dir  = (inp_dir + 'download/')
mseed_path    = (download_dir+'*.mseed')
csv_path      = (inp_dir+'info.csv')
eq_info_path  = (inp_dir+'EQinfo.csv')
analysis_dir  = (inp_dir + 'analysis/')
analysis_dir1 = (analysis_dir + 'IMF_plot/')
analysis_dir2 = (analysis_dir + 'IMF_analysis/')
result_dir    = (analysis_dir + 'result/')
try:
    os.mkdir(analysis_dir)
except:
    print('analysis directory already exists')
    pass
try:
    os.mkdir(analysis_dir1)
except:
    print('analysis directory already exists')
    pass
try:
    os.mkdir(analysis_dir2)
except:
    print('analysis directory already exists')
    pass
try:
    os.mkdir(result_dir)
except:
    print('result directory already exists')
    pass
mseed_files  = []
for file in glob.iglob(mseed_path, recursive = True):
        mseed_files.append(file)
print([str(len(mseed_files)) + ' mseed files have been imported'])

print('getting eq and station info...')
eq_info  = np.loadtxt(eq_info_path, delimiter = ',')
eq_lon   = eq_info[0]
eq_lat   = eq_info[1]
eq_depth = eq_info[2]
eq_mag   = eq_info[3]
info_imp = np.loadtxt(csv_path, delimiter=',', usecols = (1,2,3)) #imported info regarding stations
lons     = info_imp[:,0]
lats     = info_imp[:,1]
dists    = info_imp[:,2]

print('calculating PGAs...')
PGAs = []
for i in range(len(mseed_files)):
    st  = read(mseed_files[i])
    tr  = st[0]
    PGA = tr.data.max()
    PGAs.append(PGA)

plt.figure(1)
plt.scatter(lons, lats, s = 20, c = PGAs, norm=colors.LogNorm(), cmap = 'seismic')
plt.colorbar()
plt.plot(eq_lon, eq_lat, '*g')
plt.title('PGA at stations')
plt.savefig(result_dir+ 'PGA_st'+'.png',format = 'png' )

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

plt.figure(2)
fig = plt.gcf()
fig.set_size_inches(20, 10)

print('applying EMD on the signals...')
for i in range(len(mseed_files)):
            st       = read(mseed_files[i])
            tr       = st[0]
            y        = tr.data
            peak_ind = y.argmax()                                     #index of PGA
            tt       = (tr.stats.npts/tr.stats.sampling_rate)         #total time of the seismogram
            npts     = tr.stats.npts                                  #number of samples in seismogram
            dt       = tt/npts                                        #t1-t2
            t        = np.linspace(0,tt,npts)                         #time axis in sec
            peakT    = t[peak_ind]
            #HHT spectograms calculation
            #period bins
            min_period = 2*dt #Nyquist
            max_period = 500
            n_bins     = 1000
            
            peakT_t_ind        = np.where(t == (find_nearest(t, value = peakT)))
            peakT_tPlus10_ind  = np.where(t == (find_nearest(t, value = peakT+10)))
            peakT_tminus10_ind = np.where(t == (find_nearest(t, value = peakT-10)))
            y_cut_aroundPeakT  = y[(peakT_tminus10_ind[0][0]):(peakT_tPlus10_ind[0][0])]
            t_cut_aroundPeakT  = t[(peakT_tminus10_ind[0][0]):(peakT_tPlus10_ind[0][0])]
            
            #Total y
            period_bins = np.logspace(np.log10(min_period), np.log10(max_period), n_bins) #logarithmically spaced period bins
            timeSpec, periodSpec, ampSpec, IMFs = spectrogram(t_cut_aroundPeakT, y_cut_aroundPeakT, (t_cut_aroundPeakT[0]-t_cut_aroundPeakT[1]), [min_period, max_period], period_bins = period_bins)
              
            time_range       = timeSpec[0]
            peakT_timeRange  = find_nearest(time_range, value = peakT)
            peakT_timeSpec   = np.where(time_range == peakT_timeRange)
            peakT_ampSpec    = ampSpec[:,peakT_timeSpec]
            peakT_periodSpec = periodSpec[:,0]
            
            #sum IMFs
            s = np.zeros(len(IMFs[0]))
            for k in range(len(IMFs)):
                IMF = IMFs[k]
                s   = s + IMF
            
            df_s = y_cut_aroundPeakT[0:len(s)] - s
            fig_no = 333+i
            plt.figure(fig_no)
            fig = plt.gcf()
            fig.set_size_inches(20, 10)
            title = ('IMF analysis for the signal no. ' + str(i+1)+ ' of distance ' + str(dists[i]) + ' km')
            plt.subplot(3,1,1)
            plt.plot(t_cut_aroundPeakT, y_cut_aroundPeakT, '-b', label = 'Seismogram around 10 sec of PGA')
            plt.legend()
            plt.title(title)
            plt.subplot(3,1,2)
            plt.plot(t_cut_aroundPeakT[0:len(s)], s, '-r', label = 'Sum of IMFs')
            plt.legend()
            plt.subplot(3,1,3)
            plt.plot(t_cut_aroundPeakT[0:len(s)], df_s, '.g', label = 'Difference between signal and sum of IMFs')
            plt.legend()
            plt.savefig(analysis_dir2+ str(i)+'_IMF'+'.png',format = 'png' )
            plt.close(fig = fig_no)
            
            #plot IMFs
            print('plotting IMFs...')
            fig_no = i+999
            plt.figure(fig_no)           
            title = ('IMFs for the signal no. ' + str(i+1)+ ' of distance ' + str(dists[i]) + ' km')
            plt.title(title)
            for j in range(len(IMFs)):
                IMF    = IMFs[j]
                plt.subplot(len(IMFs),1,(j+1))
                plt.plot(t_cut_aroundPeakT[0:len(IMF)], IMF,'-b')
            plt.savefig(analysis_dir1 +str(fig_no)+'_IMF'+'.png',format = 'png' )
            plt.close(fig = fig_no)
            #plot FQ
            try:
                worthy_index     = [i for i,x in enumerate(peakT_ampSpec) if x>(peakT_ampSpec.max()*.1)]
                dominant_index   = np.where(peakT_ampSpec == peakT_ampSpec.max())
                fqs              = 1./peakT_periodSpec[worthy_index]
                dominant_fq      = 1./peakT_periodSpec[dominant_index[0][0]]
                plt.figure(2)
                stdv = np.std(fqs)
                #plt.errorbar(dists[i], np.mean(fqs), yerr = (fqs.max()-fqs.min()), ecolor = 'red', capsize=3)
                rng = [fqs.min(),fqs.max()]
                y_axis = [np.zeros(len(rng))+dists[i]]
                plt.plot(y_axis[0], np.array(rng), '.-b' )
                print('EMD applied on signal no, ' + str(i+1) + ' of distance ' + str(dists[i]) + ' km' )
                plt.title('Frequencies at PGA with distance')
                plt.xlabel('Distance between source and epicenter (km)')
                plt.ylabel('Frequency (Hz)')
                plt.plot(dists[i], dominant_fq, 'vr')
            except:
                pass
            plt.figure(2)
            plt.savefig(result_dir+ 'Fq_dist'+'.png',format = 'png' )
print('finally done!')




# plt.figure(3)
# plt.gca().set_yscale('log')
# plt.pcolormesh(timeSpec,periodSpec,ampSpec,cmap='hot')
# #for plotting into log scale
# fq_prcnt     = .5
# pr_fr        = fq_prcnt/100
# dsp_log_emd  = ampSpec.max()*pr_fr
# plt.ylabel('Period (s)')
# plt.pcolormesh(timeSpec,periodSpec,ampSpec,norm = colors.LogNorm(vmin = dsp_log_emd, vmax = ampSpec.max()), cmap='hot')
# plt.plot(peakT, max_period-100, 'vg', label = 'PGA')
# plt.xlabel('Time (s)')
# plt.title('frequency plot')
# plt.legend()