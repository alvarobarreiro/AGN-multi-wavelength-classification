# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 21:18:52 2022

@author: alvar
"""

import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import gridspec

# Reading table and dropping useless columns:
    
sdss_xmatch = pd.read_csv('final_xmatch.csv', keep_default_na=True)
sdss_xmatch.drop(columns=['Unnamed: 0', 'petroMag_u', 'petroMag_g', 'petroMag_r',
                          'sfr_tot_p50_2', 'sfr_tot_p16_2', 'sfr_tot_p84_2', 
                          'h_alpha_flux', 'h_alpha_flux_err', 'e_W1mag', 'e_W2mag',
                          'e_W3mag', 'e_W4mag', 'ALLWISE_sep', 'Flux4', 'e_Flux4', 
                          'Flux5', 'e_Flux5', 'Total_flux', 'Total_flux_error',
                          'XMM_sep', 'NVSS_sep', 'Flux1', 'e_Flux1', 'Flux2', 'e_Flux2',
                          'Flux3', 'e_Flux3', 'Total_luminosity_(erg/s)'], inplace=True)

#%% Same units: erg/s/cm2 for fluxes, erg/s for luminosities

freq_IR1 = (3*1e8)/(3.4*1e-6)
freq_IR2 = (3*1e8)/(4.6*1e-6)
freq_IR4 = (3*1e8)/(22.1*1e-6) 
freq_Radio = 1.4*1e9
mJy_conv = 1e-26 # constant to convert mJy to erg/s/cm2/Hz

sdss_xmatch['S1_4'] = sdss_xmatch['S1_4'].map(lambda x: mJy_conv*freq_Radio*x)
sdss_xmatch['e_S1_4'] = sdss_xmatch['e_S1_4'].map(lambda x: mJy_conv*freq_Radio*x)
sdss_xmatch['Radio_luminosity..erg.s.Hz.'] = sdss_xmatch['Radio_luminosity..erg.s.Hz.'].map(lambda x: freq_Radio*x) 

sdss_xmatch['IR_flux..erg.s.cm.2.Hz.'] = sdss_xmatch['IR_flux..erg.s.cm.2.Hz.'].map(lambda x: freq_IR4*x)
sdss_xmatch['IR_flux_error..erg.s.cm.2.Hz.'] = sdss_xmatch['IR_flux_error..erg.s.cm.2.Hz.'].map(lambda x: freq_IR4*x)
sdss_xmatch['IR_luminosity..erg.s.Hz.'] = sdss_xmatch['IR_luminosity..erg.s.Hz.'].map(lambda x: freq_IR4*x)

sdss_xmatch['IR1_flux_(erg/s/cm2/Hz)'] = sdss_xmatch['IR1_flux_(erg/s/cm2/Hz)'].map(lambda x: freq_IR1*x)
sdss_xmatch['IR1_flux_error_(erg/s/cm2/Hz)'] = sdss_xmatch['IR1_flux_error_(erg/s/cm2/Hz)'].map(lambda x: freq_IR1*x)
sdss_xmatch['IR1_luminosity_(erg/s/Hz)'] = sdss_xmatch['IR1_luminosity_(erg/s/Hz)'].map(lambda x: freq_IR1*x)

sdss_xmatch['IR2_flux_(erg/s/cm2/Hz)'] = sdss_xmatch['IR2_flux_(erg/s/cm2/Hz)'].map(lambda x: freq_IR2*x)
sdss_xmatch['IR2_flux_error_(erg/s/cm2/Hz)'] = sdss_xmatch['IR2_flux_error_(erg/s/cm2/Hz)'].map(lambda x: freq_IR2*x)
sdss_xmatch['IR2_luminosity_(erg/s/Hz)'] = sdss_xmatch['IR2_luminosity_(erg/s/Hz)'].map(lambda x: freq_IR2*x)

#%% Changing original 'naive' Hubble constant (100 km/s/Mpc) to (70/km/s/Mpc)

Hubble_ratio = (100/70)**2

sdss_xmatch['Radio_luminosity..erg.s.Hz.'] = sdss_xmatch['Radio_luminosity..erg.s.Hz.'].map(lambda x: Hubble_ratio*x)
sdss_xmatch['IR_luminosity..erg.s.Hz.'] = sdss_xmatch['IR_luminosity..erg.s.Hz.'].map(lambda x: Hubble_ratio*x)
sdss_xmatch['IR1_luminosity_(erg/s/Hz)'] = sdss_xmatch['IR1_luminosity_(erg/s/Hz)'].map(lambda x: Hubble_ratio*x)
sdss_xmatch['IR2_luminosity_(erg/s/Hz)'] = sdss_xmatch['IR2_luminosity_(erg/s/Hz)'].map(lambda x: Hubble_ratio*x)
sdss_xmatch['Hard_luminosity_(erg/s)'] = sdss_xmatch['Hard_luminosity_(erg/s)'].map(lambda x: Hubble_ratio*x)
sdss_xmatch['Soft_luminosity_(erg/s)'] = sdss_xmatch['Soft_luminosity_(erg/s)'].map(lambda x: Hubble_ratio*x)


#%% Changing order of columns

a = sdss_xmatch.iloc[:, :18]
b = sdss_xmatch.iloc[:, 21]
c = sdss_xmatch.iloc[:, 18:21]
d = sdss_xmatch.iloc[:, 22:]
sdss_xmatch = pd.concat([a,b,c,d], axis=1)

#%% Changing column names
    
sdss_xmatch.rename(columns={'sfr_tot_p50_1':'sfr_tot_p50', 'sfr_tot_p16_1':'sfr_tot_p16',
                            'sfr_tot_p84_1':'sfr_tot_p84', 'S1_4':'Radio_flux_(erg/s/cm2)',
                            'e_S1_4':'Radio_flux_error_(erg/s/cm2)',
                            'Radio_luminosity..erg.s.Hz.':'Radio_luminosity_(erg/s)',
                            'Hard_flux':'Hard_flux_(erg/s/cm2)',
                            'Hard_flux_error':'Hard_flux_error_(erg/s/cm2)',
                            'IR_flux..erg.s.cm.2.Hz.':'IR4_flux_(erg/s/cm2)',
                            'IR_flux_error..erg.s.cm.2.Hz.':'IR4_flux_error_(erg/s/cm2)',
                            'IR_luminosity..erg.s.Hz.':'IR4_luminosity_(erg/s)',
                            'Soft_flux':'Soft_flux_(erg/s/cm2)',
                            'Soft_flux_error':'Soft_flux_error_(erg/s/cm2)',
                            'IR1_flux_(erg/s/cm2/Hz)':'IR1_flux_(erg/s/cm2)',
                            'IR1_flux_error_(erg/s/cm2/Hz)':'IR1_flux_error_(erg/s/cm2)',
                            'IR1_luminosity_(erg/s/Hz)':'IR1_luminosity_(erg/s)',
                            'IR2_flux_(erg/s/cm2/Hz)':'IR2_flux_(erg/s/cm2)',
                            'IR2_flux_error_(erg/s/cm2/Hz)':'IR2_flux_error_(erg/s/cm2)',
                            'IR2_luminosity_(erg/s/Hz)':'IR2_luminosity_(erg/s)'}, inplace=True)

#%% Save as csv

sdss_xmatch.to_csv('final_xmatch_norm.csv', header=True, index=False)
#%%
'''
b = np.genfromtxt('excesses.csv', delimiter=',', dtype={'names':('IR1', 'IR2', 'IR4', 'Radio', 'Soft', 'Hard', 'IR_color', 'Radio_color', 'X_color'),
                                                        'formats': ('f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f')}, filling_values=np.nan, skip_header=1)

c = np.genfromtxt('excesses_blank.csv', delimiter=',', usecols = (0, 1), names=('a', 'b'), filling_values=np.nan, skip_header=True)
'''