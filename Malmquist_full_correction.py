# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:33:24 2021

@author: alvar
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from matplotlib import gridspec


sdss_xmatch = pd.read_csv('final_xmatch.csv', keep_default_na=True)
print(sdss_xmatch.columns)
data_length = len(sdss_xmatch['col1'])
redshift = sdss_xmatch['Z']
log_SFR = sdss_xmatch['sfr_tot_p50_1'].values
good_SFR = np.isfinite(log_SFR)
log_M = sdss_xmatch['lgm_tot_p50'].values
good_M = np.isfinite(log_M)
sigma_SFR = (sdss_xmatch['sfr_tot_p84_1'].values - sdss_xmatch['sfr_tot_p16_1'].values)/2
sigma_M = (sdss_xmatch['lgm_tot_p84'].values - sdss_xmatch['lgm_tot_p16'].values)/2


log_L_Soft = sdss_xmatch['Soft_luminosity_(erg/s)'].values
log_L_Soft = log_L_Soft*(100/70)**2
log_L_Soft = np.log10(log_L_Soft)
good_Soft = np.isfinite(log_L_Soft)
log_L_Soft[~good_Soft] = np.nan

log_L_Hard = sdss_xmatch['Hard_luminosity_(erg/s)'].values
log_L_Hard = log_L_Hard*(100/70)**2
log_L_Hard = np.log10(log_L_Hard)
good_Hard = np.isfinite(log_L_Hard)
log_L_Hard[~good_Hard] = np.nan

log_L_Radio = sdss_xmatch['Radio_luminosity..erg.s.Hz.'].values
log_L_Radio = log_L_Radio*(1.4*1e9)
log_L_Radio = log_L_Radio*(100/70)**2
log_L_Radio = np.log10(log_L_Radio)
good_Radio = np.isfinite(log_L_Radio)
log_L_Radio[~good_Radio] = np.nan

log_L_IR4 = sdss_xmatch['IR_luminosity..erg.s.Hz.'].values
log_L_IR4 = log_L_IR4*((3*1e8)/(22.1*1e-6))
log_L_IR4 = log_L_IR4*(100/70)**2
log_L_IR4 = np.log10(log_L_IR4)
good_IR4 = np.isfinite(log_L_IR4)
log_L_IR4[~good_IR4] = np.nan

log_L_IR1 = sdss_xmatch['IR1_luminosity_(erg/s/Hz)'].values
log_L_IR1 = log_L_IR1*((3*1e8)/(3.4*1e-6))
log_L_IR1 = log_L_IR1*(100/70)**2
log_L_IR1 = np.log10(log_L_IR1)
good_IR1 = np.isfinite(log_L_IR1)
log_L_IR1[~good_IR1] = np.nan

log_L_IR2 = sdss_xmatch['IR2_luminosity_(erg/s/Hz)'].values
log_L_IR2 = log_L_IR2*((3*1e8)/(4.6*1e-6))
log_L_IR2 = log_L_IR2*(100/70)**2
log_L_IR2 = np.log10(log_L_IR2)
good_IR2 = np.isfinite(log_L_IR2)
log_L_IR2[~good_IR2] = np.nan


W1 = sdss_xmatch['W1mag'].values
W2 = sdss_xmatch['W2mag'].values
IR_color_mag = W1 - W2
color_IR = log_L_IR2 - log_L_IR1
q_24 = np.log10(sdss_xmatch['IR_flux..erg.s.cm.2.Hz.'].values/(sdss_xmatch['S1_4'].values*(1e-26)))
color_Radio_IR = log_L_Radio - log_L_IR4
color_X = log_L_Hard - log_L_Soft


#%% LIMPIEZA DE DATOS Y FILTROS EN SEÑAL RUIDO

SN_SFR = log_SFR/sigma_SFR
SN_M = log_M/sigma_M

for i in range(data_length):
    
    if ((np.isnan(log_SFR[i]) == True) or (np.isnan(log_M[i]) == True)):
        
        log_SFR[i] = np.nan
        log_M[i] = np.nan
        sigma_SFR[i] = np.nan
        sigma_M[i] = np.nan
        SN_SFR[i] = np.nan
        SN_M[i] = np.nan
        
        log_L_Radio[i] = np.nan
        log_L_IR4[i] = np.nan
        log_L_IR1[i] = np.nan
        log_L_IR2[i] = np.nan
        log_L_Soft[i] = np.nan
        log_L_Hard[i] = np.nan
        color_IR[i] = np.nan
        color_Radio_IR[i] = np.nan
        color_X[i] = np.nan
        
        IR_color_mag[i] = np.nan
        q_24[i] = np.nan
        
    elif ((log_SFR[i] == -9999) or (log_M[i] == -9999)):
        
        log_SFR[i] = np.nan
        log_M[i] = np.nan
        sigma_SFR[i] = np.nan
        sigma_M[i] = np.nan
        SN_SFR[i] = np.nan
        SN_M[i] = np.nan
        
        log_L_Radio[i] = np.nan
        log_L_IR4[i] = np.nan
        log_L_IR1[i] = np.nan
        log_L_IR2[i] = np.nan
        log_L_Soft[i] = np.nan
        log_L_Hard[i] = np.nan
        color_IR[i] = np.nan
        color_Radio_IR[i] = np.nan
        color_X[i] = np.nan
        
        IR_color_mag[i] = np.nan
        q_24[i] = np.nan
        

for i in range(data_length):
        
    if ((np.isnan(SN_M[i]) == True) or (np.isnan(SN_SFR[i]) == True)):

        log_SFR[i] = np.nan
        log_M[i] = np.nan
        sigma_SFR[i] = np.nan
        sigma_M[i] = np.nan
        
        log_L_Radio[i] = np.nan
        log_L_IR4[i] = np.nan
        log_L_IR1[i] = np.nan
        log_L_IR2[i] = np.nan
        log_L_Soft[i] = np.nan
        log_L_Hard[i] = np.nan
        color_IR[i] = np.nan
        color_Radio_IR[i] = np.nan
        color_X[i] = np.nan
        
        IR_color_mag[i] = np.nan
        q_24[i] = np.nan
        
    '''elif ((SN_M[i] < 2) or (SN_SFR[i] < 2)):
        
        log_SFR[i] = np.nan
        log_M[i] = np.nan
        sigma_SFR[i] = np.nan
        sigma_M[i] = np.nan
        
        log_L_Radio[i] = np.nan
        log_L_IR4[i] = np.nan
        log_L_IR1[i] = np.nan
        log_L_IR2[i] = np.nan
        log_L_Soft[i] = np.nan
        log_L_Hard[i] = np.nan
        color_IR[i] = np.nan
        color_Radio_IR[i] = np.nan
        color_X[i] = np.nan
        
        IR_color_mag[i] = np.nan
        q_24[i] = np.nan'''
        
'''        
for i in range(data_length):
        
    if redshift[i] < 0.01:
        
        log_SFR[i] = np.nan
        log_M[i] = np.nan
        sigma_SFR[i] = np.nan
        sigma_M[i] = np.nan
        SN_SFR[i] = np.nan
        SN_M[i] = np.nan
        
        log_L_Radio[i] = np.nan
        log_L_IR4[i] = np.nan
        log_L_IR1[i] = np.nan
        log_L_IR2[i] = np.nan
        log_L_Soft[i] = np.nan
        log_L_Hard[i] = np.nan
        color_IR[i] = np.nan
        color_Radio_IR[i] = np.nan
        color_X[i] = np.nan
        
        IR_color_mag[i] = np.nan
        q_24[i] = np.nan
        '''

sigma_Soft = sdss_xmatch['Soft_flux_error'].values
SN_Soft = sdss_xmatch['Soft_flux'].values/sigma_Soft

sigma_Hard = sdss_xmatch['Hard_flux_error'].values
SN_Hard = sdss_xmatch['Hard_flux'].values/sigma_Hard

sigma_Radio = sdss_xmatch['e_S1_4'].values
SN_Radio = sdss_xmatch['S1_4'].values/sigma_Radio

sigma_IR4 = sdss_xmatch['IR_flux_error..erg.s.cm.2.Hz.'].values
SN_IR4 = sdss_xmatch['IR_flux..erg.s.cm.2.Hz.'].values/sigma_IR4

sigma_IR1 = sdss_xmatch['IR1_flux_error_(erg/s/cm2/Hz)'].values
SN_IR1 = sdss_xmatch['IR1_flux_(erg/s/cm2/Hz)'].values/sigma_IR1

sigma_IR2 = sdss_xmatch['IR2_flux_error_(erg/s/cm2/Hz)'].values
SN_IR2 = sdss_xmatch['IR2_flux_(erg/s/cm2/Hz)']/sigma_IR2


# Filtro en flujo en rayos X blandos:
    
for i in range(data_length):
    if np.isnan(SN_Soft[i]) == True:
        log_L_Soft[i] = np.nan
        
    elif SN_Soft[i] < 2:
        log_L_Soft[i] = np.nan
        
# Filtro en flujo en rayos X duros:
    
for i in range(data_length):
    if np.isnan(SN_Hard[i]) == True:
        log_L_Hard[i] = np.nan
        
    elif SN_Hard[i] < 2:
        log_L_Hard[i] = np.nan
        
# Filtro en flujo en Radio:
    
for i in range(data_length):
    if np.isnan(SN_Radio[i]) == True:
        log_L_Radio[i] = np.nan
        
    elif SN_Radio[i] < 2:
        log_L_Radio[i] = np.nan
        
# Filtro en flujo en IR (banda W4):
    
for i in range(data_length):
    if np.isnan(SN_IR4[i]) == True:
        log_L_IR4[i] = np.nan
        
    elif SN_IR4[i] < 2:
        log_L_IR4[i] = np.nan
        
# Filtro en flujo en IR (banda W1):
    
for i in range(data_length):
    if np.isnan(SN_IR1[i]) == True:
        log_L_IR1[i] = np.nan
        
    elif SN_IR1[i] < 2:
        log_L_IR1[i] = np.nan
        
# Filtro en flujo en IR (banda W2):
    
for i in range(data_length):
    if np.isnan(SN_IR2[i]) == True:
        log_L_IR2[i] = np.nan
        
    elif SN_IR2[i] < 2:
        log_L_IR2[i] = np.nan
        
        
# Filtro en color IR:
    
for i in range(data_length):
    if ((np.isnan(SN_IR1[i]) == True) or (np.isnan(SN_IR2[i]) == True)):
        color_IR[i] = np.nan
        IR_color_mag[i] = np.nan
        
    elif ((SN_IR1[i] < 2) or (SN_IR2[i] < 2)):
        color_IR[i] = np.nan
        IR_color_mag[i] = np.nan
        
# Filtro en color Radio-IR:
    
for i in range(data_length):
    if ((np.isnan(SN_Radio[i]) == True) or (np.isnan(SN_IR4[i]) == True)):
        color_Radio_IR[i] = np.nan
        q_24[i] = np.nan
        
    elif ((SN_Radio[i] < 2) or (SN_IR4[i] < 2)):
        color_Radio_IR[i] = np.nan
        q_24[i] = np.nan
        
# Filtro en color X:
    
for i in range(data_length):
    if ((np.isnan(SN_Hard[i]) == True) or (np.isnan(SN_Soft[i]) == True)):
        color_X[i] = np.nan
        
    elif ((SN_Hard[i] < 2) or (SN_Soft[i] < 2)):
        color_X[i] = np.nan
        
        
# REMOVING EXTREME MASSES AND SFRs:
    
for i in range(data_length):
    
    if (np.isnan(log_M[i]) == False and np.isnan(log_SFR[i]) == False):
        if (log_M[i] < 8.5 or log_M[i] > 11.5 or log_SFR[i] < -2 or log_SFR[i] > 1):
            
            log_L_Radio[i] = np.nan
            log_L_IR4[i] = np.nan
            log_L_IR1[i] = np.nan
            log_L_IR2[i] = np.nan
            log_L_Soft[i] = np.nan
            log_L_Hard[i] = np.nan
            color_IR[i] = np.nan
            color_Radio_IR[i] = np.nan
            color_X[i] = np.nan
            
            IR_color_mag[i] = np.nan
            q_24[i] = np.nan
            
          

#%% GETTING V MAX INVERSE FOR MALMQUIST CORRECTION


# Sensitivity (i.e. flux limit above which we have data for each band): it's been visually estimated from
# plot(redshift, L - 2*log10(redshift)) for each frequency

sensitivity_Radio = 40.85
sensitivity_IR4 = 44.8
sensitivity_IR1 = 44.6
sensitivity_IR2 = 44.2
sensitivity_Soft = 42.8
sensitivity_Hard = 43.5
sensitivity_SDSS = 12 # plot(redshift, M - 2*log10(redshift)) in this case


def get_V_max_inverse(L, luminosities, sensitivity):
    
    z_max = 10**((L - sensitivity)/2) # redshift computed from the equation above
    # z_max = z_max.min(0).clip(0.01,0.07) 
    z_max = z_max.min(0) # min of each column for the above matrix
    faint_obj = np.where(z_max < 0.015)[0]
    z_max[faint_obj] = np.nan
    
    for i in range(len(luminosities)): # filter
        luminosities[i][faint_obj] = np.nan # we exclude every object with z_max < 0.015
        
    bad = np.isnan(z_max) # finding nans
    # z_max = np.where(z_max > 0.01, z_max, 0.01) # if z < 0.01: z = 0.01, else z = z
    z_max = np.where(z_max < 0.07, z_max, 0.07) # if z > 0.07: z = 0.07, else z = z
    z_max = np.where(bad == False, z_max, np.nan) # nans still being nans

    return 1/(z_max**3 - 1e-9) # rho, or inverse max volume


def detection_filter(L, sensitivity):
    
    for p in range(data_length):
        if ((np.isnan(L[p]) == True) or (np.isnan(redshift[p]) == True)):
            L[p] = np.nan
            
        elif ((L[p] - 2*np.log10(redshift[p])) < sensitivity):
            L[p] = np.nan
            
    return L


def detection_filter_color(L1, L2, color, sensitivity1, sensitivity2):
    
    for p in range(data_length):
        if ((np.isnan(L1[p]) == True) or (np.isnan(L2[p]) == True) or (np.isnan(redshift[p]) == True)):
            color[p] = np.nan
            
        elif (((L1[p] - 2*np.log10(redshift[p])) < sensitivity1) or 
            ((L2[p] - 2*np.log10(redshift[p])) < sensitivity2)):
            color[p] = np.nan
            
    return color


log_M = detection_filter(log_M, sensitivity_SDSS)
# No sensitivity for SFR

log_L_Soft = detection_filter(log_L_Soft, sensitivity_Soft)
log_L_Hard = detection_filter(log_L_Hard, sensitivity_Hard)
log_L_Radio = detection_filter(log_L_Radio, sensitivity_Radio)
log_L_IR4 = detection_filter(log_L_IR4, sensitivity_IR4)
log_L_IR1 = detection_filter(log_L_IR1, sensitivity_IR1)
log_L_IR2 = detection_filter(log_L_IR2, sensitivity_IR2)

color_X = detection_filter_color(log_L_Soft, log_L_Hard, color_X, sensitivity_Soft, sensitivity_Hard)
color_Radio_IR = detection_filter_color(log_L_IR4, log_L_Radio, color_Radio_IR, sensitivity_IR4, sensitivity_Radio)
color_IR = detection_filter_color(log_L_IR1, log_L_IR2, color_IR, sensitivity_IR1, sensitivity_IR2)



# Updating good values:

good_M = np.isfinite(log_M)
good_SFR = np.isfinite(log_SFR)

good_Soft = np.isfinite(log_L_Soft)
good_Hard = np.isfinite(log_L_Hard)
good_Radio = np.isfinite(log_L_Radio)
good_IR4 = np.isfinite(log_L_IR4)
good_IR1 = np.isfinite(log_L_IR1)
good_IR2 = np.isfinite(log_L_IR2)

good_color_X = np.isfinite(color_X)
good_color_Radio_IR = np.isfinite(color_Radio_IR)
good_color_IR = np.isfinite(color_IR)

# Getting V max inverse for luminosities:

matrix_Soft = np.concatenate((log_L_Soft.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_Soft_vector = np.array([sensitivity_Soft, sensitivity_SDSS]).reshape(matrix_Soft.shape[0], 1)
V_max_inverse_Soft = get_V_max_inverse(matrix_Soft, (log_L_Soft, log_M), sens_Soft_vector)

matrix_Hard = np.concatenate((log_L_Hard.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_Hard_vector = np.array([sensitivity_Hard, sensitivity_SDSS]).reshape(matrix_Hard.shape[0], 1)
V_max_inverse_Hard = get_V_max_inverse(matrix_Hard, (log_L_Hard, log_M), sens_Hard_vector)

matrix_Radio = np.concatenate((log_L_Radio.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_Radio_vector = np.array([sensitivity_Radio, sensitivity_SDSS]).reshape(matrix_Radio.shape[0], 1)
V_max_inverse_Radio = get_V_max_inverse(matrix_Radio, (log_L_Radio, log_M), sens_Radio_vector)

matrix_IR4 = np.concatenate((log_L_IR4.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_IR4_vector = np.array([sensitivity_IR4, sensitivity_SDSS]).reshape(matrix_IR4.shape[0], 1)
V_max_inverse_IR4 = get_V_max_inverse(matrix_IR4, (log_L_IR4, log_M), sens_IR4_vector)

matrix_IR1 = np.concatenate((log_L_IR1.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_IR1_vector = np.array([sensitivity_IR1, sensitivity_SDSS]).reshape(matrix_IR1.shape[0], 1)
V_max_inverse_IR1 = get_V_max_inverse(matrix_IR1, (log_L_IR1, log_M), sens_IR1_vector)

matrix_IR2 = np.concatenate((log_L_IR2.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_IR2_vector = np.array([sensitivity_IR2, sensitivity_SDSS]).reshape(matrix_IR2.shape[0], 1)
V_max_inverse_IR2 = get_V_max_inverse(matrix_IR2, (log_L_IR2, log_M), sens_IR2_vector)

# Getting V max inverse for colours:

matrix_color_X = np.concatenate((log_L_Soft.reshape(1, data_length), log_L_Hard.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_color_X_vector = np.array([sensitivity_Soft, sensitivity_Hard, sensitivity_SDSS]).reshape(matrix_color_X.shape[0], 1)
V_max_inverse_color_X = get_V_max_inverse(matrix_color_X, (log_L_Soft, log_L_Hard, log_M), sens_color_X_vector)

matrix_color_Radio_IR = np.concatenate((log_L_IR4.reshape(1, data_length), log_L_Radio.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_color_Radio_IR_vector = np.array([sensitivity_IR4, sensitivity_Radio, sensitivity_SDSS]).reshape(matrix_color_Radio_IR.shape[0], 1)
V_max_inverse_color_Radio_IR = get_V_max_inverse(matrix_color_Radio_IR, (log_L_Radio, log_L_IR4, log_M), sens_color_Radio_IR_vector)

matrix_color_IR = np.concatenate((log_L_IR1.reshape(1, data_length), log_L_IR2.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_color_IR_vector = np.array([sensitivity_IR1, sensitivity_IR2, sensitivity_SDSS]).reshape(matrix_color_IR.shape[0], 1)
V_max_inverse_color_IR = get_V_max_inverse(matrix_color_IR, (log_L_IR1, log_L_IR2, log_M), sens_color_IR_vector)


#%% GRIDS Y FUNCIONES PARA EL CÁLCULO DE LOS KERNEL ESTIMATORS:

M_x = np.linspace(min(log_M[good_M]), max(log_M[good_M]), 250)
SFR_y = np.linspace(min(log_SFR[good_SFR]), max(log_SFR[good_SFR]), 250)

d_M = M_x[1] - M_x[0]
d_SFR = SFR_y[1] - SFR_y[0]

M_grid, SFR_grid = np.meshgrid(M_x, SFR_y)

def densidad(M_i, SFR_j, M_p, SFR_p, sigma_M_p, sigma_SFR_p):
    
    dangerous_terms = np.array([M_p, SFR_p, sigma_M_p, sigma_SFR_p])
    term1 = (M_i - M_p)/sigma_M_p
    term2 = (SFR_j - SFR_p)/sigma_SFR_p
    term = -term1**2 - term2**2
    term = term/2
    term = np.exp(term)
    term = term/(np.pi*sigma_M_p*sigma_SFR_p)
    
    if np.any(np.isnan(dangerous_terms)) == True:
        term = np.zeros((250, 250))
        
    return term

#%%  GETTING DENSITIES, AVERAGES AND DISPERSIONS FOR EACH LUMINOSITY AND COLOR:

def get_outputs(L, V_max_inverse):
    
    min_smoothing = 0.1
    total_dens = np.zeros((250, 250))
    weighted_dens = np.zeros((250, 250))
    weighted_dens_2 = np.zeros((250, 250))
    
    for p in range(data_length):
        L_p = L[p]
        V_max_inverse_p = V_max_inverse[p]
        if ((np.isnan(L_p) == False) and (np.isnan(V_max_inverse_p) == False)):
            
            densidad_p = densidad(M_i = M_grid, SFR_j = SFR_grid, M_p = log_M[p], SFR_p = log_SFR[p],
                                  sigma_M_p = min_smoothing + sigma_M[p], sigma_SFR_p = min_smoothing + sigma_SFR[p])
            
            w_p = densidad_p*V_max_inverse_p
            total_dens = total_dens + w_p
            weighted_dens = weighted_dens + L_p*w_p
            weighted_dens_2 = weighted_dens_2 + (L_p**2)*w_p
            
    L_mean = weighted_dens/total_dens
    L_mean_squared = weighted_dens_2/total_dens
    Dispersion = np.sqrt(L_mean_squared - L_mean**2)
    
    return total_dens, L_mean, Dispersion


def get_color_outputs(color, V_max_inverse):
    
    min_smoothing = 0.1
    total_dens = np.zeros((250, 250))
    weighted_dens = np.zeros((250, 250))
    weighted_dens_2 = np.zeros((250, 250))
    
    for p in range(data_length):
        color_p = color[p]
        V_max_inverse_p = V_max_inverse[p]
        if ((np.isnan(color_p) == False) and (np.isnan(V_max_inverse_p) == False)):
            
            densidad_p = densidad(M_i = M_grid, SFR_j = SFR_grid, M_p = log_M[p], SFR_p = log_SFR[p],
                                  sigma_M_p = min_smoothing + sigma_M[p], sigma_SFR_p = min_smoothing + sigma_SFR[p])
            
            w_p = densidad_p*V_max_inverse_p
            total_dens = total_dens + w_p
            weighted_dens = weighted_dens + color_p*w_p
            weighted_dens_2 = weighted_dens_2 + (color_p**2)*w_p
            
    color_mean = weighted_dens/total_dens
    color_mean_squared = weighted_dens_2/total_dens
    Dispersion = np.sqrt(color_mean_squared - color_mean**2)
    
    return total_dens, color_mean, Dispersion


dens_Soft, mean_Soft, disp_Soft = get_outputs(log_L_Soft, V_max_inverse_Soft)
dens_Hard, mean_Hard, disp_Hard = get_outputs(log_L_Hard, V_max_inverse_Hard)
dens_Radio, mean_Radio, disp_Radio = get_outputs(log_L_Radio, V_max_inverse_Radio)
dens_IR4, mean_IR4, disp_IR4 = get_outputs(log_L_IR4, V_max_inverse_IR4)
dens_IR1, mean_IR1, disp_IR1 = get_outputs(log_L_IR1, V_max_inverse_IR1)
dens_IR2, mean_IR2, disp_IR2 = get_outputs(log_L_IR2, V_max_inverse_IR2)

dens_color_X, mean_color_X, disp_color_X = get_color_outputs(color_X, V_max_inverse_color_X)
dens_color_Radio_IR, mean_color_Radio_IR, disp_color_Radio_IR = get_color_outputs(color_Radio_IR, V_max_inverse_color_Radio_IR)
dens_color_IR, mean_color_IR, disp_color_IR = get_color_outputs(color_IR, V_max_inverse_color_IR)

'''plt.contourf(M_grid, SFR_grid, mean_IR1, cmap='jet', extent=(9,11,0,1))
plt.xlim(left=9, right=11.5)
plt.ylim(bottom=0.15, top=1)
plt.colorbar()
plt.show()


plt.scatter(log_M, log_SFR, alpha=0.01)
plt.show()'''

def get_isocontour_levels(Density, dx, dy):
    
    mass = Density*dx*dy
    mass_vector = mass.flatten()
    dens_vector = Density.flatten()
    sorted_dens = np.argsort(dens_vector)
    density_sorted = dens_vector[sorted_dens]
    cumulative_mass = np.cumsum(mass_vector[sorted_dens])
    
    plt.plot(density_sorted, cumulative_mass)
    plt.xlabel('density')
    plt.ylabel('cumulative mass')
    plt.show()
    
    fraction_sorted = cumulative_mass/cumulative_mass[-1]
    fraction = np.interp(Density, density_sorted, fraction_sorted)
    fraction_inside = 1 - fraction
    
    plt.contourf(M_grid, SFR_grid, Density, cmap='jet', extent=(8.5,11.5,-2,1))
    plt.xlim(left=8.5, right=11.5)
    plt.xlabel('log10(M)')
    plt.ylim(bottom=-2, top=1)
    plt.ylabel('log10(SFR)')
    plt.colorbar()
    plt.show()
    
    plt.contourf(M_grid, SFR_grid, fraction, cmap='jet', extent=(8.5,11.5,-2,1))
    plt.xlim(left=8.5, right=11.5)
    plt.xlabel('log10(M)')
    plt.ylim(bottom=-2, top=1)
    plt.ylabel('log10(SFR)')
    plt.colorbar()
    plt.show()
    
    return fraction_inside
    

contours_Soft = get_isocontour_levels(dens_Soft, d_M, d_SFR)
contours_Hard = get_isocontour_levels(dens_Hard, d_M, d_SFR)
contours_Radio = get_isocontour_levels(dens_Radio, d_M, d_SFR)
contours_IR4 = get_isocontour_levels(dens_IR4, d_M, d_SFR)
contours_IR1 = get_isocontour_levels(dens_IR1, d_M, d_SFR)
contours_IR2 = get_isocontour_levels(dens_IR2, d_M, d_SFR)

contours_color_X = get_isocontour_levels(dens_color_X, d_M, d_SFR)
contours_color_Radio_IR = get_isocontour_levels(dens_color_Radio_IR, d_M, d_SFR)
contours_color_IR = get_isocontour_levels(dens_color_IR, d_M, d_SFR)


#%% KERNEL ESTIMATORS OF <L> AND <COLOR> AND THEIR RESPECTIVE DISPERSIONS:

    
def get_plots(L, contours, title):
    
    plt.contourf(M_grid, SFR_grid, L, cmap='jet', extent=(8.5,11.5,-2,1))
    plt.colorbar()
    plt.title(title)
    cont = plt.contour(M_grid, SFR_grid, contours,  levels = [0.10, 0.50, 0.90], colors = 'black')
    plt.xlim(left=8.5, right=11.5)
    plt.xlabel('log10(M)')
    plt.ylim(bottom=-2, top=1)
    plt.ylabel('log10(SFR)')
    plt.clabel(cont, inline=True, fontsize = 8)
    plt.show()
    plotted = 'Plot made!'
    
    return plotted


plot_Soft = get_plots(mean_Soft, contours_Soft, '<log L> Soft X-Rays')
plot_Soft_disp = get_plots(disp_Soft, contours_Soft, 'Dispersion Soft X-Rays')

plot_Hard = get_plots(mean_Hard, contours_Hard, '<log L> Hard X-Rays')
plot_Hard_disp = get_plots(disp_Hard, contours_Hard, 'Dispersion Hard X-Rays')

plot_Radio = get_plots(mean_Radio, contours_Radio, '<log L> Radio')
plot_Radio_disp = get_plots(disp_Radio, contours_Radio, 'Dispersion Radio')

plot_IR4 = get_plots(mean_IR4, contours_IR4, '<log L> IR (W4 band)')
plot_IR4_disp = get_plots(disp_IR4, contours_IR4, 'Dispersion IR (W4 band)')

plot_IR1 = get_plots(mean_IR1, contours_IR1, '<log L> IR (W1 band)')
plot_IR1_disp = get_plots(disp_IR1, contours_IR1, 'Dispersion IR (W1 band)')

plot_IR2 = get_plots(mean_IR2, contours_IR2, '<log L> IR (W2 band)')
plot_IR2_disp = get_plots(disp_IR2, contours_IR2, 'Dispersion IR (W2 band)')



plot_color_IR = get_plots(mean_color_IR, contours_color_IR, '<IR color>')
plot_color_IR_disp = get_plots(disp_color_IR, contours_color_IR, 'Dispersion IR color')

plot_color_Radio_IR = get_plots(mean_color_Radio_IR, contours_color_Radio_IR, '<Radio-IR (W4 band) color>')
plot_color_Radio_IR_disp = get_plots(disp_color_Radio_IR, contours_color_Radio_IR, 'Dispersion Radio-IR (W4 band) color')

plot_color_X = get_plots(mean_color_X, contours_color_X, '<X-Rays color>')
plot_color_X_disp = get_plots(disp_color_X, contours_color_X, 'Dispersion X-Rays color')


#%% CÁLCULO DE LAS LUMINOSIDADES Y COLORES MEDIOS Y DE SUS EXCESOS:
    
    
log_L_Soft_mean = np.zeros(data_length)
log_L_Hard_mean = np.zeros(data_length)
log_L_Radio_mean = np.zeros(data_length)
log_L_IR4_mean = np.zeros(data_length)
log_L_IR1_mean = np.zeros(data_length)
log_L_IR2_mean = np.zeros(data_length)
color_IR_mean = np.zeros(data_length)
color_Radio_IR_mean = np.zeros(data_length)
color_X_mean = np.zeros(data_length)

binsize_M = d_M/2
binsize_SFR = d_SFR/2

def L_mean_vector(L, mean_L_array, mean_L):
    
    for i in range(data_length):
        if ((np.isnan(L[i])==False) and (np.isnan(log_M[i])==False) and (np.isnan(log_SFR[i])==False)):
            row = 0
            column = 0
            for j in range(len(M_x)):
                if np.abs(log_M[i] - M_x[j]) < binsize_M:
                    column = j
                    break
                
            for k in range(len(SFR_y)):
                if np.abs(log_SFR[i] - SFR_y[k]) < binsize_SFR:
                    row = k
                    break
            
            if ((row != 0) and (column != 0)):
                mean_L_array[i] = mean_L[row, column]
                
            else:
                mean_L_array[i] = np.nan
                
        else:
            mean_L_array[i] = np.nan
            
    return mean_L_array


L_mean_Soft = L_mean_vector(L = log_L_Soft, mean_L_array = log_L_Soft_mean, mean_L = mean_Soft)
L_mean_Hard = L_mean_vector(L = log_L_Hard, mean_L_array = log_L_Hard_mean, mean_L = mean_Hard)
L_mean_Radio = L_mean_vector(L = log_L_Radio, mean_L_array = log_L_Radio_mean, mean_L = mean_Radio)
L_mean_IR4 = L_mean_vector(L = log_L_IR4, mean_L_array = log_L_IR4_mean, mean_L = mean_IR4)
L_mean_IR1 = L_mean_vector(L = log_L_IR1, mean_L_array = log_L_IR1_mean, mean_L = mean_IR1)
L_mean_IR2 = L_mean_vector(L = log_L_IR2, mean_L_array = log_L_IR2_mean, mean_L = mean_IR2)

color_IR_mean = L_mean_vector(L = color_IR, mean_L_array = color_IR_mean, mean_L = mean_color_IR)
color_Radio_IR_mean = L_mean_vector(L = color_Radio_IR, mean_L_array = color_Radio_IR_mean, mean_L = mean_color_Radio_IR)
color_X_mean = L_mean_vector(L = color_X, mean_L_array = color_X_mean, mean_L = mean_color_X)


# Obtaining excessses:
    
L_excess_Soft = log_L_Soft - L_mean_Soft
L_excess_Hard = log_L_Hard - L_mean_Hard
L_excess_Radio = log_L_Radio - L_mean_Radio
L_excess_IR4 = log_L_IR4 - L_mean_IR4
L_excess_IR1 = log_L_IR1 - L_mean_IR1
L_excess_IR2 = log_L_IR2 - L_mean_IR2
color_IR_excess = color_IR - color_IR_mean
color_Radio_IR_excess = color_Radio_IR - color_Radio_IR_mean
color_X_excess = color_X - color_X_mean


#%% ACTIVENESS FUNCTION AND PLOTS


def activeness(lum, xlabel, V_max_inverse):
    
    good_length = len(lum[np.isfinite(lum)])
    sorted_indices = np.argsort(lum)[0:good_length]
    sorted_lum = lum[sorted_indices]
    sorted_density = V_max_inverse[sorted_indices]
    cum_dens = np.cumsum(sorted_density)
    cum_dens = cum_dens/cum_dens[-1]
    
    total_dens = np.max(cum_dens)
    mid_dens = total_dens/2
    pivot_point = np.interp(mid_dens, cum_dens, sorted_lum)
    print(pivot_point)
    rho = total_dens
    z = np.linspace(np.min(sorted_lum), np.max(sorted_lum), 300)
    interpolation = np.interp(z, sorted_lum, cum_dens, left = 0, right = 1)
    rho_greater = rho - interpolation
    rho_normal = np.zeros(len(rho_greater))
    
    for i in range(len(z)):
        if z[i] <= pivot_point:
            rho_symmetric = 2*np.interp(pivot_point, sorted_lum, cum_dens, left = 0, right = 1) - np.interp(z[i], sorted_lum, cum_dens, left = 0, right=1)
            
        else:
            rho_symmetric = np.interp((2*pivot_point - z[i]), sorted_lum, cum_dens, left = 0, right = 1)
        
        rho_normal[i] = np.minimum(rho_symmetric, rho_greater[i])
        
        
    rho_active = rho_greater - rho_normal
    f_active = rho_active/rho_greater
    position = np.where(f_active >= 0.5)[0]
    position = np.min(position)
    z_active = z[position]
    
    plt.figure()
    plt.yscale('linear') 
    plt.plot(z, rho_greater, color = 'black')
    plt.plot(z, rho_normal, color = 'blue')
    plt.plot(z, rho_active, color = 'red')
    plt.axvline(x = z_active, linestyle='--', color='black')
    plt.xlabel(xlabel)
    plt.ylabel('normalized density')
    plt.show()
    
    return z_active


#%% COMPUTING CONVERSION PARAMETERS BETWEEN OUR COLORS AND THOSE IN LITERATURE:
    
    
W1_W2_norm = np.log10((171*3.4*1e-6)/(309*4.6*1e-6)) # cociente entre flujos cero de WISE
nuradio_nu24 = np.log10((1.4*1e9)/((3*1e8)/(22.1*1e-6))) # cociente de frecuencias (1.4 GHz/22.1 micras)

def ours_from_lit_IR(x):
    new_x = 0.4*x + W1_W2_norm
    return new_x

def lit_from_ours_IR(x):
    new_x = (x - W1_W2_norm)/0.4
    return new_x

def ours_from_lit_Radio(x):
    new_x = x + nuradio_nu24
    return new_x

def lit_from_ours_Radio(x):
    new_x = x - nuradio_nu24
    return new_x

#%% COMPUTING THRESHOLDS FOR EACH LUMINOSITY AND COLOR:
    
    
threshold_color_IR_Stern = 0.4*0.8 + W1_W2_norm # Stern threshold ---> W1-W2 > 0.8
threshold_color_IR_Assef = 0.4*0.5 + W1_W2_norm # Assef threshold ---> W1-W2 > 0.5
threshold_color_Radio_IR_Ibar = 0.23 + nuradio_nu24 # Ibar threshold ---> q_24 < -0.23
threshold_Radio_loud = 40 # L > 10^40 erg/s ---> Radio-loud AGN
threshold_Hard_lit = 42 # L > 10^42 erg/s ---> AGN
    
threshold_Hard_excess = activeness(L_excess_Hard, 'Luminosity excess (Hard X-Rays)', V_max_inverse_Hard)
threshold_Radio_excess = activeness(L_excess_Radio, 'Luminosity excess (Radio)', V_max_inverse_Radio)
threshold_IR2_excess = activeness(L_excess_IR2, 'Luminosity excess (IR, W2 band)', V_max_inverse_IR2)

threshold_color_X_excess = activeness(color_X_excess, 'color excess (X-Rays)', V_max_inverse_color_X)
threshold_color_Radio_IR_excess = activeness(color_Radio_IR_excess, 'color excess (Radio-IR W4 band)', V_max_inverse_color_Radio_IR)
threshold_color_IR_excess = activeness(color_IR_excess, 'color excess (IR)', V_max_inverse_color_IR)


#%% GALAXY CLASSIFICATION:
  
def galaxy_classification(threshold_L, threshold_color, lum1, lum2, color):
    
    good = np.where((np.isfinite(lum1)==True) & (np.isfinite(color)==True))[0]
    active = np.where((lum1 > threshold_L) & (color > threshold_color))[0]
    lum_ex = np.where((lum1 > threshold_L) & (color < threshold_color))[0]
    col_ex = np.where((lum1 < threshold_L) & (color > threshold_color))[0]
    normal = np.where((lum1 < threshold_L) & (color < threshold_color))[0]
    non_det_lum1 = np.where((np.isnan(lum1)) & (np.isfinite(lum2)))[0] #cyan
    non_det_lum2_normal = np.where((lum1 < threshold_L) & (np.isnan(lum2)))[0] # grey
    non_det_lum2_exc = np.where((lum1 > threshold_L) & (np.isnan(lum2)))[0] # black
    
    return good, active, lum_ex, col_ex, normal, non_det_lum1, non_det_lum2_normal, non_det_lum2_exc


IR_good, IR_active, IR_lum_ex, IR_col_ex, IR_normal, IR2_non_det, IR1_non_det_normal, IR1_non_det_lum_ex \
    = galaxy_classification(threshold_IR2_excess, threshold_color_IR_excess, L_excess_IR2, L_excess_IR1, color_IR_excess)
    
    
Radio_good, Radio_active, Radio_lum_ex, Radio_col_ex, Radio_normal, Radio_non_det, IR4_non_det_normal, IR4_non_det_lum_ex \
    = galaxy_classification(threshold_Radio_excess, threshold_color_Radio_IR_excess, L_excess_Radio, L_excess_IR4, color_Radio_IR_excess)
       
X_good, X_active, X_lum_ex, X_col_ex, X_normal, Hard_non_det, Soft_non_det_normal, Soft_non_det_lum_ex \
    = galaxy_classification(threshold_Hard_excess, threshold_color_X_excess, L_excess_Hard, L_excess_Soft, color_X_excess)
    
    
def weighted_hists(L, V_max_inverse, good, active, lum_ex, col_ex, normal, alpha_value, bin_number):
    plt.figure()
    n_bins = np.linspace(np.nanmin(L),np.nanmax(L), bin_number)
    L1 = L[good]
    V_max_inverse1 = V_max_inverse[good]
    plt.hist(L1, bins = n_bins, weights = V_max_inverse1, log = True, color = 'red', alpha = alpha_value)
    
    set2 = np.setdiff1d(good, active)
    L2 = L[set2]
    V_max_inverse2 = V_max_inverse[set2]
    plt.hist(L2, bins = n_bins, weights = V_max_inverse2, log = True, color = 'gold', alpha = alpha_value)
    
    set3 = np.setdiff1d(good, np.concatenate((active, col_ex)))
    L3 = L[set3]
    V_max_inverse3 = V_max_inverse[set3]
    plt.hist(L3, bins = n_bins, weights = V_max_inverse3, log = True, color = 'green', alpha = alpha_value)
    
    set4 = np.setdiff1d(good, np.concatenate((active, col_ex, lum_ex)))
    L4 = L[set4]
    V_max_inverse4 = V_max_inverse[set4]
    plt.hist(L4, bins = n_bins, weights = V_max_inverse4, log = True, color = 'blue', alpha = alpha_value)
    
    plt.show()
    printed = 'Printed!'

    return printed


IR2_weighted = weighted_hists(L_excess_IR2, V_max_inverse_IR2, IR_good, IR_active, IR_lum_ex, 
                   IR_col_ex, IR_normal, alpha_value = 1, bin_number = 51)

Radio_weighted = weighted_hists(L_excess_Radio, V_max_inverse_Radio, Radio_good, Radio_active, Radio_lum_ex, 
                   Radio_col_ex, Radio_normal, alpha_value = 1, bin_number = 31)

Hard_weighted = weighted_hists(L_excess_Hard, V_max_inverse_Hard, X_good, X_active, X_lum_ex, 
                   X_col_ex, X_normal, alpha_value = 1, bin_number = 21)


color_IR_weighted = weighted_hists(color_IR_excess, V_max_inverse_color_IR, IR_good, IR_active, IR_lum_ex, 
                   IR_col_ex, IR_normal, alpha_value = 1, bin_number = 51)

color_Radio_IR_weighted = weighted_hists(color_Radio_IR_excess, V_max_inverse_color_Radio_IR, Radio_good, Radio_active, Radio_lum_ex, 
                   Radio_col_ex, Radio_normal, alpha_value = 1, bin_number = 31)

color_X_weighted = weighted_hists(color_X_excess, V_max_inverse_color_X, X_good, X_active, X_lum_ex, 
                   X_col_ex, X_normal, alpha_value = 1, bin_number = 21)

#%% LITERATURE COMPARISON

def literature_comp(L, active, col_ex, lum_ex, normal, V_max_inverse, xlabel, mode):
    
    if mode == 'lum':
        active_gal = np.union1d(active, lum_ex)
        
    if mode == 'col':
        active_gal = np.union1d(active, col_ex)
        
    good_total = len(L[np.isfinite(L)]) # number of normal galaxies
    sorted_indices = np.argsort(L)[0:good_total] # indices sorted by L
    sorted_L = L[sorted_indices] # L sorted for normal galaxies
    sorted_density = V_max_inverse[sorted_indices]
    cum_dens = np.cumsum(sorted_density)
    total_dens = np.max(cum_dens)
    rho = total_dens
    z = np.linspace(np.min(sorted_L), np.max(sorted_L), 300)
    interpolation = np.interp(z, sorted_L, cum_dens)
    rho_greater = rho - interpolation

    
    good_active = len(L[active_gal]) # number of normal galaxies
    sorted_indices_active = np.argsort(L[active_gal])[0:good_active] # indices sorted by L
    sorted_L_active = L[active_gal][sorted_indices_active] # L sorted for normal galaxies
    sorted_density_active = V_max_inverse[active_gal][sorted_indices_active]
    cum_dens_active = np.cumsum(sorted_density_active)
    total_dens_active = np.max(cum_dens_active)
    rho_active = total_dens_active
    interpolation_active = np.interp(z, sorted_L_active, cum_dens_active)
    rho_greater_active = rho_active - interpolation_active

    f_active = rho_greater_active/rho_greater
    position_50 = np.where(f_active >= 0.5)[0]
    position_50 = np.min(position_50)
    z_active_50 = z[position_50]
    
    position_90 = np.where(f_active >= 0.9)[0]
    if len(position_90)>0:
        
        position_90 = np.min(position_90)
        z_active_90 = z[position_90]
    else:
        z_active_90 = np.nan
        
    plt.figure()
    plt.yscale('log') 
    plt.plot(z, rho_greater, color = 'black')
    plt.plot(z, rho_greater_active, color = 'red')
    plt.xlabel(xlabel)
    plt.ylabel('normalized density')
    plt.grid()
    plt.show()
    return z_active_50, z_active_90




color_IR[np.isnan(V_max_inverse_color_IR)==True]=np.nan
l_ir_col_50, l_ir_col_90 = literature_comp(color_IR, IR_active, IR_col_ex, IR_lum_ex, IR_normal,
                                           V_max_inverse_color_IR, 'color IR', mode='col')

color_Radio_IR[np.isnan(V_max_inverse_color_Radio_IR)==True]=np.nan
l_radio_col_50, l_radio_col_90 = literature_comp(color_Radio_IR, Radio_active, Radio_col_ex, Radio_lum_ex, Radio_normal,
                              V_max_inverse_color_Radio_IR, 'color Radio-IR', mode='col')

log_L_Radio[np.isnan(V_max_inverse_Radio)==True]=np.nan
l_radio_lum_50, l_radio_lum_90 = literature_comp(log_L_Radio, Radio_active, Radio_col_ex, Radio_lum_ex, Radio_normal,
                              V_max_inverse_Radio, 'log10(L Radio)', mode='lum')

log_L_Hard[np.isnan(V_max_inverse_Hard)==True]=np.nan
l_hard_lum_50, l_hard_lum_90 = literature_comp(log_L_Hard, X_active, X_col_ex, X_lum_ex, X_normal,
                             V_max_inverse_Hard, 'log10(L Hard X-Rays)', mode='lum')
#%% SCATTER PLOTS:
    
def plotting(L, color, active, lum, col, normal, alpha_value, left, right, bottom, top, 
             right_axis, f1, f2, xlab, ylab, ylab2, is_excess, v_ex, h_ex, v_ex2, h_ex2,
             v_ex3, h_ex3):
    
    fig, ax = plt.subplots()
    
    ax.scatter(L[normal], color[normal], color = 'blue', alpha = alpha_value)
    ax.scatter(L[lum], color[lum], color ='green', alpha = alpha_value)
    ax.scatter(L[col], color[col], color = 'gold', alpha = alpha_value)
    ax.scatter(L[active], color[active], color = 'red', alpha = alpha_value)
    
    plt.xlim(left = left, right = right)
    plt.ylim(bottom = bottom, top = top)
    
    if is_excess == True:
        
        #Literature
        plt.axvline(x = v_ex, linestyle='--', color='black')
        plt.axhline(y = h_ex, linestyle='--', color='black')
        
        # This work f_active = 50%
        plt.axvline(x = v_ex2, linestyle='-', color='cyan')
        plt.axhline(y = h_ex2, linestyle='-', color='cyan')
        
        # This work f_active = 90%
        plt.axvline(x = v_ex3, linestyle='-', color='blue')
        plt.axhline(y = h_ex3, linestyle='-', color='blue')
    
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    
    if right_axis == True:
        
        secax = ax.secondary_yaxis('right', functions=(f1, f2))
        secax.set_ylabel(ylab2)
    
    printed= 'printed!'
    
    return printed

# Plotting scattered excesses:

scatter_IR_excess = plotting(L_excess_IR2, color_IR_excess, IR_active, IR_lum_ex,
                             IR_col_ex, IR_normal, alpha_value = .1, left = -1, 
                             right = 2, bottom = -0.25, top = 0.75, right_axis = False, 
                             f1 = None, f2 = None,
                             xlab = 'IR (W2 band) luminosity excess', ylab = 'IR color excess', 
                             ylab2 = '', is_excess = True, v_ex = threshold_IR2_excess,
                             h_ex = threshold_color_IR_excess, v_ex2 = np.inf, h_ex2 = np.inf,
                             v_ex3 = np.inf, h_ex3 = np.inf)

scatter_Radio_excess = plotting(L_excess_Radio, color_Radio_IR_excess, Radio_active, Radio_lum_ex,
                                Radio_col_ex, Radio_normal, alpha_value = .1, left = -1.5,
                                right = 3.5, bottom = -2, top = 4, right_axis = False, 
                                f1 = None, f2 = None,
                                xlab = 'Radio luminosity excess', ylab = 'Radio-IR (W4 band) color excess',
                                ylab2 = '', is_excess = True, v_ex = threshold_Radio_excess, 
                                h_ex = threshold_color_Radio_IR_excess, v_ex2 = np.inf, h_ex2 = np.inf,
                                v_ex3 = np.inf, h_ex3 = np.inf)

scatter_X_excess = plotting(L_excess_Hard, color_X_excess, X_active, X_lum_ex,
                            X_col_ex, X_normal, alpha_value = .3, left = -1.5,
                            right = 3.5, bottom = -2, top = 4, right_axis = False,
                            f1 = None, f2 = None,
                            xlab = 'Hard luminosity excess', ylab = 'X-Rays color excess', 
                            ylab2 = '', is_excess = True, v_ex = threshold_Hard_excess,
                            h_ex = threshold_color_X_excess, v_ex2 = np.inf, h_ex2 = np.inf,
                            v_ex3 = np.inf, h_ex3 = np.inf)

#%% Plotting scattered luminosities and colors:
    
scatter_IR = plotting(log_L_IR2, color_IR, IR_active, IR_lum_ex,
                      IR_col_ex, IR_normal, alpha_value = .1, left = 40, 
                      right = 44.5, bottom = -.6 , top = .5, right_axis = True, 
                      f1 = lit_from_ours_IR, f2 = ours_from_lit_IR,
                      xlab = 'log10(IR (W2 band) luminosity)', ylab = 'IR color', 
                      ylab2 = 'W1-W2', is_excess = True, v_ex = 0,
                      h_ex = threshold_color_IR_Stern,
                      v_ex2 = np.inf, h_ex2 = l_ir_col_50,
                      v_ex3 = np.inf, h_ex3 = l_ir_col_90)

scatter_Radio = plotting(log_L_Radio, color_Radio_IR, Radio_active, Radio_lum_ex,
                         Radio_col_ex, Radio_normal, alpha_value = .1, left = 36.5,
                         right = 42, bottom = -7, top = -1, right_axis = True, 
                         f1 = lit_from_ours_Radio, f2 = ours_from_lit_Radio,
                         xlab = 'log10(Radio luminosity)', ylab = 'Radio-IR (W4 band) color',
                         ylab2 = 'q24', is_excess = True, v_ex = threshold_Radio_loud,
                         h_ex = threshold_color_Radio_IR_Ibar,
                         v_ex2 = l_radio_lum_50, h_ex2 = l_radio_col_50,
                         v_ex3 = l_radio_lum_90, h_ex3 = l_radio_col_90)

scatter_X = plotting(log_L_Hard, color_X, X_active, X_lum_ex,
                     X_col_ex, X_normal, alpha_value = .3, left = 39.5,
                     right = 44, bottom = -1.25, top = 2.5, right_axis = False,
                     f1 = None, f2 = None,
                     xlab = 'log10(Hard X-Rays luminosity)', ylab = 'X-Rays color', 
                     ylab2 = '', is_excess = True, v_ex = threshold_Hard_lit, h_ex = np.inf,
                     v_ex2 = l_hard_lum_50, h_ex2 = np.inf,
                     v_ex3 = l_hard_lum_90, h_ex3 = np.inf)

#%% STELLAR MASS vs SFR

# IR:
plt.scatter(log_M[IR_active], log_SFR[IR_active], color = 'red', alpha=0.5)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[IR_lum_ex], log_SFR[IR_lum_ex], color = 'green', alpha=0.5)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[IR_col_ex], log_SFR[IR_col_ex], color = 'gold', alpha=0.2)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[IR_normal], log_SFR[IR_normal], color = 'blue', alpha=0.005)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

# Radio:
    
plt.scatter(log_M[Radio_active], log_SFR[Radio_active], color = 'red', alpha=0.2)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[Radio_lum_ex], log_SFR[Radio_lum_ex], color = 'green', alpha=0.2)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[Radio_col_ex], log_SFR[Radio_col_ex], color = 'gold', alpha=0.3)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[Radio_normal], log_SFR[Radio_normal], color = 'blue', alpha=0.1)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

# X-Rays:
    
plt.scatter(log_M[X_active], log_SFR[X_active], color = 'red', alpha=0.5)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[X_lum_ex], log_SFR[X_lum_ex], color = 'green', alpha=0.5)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[X_col_ex], log_SFR[X_col_ex], color = 'gold', alpha=0.5)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

plt.scatter(log_M[X_normal], log_SFR[X_normal], color = 'blue', alpha=0.5)
plt.xlabel('log10(M*)')
plt.ylabel('log10(SFR)')

