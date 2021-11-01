# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:33:24 2021

@author: alvar
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec


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


def get_V_max_inverse(L, sensitivity):
    
    z_max = 10**((L - sensitivity)/2) # redshift computed from the equation above
    z_max = z_max.min(0) # min of each column for the above matrix
    bad = np.isnan(z_max) # finding nans
    z_max = np.where(z_max > 0.01, z_max, 0.01) # if z < 0.01: z = 0.01, else z = z
    z_max = np.where(z_max < 0.07, z_max, 0.07) # if z > 0.07: z = 0.07, else z = z
    z_max = np.where(bad == False, z_max, np.nan) # nans still being nans

    return 1/(z_max**3) # rho, or inverse max volume


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
V_max_inverse_Soft = get_V_max_inverse(matrix_Soft, sens_Soft_vector)

matrix_Hard = np.concatenate((log_L_Hard.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_Hard_vector = np.array([sensitivity_Hard, sensitivity_SDSS]).reshape(matrix_Hard.shape[0], 1)
V_max_inverse_Hard = get_V_max_inverse(matrix_Hard, sens_Hard_vector)

matrix_Radio = np.concatenate((log_L_Radio.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_Radio_vector = np.array([sensitivity_Radio, sensitivity_SDSS]).reshape(matrix_Radio.shape[0], 1)
V_max_inverse_Radio = get_V_max_inverse(matrix_Radio, sens_Radio_vector)

matrix_IR4 = np.concatenate((log_L_IR4.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_IR4_vector = np.array([sensitivity_IR4, sensitivity_SDSS]).reshape(matrix_IR4.shape[0], 1)
V_max_inverse_IR4 = get_V_max_inverse(matrix_IR4, sens_IR4_vector)

matrix_IR1 = np.concatenate((log_L_IR1.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_IR1_vector = np.array([sensitivity_IR1, sensitivity_SDSS]).reshape(matrix_IR1.shape[0], 1)
V_max_inverse_IR1 = get_V_max_inverse(matrix_IR1, sens_IR1_vector)

matrix_IR2 = np.concatenate((log_L_IR2.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_IR2_vector = np.array([sensitivity_IR2, sensitivity_SDSS]).reshape(matrix_IR2.shape[0], 1)
V_max_inverse_IR2 = get_V_max_inverse(matrix_IR2, sens_IR2_vector)

# Getting V max inverse for colours:

matrix_color_X = np.concatenate((log_L_Soft.reshape(1, data_length), log_L_Hard.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_color_X_vector = np.array([sensitivity_Soft, sensitivity_Hard, sensitivity_SDSS]).reshape(matrix_color_X.shape[0], 1)
V_max_inverse_color_X = get_V_max_inverse(matrix_color_X, sens_color_X_vector)

matrix_color_Radio_IR = np.concatenate((log_L_IR4.reshape(1, data_length), log_L_Radio.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_color_Radio_IR_vector = np.array([sensitivity_IR4, sensitivity_Radio, sensitivity_SDSS]).reshape(matrix_color_Radio_IR.shape[0], 1)
V_max_inverse_color_Radio_IR = get_V_max_inverse(matrix_color_Radio_IR, sens_color_Radio_IR_vector)

matrix_color_IR = np.concatenate((log_L_IR1.reshape(1, data_length), log_L_IR2.reshape(1, data_length), log_M.reshape(1, data_length)), axis=0)
sens_color_IR_vector = np.array([sensitivity_IR1, sensitivity_IR2, sensitivity_SDSS]).reshape(matrix_color_IR.shape[0], 1)
V_max_inverse_color_IR = get_V_max_inverse(matrix_color_IR, sens_color_IR_vector)


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


def activeness(lum, xlabel, colour, V_max_inverse):
    
    good_length = len(lum[np.isfinite(lum)])
    sorted_indices = np.argsort(lum)[0:good_length]
    sorted_lum = lum[sorted_indices]
    sorted_density = V_max_inverse[sorted_indices]
    cum_dens = np.cumsum(sorted_density)
    cum_dens = cum_dens/cum_dens[-1]
    
    total_dens = np.max(cum_dens)
    mid_dens = total_dens/2
    pivot_point = np.interp(mid_dens, cum_dens, sorted_lum)
    
    rho = total_dens
    z = np.linspace(np.min(sorted_lum), np.max(sorted_lum), 300)
    interpolation = np.interp(z, sorted_lum, cum_dens)
    rho_greater = rho - interpolation
    rho_normal = np.zeros(len(rho_greater))
    
    for i in range(len(z)):
        if z[i] <= pivot_point:
            rho_symmetric = 2*np.interp(pivot_point, sorted_lum, cum_dens) - np.interp(z[i], sorted_lum, cum_dens)
            
        else:
            rho_symmetric = np.interp((2*pivot_point - z[i]), sorted_lum, cum_dens)
        
        rho_normal[i] = np.minimum(rho_symmetric, rho_greater[i])
        
        
    rho_active = rho_greater - rho_normal
    f_active = rho_active/rho_greater
    position = np.where(f_active >= 0.5)
    position = np.min(position)
    z_active = f_active[position]
    
    plt.yscale('log') 
    plt.plot(z, rho_greater, color = 'black')
    plt.plot(z, rho_normal, color = 'blue')
    plt.plot(z, rho_active, color = 'red')
    plt.axvline(x = z_active, linestyle='--', color='black')
    
    return z_active

activeness(L_excess_Radio, 'Luminosity excess (Radio)', 'green', V_max_inverse_Radio)
activeness(L_excess_Hard, 'Luminosity excess (Hard X-Rays)', 'blue', V_max_inverse_Hard)