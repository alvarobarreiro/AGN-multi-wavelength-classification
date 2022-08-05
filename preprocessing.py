# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:33:24 2021

@author: alvar
"""

import numpy as np
import pandas as pd

#%% Sensitivities

# Sensitivity (i.e. flux limit above which we have data for each band): it's been visually estimated from
# plot(redshift, log_L - 2*log10(redshift)) for each frequency

sensitivity_Radio = 40.85
sensitivity_IR4 = 44.8
sensitivity_IR1 = 44.6
sensitivity_IR2 = 44.2
sensitivity_Soft = 42.8
sensitivity_Hard = 43.5
sensitivity_SDSS = 12 # plot(redshift, log_M - 2*log10(redshift)) in this case

#%% Reading cross-match catalog

sdss_xmatch = pd.read_csv('final_xmatch_norm.csv', keep_default_na=True)

#%% Preprocessing class
    
class preprocessing():
    def __init__(self, data, band, sensitivity):
        
        # Main variables
        self.log_M = data['lgm_tot_p50'].values
        self.log_M_16 = data['lgm_tot_p16'].values
        self.log_M_84 = data['lgm_tot_p84'].values
        
        self.log_SFR = data['sfr_tot_p50'].values
        self.log_SFR_16 = data['sfr_tot_p16'].values
        self.log_SFR_84 = data['sfr_tot_p84'].values
        
        self.luminosity = np.log10(data[band+'_luminosity_(erg/s)'].values)
        
        self.flux_err = data[band+'_flux_error_(erg/s/cm2)'].values
        self.signal_to_noise = data[band+'_flux_(erg/s/cm2)'].values \
             / data[band+'_flux_error_(erg/s/cm2)'].values
        
        self.sensitivity = sensitivity


    def missing_M_SFR(self, missing_value = -9999):
       
        log_L = self.luminosity
        
        M_SFR_array = np.array([self.log_M, self.log_SFR, 
                                self.log_M_16, self.log_M_84, 
                                self.log_SFR_16, self.log_SFR_84])
    
        missing_data = np.any(M_SFR_array == missing_value, axis=0)
        
        for i in range(len(M_SFR_array)):
            M_SFR_array[i][missing_data] = np.nan
        
        log_L[missing_data] = np.nan
        
        log_M, log_SFR, log_M_16, log_M_84, log_SFR_16, log_SFR_84 = M_SFR_array
        sigma_M = (log_M_84 - log_M_16) / 2
        sigma_SFR = (log_SFR_84 - log_SFR_16) / 2
        
        return log_L, log_M, log_SFR, sigma_M, sigma_SFR
    
    
    def extreme_M_SFR(self, log_L, log_M, log_SFR, sigma_M, sigma_SFR,
                      lower_M=8.5, upper_M=11.5, lower_SFR=-2, upper_SFR=1):
        
        extreme = np.where((log_M < lower_M) & (log_M > upper_M) & \
                           (log_SFR < lower_SFR) & (log_SFR > upper_SFR))
        
        log_L[extreme] = np.nan
        log_M[extreme] = np.nan
        log_SFR[extreme] = np.nan
        sigma_M[extreme] = np.nan
        sigma_SFR[extreme] = np.nan
        
        return log_L, log_M, log_SFR, sigma_M, sigma_SFR 
    
    
    def signal_to_noise_filter(self, log_L, limit=2):
        
        SN = self.signal_to_noise
        #SN[np.isnan(log_L)] = np.nan
        
        # We reject data without a proper estimate of the error or SN:
        non_valid = np.where(np.isfinite(SN) == False)[0]   
        cut_off = np.where(SN < limit)[0] # SN threshold
        rejected = np.unique(np.concatenate((non_valid, cut_off)))
        log_L[rejected] = np.nan
        
        return log_L
    
    def flux_cut_off(self, log_L, SN0=2, const=57.34):
        
        sensitivity = self.sensitivity
        flux_err = self.flux_err
        flux_err[np.isnan(log_L)] = np.nan

        SN0_dF = np.log10(SN0*flux_err) + const #cut-off due to signal-to-noise
        # const --> to proper units (maybe using astropy?)
        F0 = sensitivity*np.ones_like(SN0_dF) # cut-off due to sensitivity
        F_cut_off = np.array([SN0_dF, F0]).max(0) # F > max(F0, SN0*dF)

        return F_cut_off
    
    
    def get_V_max_inverse(self, log_L, log_M, F_cut_off, sensitivity_SDSS=12):

        z_max_lum = 10**((log_L - F_cut_off)/2)
        z_max_mass = 10**((log_M - sensitivity_SDSS)/2)
        # z_max = z_max.min(0).clip(0.01,0.07) 
        z_max = np.array([z_max_lum, z_max_mass])
        z_max = z_max.min(0) # min of each column for the above matrix
        faint_obj = np.where(z_max < 0.015)[0]
        z_max[faint_obj] = np.nan
        
        log_L[faint_obj] = np.nan # we exclude every object with z_max < 0.015

        # z_max = np.where(z_max > 0.01, z_max, 0.01) # if z < 0.01: z = 0.01, else z = z
        z_max = np.where(z_max < 0.07, z_max, 0.07) # if z > 0.07: z = 0.07, else z = z
        V_max = z_max**3 - 1e-9
        V_max[np.isnan(log_L)] = np.nan
        
        return 1/V_max, log_L # inverse max volume
    
    
    def forward_pass(self):
        
        log_L, log_M, log_SFR, sigma_M, sigma_SFR = self.missing_M_SFR()
        log_L, log_M, log_SFR, sigma_M, sigma_SFR = self.extreme_M_SFR(log_L, log_M, log_SFR, sigma_M, sigma_SFR)
        log_L = self.signal_to_noise_filter(log_L)
        F_cut_off = self.flux_cut_off(log_L)
        Vmax_inv, log_L = self.get_V_max_inverse(log_L, log_M, F_cut_off)
        
        return log_L, Vmax_inv, log_M, log_SFR, sigma_M, sigma_SFR
        
    

#%% Performing preprocessing on each band

prep = preprocessing(sdss_xmatch, 'IR1', sensitivity_IR1)
log_L_IR1, Vmax_inv_IR1, log_M, log_SFR, sigma_M, sigma_SFR = prep.forward_pass()

prep = preprocessing(sdss_xmatch, 'IR2', sensitivity_IR2)
log_L_IR2, Vmax_inv_IR2, _, _, _, _ = prep.forward_pass()

prep = preprocessing(sdss_xmatch, 'IR4', sensitivity_IR4)
log_L_IR4, Vmax_inv_IR4, _, _, _, _ = prep.forward_pass()

prep = preprocessing(sdss_xmatch, 'Radio', sensitivity_Radio)
log_L_Radio, Vmax_inv_Radio, _, _, _, _ = prep.forward_pass()

prep = preprocessing(sdss_xmatch, 'Soft', sensitivity_Soft)
log_L_Soft, Vmax_inv_Soft, _, _, _, _ = prep.forward_pass()

prep = preprocessing(sdss_xmatch, 'Hard', sensitivity_Hard)
log_L_Hard, Vmax_inv_Hard, _, _, _, _ = prep.forward_pass()

#%% Updating M, SFR and luminosities in catalog

updated = {'lgm': log_M, 'sfr': log_SFR, 'IR1': log_L_IR1, 'IR2': log_L_IR2,
           'IR4': log_L_IR4, 'Radio': log_L_Radio, 'Soft': log_L_Soft, 'Hard': log_L_Hard}

def update_data(data=sdss_xmatch, updated=updated):
    for band, array in updated.items():
        for name in data.columns:
            if band+'_luminosity' in name or band+'_tot_p50' in name:
                data[name] = array
    return data

sdss_xmatch = update_data()

sdss_xmatch.drop(columns=['lgm_tot_p16','lgm_tot_p84', 'sfr_tot_p16', 'sfr_tot_p84'],
                 inplace=True) # dropping the (now) useless columns


#%% Computing colors and their respective Vmax

IR_color = log_L_IR2 - log_L_IR1
Radio_color = log_L_Radio - log_L_IR4
X_color = log_L_Hard - log_L_Soft

IR_color_lit = sdss_xmatch['W1mag'].values - sdss_xmatch['W2mag'].values
IR_color_lit[np.isnan(IR_color)] = np.nan

Radio_color_lit = np.log10((sdss_xmatch['IR4_flux_(erg/s/cm2)'].values/((3*1e8)/(22.1*1e-6))) / \
    (sdss_xmatch['Radio_flux_(erg/s/cm2)'].values/(1.4*1e9))) #q24 (ratio of fluxes/Hz)
Radio_color_lit[np.isnan(Radio_color)] = np.nan
    
# X_color_lit = ?

Vmax_inv_IR_col = np.array([Vmax_inv_IR1, Vmax_inv_IR2]).max(0)
Vmax_inv_Radio_col = np.array([Vmax_inv_Radio, Vmax_inv_IR4]).max(0)
Vmax_inv_X_col = np.array([Vmax_inv_Soft, Vmax_inv_Hard]).max(0)

#%% Putting all together

new_data = np.column_stack((sigma_M, sigma_SFR, IR_color, IR_color_lit, Radio_color,
                            Radio_color_lit, X_color, Vmax_inv_IR1, Vmax_inv_IR2,
                            Vmax_inv_IR4, Vmax_inv_Radio, Vmax_inv_Soft, Vmax_inv_Hard,
                            Vmax_inv_IR_col, Vmax_inv_Radio_col, Vmax_inv_X_col))

columns = ['sigma_M', 'sigma_SFR', 'IR_color', 'IR_color_lit', 'Radio_color', 
           'Radio_color_lit', 'X_color', 'Vmax_inv_IR1', 'Vmax_inv_IR2',
           'Vmax_inv_IR4', 'Vmax_inv_Radio', 'Vmax_inv_Soft', 'Vmax_inv_Hard',
           'Vmax_inv_IR_col', 'Vmax_inv_Radio_col', 'Vmax_inv_X_col']

sdss_xmatch_add = pd.DataFrame(data=new_data, columns=columns)
sdss_xmatch = pd.concat([sdss_xmatch, sdss_xmatch_add], axis=1).fillna("")

#%% Saving data in new csv

sdss_xmatch.to_csv('final_xmatch_preprocessed.csv', header=True, index=False)

