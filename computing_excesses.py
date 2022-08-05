# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 03:02:10 2022

@author: alvar
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sdss_xmatch = pd.read_csv('final_xmatch_preprocessed.csv', keep_default_na=True)

#%%

class computing_excesses():
    def __init__(self, data, band, is_color=False, h=250, w=250):
        
        # Main variables
        self.log_M = data['lgm_tot_p50'].values 
        self.sigma_M = sdss_xmatch['sigma_M'].values
        self.log_SFR = data['sfr_tot_p50'].values
        self.sigma_SFR = sdss_xmatch['sigma_SFR'].values
        
        if is_color:
            self.log_L = data[band+'_color'].values
            self.Vmax_inv = data['Vmax_inv_'+band+'_col'].values
            
        else:
            self.log_L = data[band+'_luminosity_(erg/s)'].values
            self.Vmax_inv = data['Vmax_inv_'+band].values
        
        # Grid variables
        self.M_x = np.linspace(np.nanmin(self.log_M), np.nanmax(self.log_M), w)
        self.SFR_y = np.linspace(np.nanmin(self.log_SFR), np.nanmax(self.log_SFR), h)
        self.M_grid, self.SFR_grid = np.meshgrid(self.M_x, self.SFR_y)
        self.d_M = self.M_x[1] - self.M_x[0]
        self.d_SFR = self.SFR_y[1] - self.SFR_y[0]
        
        # Plots variables
        if is_color:
            self.title = '<'+band+' color>'
        else:
            self.title = '<log L> '+band+' band'
    
    
    def density(self, M_p, SFR_p, sigma_M_p, sigma_SFR_p, h=250, w=250, min_smoothing=0.1):
        
        M_i, SFR_j = self.M_grid, self.SFR_grid
        sigma_M_p += min_smoothing
        sigma_SFR_p += min_smoothing
        
        object_data = np.array([M_p, SFR_p, sigma_M_p, sigma_SFR_p])
        
        if not np.any(np.isnan(object_data)):
            term1 = (M_i - M_p)/sigma_M_p
            term2 = (SFR_j - SFR_p)/sigma_SFR_p
            term = -term1**2 - term2**2
            term = term/2
            term = np.exp(term)
            term = term/(np.pi*sigma_M_p*sigma_SFR_p)
                    
        else:
            term = np.zeros((h, w))
            
        return term
    
    
    def get_outputs(self, h=250, w=250):
        
        total_dens = np.zeros((h, w))
        weighted_dens = np.zeros((h, w))
        weighted_dens_2 = np.zeros((h, w))
        
        for p in range(len(self.log_L)):
            L_p, Vmax_inv_p = self.log_L[p], self.Vmax_inv[p]
            
            if not np.isnan(L_p):
                densidad_p = self.density(M_p=self.log_M[p], SFR_p=self.log_SFR[p],
                                          sigma_M_p =self.sigma_M[p], 
                                          sigma_SFR_p=self.sigma_SFR[p])
                
                w_p = densidad_p*Vmax_inv_p
                total_dens = total_dens + w_p
                weighted_dens = weighted_dens + L_p*w_p
                weighted_dens_2 = weighted_dens_2 + (L_p**2)*w_p
                
        L_mean = weighted_dens/total_dens
        L_mean_squared = weighted_dens_2/total_dens
        Dispersion = np.sqrt(L_mean_squared - L_mean**2)
        
        return total_dens, L_mean, Dispersion
    
    
    def get_isocontour_levels(self, Density, left=8.5, right=11.5, bottom=-2, top=1):
        
        dx, dy = self.d_M, self.d_SFR
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
        
        plt.contourf(self.M_grid, self.SFR_grid, Density, cmap='jet',
                     extent=(left,right,bottom,top))
        plt.xlim(left=left, right=right)
        plt.xlabel('log10(M)')
        plt.ylim(bottom=bottom, top=top)
        plt.ylabel('log10(SFR)')
        plt.colorbar()
        plt.show()
        
        plt.contourf(self.M_grid, self.SFR_grid, fraction, cmap='jet', 
                     extent=(left,right,bottom,top))
        plt.xlim(left=left, right=right)
        plt.xlabel('log10(M)')
        plt.ylim(bottom=bottom, top=top)
        plt.ylabel('log10(SFR)')
        plt.colorbar()
        plt.show()
        
        return fraction_inside
    
    
    def get_plots(self, L, contours, title, left=8.5, right=11.5, bottom=-2, top=1):
        
        plt.contourf(self.M_grid, self.SFR_grid, L, cmap='jet', 
                     extent=(left,right,bottom,top))
        plt.colorbar()
        plt.title(title)
        cont = plt.contour(self.M_grid, self.SFR_grid, contours,  
                           levels = [0.10, 0.50, 0.90], colors = 'black')
        plt.xlim(left=left, right=right)
        plt.xlabel('log10(M)')
        plt.ylim(bottom=bottom, top=top)
        plt.ylabel('log10(SFR)')
        plt.clabel(cont, inline=True, fontsize = 8)
        plt.show()
        plotted = 'Plot made!'
        
        return plotted
    
    
    def L_mean_vector(self, L, mean_L):
        
        mean_L_array = np.zeros(len(L))
        binsize_M, binsize_SFR = self.d_M / 2, self.d_SFR / 2
        
        for i in range(len(L)):
            if ((np.isnan(L[i])==False) and (np.isnan(self.log_M[i])==False) and \
                (np.isnan(self.log_SFR[i])==False)):
                row = 0
                column = 0
                for j in range(len(self.M_x)):
                    if np.abs(self.log_M[i] - self.M_x[j]) < binsize_M:
                        column = j
                        break
                    
                for k in range(len(self.SFR_y)):
                    if np.abs(self.log_SFR[i] - self.SFR_y[k]) < binsize_SFR:
                        row = k
                        break
                
                if ((row != 0) and (column != 0)):
                    mean_L_array[i] = mean_L[row, column]
                    
                else:
                    mean_L_array[i] = np.nan
                    
            else:
                mean_L_array[i] = np.nan
                
        return mean_L_array
    
    
    def forward_pass(self):
        
        total_dens, L_mean, _ = self.get_outputs()
        #contours = self.get_isocontour_levels(total_dens)
        #_ = self.get_plots(L_mean, contours, self.title)
        mean_L_array = self.L_mean_vector(self.log_L, L_mean)
        log_L_excess = self.log_L - mean_L_array
        
        return log_L_excess
    
            
#%% Getting the excesses

comp = computing_excesses(sdss_xmatch, 'IR1')
IR1_excess = comp.forward_pass()

comp = computing_excesses(sdss_xmatch, 'IR2')
IR2_excess = comp.forward_pass()

comp = computing_excesses(sdss_xmatch, 'IR4')
IR4_excess = comp.forward_pass()

comp = computing_excesses(sdss_xmatch, 'Radio')
Radio_excess = comp.forward_pass()

comp = computing_excesses(sdss_xmatch, 'Soft')
Soft_excess = comp.forward_pass()

comp = computing_excesses(sdss_xmatch, 'Hard')
Hard_excess = comp.forward_pass()



comp = computing_excesses(sdss_xmatch, 'IR', is_color=True)
IR_color_excess = comp.forward_pass()

comp = computing_excesses(sdss_xmatch, 'Radio', is_color=True)
Radio_color_excess = comp.forward_pass()

comp = computing_excesses(sdss_xmatch, 'X', is_color=True)
X_color_excess = comp.forward_pass()


#%% Writing csv with excesses and Vmax

excesses = np.column_stack((IR1_excess, IR2_excess, IR4_excess, Radio_excess,
                            Soft_excess, Hard_excess, IR_color_excess,
                            Radio_color_excess, X_color_excess))

columns = ['IR1_excess', 'IR2_excess', 'IR4_excess', 'Radio_excess',
           'Soft_excess', 'Hard_excess', 'IR_color_excess',
           'Radio_color_excess', 'X_color_excess']

excesses_df = pd.DataFrame(data=excesses, columns=columns)

Vmax = sdss_xmatch.iloc[:, -9:]

excesses_Vmax = pd.concat([excesses_df, Vmax], axis=1).fillna("")

excesses_Vmax.to_csv('excesses_Vmax.csv', header=True, index=False)
