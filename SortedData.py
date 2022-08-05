#!/usr/bin/env python
# coding: utf-8

# # SortedData class
# 
# ## Definition

#%% Importing libraries


get_ipython().run_line_magic('matplotlib', 'ipympl')
from matplotlib import pyplot as plt
import numpy as np


#%% Defining SortedData class


class SortedData(object):
    
    def __init__(self, data, weights=None, figname=None):
        
        self.original = data
        sorted_by_data = np.argsort(data.flatten())
        sorted_data = data.flatten()[sorted_by_data]
        valid = np.isfinite(sorted_data)
        self.sorted_by_data = sorted_by_data[valid]
        self.sorted_data = sorted_data[valid]
    
        if weights is None:
            weights = np.ones_like(data)
        self.weights = weights
        self.cumulative_mass = np.nancumsum(self.weights.flatten()[self.sorted_by_data])
        
        self.x0 = self.find_peak()
        self.p_normal = self.compute_p_normal()

        
    def find_peak(self, delta=.5):

        m = np.linspace(delta**2,  # somewhat arbitrary :^(
                    1/(1+delta),  # do not overshoot
                    int(np.sqrt(self.sorted_data.size)) # a reasonable number of trials
                   )
        cumulative_fraction = self.cumulative_mass / self.cumulative_mass[-1]
        x_top = np.interp((1+delta)*m, cumulative_fraction, self.sorted_data)
        x_tm = np.interp((1+delta/2)*m, cumulative_fraction, self.sorted_data)
        x_mid = np.interp(m, cumulative_fraction, self.sorted_data)
        x_bm = np.interp((1-delta/2)*m, cumulative_fraction, self.sorted_data)
        x_bot = np.interp((1-delta)*m, cumulative_fraction, self.sorted_data)
        rho_top = delta/2 * m / (x_top - x_tm)
        rho_bot = delta/2 * m / (x_bm - x_bot)
        peak = np.nanargmin((rho_top - rho_bot) ** 2)

        return min(x_mid[peak], (x_bot[peak]+x_top[peak])/2, (x_bm[peak]+x_tm[peak])/2) 


    def probability_density(self, nbins=None):
        
        if nbins is None:
            nbins = min(int(np.sqrt(self.sorted_data.size)), 100)

        m_bins_m = np.linspace(0, self.cumulative_mass[-1], nbins)
        m_bins_x = np.interp(m_bins_m, self.cumulative_mass, self.sorted_data)

        x_bins_x = np.linspace(self.sorted_data[0], self.sorted_data[-1], nbins)
        bins_x = (x_bins_x + m_bins_x) / 2
        bins_m = np.interp(bins_x, self.sorted_data, self.cumulative_mass)

        x = (bins_x[1:] + bins_x[:-1]) / 2
        dm = bins_m[1:] - bins_m[:-1]
        dx = bins_x[1:] - bins_x[:-1]

        return x, dm/dx/self.cumulative_mass[-1]

    
    def compute_p_normal(self):

        a_x, a_hist = self.probability_density()
        left = np.where(a_x < self.x0)
        x_sym = 2*self.x0-a_x[left]
        f_normal = np.clip(a_hist[left]/np.interp(x_sym, a_x, a_hist), 0, 1)
        p_normal = np.interp(self.original, x_sym[::-1], f_normal[::-1], left=1, right=0)

        return p_normal
    


#%% Loading data


#excess = np.genfromtxt('excesses.csv', delimiter=',',
#                      dtype={'names': ('IR1', 'IR2', 'IR4', 'Radio', 'Soft X', 'Hard X', 'IR color','Radio color', 'X color'),
#                            'formats': ('f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f')}
#                      )
excess = np.genfromtxt('excesses_Vmax.csv', delimiter=',', usecols = np.array(range(0, 9)),
                        names=('IR1', 'IR2', 'IR4', 'Radio', 'Soft X', 'Hard X', 'IR color', 'Radio color', 'X color'),
                        filling_values=np.nan, skip_header=True)

Vmax = np.genfromtxt('excesses_Vmax.csv', delimiter=',', usecols = np.array(range(9, 18)),
                      names=('IR1', 'IR2', 'IR4', 'Radio', 'Soft X', 'Hard X', 'IR color', 'Radio color', 'X color'),
                      filling_values=np.nan, skip_header=True)

excess_sorted = {}
for band in excess.dtype.names:    
    excess_sorted[band] = SortedData(data=excess[band], weights=None)



#%% Plots


def plot_histograms(data, ax):

    x = data.sorted_data
    a_x, a_hist = data.probability_density()

    w = data.weights[data.sorted_by_data]
    p = data.p_normal[data.sorted_by_data]
    normal_hist, bins = np.histogram(x, weights=w*p, bins=a_x, density=True)
    active_hist, bins = np.histogram(x, weights=w*(1-p), bins=a_x, density=True)
    a_x_mid = (a_x[:-1] + a_x[1:])/2
    normal_fraction = np.nansum(p)/np.nansum(w)
    normal_hist *= normal_fraction
    active_hist *= 1-normal_fraction

    ax.plot(a_x, a_hist, 'k-', alpha=.2, label=f'{x.size} galaxies')
    ax.plot(a_x_mid, normal_hist, 'b:',
            label=f'{100*np.nansum(data.p_normal)/x.size:.1f}({100*normal_fraction:.1f})% normal')
    ax.plot(a_x_mid, active_hist, 'r:',
            label=f'{100*np.nansum(1-data.p_normal)/x.size:.1f}({100*(1-normal_fraction):.1f})% active')
    ax.axvline(data.x0, c='k', ls='--')
    
    ax.grid(alpha=.25)
    ax.set_yscale('log')
    #ax.legend(title=data.name)
    ax.legend(loc='upper right')

    
def plot_correlation(data_x, data_y, ax):

    x = data_x.original
    y = data_y.original
    active_x = np.where(data_x.p_normal < 0.5)
    active_y = np.where(data_y.p_normal < 0.5)
    active_both = np.where((data_x.p_normal < 0.5) & (data_y.p_normal < 0.5))
    
    color = np.full(x.shape, 'c')
    color[active_x] = 'g'
    color[active_y] = 'y'
    color[active_both] = 'r'
    
    ax.scatter(x, y, c=color, s=1)

def fig_cornerplot(figname, sorted_data, bands):
    plt.close(figname)
    fig = plt.figure(figname, figsize=(10, 10))

    nbands = len(bands)
    ax = fig.subplots(nrows=nbands, sharey='row',
                      ncols=nbands, sharex='col',
                      gridspec_kw={'hspace':0, 'wspace':0}, squeeze=False)

    # First row: probability density

    for col in range(nbands):
        band = bands[col]
        xs = sorted_data[band]
        xmin = xs.sorted_data[0]
        xmax = xs.sorted_data[-1]

        plot_histograms(xs, ax[0, col])

        ax[0, col].set_title(f'{band} excess')
        #ax[0, col].set_xlim(xmin, xmin+1.7*(xmax-xmin))
        #ax[0, col].set_ylim(xmin, xmin+1.7*(xmax-xmin))

    ax[0, 0].set_ylabel('probability density')
    ax[0, -1].set_xlabel(f'{bands[-1]} excess')

    # Additional rows: correlation

    for row in range(1, nbands):
        band_y = bands[nbands-row]
        for col in range(nbands):
            band_x = bands[col]

            if col >= nbands-row:
                ax[row, col].axis('off')
            else:
                if col == nbands-row-1:
                    ax[row, col].set_xlabel(f'{band_x} excess')
                plot_correlation(sorted_data[band_x], sorted_data[band_y], ax[row, col])

        ax[row, 0].set_ylabel(f'{band_y} excess')


#%% Example of usage: Radio


fig_cornerplot('Radio', excess_sorted, ['IR4', 'Radio', 'Radio_color'])


#%% Example of usage: IR


fig_cornerplot('IR', excess_sorted, ['IR1', 'IR2', 'IR_color'])


#%% Example of usage: X


fig_cornerplot('X', excess_sorted, ['Soft_X', 'Hard_X', 'X_color'])

