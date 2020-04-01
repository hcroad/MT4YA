# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 07:42:44 2016


This code is used to drive the dao model (you need to make sure you have both 
files present)

There are three options here to change (you could explore the space by looping)
over this function for example
alpha - efficieny of wave damping (between 0 and 1)
trans - transit time for wave to cross Pacific in days 
N - length of integration (in units of half a day)
if you run with no options you get alpha=0.6, trans=400, N=20000


There is one additional parameters used in the MT4YA problem following 
Bouttle et al. (2006; doi: 10.1119/1.2358155)
beta - a crude way of representing warming SST as a constant increase over time
with units of degrees per year

This is set to zero as standard

alpha (unitless)
beta (K / year)
delta = trans (days)
N = t = (half days)


@author: andrewcharlton-perez
"""

import matplotlib.pylab as plt
import dao
from matplotlib import mlab
import numpy as np
from matplotlib import colors
import scipy.signal as signal
from matplotlib import cm


#1 ----------------------------------------------------------------------------
#Change alpha in the DAO where beta = 0, trans = 400.
#Lowest alpha for oscillatory motion: 0.59

#Pick a range of alphas:

a_range = [0.4, 0.5, 0.6, 0.7, 0.8]
label = [r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.7$', r'$\alpha=0.8$']

fig, ax = plt.subplots(3 , 1, figsize=(10,8))

for i in range(5):
    
    time,T,dt=dao.calc_dao(alpha=a_range[i], N = 10000, trans=400, beta=0)
    
    # computer the power spectrum using Welch's method
    time_step=(time[1]-time[0])
    power, freq = mlab.psd(T, NFFT=10001, Fs=1. / time_step,
                      window=mlab.window_hanning, noverlap=2000)

    period=1./freq
    
    ax[0].plot(time, T, label=label[i])
    ax[0].set_xlabel('Time (years)', fontsize=17)
    ax[0].set_ylabel('T anomaly (K)', fontsize=17)
    ax[0].tick_params(labelsize=14)
    ax[0].set_title('a)', fontsize=17)
    
    ax[1].plot(T, dt, label=label[i])
    ax[1].set_xlabel('T anomaly (K)', fontsize=17)
    ax[1].set_ylabel(r'$\frac{\partial T}{\partial t}$ (K yr$^{-1}$)', fontsize=17)
    ax[1].tick_params(labelsize=14)
    ax[1].set_title('b)', fontsize=17)
    
    ax[2].plot(period, power, label=label[i])
    ax[2].set_xlabel('Time period (years)', fontsize=17)
    ax[2].set_ylabel('Power', fontsize=17)
    ax[2].set_xlim([0,10])
    ax[2].set_title('c)', fontsize=17)
    ax[2].tick_params(labelsize=14)
    ax[2].legend(fontsize=14, loc=6)
     
plt.subplots_adjust(hspace=0.65, top=0.95, bottom=0.1)
plt.savefig('Fig1.png', dpi=200)

#2 ----------------------------------------------------------------------------
#Change delta in the DAO where beta = 0, alpha = 0.6.
#Lowest delta for oscillatory motion: 347 days

#Pick a range of deltas:
d_range = [300, 350, 400, 450, 500]
label = [r'$\delta=300$', r'$\delta=350$', r'$\delta=400$', r'$\delta=450$', r'$\delta=500$']

fig, ax = plt.subplots(3 , 1, figsize=(10,8))

for i in range(5):
    
    time,T,dt=dao.calc_dao(alpha=0.6, N=10000, trans=d_range[i], beta=0)
    
    # computer the power spectrum using Welch's method
    time_step=(time[1]-time[0])
    power, freq = mlab.psd(T, NFFT=10001, Fs=1. / time_step,
                      window=mlab.window_hanning, noverlap=2000)

    period=1./freq
    
    ax[0].plot(time, T, label=label[i])
    ax[0].set_xlabel('Time (years)', fontsize=17)
    ax[0].set_ylabel('T anomaly (K)', fontsize=17)
    ax[0].tick_params(labelsize=14)
    ax[0].set_title('a)', fontsize=17)
    
    ax[1].plot(T, dt, label=label[i])
    ax[1].set_xlabel('T anomaly (K)', fontsize=17)
    ax[1].set_ylabel(r'$\frac{\partial T}{\partial t}$ (K yr$^{-1}$)', fontsize=17)
    ax[1].tick_params(labelsize=14)
    ax[1].set_title('b)', fontsize=17)
    
    ax[2].plot(period, power, label=label[i])
    ax[2].set_xlabel('Time period (years)', fontsize=17)
    ax[2].set_ylabel('Power', fontsize=17)
    ax[2].set_xlim([0,10])
    ax[2].set_title('c)', fontsize=17)
    ax[2].tick_params(labelsize=14)
    ax[2].legend(fontsize=14, loc=6)
     
plt.subplots_adjust(hspace=0.65, top=0.95, bottom=0.1)
plt.savefig('Fig2.png', dpi=200)

#3 ----------------------------------------------------------------------------
#Change beta in the DAO where alpha = 0.6, delta = 400
#Lowest delta for oscillatory motion: 0.17

b_range = [0, 0.05, 0.1, 0.15, 0.2]
label = [r'$\beta=0.00$', r'$\beta=0.05$', r'$\beta=0.10$', r'$\beta=0.15$', r'$\beta=0.20$']

fig, ax = plt.subplots(3 , 1, figsize=(10,8))

for i in range(5):
    
    time,T,dt=dao.calc_dao(alpha=0.6, N=10000, trans=400, beta=b_range[i])
    
    # computer the power spectrum using Welch's method
    time_step=(time[1]-time[0])
    power, freq = mlab.psd(T, NFFT=10001, Fs=1. / time_step,
                      window=mlab.window_hanning, noverlap=2000)

    period=1./freq
    
    ax[0].plot(time, T, label=label[i])
    ax[0].set_xlabel('Time (years)', fontsize=17)
    ax[0].set_ylabel('T anomaly (K)', fontsize=17)
    ax[0].tick_params(labelsize=14)
    ax[0].set_title('a)', fontsize=17)
    
    ax[1].plot(T, dt, label=label[i])
    ax[1].set_xlabel('T anomaly (K)', fontsize=17)
    ax[1].set_ylabel(r'$\frac{\partial T}{\partial t}$ (K yr$^{-1}$)', fontsize=17)
    ax[1].tick_params(labelsize=14)
    ax[1].set_title('b)', fontsize=17)
    
    ax[2].plot(period, power, label=label[i])
    ax[2].set_xlabel('Time period (years)', fontsize=17)
    ax[2].set_ylabel('Power', fontsize=17)
    ax[2].set_xlim([0,10])
    ax[2].set_title('c)', fontsize=17)
    ax[2].tick_params(labelsize=14)
    ax[2].legend(fontsize=14, loc=6)

plt.subplots_adjust(hspace=0.65, top=0.95, bottom=0.1)
plt.savefig('Fig3.png', dpi=200)


#4 --------------------------------------------------------------------------

beta_range = np.arange(0, 0.205, 0.005)
beta_periods = np.zeros(41)

for i in range(41):
    
    time,T,dt=dao.calc_dao(alpha=0.6, N=20000, trans=400, beta=beta_range[i])
    
    peakind = signal.find_peaks(T)

    Tn = np.zeros(len(peakind[0]))
    tn = Tn.copy()
    period = np.zeros(len(peakind[0])-1)

    for j in range(len(peakind[0])):
        Tn[j] = T[peakind[0][j]]
        tn[j] = time[peakind[0][j]]

    for k in range(len(peakind[0])-1):
        period[k] = tn[k+1] - tn[k]
        mean_period = np.mean(period)
        
    beta_periods[i] = mean_period
    
fig, ax = plt.subplots(1,1, figsize=(7,6))

ax.plot(beta_range, beta_periods, color='royalblue', linewidth=2.5, zorder=2)
ax.set_xlabel(r'$\beta$ (K yr$^{-1}$)', fontsize=17)
ax.set_ylabel('Time period (years)', fontsize=17)
ax.tick_params(labelsize=14)
ax.set_xticks([0, 0.04, 0.08, 0.12, 0.16, 0.2])
ax.set_ylim([2,4])
ax.vlines(0.16, 2, 4, linestyle='--', linewidth=2.5, zorder=1)

plt.savefig('Fig4.png', dpi=200)


#5 ----------------------------------------------------------------------------
#Change beta and alpha, where delta=400

a_range2 = np.arange(0.5, 0.72, 0.02)
b_range2 = np.arange(0.0, 0.22, 0.02)

xticklabels = ['0.50', '0.52', '0.54', '0.56', '0.58', '0.60', '0.62', '0.64', '0.66', '0.68', '0.70']
yticklabels = ['0.20', '0.18', '0.16', '0.14', '0.12', '0.10', '0.08', '0.06', '0.04', '0.02', '0.00']


ab_grid = np.zeros([11,11])

for i in range(11):
    for j in range(11):
        # drive calculation and plot output
        time,T,dt=dao.calc_dao(alpha=a_range2[j], trans=400, N=70000, beta=b_range2[(10 - i)])

        # computer the power spectrum using Welch's method
        time_step=(time[1]-time[0])
        power, freq = mlab.psd(T, NFFT=10001, Fs=1. / time_step,
                      window=mlab.window_hanning, noverlap=2000)
        #period
        period=1./freq
        
        if np.argmax(power) == 0:
            ab_grid[i,j] = 5
        else:
            ab_grid[i,j] = 15
            
cmap = colors.ListedColormap(['red','blue'])
bounds = [0, 10, 20]
norm = colors.BoundaryNorm(bounds, cmap.N)


fig, ax = plt.subplots(figsize=(10,10))
ax.imshow(ab_grid, cmap=cmap, norm=norm)

ax.grid(which='major', axis='both', linestyle='--', color='k', linewidth=0.5)
ax.set_xticks(np.arange(0, 11, 1));
ax.set_yticks(np.arange(0, 11, 1));
ax.set_xticklabels(xticklabels, fontsize=20, rotation=20)
ax.set_yticklabels(yticklabels, fontsize=20)
ax.hlines(0.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(1.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(2.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(3.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(4.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(5.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(6.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(7.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(8.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(9.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(10.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(0.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(1.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(2.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(3.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(4.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(5.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(6.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(7.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(8.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(9.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(10.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.set_ylabel(r'$\beta$ (K yr$^{-1}$)', fontsize=28)
ax.set_xlabel(r'$\alpha$', fontsize=28)

plt.savefig('Fig5.png', dpi=200)

#6 ----------------------------------------------------------------------------
#Change beta and delta, where alpha=0.6

d_range2 = np.arange(325, 407.5, 7.5)
b_range2 = np.arange(0.0, 0.22, 0.02)

xticklabels = ['325', '332.5', '340', '347.5', '355', '362.5', '370', '377.5', '385', '392.5', '400']
yticklabels = ['0.20', '0.18', '0.16', '0.14', '0.12', '0.10', '0.08', '0.06', '0.04', '0.02', '0.00']


ab_grid = np.zeros([11,11])

for i in range(11):
    for j in range(11):
        # drive calculation and plot output
        time,T,dt=dao.calc_dao(alpha=0.6, N=70000, trans=d_range2[j], beta=b_range2[(10 - i)])

        # computer the power spectrum using Welch's method
        time_step=(time[1]-time[0])
        power, freq = mlab.psd(T, NFFT=10001, Fs=1. / time_step,
                      window=mlab.window_hanning, noverlap=2000)
        #period
        period=1./freq
        
        if np.argmax(power) == 0:
            ab_grid[i,j] = 5
        else:
            ab_grid[i,j] = 15
            
cmap = colors.ListedColormap(['red','blue'])
bounds = [0, 10, 20]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(10,10))
ax.imshow(ab_grid, cmap=cmap, norm=norm)

ax.grid(which='major', axis='both', linestyle='--', color='k', linewidth=0.5)
ax.set_xticks(np.arange(0, 11, 1));
ax.set_yticks(np.arange(0, 11, 1));
ax.set_xticklabels(xticklabels, fontsize=20, rotation=20)
ax.set_yticklabels(yticklabels, fontsize=20)
ax.hlines(0.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(1.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(2.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(3.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(4.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(5.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(6.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(7.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(8.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(9.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.hlines(10.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(0.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(1.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(2.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(3.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(4.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(5.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(6.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(7.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(8.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(9.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.vlines(10.5, -0.5, 10.5, linestyle='-', linewidth=3)
ax.set_ylabel(r'$\beta$ (K yr$^{-1}$)', fontsize=28)
ax.set_xlabel(r'$\delta$ (days)', fontsize=28)

plt.savefig('Fig5.png', dpi=200)

#beta = 0.05 is stable for alpha = 0.59 (delta=400)
#beta = 0.05 is unstable for delta = 354 (alpha=0.6)
#Dominant timescale = 3.43 years

#6 ----------------------------------------------------------------------------
beta_range = np.arange(0, 0.21, 0.01)
delta_range = np.arange(325, 403.75, 3.75)
levs1 = [2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4]

beta_u = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18]
delta_u = [347.5, 347.5, 347.5, 355, 362.5, 370, 377.5, 385, 392.5, 400]


array = np.zeros([21,21])

for i in range(21):
    for j in range(21):
        time,T,dt=dao.calc_dao(alpha=0.6, N=20000, trans=delta_range[i], beta=beta_range[j])
        
        peakind = signal.find_peaks(T)

        Tn = np.zeros(len(peakind[0]))
        tn = Tn.copy()
        period = np.zeros(len(peakind[0])-1)

        for k in range(len(peakind[0])):
            Tn[k] = T[peakind[0][k]]
            tn[k] = time[peakind[0][k]]

        for m in range(len(peakind[0])-1):
            period[m] = tn[m+1] - tn[m]
            mean_period = np.mean(period)
        
        array[j,i] = mean_period
        

 
plt.figure(figsize=(8.5,6))    
cs = plt.contourf(delta_range, beta_range, array, levs1, cmap=cm.viridis, extend='both')
cbar = plt.colorbar(cs, label='Time period (years)')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='Time period (years)', size=17)
plt.tick_params(labelsize=14)
plt.xlabel(r'$\delta$ (days)', fontsize=17)
plt.ylabel(r'$\beta$ (K yr$^{-1}$)', fontsize=17)
plt.plot(delta_u, beta_u, color='black', linewidth=3)
plt.savefig('Fig6.png', dpi=200)


        

