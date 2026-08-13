#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:09:48 2026

@author: madeline.evenson
"""

#exploring cigale results
#these results are based on the vf_data.txt file Kim sent

import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import unumpy as unp

from astropy.io import fits
from astropy.table import Table

#to calculate and plot line of best fit in loglog space, taking each point's error into
#account, going to use Orthogonal Distance Regression from SciPy ODR Module
from scipy import odr


#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/Virgo/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are
cigaledir = homedir+'/cigale/'

n_ephot_dir = cigaledir + '/vf_ephot_north/out/'
s_ephot_dir = cigaledir + 'vf_ephot_south/out/'

n_herschel_dir = cigaledir + '/vf_north/out/'
s_herschel_dir = cigaledir + '/vf_south/out/'

# %%


n_ephot_file = n_ephot_dir + 'results.fits'
n_ephot_table = Table.read(n_ephot_file, format='fits')
n_ephot_df = n_ephot_table.to_pandas()

s_ephot_file = s_ephot_dir + 'results.fits'
s_ephot_table = Table.read(s_ephot_file, format='fits')
s_ephot_df = s_ephot_table.to_pandas()

n_herschel_file = n_herschel_dir + 'results.fits'
n_herschel_table = Table.read(n_herschel_file, format='fits')
n_herschel_df = n_herschel_table.to_pandas()

s_herschel_file = s_herschel_dir + 'results.fits'
s_herschel_table = Table.read(s_herschel_file, format='fits')
s_herschel_df = s_herschel_table.to_pandas()

kim_file = cigaledir + 'kim_results.txt'
kim_df = pd.read_csv(kim_file, sep=r'\s+')

ephot_df = pd.concat([n_ephot_df, s_ephot_df], ignore_index=True)
herschel_df = pd.concat([n_herschel_df, s_herschel_df], ignore_index=True)

# %%

#define linear model for ODR: y = m*x + c
def linear_func(p, x):
    return p[0] * x + p[1]



def odr_fit(df_x, df_y, value_col, err_col):
    
    x = np.asarray(df_x[value_col])
    y = np.asarray(df_y[value_col])
    xerr = np.asarray(df_x[err_col])
    yerr = np.asarray(df_y[err_col])
    
    mask = (
        np.isfinite(x) &
        np.isfinite(y) &
        np.isfinite(xerr) &
        np.isfinite(yerr) &
        (x > 0) &
        (y > 0))
    
    print(f'Keeping {mask.sum()} of {len(mask)} galaxies')
    
    x = x[mask]
    y = y[mask]
    xerr = xerr[mask]
    yerr = yerr[mask]
    
    #convert to log space
    log_x = np.log10(x)
    log_y = np.log10(y)
    
    #error propagation for log10 --> sigma_log(v) = sigma_v / (v*ln(10))
    log_xerr = xerr / (x * np.log(10))
    log_yerr = yerr / (y * np.log(10))
    
    model = odr.Model(linear_func)
    data = odr.RealData(log_x, log_y, sx=log_xerr, sy=log_yerr)
    fit = odr.ODR(data, model, beta0=[1.0, 1.0]).run()
    
    return x, y, xerr, yerr, fit



def plot_odr(x, y, xerr, yerr, fit, xlabel, ylabel, title, output_file, x_lowlim=None, x_highlim=None, y_lowlim=None, y_highlim=None):
    
    m, c = fit.beta
    
    xfit_log = np.linspace(np.log10(x.min()), np.log10(x.max()), 100)
    yfit_log = linear_func(fit.beta, xfit_log)
    
    xfit = 10**xfit_log
    yfit = 10**yfit_log
    
    xmin = min(x.min(), y.min())
    xmax = max(x.max(), y.max())
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', alpha=0.5, ecolor='gray', elinewidth=1.5, capsize=2)
    
    ax.plot(xfit, yfit, color='red', label=rf'Fit: $y=10^{{{c:.2f}}}x^{{{m:.2f}}}$')
    
    ax.plot([xmin, xmax], [xmin, xmax], 'k--', label=r'$y=x$')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    if x_lowlim is not None:
        xmin = x_lowlim
    if x_highlim is not None:
        xmax = x_highlim
    ax.set_xlim(xmin, xmax)
    
    if y_lowlim is not None:
        ymin = y_lowlim
    if y_highlim is not None:
        ymax = y_highlim
    ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()
    
    output_path = plotdir + output_file
    
    plt.savefig(output_path, dpi=150)
    plt.show()
    
    
    
def err_prop(df_1, df_2, value_col_1, value_col_2, err_col_1, err_col_2):
    
    #define data
    a_nom, a_err = np.array(df_1[value_col_1]), np.array(df_1[err_col_1])
    b_nom, b_err = np.array(df_1[value_col_2]), np.array(df_1[err_col_2])
    c_nom, c_err = np.array(df_2[value_col_1]), np.array(df_2[err_col_1])
    d_nom, d_err = np.array(df_2[value_col_2]), np.array(df_2[err_col_2])
    
    #package into UFloat arrays for propagation
    A = unp.uarray(a_nom, a_err)
    B = unp.uarray(b_nom, b_err)
    C = unp.uarray(c_nom, c_err)
    D = unp.uarray(d_nom, d_err)
    
    #calculate quotients and extract nominals and uncertainties
    #result X = a/b, result Y = c/d
    X = A / B
    Y = C / D
    
    x_vals = unp.nominal_values(X)
    x_errs = unp.std_devs(X)
    y_vals = unp.nominal_values(Y)
    y_errs = unp.std_devs(Y)
    
    return x_vals, x_errs, y_vals, y_errs
    
# %%

################ Stellar Mass Plotting ################

##### Stellar Mass: legacy ephot vs herschel #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    herschel_df,
    'bayes.stellar.m_star',
    'bayes.stellar.m_star_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy Stellar Mass ($M_\odot$)',
         ylabel=r'Herschel Stellar Mass ($M_\odot$)',
         title='Legacy vs Herschel Stellar Mass',
         output_file='Legacy-vs-Herschel-Stellar-Mass.png')



##### Stellar Mass: legacy ephot vs kim #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    kim_df,
    'bayes.stellar.m_star',
    'bayes.stellar.m_star_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy Stellar Mass ($M_\odot$)',
         ylabel=r'Kim Stellar Mass ($M_\odot$)',
         title='Legacy vs Kim Stellar Mass',
         output_file='Legacy-vs-Kim-Stellar-Mass.png')



##### Stellar Mass: herschel vs kim #####
x, y, xerr, yerr, fit = odr_fit(
    herschel_df,
    kim_df,
    'bayes.stellar.m_star',
    'bayes.stellar.m_star_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Herschel Stellar Mass ($M_\odot$)',
         ylabel=r'Kim Stellar Mass ($M_\odot$)',
         title='Herschel vs Kim Stellar Mass',
         output_file='Herschel-vs-Kim-Stellar-Mass.png')

# %%

################ SFR Plotting ################

##### SFR: legacy ephot vs herschel #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    herschel_df,
    'bayes.sfh.sfr',
    'bayes.sfh.sfr_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy SFR ($M_\odot yr^{-1}$)',
         ylabel=r'Herschel SFR ($M_\odot yr^{-1}$)',
         title='Legacy vs Herschel SFR',
         output_file='Legacy-vs-Herschel-SFR.png')



##### SFR: legacy ephot vs kim #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    kim_df,
    'bayes.sfh.sfr',
    'bayes.sfh.sfr_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy SFR ($M_\odot yr^{-1}$)',
         ylabel=r'Kim SFR ($M_\odot yr^{-1}$)',
         title='Legacy vs Kim SFR',
         output_file='Legacy-vs-Kim-SFR.png')



##### SFR: herschel vs kim #####
x, y, xerr, yerr, fit = odr_fit(
    herschel_df,
    kim_df,
    'bayes.sfh.sfr',
    'bayes.sfh.sfr_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Herschel SFR ($M_\odot yr^{-1}$)',
         ylabel=r'Kim SFR ($M_\odot yr^{-1}$)',
         title='Herschel vs Kim SFR',
         output_file='Herschel-vs-Kim-SFR.png')

# %%

################ Preparing for Dust Mass plotting ################

##### Dust mass --> preparing columns #####

factor = 1.98847e30 

ephot_df['bayes.dust.mass.msol'] = ephot_df['bayes.dust.mass'] / factor
ephot_df['bayes.dust.mass.msol_err'] = ephot_df['bayes.dust.mass_err'] / factor

herschel_df['bayes.dust.mass.msol'] = herschel_df['bayes.dust.mass'] / factor
herschel_df['bayes.dust.mass.msol_err'] = herschel_df['bayes.dust.mass_err'] / factor

kim_df['bayes.dust.mass.msol'] = kim_df['bayes.dust.mass'] / factor
kim_df['bayes.dust.mass.msol_err'] = kim_df['bayes.dust.mass_err'] / factor

# %%

################ Dust Mass Plotting ################

##### dust mass: legacy ephot vs herschel #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    herschel_df,
    'bayes.dust.mass.msol',
    'bayes.dust.mass.msol_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy Dust Mass ($M_\odot$)',
         ylabel=r'Herschel Dust Mass ($M_\odot$)',
         title='Legacy vs Herschel Dust Mass',
         output_file='Legacy-vs-Herschel-Dust-Mass.png',
         x_lowlim = 10**1,
         x_highlim = 10**9,
         y_lowlim = 10**1,
         y_highlim = 10**9)



##### dust mass: legacy ephot vs kim #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    kim_df,
    'bayes.dust.mass.msol',
    'bayes.dust.mass.msol_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy Dust Mass ($M_\odot$)',
         ylabel=r'Kim Dust Mass ($M_\odot$)',
         title='Legacy vs Kim Dust Mass',
         output_file='Legacy-vs-Kim-Dust-Mass.png',
         x_lowlim = 10**1,
         x_highlim = 10**9,
         y_lowlim = 10**1,
         y_highlim = 10**9)




##### dust mass: herschel vs kim #####
x, y, xerr, yerr, fit = odr_fit(
    herschel_df,
    kim_df,
    'bayes.dust.mass.msol',
    'bayes.dust.mass.msol_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'LHerschel Dust Mass ($M_\odot$)',
         ylabel=r'Kim Dust Mass ($M_\odot$)',
         title='Herschel vs Kim Dust Mass',
         output_file='Herschel-vs-Kim-Dust-Mass.png',
         x_lowlim = 10**1,
         x_highlim = 10**9,
         y_lowlim = 10**1,
         y_highlim = 10**9)

# %%

################ Preparing for Specific SFR plotting ################

##### Specific SFR --> prepping dataframes #####


# legacy ephot vs herschel #
x_vals, x_errs, y_vals, y_errs = err_prop(ephot_df, herschel_df, 
                                          'bayes.sfh.sfr', 'bayes.stellar.m_star', 
                                          'bayes.sfh.sfr_err', 'bayes.stellar.m_star_err')

ephot_df['ephot.herschel.dust.mass'] = x_vals
ephot_df['ephot.herschel.dust.mass_err'] = x_errs
herschel_df['ephot.herschel.dust.mass'] = y_vals
herschel_df['ephot.herschel.dust.mass_err'] = y_errs


# legacy ephot vs kim #
x_vals, x_errs, y_vals, y_errs = err_prop(ephot_df, kim_df,
                                          'bayes.sfh.sfr', 'bayes.stellar.m_star',
                                          'bayes.sfh.sfr_err', 'bayes.stellar.m_star_err')

ephot_df['ephot.kim.dust.mass'] = x_vals
ephot_df['ephot.kim.dust.mass_err'] = x_errs
kim_df['ephot.kim.dust.mass'] = y_vals
kim_df['ephot.kim.dust.mass_err'] = y_errs


# herschel vs kim #
x_vals, x_errs, y_vals, y_errs = err_prop(herschel_df, kim_df,
                                          'bayes.sfh.sfr', 'bayes.stellar.m_star',
                                          'bayes.sfh.sfr_err', 'bayes.stellar.m_star_err')

herschel_df['herschel.kim.dust.mass'] = x_vals
herschel_df['herschel.kim.dust.mass_err'] = x_errs
kim_df['herschel.kim.dust.mass'] = y_vals
kim_df['herschel.kim.dust.mass_err'] = y_errs

# %%

################ Specific SFR Plotting ################

##### Specific SFR: legacy vs Herschel #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    herschel_df,
    'ephot.herschel.dust.mass',
    'ephot.herschel.dust.mass_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy Specific SFR ($yr^{-1}$)',
         ylabel=r'Herschel Specific SFR ($yr^{-1}$)',
         title='Legacy vs Herschel Specific Star Formation Rate',
         output_file='Legacy-vs-Herschel-SSFR.png')


##### Specific SFR: legacy vs Kim #####
x, y, xerr, yerr, fit = odr_fit(
    ephot_df,
    kim_df,
    'ephot.kim.dust.mass',
    'ephot.kim.dust.mass_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Legacy Specific SFR ($yr^{-1}$)',
         ylabel=r'Kim Specific SFR ($yr^{-1}$)',
         title='Legacy vs Kim Specific Star Formation Rate',
         output_file='Legacy-vs-Kim-SSFR.png')


##### Specific SFR: herschel vs Kim #####
x, y, xerr, yerr, fit = odr_fit(
    herschel_df,
    kim_df,
    'herschel.kim.dust.mass',
    'herschel.kim.dust.mass_err')

plot_odr(x, y, xerr, yerr, fit,
         xlabel=r'Herschel Specific SFR ($yr^{-1}$)',
         ylabel=r'Kim Specific SFR ($yr^{-1}$)',
         title='Herschel vs Kim Specific Star Formation Rate',
         output_file='Herschel-vs-Kim-SSFR.png')
