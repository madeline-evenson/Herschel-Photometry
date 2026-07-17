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
import uncertainties as u
from uncertainties import unumpy as unp

from astropy.io import fits
from astropy.table import Table

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


ephot_df = pd.concat([n_ephot_df, s_ephot_df], ignore_index=True)
herschel_df = pd.concat([n_herschel_df, s_herschel_df], ignore_index=True)

# %%

fig, ax = plt.subplots(figsize=(8,6))



ax.errorbar(ephot_df['bayes.stellar.m_star'], herschel_df['bayes.stellar.m_star'],
            xerr=ephot_df['bayes.stellar.m_star_err'], yerr=herschel_df['bayes.stellar.m_star_err'],
            fmt='o', alpha=0.5, ecolor='gray', elinewidth=1.5, capsize=2)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Legacy Stellar Mass ($M_{\odot}$)')
ax.set_ylabel(r'Herschel Stellar Mass ($M_{\odot}$)')
ax.set_title('Stellar Mass')

ax.grid(True, linestyle='--', color='gray', alpha=0.7, linewidth=0.5)

plt.show()
plt.close()


# %%


fig, ax = plt.subplots(figsize=(8,6))



ax.errorbar(ephot_df['bayes.sfh.sfr'], herschel_df['bayes.sfh.sfr'],
            xerr=ephot_df['bayes.sfh.sfr_err'], yerr=herschel_df['bayes.sfh.sfr_err'],
            fmt='o', alpha=0.5, ecolor='gray', elinewidth=1.5, capsize=2)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Legacy SFR ($M_{\odot}yr^{-1}$)')
ax.set_ylabel(r'Herschel SFR ($M_{\odot}yr^{-1}$)')
ax.set_title('Star Formation Rate')

ax.grid(True, linestyle='--', color='gray', alpha=0.7, linewidth=0.5)

plt.show()
plt.close()

# %%


fig, ax = plt.subplots(figsize=(8,6))

factor = 1.98847e30


ax.scatter(ephot_df['bayes.dust.mass'] / factor, herschel_df['bayes.dust.mass'] / factor,
            marker='o', alpha=0.5, linewidths=1, edgecolors='black')

ax.errorbar(ephot_df['bayes.dust.mass'] / factor, herschel_df['bayes.dust.mass'] / factor,
            xerr=ephot_df['bayes.dust.mass_err'] / factor, yerr=herschel_df['bayes.dust.mass_err'] / factor,
            fmt='o', capsize=2, alpha=0.5, ecolor='gray', elinewidth=1.5)

ax.set_xscale('log')
ax.set_xlim(10**1, 10**9)

ax.set_yscale('log')
ax.set_ylim(10**1, 10**9)

ax.set_xlabel(r'Legacy Dust Mass ($M_{\odot}$)')
ax.set_ylabel(r'Herschel Dust Mass ($M_{\odot}$)')
ax.set_title('Dust Mass')

ax.grid(True, linestyle='--', color='gray', alpha=0.7, linewidth=0.5)

plt.show()
plt.close()


# %%

#error propagation for specific SFR

#1: define data
a_nom, a_err = np.array(ephot_df['bayes.sfh.sfr']), np.array(ephot_df['bayes.sfh.sfr_err'])
b_nom, b_err = np.array(ephot_df['bayes.stellar.m_star']), np.array(ephot_df['bayes.stellar.m_star_err'])
c_nom, c_err = np.array(herschel_df['bayes.sfh.sfr']), np.array(herschel_df['bayes.sfh.sfr_err'])
d_nom, d_err = np.array(herschel_df['bayes.stellar.m_star']), np.array(herschel_df['bayes.stellar.m_star_err'])

#2: package into UFloat arrays for propagation
A = unp.uarray(a_nom, a_err)
B = unp.uarray(b_nom, b_err)
C = unp.uarray(c_nom, c_err)
D = unp.uarray(d_nom, d_err)

#3: calculate quotients and extract nominals + uncertainties
#result X = a/b, result Y = c/d
X = A / B
Y = C / D

x_vals = unp.nominal_values(X)
x_errs = unp.std_devs(X)
y_vals = unp.nominal_values(Y)
y_errs = unp.std_devs(Y)

fig, ax = plt.subplots(figsize=(8,6))

ax.errorbar(x_vals, y_vals, xerr=x_errs, yerr=y_errs,
            fmt='o',  capsize=2, alpha=0.5, ecolor='gray', elinewidth=1.5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Legacy Specific SFR ($yr^{-1}$)')
ax.set_ylabel(r'Herschel Specific SFR ($yr^{-1}$)')
ax.set_title('Specific Star Formation Rate')

ax.grid(True, linestyle='--', color='gray', alpha=0.7, linewidth=0.5)

plt.show()
plt.close()

# %%



fig, ax = plt.subplots(figsize=(8,6))

ax.scatter(ephot_df['bayes.sfh.sfr'] / ephot_df['bayes.stellar.m_star'], herschel_df['bayes.sfh.sfr'] / herschel_df['bayes.stellar.m_star'],
            marker='o', alpha=0.5, linewidths=1, edgecolors='black')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'Legacy Specific SFR ($yr^{-1}$)')
ax.set_ylabel(r'Herschel Specific SFR ($yr^{-1}$)')
ax.set_title('Specific Star Formation Rate')

ax.grid(True, linestyle='--', color='gray', alpha=0.7, linewidth=0.5)

plt.show()
plt.close()







