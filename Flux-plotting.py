#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 11:05:14 2026

@author: madeline.evenson
"""

# this notebook is to start plotting fluxes :)

# %%



#import necessary libraries

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#get_ipython().run_line_magic('matplotlib', 'inline')
import astropy.units as u
import numpy as np
import pandas as pd
import glob
import sys
import os
#import wget

from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import Angle
from astropy.io.ascii import masked
from astropy.io import ascii
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch
from astropy.visualization import ImageNormalize
from astropy.visualization import LogStretch
from astropy.visualization import MinMaxInterval

from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import colors

from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture
from photutils.isophote import Ellipse
from photutils.aperture import aperture_photometry
from astropy.wcs.utils import proj_plane_pixel_scales

from scipy.stats import scoreatpercentile
from scipy import stats

import warnings
warnings.filterwarnings('ignore')

from reproject import reproject_interp
from IPython.display import clear_output


mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are

#make directory if doesn't exist
os.makedirs(plotdir, exist_ok=True)


# %%

csv_file = tabledir + 'Photometrytesting2.csv'

df = pd.read_csv(csv_file)

#print(df.columns)


#super conservative cuts --> can maybe keep 2665/2667, 4177/4179
overlap_VFID_list = [711, 712, 2665, 2667, 4177, 4179, 4961, 4962]

#df_clean = df.copy()

#removes rows where VF_ID column contains int value in overlap_VFID_list
#df_clean = df[~df['VF_ID'].isin(overlap_VFID_list)]

sma_labels = [f'SMA_AP0{i}' for i in range(1,9)]


# %%

#create and save all the 70 vs. 100 flux plots

for sma in sma_labels:

    fig, ax = plt.subplots(figsize=(8,6))

    x_axis = f'70Flux_{sma}'
    y_axis = f'100Flux_{sma}'

    ax.scatter(df[x_axis], df[y_axis], alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(f'{x_axis} (Jy)')
    ax.set_ylabel(f'{y_axis} (Jy)')
    ax.set_title(x_axis+' vs. '+ y_axis)


    savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

    print('Saving:', savepath)

    plt.savefig(savepath, dpi=150)

    plt.close()
    

#create and save all the 100 vs. 160 flux plots
for sma in sma_labels:

    fig, ax = plt.subplots(figsize=(8,6))

    x_axis = f'100Flux_{sma}'
    y_axis = f'160Flux_{sma}'

    ax.scatter(df[x_axis], df[y_axis], alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(f'{x_axis} (Jy)')
    ax.set_ylabel(f'{y_axis} (Jy)')
    ax.set_title(x_axis+' vs. '+ y_axis)


    savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

    print('Saving:', savepath)

    plt.savefig(savepath, dpi=150)

    plt.close()


#create and save all the 70 vs. 160 flux plots
for sma in sma_labels:

    fig, ax = plt.subplots(figsize=(8,6))

    x_axis = f'70Flux_{sma}'
    y_axis = f'160Flux_{sma}'

    ax.scatter(df[x_axis], df[y_axis], alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(f'{x_axis} (Jy)')
    ax.set_ylabel(f'{y_axis} (Jy)')
    ax.set_title(x_axis+' vs. '+ y_axis)


    savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

    print('Saving:', savepath)

    plt.savefig(savepath, dpi=150)

    plt.close()
    
 