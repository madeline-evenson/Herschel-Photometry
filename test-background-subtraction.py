#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 12:48:34 2026

@author: madeline.evenson
"""

#creating fake data set to use to test background subtraction pipeline

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

from scipy import stats
from matplotlib import colors
from scipy.stats import scoreatpercentile

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



#function to create square 2D array of data points pulled from gaussian

np.random.seed(42)

size = 100
mean = 1
std_dev = .2

gaussian_array = np.random.normal(loc=mean, scale=std_dev, size=(size, size))

print(np.round(gaussian_array, 2))


# %%

#function for displaying 2D np array as image

def np_image(image, v1perc=10, v2perc=95, logscale=True):
    
    #one end of color map assigned to v1perc percent lowest flux
    #other end of color map assigned to v2perc percent highest flux
    
    #print(image)
    
    fig, ax = plt.subplots(figsize=(6,6))
    
    #make sure image is np array
    nimage = np.array(image)
    
    #print(nimage)
    
    #determine the pixel values at the 10th and 95th percentile
    v1 = scoreatpercentile(nimage, v1perc)
    v2 = scoreatpercentile(nimage, v2perc)
    
    #print(v1, v2)
    
    #display using imshow
    if (logscale):
        im = ax.imshow(nimage, cmap='gray', aspect='equal', norm=colors.LogNorm(vmin=v1, vmax=v2), origin='lower')
    else:
        im = ax.imshow(nimage, cmap='gray', aspect='equal', vmin=v1, vmax=v2, origin='lower')
        
    cbar = fig.colorbar(im, ax=ax, fraction=0.04)
    cbar.ax.set_ylabel('counts')
    
    return fig, ax
        
        
        
# %%

#plot original fake data

fig, ax = np_image(gaussian_array, v1perc=10, v2perc=95, logscale=False)
ax.set_title('Original Data')

savepath = os.path.join(plotdir, "fake-100x100-array.png")
plt.savefig(savepath, dpi=150)
        
#calcualte median and mode of original fake data

bkgmode = stats.mode(gaussian_array, axis=None).mode
bkgmed = np.median(gaussian_array)

print(f'mode: {bkgmode}')
print(f'median: {bkgmed}')



#subtract background from fake data

median_subtracted = gaussian_array - bkgmed
mode_subtracted = gaussian_array - bkgmode


# %%



#plot background subtracted images

## when we subtract background we need to make sure that the color bar v1perc and v2perc stays consistent because otherwise will just scale with images
v1 = scoreatpercentile(gaussian_array, 10)
v2 = scoreatpercentile(gaussian_array, 95)

fig, ax = plt.subplots()
im = ax.imshow(median_subtracted, cmap='gray', vmin=v1, vmax=v2, origin='lower')
fig.colorbar(im, ax=ax, label='counts')
ax.set_title('Background (Median) Subtracted')

savepath = os.path.join(plotdir, "median-subtracted-100x100.png")
plt.savefig(savepath, dpi=150)




fig, ax = plt.subplots()
im = ax.imshow(mode_subtracted, cmap='gray', vmin=v1, vmax=v2, origin='lower')
fig.colorbar(im, ax=ax, label='counts')
ax.set_title('Background (Mode) Subtracted')

savepath = os.path.join(plotdir, "mode-subtracted-100x100.png")
plt.savefig(savepath, dpi=150)


        
        
        
        
        
        
        
        
        
        
        