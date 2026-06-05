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
import photutils

from scipy import stats
from photutils import datasets
from matplotlib import colors
from scipy.stats import scoreatpercentile
from astropy.io import fits
from astropy.visualization import ZScaleInterval

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

#print(np.round(gaussian_array, 2))


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

##### plot original fake data ######

fig, ax = np_image(gaussian_array, v1perc=10, v2perc=95, logscale=False)
ax.set_title('Original Data')

savepath = os.path.join(plotdir, "fake-100x100-array.png")
plt.savefig(savepath, dpi=150)





##### plot histogram of fake data #####

flat_gaussian_array = gaussian_array.ravel()
#print(flat_gaussian_array.round(2))

fig, ax = plt.subplots(figsize=(6,6))

plt.hist(flat_gaussian_array, bins=100, edgecolor='black')
plt.title('Fake Data Histogram')
plt.xlabel('Flux')
plt.ylabel('Counts')

savepath = os.path.join(plotdir, "Fake-data-histogram.png")
plt.savefig(savepath, dpi=150)


     

# %%
    

##### calcualte median and mode of original fake data #####

bkgmode = stats.mode(gaussian_array, axis=None).mode
bkgmed = np.median(gaussian_array)

print(f'mode: {bkgmode}')
print(f'median: {bkgmed}')



##### subtract background from fake data #####

median_subtracted = gaussian_array - bkgmed
mode_subtracted = gaussian_array - bkgmode


# %%



#####plot background subtracted images #####


#median subtracted image
fig, ax = np_image(median_subtracted, v1perc=10, v2perc=95, logscale=False)
ax.set_title('Background (Median) Subtracted')

savepath = os.path.join(plotdir, "median-subtracted-100x100.png")
plt.savefig(savepath, dpi=150)



#median subtracted histogram
flat_median_array = median_subtracted.ravel()

fig, ax = plt.subplots(figsize=(6,6))

plt.hist(flat_median_array, bins=100, edgecolor='black')
plt.title('Median-Subtracted Histogram')
plt.xlabel('Flux')
plt.ylabel('Counts')

savepath = os.path.join(plotdir, "Median-subtracted-histogram.png")
plt.savefig(savepath, dpi=150)




#mode subtracted image
fig, ax = np_image(mode_subtracted, v1perc=10, v2perc=95, logscale=False)
ax.set_title('Background (Mode) Subtracted')

savepath = os.path.join(plotdir, "mode-subtracted-100x100.png")
plt.savefig(savepath, dpi=150)



#mode subtracted histogram
flat_mode_array = mode_subtracted.ravel()

fig, ax = plt.subplots(figsize=(6,6))

plt.hist(flat_mode_array, bins=100, edgecolor='black')
plt.title('Mode-Subtracted Histogram')
plt.xlabel('Flux')
plt.ylabel('Counts')

savepath = os.path.join(plotdir, "Mode-subtracted-histogram.png")
plt.savefig(savepath, dpi=150)


# %%


#adding stars to fake data --> function from astropy tutorials
def stars(image, number, max_counts=10000, fwhm=9.4, min_separation=20, gain=1):
    """
    Add some stars to the image.
    """
    from photutils.datasets import make_model_image, make_model_params
    from photutils.psf import CircularGaussianPSF

    psf_model = CircularGaussianPSF(fwhm=fwhm)
    max_counts *= 100  # approx. peak amplitude to flux
    params = make_model_params(image.shape, 
                               n_sources=number,
                               flux=(max_counts / 10, max_counts),
                               min_separation=min_separation,
                               border_size=10)

    return make_model_image(image.shape, psf_model, params,
                            progress_bar=True)



synthetic_image = np.zeros([1000, 1000])


#add stars onto blank iamge
stars_only = stars(synthetic_image, 50, max_counts=2000)

fig, ax = plt.subplots(figsize=(6, 6))
ax.imshow(stars_only, cmap='gray', origin='lower', vmin=0, vmax=np.percentile(stars_only, 99.9))

ax.set_title('Stars only')


# %%


     

    
#creating synthetic extended source to add to fake data to later remove

size = 100

center_x, center_y = 50, 50
sigma_x, sigma_y = 5, 5  #controls the "extension" of the source

#create the coordinate grids
y, x = np.ogrid[:size, :size]

#calcualte the gaussian extended source

source_peak_flux = 1.7
extended_source = source_peak_flux * np.exp(-((x - center_x)**2 / (2*sigma_x**2) + (y - center_y)**2 / (2*sigma_y**2)))

final_image = extended_source + gaussian_array

fig, ax = np_image(final_image, v1perc=10, v2perc=95, logscale=False)
ax.set_title('Fake Data with Synthetic Source')

savepath = os.path.join(plotdir, "Fake-data-plus-synthetic-source.png")
plt.savefig(savepath, dpi=150)




#synthetic source histogram
flat_final_array = final_image.ravel()

fig, ax = plt.subplots(figsize=(6,6))

plt.hist(flat_final_array, bins=100, edgecolor='black')
plt.title('Synthetic Source Histogram')
plt.xlabel('Flux')
plt.ylabel('Counts')

savepath = os.path.join(plotdir, "Synthetic-source-histogram.png")
plt.savefig(savepath, dpi=150)

        
        

#adding stars to fake data with synthetic source
star_data = stars(final_image, number=10, max_counts=20, fwhm=1, min_separation=10)
final_image = final_image + star_data

fig, ax = np_image(final_image, v1perc=10, v2perc=95, logscale=False)
ax.set_title('Synthetic Source + Stars')

savepath = os.path.join(plotdir, "Synthetic-source-with-stars.png")
plt.savefig(savepath, dpi=150)



# %%


#this part of the notebook is dedicated to learning how to create a mask and then using routine to mask bright sources
#once making masks are understood, then will figure out how to make masks bigger/smaller

data = final_image


#define brightness threshold
threshold = 1.5

#create boolean mask
    # True = 1 = good pixels (keep)
    # False = 0 = bright pixels (mask)

mask = data < threshold

#apply mask to data
masked_data = np.ma.masked_array(data, mask=mask)

#visualize data using astropy's ZScale
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
interval = ZScaleInterval()
vmin, vmax = interval.get_limits(data)

#original image
im1 = ax[0].imshow(data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
ax[0].set_title('Original Image')
ax[0].axis('off')

#masked image
#matplotlib automatically renders masked values (False) transparent/white
im2 = ax[1].imshow(masked_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
ax[1].set_title(f'Masked Pixels (> {threshold})')
ax[1].axis('off')

fig.colorbar(im2, ax=ax[1])

plt.tight_layout()

savepath = os.path.join(plotdir, "Masked-plus-unmasked.png")
plt.savefig(savepath, dpi=150)





        
        