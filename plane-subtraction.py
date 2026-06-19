#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 17:37:33 2026

@author: madeline.evenson
"""

#the goal of this notebook is to eventually fit a plane to the background of the galaxies with the new masks

#first we'll just calculate the median background for each galaxy with the new masks, and then we will fit planes 
    #to the 2D background array
    
    
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval

#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are

#directory for original fits 
fitsdir = datadir + 'pipeline/'

#directory for new masks
new_mask_dir = datadir + '/new_masked/'

#directory for visualizations of the new masks
new_mask_plot_dir = datadir + '/new_mask_images/'

#directory for visualizations of the new masks on herschel images
new_masked_herschel_dir = datadir + '/new_masked_herschel/'

# %%

#import masking_funct.py
from masking_funct import *
    
# %%

#first try for one single FITS file

#file = new_mask_dir + '/new_IC0694_maskedG.fits'

#with fits.open(file) as hdul:
    #print structure of file
    #hdul.info()
    
    #extract data and header from primary HDU
    #data = hdul[0].data
    #header = hdul[0].header
    
    #calculate the median using np.nanmedian to avoid a median of nan
    #median_val = np.nanmedian(data)
    
#print(f"The median value is: {median_val}")

#image_data = fits.getdata(file, ext=0)



#plot mask
#zscale = ZScaleInterval()
#vmin, vmax = zscale.get_limits(image_data)

#plt.figure(figsize=(8,8))
#plt.imshow(image_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()




#plot median-subtracted mask 
#median_subtracted = image_data - median_val
#print(f"The median of the median-subtracted fits is: {np.nanmedian(median_subtracted)}")

#vmin, vmax = zscale.get_limits(median_subtracted)

#plt.figure(figsize=(8,8))
#plt.imshow(median_subtracted, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()




#plot original fits
#fits_file = fitsdir + 'IC0694/HPPUNIMAPG/hpacs_25HPPUNIMAPB_green_1129_p5834_00_v1.0_1471613232020.fits'
#fits_data = fits.getdata(fits_file, ext=1)

#zscale = ZScaleInterval()
#vmin, vmax = zscale.get_limits(fits_data)

#plt.figure(figsize=(8,8))
#plt.imshow(fits_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()



#plot median-subtracted original fits
#subtracted_fits_data = fits_data - median_val

#plt.figure(figsize=(8,8))
#plt.imshow(subtracted_fits_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()

# %%

#now subtract median from ALL fits files

csv_file = tabledir+'/Photometrytesting.csv'
#csv_file2 = tabledir+'/Herschelstuff.csv'

galaxy = Table.read(csv_file)
#galaxy2 = Table.read(csv_file2)

galaxy_length = len(galaxy)

sma_labels = [f'SMA_AP0{i}' for i in range(1, 9)]
flux_labels = [f'{wave}Flux_{sma}' for wave in (70, 100, 160) for sma in sma_labels]
good_pixel_labels = [f'{wave}GoodPixels_{sma}' for wave in (70, 100, 160) for sma in sma_labels]

#initialize new columns in the galaxy table
for label in flux_labels:
    galaxy[label] = np.full(len(galaxy), np.nan)
for label2 in good_pixel_labels:
    galaxy[label2] = np.full(len(galaxy), np.nan)


#initialize dictionary to streamline code
bands = [
    {"wave": 70, "pixscale": 1.60000001784, "suffix": "B", "folder": "HPPUNIMAPB", "partial": "hpacs_25HPPUNIMAPB"},
    {"wave": 100, "pixscale": 1.60000001784, "suffix": "G", "folder": "HPPUNIMAPG", "partial": "hpacs_25HPPUNIMAPG"},
    {"wave": 160, "pixscale": 3.20000003568, "suffix": "R", "folder": "HPPUNIMAPR", "partial": "hpacs_25HPPUNIMAPR"}
]




for i in range(len(galaxy)):
    
    galaxy_name = str(galaxy['GALAXY'][i])
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"
    
    path = os.path.join(datadir, 'pipeline', galaxy_name)
    legacy_path = os.path.join('/Users/madeline.evenson/Research/Virgo/HTML-building/galaxy/png/', f"{VFID}-{galaxy_name}-LS.jpg")
    
    if not os.path.exists(path):
        print('Missing pipeline folder:', path)
        continue
    
    print('Processing:', galaxy_name)
    
    for band in bands:
        
        wave = str(band['wave'])
        pixscale = band['pixscale']
        suffix = band['suffix']
        
        
        
        ###### find image file #######
        image_fits_path = os.path.join(path, f'HPPUNIMAP{suffix}')
        masked_path = os.path.join(datadir, 'new_masked', f'new_{galaxy_name}_masked{suffix}.fits')
        
        if os.path.exists(masked_path):
            
            #masked path
            mask_found_file = masked_path
            #print(f'new mask path = {mask_found_file}')
            
            #loop over all color bands
            if suffix in ['B', 'G']: #blue and green bands
                partial_name = f'hpacs_25HPPUNIMAPB'
            else: #red band
                partial_name = f'hpacs_25HPPUNIMAPR'
            
            #original fits file path
            image_found_files = find_files(image_fits_path, partial_name)
            #print(f'original fits file path = {image_found_files}')
            
            if not image_found_files:
                #print(f'No original FITS found for {galaxy_name} {wave}')
                continue
            
            image_file = image_found_files[0]
            
        else: 
            destination_folder = os.path.join(path, band['folder'])
            mask_found_files = find_files(destination_folder, band['partial'])
            
            if not mask_found_files:
                #print('No file found:', galaxy_name, wave)
                continue
            
            mask_found_file = mask_found_files[0]
            
            
            
        #read in new mask data and calculate median   
        mask_data, mask_header = fits.getdata(mask_found_file, header=True)
        mask_median = np.nanmedian(mask_data)
        #print(f'mask median: {mask_median}')
        
        #read in original image fits data and subtract mask median
        image_data, image_header = fits.getdata(image_file, header=True)
        subtracted_image_data = image_data - mask_median
        #print(f'median of subtracted image: {np.nanmedian(subtracted_image_data)}')
        
        #output FITS file path
        output_fits = os.path.join(datadir, 'median-subtracted-fits', f'{galaxy_name}_{suffix}_median_subtracted.fits') 
        #print(f'output path: {output_fits}')
        
        #write median-subtracted fits files to output FITS file path
        fits.writeto(output_fits, subtracted_image_data, header=image_header, overwrite=True)


# %%

#now we want to fit a plane to the background of one example galaxy, then adjust the code to do it for all galaxies


#new numpy functions used:
    # numpy.column_stack() --> stacks 1D arrays as columns into a single 2D array
    # numpy.linalg.lstsq --> returns the least-squares solution to a linear matrix equation Ax = b
            # A - coefficent matrix
            # b - solution/"dependent" matrix
    # numpy.ones_like() --> returns a new array filled with ones that shares the exact shape and data type of an existing array

#file = new_mask_dir + '/new_NGC3622_maskedG.fits'

#data, header = fits.getdata(file, header=True)
#print(f'median of data: {np.nanmedian(data)}')

#plot original mask
#zscale = ZScaleInterval()
#vmin, vmax = zscale.get_limits(data)

#plt.figure(figsize=(8,8))
#plt.imshow(data, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()



#create a boolean mask of valid (unmasked) pixels
#mask = ~np.isnan(data)

#generate x and y coordinate grids for the valid pixels
#y, x = np.indices(data.shape)

#flatten arrays and filter out masked pixels
#x_flat = x[mask].flatten()
#y_flat = y[mask].flatten()
#z_flat = data[mask].flatten()

#set up the coefficient matrix A and solve using lstsq
#A = np.column_stack((x_flat, y_flat, np.ones_like(x_flat)))
#c0, c1, c2 = np.linalg.lstsq(A, z_flat, rcond=None)[0]

#c0 --> x slope
#c1 --> y slope
#c2 --> intercept

#create full background plane model
#background_plane = c0*x + c1*y + c2

#subtract plane from original image
#subtracted_data = data - background_plane

#print(f'median of plane-subtracted data: {np.nanmedian(subtracted_data)}')

#plt.figure(figsize=(8,8))
#plt.imshow(subtracted_data, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
#plt.colorbar()


# %%

#now apply plane-subtraction pipeline to all of the galaxies in all bands



for i in range(len(galaxy)):
    
    galaxy_name = str(galaxy['GALAXY'][i])
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"
    
    path = os.path.join(datadir, 'pipeline', galaxy_name)
    legacy_path = os.path.join('/Users/madeline.evenson/Research/Virgo/HTML-building/galaxy/png/', f"{VFID}-{galaxy_name}-LS.jpg")
    
    if not os.path.exists(path):
        print('Missing pipeline folder:', path)
        continue
    
    print('Processing:', galaxy_name)
    
    for band in bands:
        
        wave = str(band['wave'])
        pixscale = band['pixscale']
        suffix = band['suffix']
        
        
        
        ###### find image file #######
        image_fits_path = os.path.join(path, f'HPPUNIMAP{suffix}')
        masked_path = os.path.join(datadir, 'new_masked', f'new_{galaxy_name}_masked{suffix}.fits')
        
        if os.path.exists(masked_path):
            
            #masked path
            mask_found_file = masked_path
            #print(f'new mask path = {mask_found_file}')
            
            #loop over all color bands
            if suffix in ['B', 'G']: #blue and green bands
                partial_name = f'hpacs_25HPPUNIMAPB'
            else: #red band
                partial_name = f'hpacs_25HPPUNIMAPR'
            
            #original fits file path
            image_found_files = find_files(image_fits_path, partial_name)
            #print(f'original fits file path = {image_found_files}')
            
            if not image_found_files:
                #print(f'No original FITS found for {galaxy_name} {wave}')
                continue
            
            image_file = image_found_files[0]
            
        else: 
            destination_folder = os.path.join(path, band['folder'])
            mask_found_files = find_files(destination_folder, band['partial'])
            
            if not mask_found_files:
                #print('No file found:', galaxy_name, wave)
                continue
            
            mask_found_file = mask_found_files[0]
            
            
            
        #read in new mask data and calculate median   
        mask_data, mask_header = fits.getdata(mask_found_file, header=True)
        #print(f'median of mask data: {np.nanmedian(mask_data)}')
        
        #create boolean mask of valid (unmasked pixels)
        mask = ~np.isnan(mask_data)
        
        #generate x and y coordinate grids for the valid pixels
        y, x = np.indices(mask_data.shape)
        
        #flatten arrays and filter out masked pixels
        x_flat = x[mask].flatten()
        y_flat = y[mask].flatten()
        z_flat = mask_data[mask].flatten()
        
        #set up the coefficient matrix A and solve using lstsq
        A = np.column_stack((x_flat, y_flat, np.ones_like(x_flat)))
        c0, c1, c2 = np.linalg.lstsq(A, z_flat, rcond=None)[0]
        
        #create full background plane model
        background_plane = c0*x + c1*y + c2
        
        #create new background
        plane_subtracted_background_fits = mask_data - background_plane
        
        plane_subtracted_path = os.path.join(datadir, 'plane-subtracted-background-fits', f'{galaxy_name}_{suffix}_plane_subtracted_background.fits')
        fits.writeto(plane_subtracted_path, plane_subtracted_background_fits, header=mask_header, overwrite=True)
        
        #read in original image fits data and subtract plane
        image_data, image_header = fits.getdata(image_file, header=True)
        subtracted_image_data = image_data - background_plane
        #print(f'median of subtracted image: {np.nanmedian(subtracted_image_data)}')
        
        #output FITS file path
        output_fits = os.path.join(datadir, 'plane-subtracted-fits', f'{galaxy_name}_{suffix}_plane_subtracted.fits') 
        #print(f'output path: {output_fits}')
        
        #write median-subtracted fits files to output FITS file path
        fits.writeto(output_fits, subtracted_image_data, header=image_header, overwrite=True)




