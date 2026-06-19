#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 14:14:47 2026

@author: madeline.evenson
"""


#the purpose of this notebook is to determine whether a plane is a good model for the background of the Herschel images

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval

from masking_funct import *

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

#directory for median-subtracted fits
median_subtracted_fits = datadir + '/median-subtracted-fits/'

#directory for plane-subtracted background fits
plane_subtracted_background_fits = datadir + '/plane-subtracted-background-fits'

#directory for plane-subtracted image fits
plane_subtracted_fits = datadir + '/plane-subtracted-fits/'

#directory for plane-subtracted background flux histograms
psub_bg_flux_hists = datadir + '/psub-bg-flux-hists/'
    #psub = plane subtracted
    #bg = background
    
#directory for horizontal and vertical cut flux profiles
psub_bg_h_cut_profiles = datadir + '/psub-bg-h-cut-profiles'
psub_bg_v_cut_profiles = datadir + '/psub-bg-v-cut-profiles'


# %%

#first create flux histogram for one example galaxy --> NGC5290

file = plane_subtracted_background_fits + '/NGC5290_R_plane_subtracted_background.fits'

data, header = fits.getdata(file, header=True)
#print(header)

plt.figure(figsize=(6,6))
plt.imshow(data, cmap='viridis', origin='lower')
plt.colorbar()


fig, ax = plt.subplots(figsize=(6,6))
ax.hist(data, bins=20, edgecolor='black')
ax.set_title(f'NGC5290 R Background Fluxes')
ax.set_ylabel('Counts')
ax.set_xlabel('Flux (Jy)')


# %%

#now create directory of background flux histograms

csv_file = tabledir+'/Photometrytesting.csv'
#csv_file2 = tabledir+'/Herschelstuff.csv'

galaxy = Table.read(csv_file)
#galaxy2 = Table.read(csv_file2)

galaxy_length = len(galaxy)

sma_labels = [f'SMA_AP0{i}' for i in range(1, 9)]

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
        #print('Missing pipeline folder:', path)
        continue
    
    print('Processing:', galaxy_name)
    
    for band in bands:
        
        suffix = band['suffix']
        
        
        ###### find plane-subtracted background file #######
        background_path = os.path.join(plane_subtracted_background_fits, f'{galaxy_name}_{suffix}_plane_subtracted_background.fits')
        
        if not os.path.exists(background_path):
            print('No file found:', galaxy_name, suffix)
            continue

            
        #read in new mask data and calculate median   
        background_data, background_header = fits.getdata(background_path, header=True)
        
        
        fig, ax = plt.subplots(figsize=(6,6))
        ax.hist(background_data.ravel(), bins=30, edgecolor='black')
        ax.set_xlim(-0.07, 0.07)
        ax.set_title(f'{galaxy_name} {suffix} Background Fluxes')
        ax.set_ylabel('Counts')
        ax.set_xlabel('Flux (Jy)')
        
        output_path = psub_bg_flux_hists + f'{galaxy_name}_{suffix}_psub_bg_flux_hist.png'
        fig.savefig(output_path, dpi=100)
        plt.close(fig)
        
        
# %%

#now we need to make horizontal and vertial cuts along plane-subtracted background fits to 
    #see if any trends in residual
    #like usual first start with one singular galaxy in one band and then apply to all galaxies
    
    
#example galaxy: ICO902_G --> one of the ones that looked weird in the psub_bg_flux_hists

file = plane_subtracted_background_fits + '/IC0902_G_plane_subtracted_background.fits'

data, header = fits.getdata(file, header=True)

median = np.nanmedian(data)
print(f'median: {median}')
print(data.shape)
#plt.figure(figsize=(6,6))
#plt.imshow(data, cmap='viridis', origin='lower')
#plt.colorbar()


#make horizontal cut --> extracting row index 50, all columns
    # FITS uses (y,x) indexing!!
    
h_slice_50 = data[50, :]
h_slice_150 = data[150, :]
h_slice_250 = data[250, :]

plt.figure(figsize=(8,6))
plt.plot(h_slice_50, label = 'Row 50')
plt.plot(h_slice_150, label='Row 150')
plt.plot(h_slice_250, label='Row 250')

plt.title('IC0902 G Horizontal Cut Flux Profiles')
plt.xlabel('Pixel Column (x-axis)')
plt.ylabel('Flux (Jy)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()

plt.show()

# %%

#now create horizontal and vertical cut flux profiles for every galaxy in every band
#want to choose three h/v cuts from each galaxy that aren't completely full of nans and not on edges of image

for i in range(len(galaxy)):

    galaxy_name = str(galaxy['GALAXY'][i])
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"

    path = os.path.join(datadir, 'pipeline', galaxy_name)

    if not os.path.exists(path):
        continue

    print(f'Processing: {galaxy_name}')

    for band in bands:

        suffix = band['suffix']

        background_path = os.path.join(
            plane_subtracted_background_fits,
            f'{galaxy_name}_{suffix}_plane_subtracted_background.fits'
        )

        if not os.path.exists(background_path):
            print(f'No file found: {galaxy_name} {suffix}')
            continue

        
        #read in image and determine shape
        data, header = fits.getdata(background_path, header=True)

        ny, nx = data.shape

        #skip completely empty images
        if np.sum(np.isfinite(data)) < 100:
            print(f'   Skipping {galaxy_name} {suffix}: almost no valid pixels')
            continue

        
        #choose rows and columns near 25%, 50%, and 75%
        #avoid the outer 15% of the image

        edge_frac = 0.15

        #find rows with at least 20 finite pixels
        valid_rows = np.where(np.sum(np.isfinite(data), axis=1) > 0.5 * nx)[0]
    
        #keep only rows away from the edges
        valid_rows = valid_rows[
            (valid_rows > edge_frac * ny) &
            (valid_rows < (1 - edge_frac) * ny)
            ]

        #find columns with at least 20 finite pixels
        valid_cols = np.where(np.sum(np.isfinite(data), axis=0) > 0.5 * ny)[0]

        #keep only columns away from edges
        valid_cols = valid_cols[
            (valid_cols > edge_frac * nx) &
            (valid_cols < (1 - edge_frac) * nx)
            ]

        #if there are no "valid" rows or columns then continue and print message (tab included to be seen easier when running)
        if len(valid_rows) == 0 or len(valid_cols) == 0:
            print(f"   No good rows/columns for {galaxy_name} {suffix}")
            continue

        #create arrays near 25%, 50%, and 75% width/height of image
        target_rows = np.array([0.25*ny, 0.50*ny, 0.75*ny])
        target_cols = np.array([0.25*nx, 0.50*nx, 0.75*nx])

        #find nearest valid row/column to each target %
        rows = [valid_rows[np.argmin(np.abs(valid_rows - r))] for r in target_rows]
        cols = [valid_cols[np.argmin(np.abs(valid_cols - c))] for c in target_cols]
        
        #remove duplicates if two targets map to the same row --> unlikely, but just in case!!
        rows = np.unique(rows)
        cols = np.unique(cols)

        #if duplicates occurred, fill in with additional valid rows/cols --> again unlikely, but just in case!!
        while len(rows) < 3:
            remaining = np.setdiff1d(valid_rows, rows)
            if len(remaining) == 0:
                break
            rows = np.append(rows, remaining[len(remaining)//2])

        while len(cols) < 3:
            remaining = np.setdiff1d(valid_cols, cols)
            if len(remaining) == 0:
                break
            cols = np.append(cols, remaining[len(remaining)//2])

        rows = np.sort(rows.astype(int))
        cols = np.sort(cols.astype(int))

        


        #make the horizontal cuts and plot with lines of best fit

        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

        x = np.arange(nx)

        for ax, row in zip(axes, rows):

            profile = data[row, :]
            good = np.isfinite(profile)

            ax.plot(x, profile, color='C0', label='Data')

            #only find fit if there are more than 10 finite points on the profile
            if np.sum(good) >= 10:

                try:
                    #find line of best fit
                    coeffs = np.polyfit(x[good], profile[good], 1)
                    fit = np.polyval(coeffs, x)

                    ax.plot(
                        x,
                        fit,
                        'r--',
                        linewidth=2,
                        label=f'Slope = {coeffs[0]:.2e}' #label the LOBT on the image
                    )

                #make sure if fit fails the code doesn't crash
                except np.linalg.LinAlgError:
                    ax.text(
                        0.5,
                        0.5,
                        'Fit failed',
                        transform=ax.transAxes,
                        ha='center'
                    )

            else:
                ax.text(
                    0.5,
                    0.5,
                    'Too few valid pixels',
                    transform=ax.transAxes,
                    ha='center'
                )

            ax.set_title(f'Row {row}')
            ax.set_ylabel('Flux (Jy)')
            ax.grid(alpha=0.4)
            ax.legend(fontsize=8)

        axes[-1].set_xlabel('Pixel Column')

        fig.suptitle(f'{galaxy_name} {suffix} Horizontal Cuts')

        plt.tight_layout()

        #save figure to directory
        output_path = os.path.join(
            psub_bg_h_cut_profiles,
            f'{galaxy_name}_{suffix}_horizontal_profiles.png'
        )

        fig.savefig(output_path, dpi=100)
        plt.close(fig)



        
        #make vertical cuts and plot/save fits
        #code the exact same as above except replacing all the x's with y's

        fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

        y = np.arange(ny)

        for ax, col in zip(axes, cols):

            profile = data[:, col]
            good = np.isfinite(profile)

            ax.plot(y, profile, color='C0', label='Data')

            if np.sum(good) >= 10:

                try:
                    coeffs = np.polyfit(y[good], profile[good], 1)
                    fit = np.polyval(coeffs, y)

                    ax.plot(
                        y,
                        fit,
                        'r--',
                        linewidth=2,
                        label=f'Slope = {coeffs[0]:.2e}'
                    )

                except np.linalg.LinAlgError:
                    ax.text(
                        0.5,
                        0.5,
                        'Fit failed',
                        transform=ax.transAxes,
                        ha='center'
                    )

            else:
                ax.text(
                    0.5,
                    0.5,
                    'Too few valid pixels',
                    transform=ax.transAxes,
                    ha='center'
                )

            ax.set_title(f'Column {col}')
            ax.set_ylabel('Flux (Jy)')
            ax.grid(alpha=0.4)
            ax.legend(fontsize=8)

        axes[-1].set_xlabel('Pixel Row')

        fig.suptitle(f'{galaxy_name} {suffix} Vertical Cuts')

        plt.tight_layout()

        output_path = os.path.join(
            psub_bg_v_cut_profiles,
            f'{galaxy_name}_{suffix}_vertical_profiles.png'
        )

        fig.savefig(output_path, dpi=100)
        plt.close(fig)




