#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 17:40:18 2026

@author: madeline.evenson
"""

#the point of this notebook is to make the background comparison plots for the html website

import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval

from masking_funct import *


os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')

tabledir = homedir + '/Virgo/tables/'
htmldir = homedir + '/Virgo/HTML-building/galaxy/'
datadir = homedir + '/masking/'

#original Herschel pipeline directory
pipeline_dir = os.path.join(datadir, 'pipeline')

#plane-subtracted background FITS
    #using plane_subtracted.fits instead of plane-subtracted-background-fits because don't want to show masked ellipses
plane_subtracted_fits = os.path.join(
    datadir,
    'plane-subtracted-fits'
)

#output directory
output_dir = os.path.join(htmldir, 'bg-comparison-plots')
os.makedirs(output_dir, exist_ok=True)



# %%


#load photometry data

csv_file = os.path.join(tabledir, 'Photometrytesting2.csv')
galaxy = Table.read(csv_file)

sma_labels = [f'SMA_AP0{i}' for i in range(1, 9)]

#create dictionary for bands --> includes prefixes for /masking/pipeline/...

bands = {
    70: {
        "suffix": "B",
        "prefix": "hpacs_25HPPUNIMAPB"
    },
    100: {
        "suffix": "G",
        "prefix": "hpacs_25HPPUNIMAPB"  #this is supposed to be like this
    },
    160: {
        "suffix": "R",
        "prefix": "hpacs_25HPPUNIMAPR"
    }
}

mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


#initialize zscale to make sure everything plotted in Zscale

interval = ZScaleInterval()



# %%

#loop over all the galaxies


for i in range(len(galaxy)):

    galaxy_name = str(galaxy['GALAXY'][i])
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"

    print(f"\nProcessing {VFID} : {galaxy_name}")


    #first determine which bands have valid data

    valid_bands = []

    for j, wave in enumerate([70, 100, 160]):

        has_valid = False

        for sma_label in sma_labels:

            flux_col = f"{wave}Flux_{sma_label}"

            if flux_col in galaxy.colnames:

                value = galaxy[flux_col][i]

                if not np.isnan(value):
                    has_valid = True
                    break

        if has_valid:
            valid_bands.append((wave, mycolors[j]))

    if len(valid_bands) == 0:
        print("   No valid bands.")
        continue

    

    #create figure with dynamic number of columns, depending on valid_bands

    ncols = len(valid_bands)

    fig, axes = plt.subplots(
        2,
        ncols,
        figsize=(5 * ncols, 10),
        squeeze=False
    )

    
    #plot each band that has valid data

    for col_idx, (wave, color) in enumerate(valid_bands):

        suffix = bands[wave]["suffix"]
        prefix = bands[wave]["prefix"]

        

        #locate original data fits from /masking/pipeline/{galaxy_name}...
        band_dir = os.path.join(
            pipeline_dir,
            galaxy_name,
            f"HPPUNIMAP{suffix}"
        )

        if not os.path.isdir(band_dir):
            print(f"   Missing directory: {band_dir}")
            continue

        fits_files = sorted([
            f for f in os.listdir(band_dir)
            if f.startswith(prefix)
        ])

        if len(fits_files) == 0:
            print(f"   No science FITS found in {band_dir}")
            continue

        original_path = os.path.join(
            band_dir,
            fits_files[0]
        )

        
        #locate plane-subtracted fits
        psub_path = os.path.join(
            plane_subtracted_fits,
            f"{galaxy_name}_{suffix}_plane_subtracted.fits"
        )

        if not os.path.exists(psub_path):
            print(f"   Missing: {psub_path}")
            continue

        
        #read in both sets of images
        original_data = fits.getdata(original_path)
        psub_data = fits.getdata(psub_path)

        #calculate their respective images ranges
        vmin1, vmax1 = interval.get_limits(original_data)
        vmin2, vmax2 = interval.get_limits(psub_data)

        

        #plot original fits image
        im1 = axes[0, col_idx].imshow(
            original_data,
            origin='lower',
            cmap='viridis',
            vmin=vmin1,
            vmax=vmax1
        )

        axes[0, col_idx].set_title(f'{wave} μm Original')

        fig.colorbar(
            im1,
            ax=axes[0, col_idx],
            fraction=0.046,
            pad=0.04
        )

        
        #plot plane subtracted fits image
        im2 = axes[1, col_idx].imshow(
            psub_data,
            origin='lower',
            cmap='viridis',
            vmin=vmin2,
            vmax=vmax2
        )

        axes[1, col_idx].set_title(
            f'{wave} μm Plane-Subtracted'
        )

        fig.colorbar(
            im2,
            ax=axes[1, col_idx],
            fraction=0.046,
            pad=0.04
        )

    #create suptitle
    fig.suptitle(
        f'{VFID} — {galaxy_name}',
        fontsize=18,
        fontweight='bold'
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    output_path = os.path.join(
        output_dir,
        f'{VFID}-{galaxy_name}_bg_comparison.png'
    )

    fig.savefig(output_path, dpi=150)
    plt.close(fig)

print("\nAll background comparison plots have been generated.")





