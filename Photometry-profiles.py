#!/usr/bin/env python
# coding: utf-8

# **Tom's notes:** One of the two most important notebooks of this project. This notebook covers both calculating photometry data from the Herschel .fits files, generating images from those .fits files, as well as generating individual galaxies' HTML files that contain these images amongst other data. 
# 
# The photometry part of the program is attached alongside the generating images part because during troubleshooting this was the easiest. But, because the image generation section takes by far the most time, feel free to comment the sections out in order to run photometry more quickly.
# 
# There is also an unfinished profile graphing code block in progress. The idea is to calculate both the flux change from ellipse to ellipse of these galaxies, as well as the total flux as you get to bigger and bigger ellipses, and then graph both of these for all galaxies onto a postage stamp of 4-6 plots. 
# 
# Finally, the HTML generation outputs the galaxies and links them together, ordered by their VFID.
# 
# The first two blocks are to import the necessary packages. 

# In[9]:


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

#define an empty dictionary that will contain the EllipseGeometry instance
geometry = {}
initparams = {}

#initialize dictionary for half-light radii
rhalfpix = {}
rhalfasec = {}

#initialize dictionary for ellipse fitting
ellipse = {}
isolist = {}

#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are


#double check what the directories should be named --> 
    # 'maskdir' in Tom's version of this notebook is 'datadir' in this notebook to stay consistent with previous notebooks
    # 'datadir' in Tom's version of this notebook is 'htmldir' in this notebook

#also double check exactly the referenced html resources are


# In[10]:


#import masking_funct.py
from masking_funct import *
from photometry_funct import *


# **Tom's notes:** This is the code to both calculate the photometry AND generating .png files to visualize the galaxies. I do think it's probably best to separate them, but I never managed to get around to it. The code reads in the photometry table's data regarding namely galaxy, RA, DEC, as well as their ellipse data, namely BA, PA, SMA. The program does do photometry over all 8 SMA values, so it calculates the Herschel photometry 8 times for each image. 
# 
# Because Herschel's naming isn't straightforward, I find individual files inside the code itself, and then perform photometry on them using photutils. To disable the image generation, comment out both the v1 and v2 calculations, as well as the imdisplay2 functions. There are thesebits of code for every band, so will need to comment out 3 times.

# In[11]:


#calculating the photometry and generating ellipse images (masked)
#ellipse generation for the herschel matches
debugflag= True 
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

#make directory if doesn't exist
os.makedirs(datadir + "/mask", exist_ok=True)



for i in range(len(galaxy)):
    
    galaxy_name = str(galaxy['GALAXY'][i])
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"
    RA = galaxy['RA_MOMENT'][i]
    DEC = galaxy['DEC_MOMENT'][i]
    EPLI = galaxy['BA_MOMENT'][i]
    
    #two different position angles
    PAN_sky = (galaxy['PA_MOMENT'][i] + 90) * np.pi/180
    
    path = os.path.join(datadir, 'pipeline', galaxy_name)
    legacy_path = os.path.join('/Users/madeline.evenson/Research/HTML-building/galaxy/png/', f"{VFID}-{galaxy_name}-LS.jpg")
    
    if not os.path.exists(path):
        print('Missing pipeline folder:', path)
        continue
    
    print('Processing:', galaxy_name)
    
    for band in bands:
        
        wave = str(band['wave'])
        pixscale = band['pixscale']
        suffix = band['suffix']
        
        
        
        ###### find image file #######
        
        masked_path = os.path.join(datadir, 'masked', f'{galaxy_name}_masked{suffix}.fits')
        
        if os.path.exists(masked_path):
            found_file = masked_path
        else: 
            destination_folder = os.path.join(path, band['folder'])
            found_files = find_files(destination_folder, band['partial'])
            
            if not found_files:
                print('No file found:', galaxy_name, wave)
                continue
            
            found_file = found_files[0]
            
            
            
        ###### read image ######
        
        data, header = fits.getdata(found_file, header=True)
        wcs = WCS(header)
        
        #get CD matrix
        cd = wcs.pixel_scale_matrix
        
        #compute rotation angle (radians)
        theta = np.arctan2(-cd[0, 1], cd[1,1])
        
        #convert to pixel frame of Herschel
        PAN_herschel = PAN_sky - theta
        
        #for legacy (approx --> assume north-up)
        PAN_legacy = - PAN_sky
        
        x0, y0 = wcs.all_world2pix(RA, DEC, 0)
        x0_plot, y0_plot = x0, y0
        
        ###### display LEGACY image ######
        
        fig, ax = plt.subplots(figsize=(10, 10))
        nimage = np.array(data)
        
        
        #prevents code from crashing if Legacy image isn't found
        scale_x = 1.0
        scale_y = 1.0
    
        
        if os.path.exists(legacy_path):
            legacy_img = mpimg.imread(legacy_path)
            ax.imshow(legacy_img, origin='upper')
            
            #rescale coordinates to match image size
            ny, nx = legacy_img.shape[0], legacy_img.shape[1]
            hy, hx = nimage.shape
            
            scale_x = nx / hx
            scale_y = ny / hy
            
            x0_plot = x0 * scale_x
            y0_plot = y0 * scale_y
            
        else:
            print("Missing Legacy image:", legacy_path)
            ax.imshow(nimage, origin='lower', cmap='gray')
        

        ###### photometry ######
        
        for j, sma_label in enumerate(sma_labels):
            SMAO = galaxy[sma_label][i]
            SMAN = (SMAO * 0.262) / pixscale 
            SMAN_display = SMAN * scale_x


            # photometry aperture (Herschel space)
            aper_phot = EllipticalAperture(
                (x0, y0),
                SMAN,
                SMAN * EPLI,
                PAN_herschel
            )
            
            # plotting aperture (Legacy display space)
            aper_legacy = EllipticalAperture(
                (x0_plot, y0_plot),
                SMAN_display,
                SMAN_display * EPLI,
                PAN_legacy
            )
            
            aperture_mask = aper_phot.to_mask(method='center')
            aperture_data = aperture_mask.multiply(nimage)
            
            total_pixels = np.isfinite(aperture_data).sum()
            empty_pixels = np.isnan(aperture_data).sum()
            
            if total_pixels + empty_pixels > 0:
                good_percentage = 1 - (empty_pixels / (total_pixels + empty_pixels))
            else:
                good_percentage = 1
                
            data_nonan = np.nan_to_num(nimage, nan=0)
            
            phot_table = aperture_photometry(data_nonan, aper_phot)
            flux_value = phot_table['aperture_sum'][0]
            
            galaxy[f"{wave}Flux_{sma_label}"][i] = flux_value
            galaxy[f"{wave}GoodPixels_{sma_label}"][i] = good_percentage
            
            print(f'{galaxy_name} {wave} µm SMA{j+1}: {flux_value:.3f}')
            
            aper_legacy.plot(ax=ax, color='red', lw=1, alpha=0.7)
            
        
        ###### save image ######
        
        savepath = os.path.join(htmldir, "AP", f"{VFID}-{galaxy_name}-{wave}-AP.png")
        
        print('Saving:', savepath)
        
        plt.savefig(savepath, dpi=150)
        plt.close(fig)
        
galaxy.write(tabledir + '/Photometrytesting2.csv', format='csv', overwrite=True)

print('done')
# In[12]:


csv_file = tabledir+'/Photometrytesting2.csv'

galaxy = Table.read(csv_file)
galaxy

# **Tom's notes:** This is theunfinished profile graphing code block. 
#The desired output is a series of 4 or 6 plots depending on number of Herschel bands a particular galaxy has. 
#One set of plots is the graph of the ellipse number vs total flux contained within that ellipse for a given galaxy in a given band. 
#The other set of plots should plot the ellipse number vs the difference in flux of the respective ellipses (effectively the rate of change. 
# was having trouble with setting up the rate of change correctly. 

# In[13]:


#profile graphing (unfinished)

#load the photometry data
csv_file = tabledir + '/Photometrytesting2.csv'
galaxy = Table.read(csv_file)

#define semi-major axis (SMA) and flux labels
sma_labels = [f'SMA_AP0{i}' for i in range(1, 9)]
wavelengths = [70, 100, 160]
colors = ['blue', 'green', 'red'] #corresponding colors for wavelengths


#loop over each galaxy
for i in range(len(galaxy)):
    galaxy_name = str(galaxy['GALAXY'][i])
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"
    EPLI = galaxy['BA_MOMENT'][i]

    #determine which bands have valid data
    valid_bands = []


    for j, wave in enumerate(wavelengths):
        has_valid = False
    
        for sma_label in sma_labels:
            flux_col = f'{wave}Flux_{sma_label}'
            if flux_col in galaxy.colnames:
                val = galaxy[flux_col][i]
                if not np.isnan(val):
                    has_valid = True
                    break
            
        if has_valid:
            valid_bands.append((wave, colors[j]))

    #skip galaxy if no valid bands:
    if len(valid_bands) == 0:
        continue


    #create dynamic subplot grid
    ncols = len(valid_bands)
    fig, axes = plt.subplots(2, ncols, figsize = (5*ncols, 10))


    #handle case where ncols==1
    if ncols == 1:
        axes = np.array([[axes[0]], [axes[1]]])
    
    
    #loop over only valid bands
    for col_idx, (wave, color) in enumerate(valid_bands):
    
        sma_values = []
        sma_values_2 = []
        flux_values = []
        sb_values = []
    
        for k in range(len(sma_labels) - 1):
            sma_label_1 = sma_labels[k]
            sma_label_2 = sma_labels[k+1]
        
            flux_col_1 = f'{wave}Flux_{sma_label_1}'
            flux_col_2 = f'{wave}Flux_{sma_label_2}'
        
            if flux_col_1 in galaxy.colnames and flux_col_2 in galaxy.colnames:
                flux_1 = galaxy[flux_col_1][i]
                flux_2 = galaxy[flux_col_2][i]
            
                if np.isnan(flux_1) or np.isnan(flux_2):
                    continue
            
                sma_pixels_1 = galaxy[sma_label_1][i]
                sma_pixels_2 = galaxy[sma_label_2][i]
            
                sma_arcsec_1 = sma_pixels_1 * 0.262 #converting to arcsec
                sma_arcsec_2 = sma_pixels_2 * 0.262 
            
                sma_values.append(sma_arcsec_1)
                sma_values_2.append(sma_arcsec_2)
            
                #scale by pixel size
                if wave in [70, 100]:
                    scale = 1.60000001784
                else:
                    scale = 3.20000003568
                
                flux_1 *= scale
                flux_2 *= scale
            
                flux_values.append(flux_2)
            
                ### old area_annulus --___> area_annulus = np.pi * (sma_arcsec_2**2 - sma_arcsec_1**2)
                
                #fix area_annulus to calculate area of **ellipse**
                a1 = sma_arcsec_1
                a2 = sma_arcsec_2
                
                b1 = a1 * EPLI
                b2 = a2 * EPLI
                
                area_annulus = np.pi * (a2*b2 - a1*b1)
                
                
                if area_annulus > 0:
                    sb_values.append((flux_2 - flux_1) / area_annulus)
                
                
        ##### plotting #####
        if sma_values:
            axes[0, col_idx].plot(sma_values_2, flux_values, marker='o', color=color)
            axes[0, col_idx].set_title(f'{wave} µm Flux Profile')
            axes[0, col_idx].set_xlabel('Semi-Major Axis (arcsec)')
            axes[0, col_idx].set_ylabel('Flux')
        
        if sb_values:
            mid_radii = 0.5 * (np.array(sma_values) + np.array(sma_values_2))
            axes[1, col_idx].plot(mid_radii, sb_values, marker='s', color=color)
            #axes[1, col_idx].plot(sma_values, sb_values, marker='s', color=color)
            axes[1, col_idx].set_title(f'{wave} µm Surface Brightness')
            axes[1, col_idx].set_xlabel('Semi-Major Axis (arcsec)')
            axes[1, col_idx].set_ylabel('Surface Brightness (Jy/arcsec²)')
        
    #adjust layout
    plt.tight_layout()

    #save figure
    output_path = f'{htmldir}/profiles/{VFID}-{galaxy_name}_flux_sb_profile.png'
    plt.savefig(output_path, dpi=150)
    plt.close(fig)

print('All flux and surface brightness profiles have been generated.')

print(galaxy.colnames)

# %%
