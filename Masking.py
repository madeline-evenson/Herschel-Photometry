#!/usr/bin/env python
# coding: utf-8

# **Tom's notes:** In order to measure photometry in the Herschel bands in the same areas of the other 
#photometry measurements, we will also need to apply the same masks that are applied within the other 
#measurements. For the most parts these masks mask out problematic stars and things that are not part 
#of the galaxies, but even if there are no star emissions within the Herschel bands, we still want to 
#keep a consistent area of measurements within the ellipses. This notebook is purely to get the masking 
#fits and impose it onto the Herschel images correctly, as well as output those as FITS. It does not 
#calculate the photometry of these galaxies.

# In[12]:


#import necessary libraries

import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import pandas as pd
import os
#import wget

from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits


import warnings
warnings.filterwarnings('ignore')

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
datadir = homedir+'/masking/' #set to where to find the masking data and images -> all the contents of the masking here should be in the zip file


# In[13]:


#import masking_funct.py
from masking_funct import *


# **Tom's notes:** This next codeblock is used to find the central pixels that corresponds to the galaxy's 
    #center. We need to do this in order to impose the masks on correctly, because the masking 
    #files always have the galaxy's center as it central pixel, and so we want to impose those 
    #masks by lining up the central pixels on top of the galaxy's center pixel in the Herschel 
    #images. The program does it by first reading in the respective .fits images, taking the galaxy's 
    #RA and DEC and converting it into pixels, and then recording the pixel value into the spreadsheet. 
    #There is also a check on the WISE mask images to make sure the pixel coordinates are correct 
    #there too. 

# In[14]:


#finding the central pixels of the galaxies

#input CSV file with galaxy data
csv_file = tabledir + 'Photometrytesting.csv'
galaxy_data = pd.read_csv(csv_file)
print(galaxy_data.columns)

#prepare output CSV
output_data = []

#process each galaxy
for i in range(len(galaxy_data)):
    galaxy_name = str(galaxy_data['GALAXY'][i])
    RA = galaxy_data['RA_MOMENT'][i]
    DEC = galaxy_data['DEC_MOMENT'][i]
    VFID = f"VFID{int(galaxy_data['VF_ID'][i]):04d}"

    #paths for the images
    base_path = os.path.join(datadir, 'pipeline', galaxy_name)
    mask_path = os.path.join(datadir, 'masks')

    #search for HPPUNIMAPB FITS file dynamically
    destination_folder = os.path.join(datadir, 'pipeline', galaxy_name, 'HPPUNIMAPB')
    partial_name = 'hpacs_25HPPUNIMAPB'
    found_files = find_files(destination_folder, partial_name)
    hppunimapb_image_path = found_files[0] if found_files else None

    #search for HPPUNIMAPG FITS file dynamically
    destination_folder = os.path.join(datadir, 'pipeline', galaxy_name, 'HPPUNIMAPG')
    partial_name = 'hpacs_25HPPUNIMAPB' #checked with Kim and Rudnick --> supposed to be this way
    found_files = find_files(destination_folder, partial_name)
    hppunimapg_image_path = found_files[0] if found_files else None

    #search for HPPUNIMAPR FITS file dynamically
    destination_folder = os.path.join(datadir, 'pipeline', galaxy_name, 'HPPUNIMAPR')
    partial_name = 'hpacs_25HPPUNIMAPR'
    found_files = find_files(destination_folder, partial_name)
    hppunimapr_image_path = found_files[0] if found_files else None

    #check and process WISE mask image
    wise_mask_image_path = os.path.join(mask_path, f'{galaxy_name}-custom-image-wise-mask.fits')
    wise_mask_x, wise_mask_y = None, None
    
    if os.path.exists(wise_mask_image_path):
        wise_mask_data, wise_mask_header = fits.getdata(wise_mask_image_path, header=True)
        wise_mask_wcs = WCS(wise_mask_header)
        wise_mask_x, wise_mask_y = wise_mask_wcs.all_world2pix(RA, DEC, 0)
    else:
        #if the original WISE mask is not found, check for the alternative '-custom-image-r-mask.fits'
        wise_mask_image_path_r = os.path.join(mask_path, f'{galaxy_name}-custom-image-r-mask.fits')

        if os.path.exists(wise_mask_image_path_r):
            wise_mask_data, wise_mask_header = fits.getdata(wise_mask_image_path_r, header=True)
            wise_mask_wcs = WCS(wise_mask_header)
            wise_mask_x, wise_mask_y = wise_mask_wcs.all_world2pix(RA, DEC, 0)
        else:
            #if neither WISE mask if found, set the pixel coordinates to None
            wise_mask_x, wise_mask_y = None, None

    #initialize pixel coordiantes
    hppunimapb_x, hppunimapb_y = None, None
    hppunimapg_x, hppunimapg_y = None, None
    hppunimapr_x, hppunimapr_y = None, None

    #check and process HPPUNIMAPB image
    if hppunimapb_image_path and os.path.exists(hppunimapb_image_path):
        hppunimapb_data, hppunimapb_header = fits.getdata(hppunimapb_image_path, header=True)
        hppunimapb_wcs = WCS(hppunimapb_header)
        hppunimapb_x, hppunimapb_y = hppunimapb_wcs.all_world2pix(RA, DEC, 0)
        
    #check and process HPPUNIMAPG image
    if hppunimapg_image_path and os.path.exists(hppunimapg_image_path):
        hppunimapg_data, hppunimapg_header = fits.getdata(hppunimapg_image_path, header=True)
        hppunimapg_wcs = WCS(hppunimapg_header)
        hppunimapg_x, hppunimapg_y = hppunimapg_wcs.all_world2pix(RA, DEC, 0)

    #check and process HPPUNIMAPR image
    if hppunimapr_image_path and os.path.exists(hppunimapr_image_path):
        hppunimapr_data, hppunimapr_header = fits.getdata(hppunimapr_image_path, header=True)
        hppunimapr_wcs = WCS(hppunimapr_header)
        hppunimapr_x, hppunimapr_y = hppunimapr_wcs.all_world2pix(RA, DEC, 0)

    #append the results to the output data
    output_data.append({
        'VFID': VFID,
        'Galaxy': galaxy_name,
        'RA': RA,
        'DEC': DEC,
        'HPPUNIMAPB_X': hppunimapb_x,
        'HPPUNIMAPB_Y': hppunimapb_y,
        'HPPUNIMAPG_X': hppunimapg_x,
        'HPPUNIMAPG_Y': hppunimapg_y,
        'HPPUNIMAPR_X': hppunimapr_x,
        'HPPUNIMAPR_Y': hppunimapr_y,
        'WISE_Mask_X': wise_mask_x,
        'WISE_Mask_Y': wise_mask_y
    })

#convert to DataFrame and save to CSV
output_df = pd.DataFrame(output_data)
output_csv_path = datadir + '/pixel_coordinates.csv'
output_df.to_csv(output_csv_path, index=False)
print('done')


# In[15]:


#input CSV file with galaxy data
csv_file = datadir + 'pixel_coordinates.csv'
pixel_data = pd.read_csv(csv_file)
pixel_data


# **Tom's note:** Once we have the coordinates of the central pixels, we can then impose the masks onto the 
    #Herschel images using the following code, which takes respective galaxies and their masks by cross-checking 
    #the galaxy's names, scales the masks according to the pixel scales of the Herschel images, and imposes them 
    #onto the image. The program will then re-output the fits images in a seperate mask folder. 

# In[16]:


#imposing the masks onto the images by lining up the cental mask pixel with the central galaxy pixels calculated from the previous code block

#file paths
csv_file = os.path.join(datadir, 'pixel_coordinates.csv')
galaxy = Table.read(csv_file)

for i in range(len(galaxy)):
    galaxy_name = str(galaxy['Galaxy'][i])
    path = os.path.join(datadir, 'pipeline', galaxy_name)
    VFID = str(galaxy['VFID'][i])

    if os.path.exists(path):
        for color in ['B', 'G', 'R']: #loop over the three color channels
            destination_folder = os.path.join(path, f'HPPUNIMAP{color}')
            if color in ['B', 'G']: #blue and green bands
                partial_name = f'hpacs_25HPPUNIMAPB'
            else: #red band
                partial_name = f'hpacs_25HPPUNIMAPR'
            found_files = find_files(destination_folder, partial_name)

            if found_files:
                found_file = found_files[0]

                #mask file paths
                wise_mask_file = os.path.join(datadir, 'masks', f'{galaxy_name}-custom-image-wise-mask.fits')
                r_mask_file = os.path.join(datadir, 'masks', f'{galaxy_name}-custom-image-r-mask.fits')

                #output FITS file path
                output_fits = os.path.join(datadir, 'masked', f'{galaxy_name}_masked{color}.fits') #where the masked file should go

                #read the CSV file for the central pixel coordinates
                coords_x = galaxy[f'HPPUNIMAP{color}_X'][i]
                coords_y = galaxy[f'HPPUNIMAP{color}_Y'][i]

                #check if the wise-mask file exists, otherwise check for the r-mask
                if os.path.exists(wise_mask_file):
                    #if wise-mask exists, overlay it
                    overlay_mask_on_fits(found_file, wise_mask_file, csv_file, output_fits, i, coords_x, coords_y)
                elif os.path.exists(r_mask_file):
                    #if r-mask exists, overlay it
                    overlay_mask_on_fits(found_file, r_mask_file, csv_file, output_fits, i, coords_x, coords_y)
                else:
                    #if no mask files exist, copy the original FITS file to the output location
                    if not os.path.exists(os.path.dirname(output_fits)):
                        os.makedirs(os.path.dirname(output_fits))
                    with fits.open(found_file) as hdul:
                        hdul.writeto(output_fits, overwrite=True)
                    print(f'Mask file not found for {galaxy_name} or {VFID}. Saved original FITS as {output_fits}.')


# In[17]:


galaxy


# In[ ]:




