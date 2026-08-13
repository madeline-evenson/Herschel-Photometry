#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:29:27 2026

@author: madeline.evenson
"""

#one of the main tasks of the project is generating an HTML website with the purpose of 
#having a visual catalog of the sample of Herschel galaxies
#we need this to more easily reference each galaxy, as well as visually inspect its 
#various properties like image quality, ellipse generation, and profile trends


# %%
#import the necessary packages (some might be unnecessary)

from astropy.table import Table
from matplotlib import pyplot as plt
import os
import numpy as np
import wget
from astropy.io import fits
from astropy.wcs import WCS
import warnings
import time
from urllib.error import ContentTooShortError

warnings.filterwarnings('ignore')

mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

#set the home path


os.environ['HOME'] = '/Users/madeline.evenson/Research' #general path for all the coding space
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are

from masking_funct import *
from photometry_funct import *


# %%

#this next code block generates cutout pngs from the Legacy Survey by first finding if there is 
#a galaxy .fits file inside the specified folders of galaxies
#then it extracts the RA and DEC of those galaxies
#using these two points, it searches the Legacy Survey to find an image of the galaxy and 
#saves it locally

# %%

#get optical cutout pngs from the Legacy Survey  
csv_file = tabledir + 'Photometrytesting2.csv'
galaxy = Table.read(csv_file)
pixscale = 1 #standard
pscale = 0
xsize = 0
max_retries = 5

for i in range(len(galaxy)):
    galaxy_name = str(galaxy['GALAXY'][i])
    path = datadir + '/pipeline/' + galaxy_name
    
    destination_folder = path + '/HPPUNIMAPR/'
    partial_name = 'hpacs_25HPPUNIMAPR'
    found_files = find_files(destination_folder, partial_name)
    
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"
    filename_LS = os.path.join(homedir, "HTML-building/galaxy/png", f"{VFID}-{galaxy_name}-LS.jpg")
    
    os.makedirs(os.path.join(homedir, "HTML-building/galaxy/png"), exist_ok=True)
    
    RA = galaxy['RA_MOMENT'][i]
    DEC = galaxy['DEC_MOMENT'][i]
    
    if os.path.exists(path):
        if found_files:
            found_file = found_files[0]
            image, head = fits.getdata(found_file, header=True)
            wcs_info = WCS(head)
            
            pscale = np.abs(float(head['CDELT1'])) #grab transformation matrix of Herschel image
            xsize = np.abs(float(head['NAXIS1'])) #grab length of Herschel image
            xsize_arcsec = pscale*3600*xsize #length convert to arcseconds
            imsize = int(xsize_arcsec / pixscale) #convert length to an integer
            imsize = str(imsize) #convert integer length to a string
            
            image_url = f'https://www.legacysurvey.org/viewer/cutout.jpg?ra={RA}&dec={DEC}&layer=ls-dr9&size={imsize}&pixscale={1}'
            print(image_url)
            print(f'filename_LS = {filename_LS:s}')
            
            for attempt in range(max_retries):
                try:
                    if os.path.exists(filename_LS):
                        os.remove(filename_LS)
                        
                    imageLS = wget.download(image_url, out=filename_LS)
                    print(f"\nDownloaded successfully: {filename_LS}")
                    break #success --> exit retry loop
                    
                except ContentTooShortError as e:
                    print(f"\nDownload failed (attempt {attempt+1}/{max_retries}): {e}")
                    time.sleep(2) #wait before retrying
                    
                except Exception as e:
                    print(f"\nUnexpected error: {e}")
                    break
            else:
                print(f"\nFAILED after {max_retries} attempts: {filename_LS}")
            #if os.path.exists(filename_LS):
                #os.remove(filename_LS)
                #imageLS = wget.download(image_url, out = filename_LS)
            #else:
                #imageLS = wget.download(image_url, out = filename_LS)
                
    #else:
       # print(f'Galaxy not found: {VFID}-{galaxy_name}')
            
print('done')

