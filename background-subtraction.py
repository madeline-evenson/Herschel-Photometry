#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:54:36 2026

@author: madeline.evenson
"""

#loading in FITS file for NGC 5713

from astropy.io import fits
import os

os.environ['HOME'] = '/Users/madeline.evenson/Research' #general path for all the coding space
homedir = os.getenv('HOME')
fitsdir = homedir+'/masking/masks/'

file = fitsdir + 'NGC5713-custom-image-wise-mask.fits'

#open FITS file
with fits.open(file) as hdul:
    
    #print summary of file
    hdul.info()
    
    #extract header metadata and pixel/table data
    header = hdul[0].header
    data = hdul[0].data