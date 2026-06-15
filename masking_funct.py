#!/usr/bin/env python
# coding: utf-8

# In[2]:


#import necessary libraries

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import astropy.units as u
import pandas as pd
import numpy as np
import glob 
import sys
import os
#import wget

import warnings
warnings.filterwarnings('ignore')

from numpy.ma import is_masked

from astropy.table import Table
from astropy.io.ascii import masked
from astropy.io import ascii
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle
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

from scipy import stats
from scipy.stats import scoreatpercentile

from reproject import reproject_interp

from IPython.display import clear_output

from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture


# In[5]:


#define empty dictionary that will contain the EllipseGeometry instance
geometry = {}
initparams = {}

#initialize dictionary for half-light radii
rhalfpix = {}
rhalfasec = {}

#initalize dictionary for ellipse fitting
ellipse = {}
isolist = {}



#display an image
def imdisplay(image, x, y, sma, ellip, pa, width=100, height=100, v1perc=10, v2perc=95, logscale=True):

    '''
    OPTIONAL KEYWORD PARAMETERS
    v1perc: one end of the colormap assigned to the v1perc percent lowest flux
    v2perc: the other end of the colormap assigned to the v2perc percent highest flux
    '''

    #make sure image is an np array
    nimage = np.array(image)

    #determine the pixel values at the 10th and 95th percentile
    v1 = scoreatpercentile(nimage, v1perc)
    v2 = scoreatpercentile(nimage, v2perc)

    #display using imshow - you can play with alternate cmaps
    #the 'gray_r' color map reverses the color-scale so that dark display pixels are the brightest in the image

    #vmin and vmax set the mix and max pixel calues that will be mapped to the extremes of the colormap
    print(v1, v2)
    
    norm = None
    if (logscale):
        norm = ImageNormalize(vmin=v1, vmax=v2, stretch=LogStretch())
    else:
        norm = ImageNormalize(vmin=v1, vmax=v2)

    im = ax.imshow(image, origin='lower', norm=norm)
    aper = EllipticalAperture((x, y), sma, (ellip)*sma, pa)

    #adjust the view to focus on (x,y) without cropping the image
    x_min = max(x - width // 2, 0)
    x_max = min(x + width // 2, nimage.shape[1])
    y_min = max(y - height // 2, 0)
    y_max = min(y + height // 2, nimage.shape[0])

    #hide axis ticks and labels
    plt.axis('off')
    aper.plot(color = 'red')

    return fig, ax
    #fig.colorbar(fraction=.08)
    


#exact same as function in photometry_funct.py
def find_files(destination_folder, partial_name):
    matching_files = []

    for root, dirs, files in os.walk(destination_folder):
        for file in files:
            if partial_name.lower() in file.lower():
                matching_files.append(os.path.join(root, file))

    return matching_files



#function to overlay mask on original FITS image and save as a new FITS file
def overlay_mask_on_fits(
    fits_file,
    mask_file,
    csv_file,
    output_fits,
    n,
    coords_x,
    coords_y,
    RA,
    DEC
):

    # load the original FITS image and its WCS
    original_fits_data, original_fits_header = fits.getdata(
        fits_file,
        header=True
    )
    original_wcs = WCS(original_fits_header)

    # load the mask FITS image and its WCS
    mask_fits_data, mask_fits_header = fits.getdata(
        mask_file,
        header=True
    )
    mask_wcs = WCS(mask_fits_header)

    # if coordinates are missing, recompute them from RA/DEC
    if (
        np.ma.is_masked(coords_x)
        or np.ma.is_masked(coords_y)
        or not np.isfinite(coords_x)
        or not np.isfinite(coords_y)
    ):

        print(
            f"Coordinates missing; recomputing from RA/DEC "
            f"(RA={RA}, DEC={DEC})"
        )

        coords_x, coords_y = original_wcs.all_world2pix(
            RA,
            DEC,
            0
        )

    # if still bad, skip
    if (
        not np.isfinite(coords_x)
        or not np.isfinite(coords_y)
    ):
        print("Could not determine valid image coordinates")
        return

    # extract central pixel coordinates
    central_x = int(coords_x)
    central_y = int(coords_y)

    # reproject mask to Herschel image WCS
    reprojected_mask, footprint = reproject_interp(
        (mask_fits_data, mask_wcs),
        original_wcs,
        shape_out=original_fits_data.shape
    )

    # clip coordinates to image boundaries
    height, width = original_fits_data.shape

    central_x = np.clip(
        central_x,
        0,
        width - 1
    )

    central_y = np.clip(
        central_y,
        0,
        height - 1
    )

    # overlay mask
    combined_data = original_fits_data.copy()

    combined_data[reprojected_mask > 0] = np.max(
        original_fits_data
    )

    # restore central pixel if it was not masked
    if reprojected_mask[central_y, central_x] == 0:
        combined_data[central_y, central_x] = (
            original_fits_data[central_y, central_x]
        )

    # save result
    fits.writeto(
        output_fits,
        combined_data,
        original_fits_header,
        overwrite=True
    )

# In[ ]:




