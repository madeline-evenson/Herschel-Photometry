#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 13:37:15 2026

@author: madeline.evenson

"""


#this notebook is an edited version of Tom's original Masking.py notebook
#it is mostly the same as original, except it adds an additional ellipse mask onto the WISE mask
#the additional ellipse mask should be 3 times larger than the original central ellipse mask
#anything labeled 'new' involves the new masks unique to this notebook (not in Tom's original .py notebook)



import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from photutils.aperture import EllipticalAperture

from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.wcs.utils import proj_plane_pixel_scales

#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are

#make directory if doesn't exist
os.makedirs(plotdir, exist_ok=True)

os.makedirs(os.path.join(datadir, 'new_masked'), exist_ok=True)

# directory for visualizations of the new masks
new_mask_plot_dir = os.path.join(datadir, 'new_mask_images')
os.makedirs(new_mask_plot_dir, exist_ok=True)


# directory for visualizations of the new masks on herschel images
new_masked_herschel_dir = os.path.join(datadir, 'new_masked_herschel')
os.makedirs(new_masked_herschel_dir, exist_ok=True)

from masking_funct import *


# %%

astropy_table = Table.read(tabledir + 'vf_v2_legacy_ephot.fits', format='fits')

#convert astropy table directly to pandas DataFrame
galaxy_data = astropy_table.to_pandas()

   
    
#prepare output CSV
output_data = []


#process each galaxy
for i in range(len(galaxy_data)):
    
    galaxy_name = str(galaxy_data['GALAXY'][i])
    VFID = f"VFID{int(galaxy_data['VF_ID'][i]):04d}"
    
    SMA = galaxy_data['SMA_MOMENT'][i]
    RA = galaxy_data['RA_MOMENT'][i]
    DEC = galaxy_data['DEC_MOMENT'][i]
    EPLI = galaxy_data['BA_MOMENT'][i] #moment ratio --> a/b
    
    #convert position angle into radians
    PAN = (galaxy_data['PA_MOMENT'][i] + 90) * np.pi/180
    
    #calculate semi-minor axis from moment ratio
    SMB = SMA * EPLI

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


    #clean galaxy_name part of csv to avoid bytes object string conversion
    gal_name = galaxy_data['GALAXY'][i]
    
    if isinstance(gal_name, bytes):
        gal_name = gal_name.decode('utf-8')
        
    galaxy_name = str(gal_name)
    
    
    #append the results to the output data
    output_data.append({
        'VFID': VFID,
        'GALAXY': galaxy_name,
        'SMA': SMA,
        'RA': RA,
        'DEC': DEC,
        'EPLI': EPLI,
        'PAN': PAN,
        'SMB': SMB,
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
output_csv_path = datadir + '/new_pixel_coordinates.csv'
output_df.to_csv(output_csv_path, index=False)
print('done')


# In[15]:


#input CSV file with galaxy data
csv_file = datadir + 'new_pixel_coordinates.csv'
pixel_data = pd.read_csv(csv_file)
print(pixel_data.columns)


# **Tom's note:** Once we have the coordinates of the central pixels, we can then impose the masks onto the 
    #Herschel images using the following code, which takes respective galaxies and their masks by cross-checking 
    #the galaxy's names, scales the masks according to the pixel scales of the Herschel images, and imposes them 
    #onto the image. The program will then re-output the fits images in a seperate mask folder. 


# %%

def convert_to_pixels(value, header):
    wcs = WCS(header)
    pixscale_arcsec = np.mean(proj_plane_pixel_scales(wcs)) * 3600.0
    return value / pixscale_arcsec


# %%

def create_combined_mask(
    image_data,
    image_header,
    mask_file,
    coords_x,
    coords_y,
    RA,
    DEC,
    SMA,
    SMB,
    PAN
):
    
    #Create a combined mask consisting of original WISE/R mask and new 3x larger ellipse


    current_shape = image_data.shape

    original_wcs = WCS(image_header)
    
    #proj_plane_pixel_scales used to calculate spatial size of pixel along each axis of image
    pixscale = proj_plane_pixel_scales(original_wcs) * 3600

    print("Pixel scale (arcsec/pixel):", pixscale)

    # --------------------------------------------------
    # Recover coordinates if masked or invalid
    # --------------------------------------------------

    bad_coords = (
        np.ma.is_masked(coords_x)
        or np.ma.is_masked(coords_y)
        or not np.isfinite(coords_x)
        or not np.isfinite(coords_y)
    )

    
    #make sure code doesn't crash because of bad pixels
    if bad_coords:
        coords_x, coords_y = original_wcs.all_world2pix(
            RA,
            DEC,
            0
        )

        #print(
            #f"Recovered coordinates from WCS: "
            #f"{coords_x:.1f}, {coords_y:.1f}"
        #)

    if (
        not np.isfinite(coords_x)
        or not np.isfinite(coords_y)
    ):
        #print("Could not recover coordinates")
        return None
    
    #print("SMA =", SMA)
    #print("SMB =", SMB)
    #print("coords =", coords_x, coords_y)
    #print("image shape =", current_shape)



    #convert SMA / SMB to pixel units 
    SMA_pix = convert_to_pixels(SMA, image_header)
    SMB_pix = convert_to_pixels(SMB, image_header)

    #prevent crazy large apertures
    ny, nx = image_data.shape
    max_axis = 0.25 * min(nx, ny)

    scale = min(1.0, max_axis / max(SMA_pix, SMB_pix))

    SMA_pix *= scale
    SMB_pix *= scale

    

    #build large ellipse --> started as 3x but that was way too big on images so reduced to 2.5
    aperture = EllipticalAperture(
        (coords_x, coords_y),
        2.5 * SMA_pix,
        2.5 * SMB_pix,
        theta=PAN
        )

    #print(f"[{galaxy_name}] SMA arcsec={SMA:.2f} → pix={SMA_pix:.2f}")
    #print(f"[{galaxy_name}] SMB arcsec={SMB:.2f} → pix={SMB_pix:.2f}")
    #print(f"[{galaxy_name}] image shape={image_data.shape}")
    
    additional_mask = aperture.to_mask(method="exact").to_image(shape=current_shape)

    if additional_mask is None:
        return None

    #make sure no negatively masked pixels
    additional_mask = additional_mask > 0

    
    #reproject original wise/R mask
    mask_data, mask_header = fits.getdata(
        mask_file,
        header=True
    )

    mask_wcs = WCS(mask_header)

    reprojected_mask, footprint = reproject_interp(
        (mask_data, mask_wcs),
        original_wcs,
        shape_out=current_shape
    )

    #any nans in the original mask converted to 0 values to match the new mask
    reprojected_mask = np.nan_to_num(reprojected_mask, nan=0)
    reprojected_mask = reprojected_mask > 0
    
    #prints the percent of image covered by ellipse 
    print(
        "Ellipse fraction:",
        np.sum(additional_mask) / additional_mask.size
    )

    print(
        "Reprojected fraction:",
        np.sum(reprojected_mask) / reprojected_mask.size
    )



    #combine masks
    combined_mask = (additional_mask |reprojected_mask)
    
    plt.imshow(combined_mask, origin='lower')
    plt.colorbar()
    #plt.show()

    return combined_mask




# In[16]:


#imposing the masks onto the images by lining up the cental mask pixel with the central galaxy pixels calculated from the previous code block

#file paths
csv_file = os.path.join(datadir, 'new_pixel_coordinates.csv')
galaxy = Table.read(csv_file)

print(galaxy['GALAXY', 'SMA', 'SMB'][:5])

#print("Number of galaxies:", len(galaxy))

#print(type(galaxy['GALAXY'][0]))
#print(repr(galaxy['GALAXY'][0]))

#galaxy_name = str(galaxy['GALAXY'][0])

#print("galaxy_name =", galaxy_name)
#print("repr =", repr(galaxy_name))

#bad = galaxy['GALAXY'] == 'NGC4589'
#print(galaxy[bad])
#print(galaxy['HPPUNIMAPB_X'][bad])
#print(galaxy['HPPUNIMAPB_Y'][bad])

#for i in range(min(3, len(galaxy))):
#    galaxy_name = galaxy['GALAXY'][i]
#    path = os.path.join(datadir, 'pipeline', galaxy_name)

#    print(galaxy_name)
#    print("path exists:", os.path.exists(path))


# %%


#all of this code is edited from Tom's original Masking.py file but then edited to add additional elliptical mask

##### NOTES FROM CODING ######
    # overlay_mask_on_fits defined in masking_funct.py
    # HDUL = Header Data Unit List
    # lazy_load_hdus --> allows astropy to open FITS files by parsing the headers as needed instead of loading all HDUs into memory at once
            # saves processing time and memory for large datasets

for i in range(len(galaxy)):
    
    galaxy_name = galaxy['GALAXY'][i]
    VFID = str(galaxy['VFID'][i])
    
    path = os.path.join(datadir, 'pipeline', galaxy_name)
    
    
    SMA = galaxy['SMA'][i]
    SMB = galaxy['SMB'][i]
    RA = galaxy['RA'][i]
    DEC = galaxy['DEC'][i]
    EPLI = galaxy['EPLI'][i]
    PAN = galaxy['PAN'][i]



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

                #make output FITS file path for 'new' masks that have additional elliptical aperture mask
                new_output_fits = os.path.join(datadir, 'new_masked', f'new_{galaxy_name}_masked{color}.fits')                
                
                #read the CSV file for the central pixel coordinates
                coords_x = galaxy[f'HPPUNIMAP{color}_X'][i]
                coords_y = galaxy[f'HPPUNIMAP{color}_Y'][i]
                
                
                
                
                #####checking wise-mask#####
                
                #check if the wise-mask file exists, otherwise check for the r-mask
                #create two seperate fits files: one for the original masked image, one for the 'new' masked image
                
                if os.path.exists(wise_mask_file):
                    
                    #if wise-mask exists, overlay it
                    overlay_mask_on_fits(found_file, wise_mask_file, csv_file, output_fits, i, coords_x, coords_y, RA, DEC)  
                    
                    with fits.open(found_file) as hdul:
                        hdul.info()
                        image_hdu = None

                        with fits.open(found_file, lazy_load_hdus=False) as hdul:
                            print("HDU count:", len(hdul))
                            image_hdu = None

                            for k in range(len(hdul)):
                                if hdul[k].data is not None:
                                    image_hdu = hdul[k]
                                    break


                        image_data = image_hdu.data.copy()
                        image_header = image_hdu.header
                                
                    
                    #call the create_combined_mask function for the wise masks
                    combined_mask = create_combined_mask(
                        image_data=image_data,
                        image_header=image_header,
                        mask_file=wise_mask_file,
                        coords_x=coords_x,
                        coords_y=coords_y,
                        RA=RA,
                        DEC=DEC,
                        SMA=SMA,
                        SMB=SMB,
                        PAN=PAN
                        )
                    
                    print("Image shape:", image_data.shape)

                    print(
                        "Masked fraction:",
                        np.sum(combined_mask) / combined_mask.size
                        )
                    

                    if combined_mask is None:
                        continue
                    
                    # save image of the new combined mask
                    plt.figure(figsize=(6, 6))
                    plt.imshow(combined_mask, origin='lower', cmap='gray')
                    plt.title(f'{galaxy_name} {color} Combined Mask')
                    plt.colorbar(label='Masked')
                    plt.tight_layout()

                    mask_png = os.path.join(new_mask_plot_dir, f'{galaxy_name}_{color}_new_mask.png')

                    plt.savefig(mask_png, dpi=150)
                    plt.close()


                    #save image of new combined mask on herschel images
                    new_image = image_data.copy()
                    new_image[combined_mask] = np.nan
                    
                    plt.figure(figsize=(6, 6))
                    plt.imshow(new_image, origin='lower')
                    plt.title(f'{galaxy_name} {color} New Masked Image')
                    plt.colorbar()
                    plt.tight_layout()

                    masked_png = os.path.join(new_masked_herschel_dir, f'{galaxy_name}_{color}_new_masked.png')

                    plt.savefig(masked_png, dpi=150)
                    plt.close()


                    #write mask to FITS file
                    fits.writeto(
                        new_output_fits,
                        new_image,
                        image_header,
                        overwrite=True
                        )
                    
                    
                    
                    
                #####checking r_mask_file#####    
                #same as previous code block except for r_mask_file instead of wise_mask_file
                
                elif os.path.exists(r_mask_file):
                    #if r-mask exists, overlay it
                        
                    overlay_mask_on_fits(found_file, r_mask_file, csv_file, output_fits, i, coords_x, coords_y, RA, DEC)
                    
                        
                    with fits.open(found_file) as hdul:
                        image_hdu = None

                        for hdu in hdul:
                            if hdu.data is not None:
                                image_hdu = hdu
                                break

                        if image_hdu is None:
                            print(f"No image data found in {found_file}")
                            continue


                        image_data = image_hdu.data.copy()
                        image_header = image_hdu.header
                    
                    
                    combined_mask = create_combined_mask(
                        image_data=image_data,
                        image_header=image_header,
                        mask_file=r_mask_file,
                        coords_x=coords_x,
                        coords_y=coords_y,
                        RA=RA,
                        DEC=DEC,
                        SMA=SMA,
                        SMB=SMB,
                        PAN=PAN
                        )

                    if combined_mask is None:
                        continue
                    
                    # save image of the new combined mask
                    plt.figure(figsize=(6, 6))
                    plt.imshow(combined_mask, origin='lower', cmap='gray')
                    plt.title(f'{galaxy_name} {color} Combined Mask')
                    plt.colorbar(label='Masked')
                    plt.tight_layout()

                    mask_png = os.path.join(new_mask_plot_dir, f'{galaxy_name}_{color}_new_mask.png')

                    plt.savefig(mask_png, dpi=150)
                    plt.close()


                    #save image of combined mask on top of herschel image
                    new_image = image_data.copy()
                    new_image[combined_mask] = np.nan
                    
                    plt.figure(figsize=(6, 6))
                    plt.imshow(new_image, origin='lower')
                    plt.title(f'{galaxy_name} {color} New Masked Image')
                    plt.colorbar()
                    plt.tight_layout()

                    masked_png = os.path.join(new_masked_herschel_dir, f'{galaxy_name}_{color}_new_masked.png')

                    plt.savefig(masked_png, dpi=150)
                    plt.close()

                    fits.writeto(
                        new_output_fits,
                        new_image,
                        image_header,
                        overwrite=True
                        )
                    
                    
                else:
                    #if no mask files exist, copy the original FITS file to the output location
                    if not os.path.exists(os.path.dirname(output_fits)):
                        os.makedirs(os.path.dirname(output_fits))
                    with fits.open(found_file) as hdul:
                        hdul.writeto(output_fits, overwrite=True)
                        hdul.writeto(new_output_fits, overwrite=True)
                    print(f'Mask file not found for {galaxy_name} or {VFID}. Saved original FITS as {output_fits} and new FITS as {new_output_fits}.')
                    
                    
                    
                    


# In[17]:



# In[ ]:




 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    