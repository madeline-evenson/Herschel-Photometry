#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:51:37 2026

@author: madeline.evenson
"""

#the objective of this notebook is to gather avaliable data from all of the tables with photometry data from FUV to WISE,
    #combine them with the measured Herschel data, and input everything in a .txt file that's labelled according to CIGALE's specifications
    
#the first 2 code blocks are used to import all the necessary packages to efficiently run the notebook



# %%

#import the necessary packages

from astropy.table import Table
from matplotlib import pyplot as plt
import os
import numpy as np
from astropy.io.ascii import masked
from astropy.io import ascii
import glob
from astropy.io import fits
import wget
import matplotlib.image as mpimg
from astropy.wcs import WCS
from scipy.stats import scoreatpercentile
from astropy.visualization import simple_nome
from reproject import reproject_interp
import sys
from IPython.display import clear_output
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization import ImageNormalize
from astropy.visualization import LogStretch
import astropy.units as u
from matplotlib import colors
import warnings
import pandas as pd


os.environ['HOME'] = '/Users/madeline.evenson/Research' #general path for all the coding space
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are


from CIGALEInputprep import *




# %%

#first thing we'll do is convert the photometry table'sinverse varianceinto uncertainty
#the function works by looking for any columns that starts with "FLUX_IVAR_AP06" in it, and transforms it by doing an inverse square of that value
#in order tomake it work for all 6 apertures though, we will need to change that code so that it will iterate through all 6


#first we need to convert the table's inverse variance into uncertainty
input_file = tabledir+'Photometrytesting2.csv'
output_file = 'outputuncert.csv'

inverse_variance_into_uncertainty(input_file,output_file)

# %%

#once we have the uncertainties of these measurements, we need to add extinction data from the VFS into it using the extinction table fits
#we can do this by finding the appropriate extinction columns, and multiplying it by the appropriate wavelengths inside the photometry tables

#then we add extinction into the table measurements
#the extinction data is fromthe extinction table fits, trimming to only have Herschel galaxies

file1 = tabledir + 'outputuncert.csv'
file2 = tabledir + 'trimmedextinction.csv'

#listof column pairs to multiply (from file1, file2)
column_pairs = [
    ('FLUX_AP06_G', 'A(G)_SFD'),
    ('FLUX_UNCERT_AP06_G', 'A(G)_SFD'),
    ('FLUX_AP06_R','A(R)_SFD'),
    ('FLUX_UNCERT_AP06_R', 'A(Z)_SFD'),
    ('FLUX_AP06_Z', 'A(Z)_SFD'),
    ('FLUX_UNCERT_AP06_Z', 'A(Z)_SFD'),
    ('FLUX_AP06_FUV', 'A(FUV)_SFD'),
    ('FLUX_UNCERT_AP06_FUV','A(FUV)_SFD'),
    ('FLUX_AP06_NUV', 'A(NUV)_SFD'),
    ('FLUX_UNCERT_AP06_NUV','A(NUV)_SFD'),
    ('FLUX_AP06_W1', 'A(W1)_SFD'),
    ('FLUX_UNCERT_AP06_W1', 'A(W1)_SFD'),
    ('FLUX_AP06_W2', 'A(W2)_SFD'),
    ('FLUX_UNCERT_AP06_W2', 'A(W2)_SFD'),
    ('FLUX_AP06_W3','A(W3)_SFD'),
    ('FLUX_UNCERT_AP06_W3','A(W3)_SFD'),
    ('FLUX_AP06_W4','A(W4)_SFD'),
    ('FLUX_UNCERT_AP06_W4', 'A(W4)_SFD')
]

output_file = 'extinctoutput.csv'  #name of the output file

multiply_columns_and_save(file1, file2, column_pairs, output_file)


# %%



# %%
#because the photometry table still has a lot of information that's not necessary for CIGALE inputs,Ithen remove any columns that I won't be using in order to make everything easier to keep track of
#in addition to this, I also convert all fo the Fluxes (that are not Herschel's) from nanomaggies to micro-Janskies, as that is what CIGALE takes in
#the factor used is after consulting with Kim about it
#as for the unit converstions on Herschel galaxies, I instead multiply it by 1000 to go from Janskies to micro-Janskies
#then, I create error columns for the Herschel fluxes by creating a column and multtiplying every value by 5% --> THIS NUMBER CAN BE FURTHER FINETUNED
#then,we do a few more checks such as an error limit cut, as well as converting negative fluxes to blanks, as well as rounding the inputs so that they can be compatible with CIGALE
#then we output only what we need into a csv file



#exporting the ephots table to only the columns I need
csv_input_file = tabledir + 'extinctoutput.csv'
csv_output_file = 'finalreal.csv'

#specify the columns you want to output
columns_to_keep = ['VF_ID', 'GALAXY','RA_MOMENT','DEC_MOMENT', 'FLUX_AP06_G', 'FLUX_UNCERT_AP06_G', 'FLUX_AP06_R', 'FLUX_UNCERT_AP06_R','FLUX_AP06_Z', 'FLUX_UNCERT_AP06_Z','FLUX_AP06_FUV', 'FLUX_UNCERT_AP06_FUV','FLUX_AP06_NUV', 'FLUX_UNCERT_AP06_NUV',
                   'FLUX_AP06_W1', 'FLUX_UNCERT_AP06_W1','FLUX_AP06_W2', 'FLUX_UNCERT_AP06_W2','FLUX_AP06_W3', 'FLUX_UNCERT_AP06_W3','FLUX_AP06_W4', 'FLUX_UNCERT_AP06_W4','70Flux_AP06','100Flux_AP06','160Flux_AP06']

#read the csv file
df = pd.read_csv(csv_input_file)


#format the VF_ID column to 'VFIDxxxx' with leading zeros
df['VF_ID'] = df['VF_ID'].apply(lambda x: f'VFID{x:04d}')


#multiply all columns that start with 'FLUX' by 3.631e-3 as the conversion for nanomaggies to Jansky
factor = 3.631e-3 
flux_columns = [col for col in df.columns if col.startswith('FLUX')]
df[flux_columns] = df[flux_columns] * factor


#multiply specific columns by 1000 before any further calculations to convert to Jansky for Herschel
df[['70Flux_AP06', '100Flux_AP06', '160Flux_AP06']] = df[['70Flux_AP06', '100Flux_AP06', '160Flux_AP06']] * 1000


#create error columns by multiplying the original flux columns by 5%
error_cols = {
    '70Flux_AP06_err': df['70Flux_AP06'] * 0.05,
    '100Flux_AP06_err': df['100Flux_AP06'] * 0.05,
    '160Flux_AP06_err': df['160Flux_AP06'] * 0.05
}


#add the new error columns using pd.concat to avoid fragmentation
df = pd.concat([df, pd.DataFrame(error_cols)], axis=1)


#add the new error columns to the list of columns
columns_to_keep.extend(['70Flux_AP06_err', '100Flux_AP06_err', '160Flux_AP06_err'])


#iterate over each pair of FLUX and FLUX_UNCERT columns and apply the check
for flux_col in flux_columns:
    uncert_col = flux_col.replace('FLUX_', 'FLUX_UNCERT_')
    if uncert_col in df.columns:
        #replace the value in the FLUX_UNCERT column if it's smaller than 5% of the corresponding FLUX column, this is to keep the errors to a minimum of 5% to keep results somewhat accurate
        df[uncert_col] = df[[flux_col, uncert_col]].apply(lambda row: max(row[uncert_col], 0.05 * row[flux_col]), axis=1)
  
        
#check for negative values and replace with NaN
# ******** this can be commented out if we want potentially negative fluxes in the CIGALE inputs **********
columns_to_check = flux_columns + ['70Flux_AP06', '100Flux_AP06', '160Flux_AP06'] + [col for col in df.columns if col.endswith('_err')]
#df[columns_to_check] = df[columns_to_check].map(lambda x: np.nan if x < 0 else x)

#round all numerical values to 3 decimal places
df[columns_to_check] = df[columns_to_check].round(3)

#select only the specified columns
df_selected = df[columns_to_keep].copy()

#rename the VF_ID column to VFID
df_selected.rename(columns={'VF_ID': 'VFID'}, inplace=True)

#write the selected columns to a new CSV file
df_selected.to_csv(csv_output_file, index=False)


print(f'Selected columns successfully written to {csv_output_file}')




# %%

#the final block is used to generate the two .txt files directly from the table generated above
#the two files are separated using the north and south designation using the galaxies' DEC
#we also need a redshift measurement of these galaxies, which is taken from the environments .fits table from the VFS
#once everything is read in, I then follow the template Kim gave to write my two input files

# %%

#write out the flux data to a text file for CIGALE
#define paths for northa nd south output files
north_path = homedir + '/Table-generating/vf_data_north.txt'
south_path = homedir + '/Table-generating/vf_data_south.txt'


#read the main CSV file
csv_input_file = homedir + '/Table-generating/finalreal.csv'
flux_tab = pd.read_csv(csv_input_file)


#replace zeros with NaN
flux_tab.replace(0, np.nan, inplace=True)


#read the second CSV file containing the redshift information
redshift_file = homedir + '/Table-generating/trimmedenvironment2.csv'
redshift_tab = pd.read_csv(redshift_file)


#merge flux_tab with redshift_tab based on 'VFID'
flux_tab = pd.merge(flux_tab, redshift_tab[['VFID', 'redshift']], on='VFID', how='left')


#create flags for north and south based on DEC_MOMENT
north_flag = flux_tab['DEC_MOMENT'] > 32
south_flag = flux_tab['DEC_MOMENT'] < 32


#function to check and create directories if they don't exist
def check_dir(*paths):
    for path in paths:
        directory = os.path.dirname(path)
        if not os.path.exists(directory):
            os.makedirs(directory)


#check and create directories for the output files
check_dir(north_path, south_path)


#write north data
with open(north_path, 'w') as file:
    #create file header
    s = '# id redshift FUV FUV_err NUV NUV_err BASS-g BASS-g_err BASS-r BASS-r_err WISE1 WISE1_err WISE2 WISE2_err WISE3 WISE3_err WISE4 WISE4_err PACS_blue PACS_blue_err PACS_green PACS_green_err PACS_red PACS_red_err\n'
    file.write(s)
    
    #write data rows
    for _, n in flux_tab[north_flag].iterrows():
        s_gal = f"{n['VFID']} {n['redshift']} {n['FLUX_AP06_FUV']} {n['FLUX_UNCERT_AP06_FUV']} {n['FLUX_AP06_NUV']} {n['FLUX_UNCERT_AP06_NUV']} " \
                f"{n['FLUX_AP06_G']} {n['FLUX_UNCERT_AP06_G']} {n['FLUX_AP06_R']} {n['FLUX_UNCERT_AP06_R']} " \
                f"{n['FLUX_AP06_W1']} {n['FLUX_UNCERT_AP06_W1']} {n['FLUX_AP06_W2']} {n['FLUX_UNCERT_AP06_W2']} " \
                f"{n['FLUX_AO06_W3']} {n['FLUX_UNCERT_AP06_W3']} {n['FLUX_AP06_W4']} {n['FLUX_UNCERT_AP06_W4']} " \
                f"{n['70Flux_AP06']} {n['70Flux_AP06_err']} {n['100Flux_AP06']} {n['100Flux_AP06_err']} {n['160Flux_AP06']} {n['160Flux_AP06_err']}\n"
        file_.write(s_gal)
        
        
#write south data
with open(south_path, 'w') as file:
    #create file header
    s = '# id redshift FUV FUV_err NUV NUV_err decamDR1-g decamDR1-g_err decamDR1-r decamDR1-r_err decamDR1-z decamDR1-z_err WISE1 WISE1_err WISE2 WISE2_err WISE3 WISE3_err WISE4 WISE4_err PACS_blue PACS_blue_err PACS_green PACS_green_err PACS_red PACS_red_err\n'
    file.write(s)
    
    #write data rows
    for _, n in flux_tab[south_flag].iterrows():
        s_gal = f"{n['VFID']} {n['redshift']} {n['FLUX_AP06_FUV']} {n['FLUX_UNCERT_AP06_FUV']} {n['FLUX_AP06_NUV']} {n['FLUX_UNCERT_AP06_NUV']} " \
                f"{n['FLUX_AP06_G']} {n['FLUX_UNCERT_AP06_G']} {n['FLUX_AP06_R']} {n['FLUX_UNCERT_AP06_R']} {n['FLUX_UNCERT_AP06_Z']} {n['FLUX_UNCERT_AP06_Z']} " \
                f"{n['FLUX_AP06_W1']} {n['FLUX_UNCERT_AP06_W1']} {n['FLUX_AP06_W2']} {n['FLUX_UNCERT_AP06_W2']} " \
                f"{n['FLUX_AP06_W3']} {n['FLUX_UNCERT_AP06_W3']} {n['FLUX_AP06_W4']} {n['FLUX_UNCERT_AP06_W4']} " \
                f"{n['70Flux_AP06']} {n['70Flux_AP06_err']} {n['100Flux_AP06']} {n['100Flux_AP06_err']} {n['160Flux_AP06']} {n['160Flux_AP06_err']}\n"
            file.write(s_gal)








