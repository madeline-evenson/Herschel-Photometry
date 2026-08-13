#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 13:40:24 2026

@author: madeline.evenson
"""

#this notebook has the goal of creating a proper .txt file for the vf ephot catalog that can be run through CIGALE
#mostly based on Tom's original CIGALE-inputs-test-file.py notebook, but had to be edited a bunch

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
cigaledir = homedir + '/cigale/'

from CIGALEInputprep import *

# %%

#convert legacy ephot catalog into csv file so can plug into CIGALE-inputs-text-file pipeline
ephot_input = tabledir + 'vf_v2_legacy_ephot.fits'
fits_table = Table.read(ephot_input, hdu=1)

#convert astropy table to pandas dataframe
df = fits_table.to_pandas()

#decode byte columns to strings
for col in df.columns:
    if df[col].dtype == object:
        try:
            df[col] = df[col].str.decode("utf-8")
        except AttributeError:
            pass
        



        
csv_path = tabledir + 'ephot_output_table.csv'
        
#save to csv
df.to_csv(csv_path, index = False)


# %%

#first we need to convert the table's inverse variance into uncertainty
#inverse_variance_into_uncertainty is defined in CIGALEInputprep.py

input_file = tabledir + 'ephot_output_table.csv'
output_file = tabledir + 'ephot_outputuncert.csv'

inverse_variance_into_uncertainty(input_file, output_file)

# %%


#Tom's notes: once we have the uncertainties of these measurements, we need to add extinction data from the 
    #VFS into it using the extinction table fits. we do this by finding the appropriate extinction columns and
    #multiplying it by the appropriate wavelengths inside the photometry tables
    
   
#now we add extinction into the table measurements --> the extinction data is from the extinction table fits,
    #trimmed to only have Herschel galaxies
    
   
#upload vf_v2_environment table --> originally in ascii format so need to convert into csv
#this table from https://ui.adsabs.harvard.edu/abs/2022ApJS..259...43C/abstract (link from vf v2 tables README doc)
table = Table.read(
    "/Users/madeline.evenson/Research/Virgo/tables/vf_v2_environment.txt",
    format="ascii.cds")

table.write(
    "/Users/madeline.evenson/Research/Virgo/tables/vf_v2_environment.csv",
    format="csv",
    overwrite=True)


file1 = tabledir + 'ephot_outputuncert.csv'
file2 = tabledir + 'trimmedextinction.csv'

#debug check for input files
df1_check = pd.read_csv(file1)
df2_check = pd.read_csv(file2)




#list of column pairs to multiply (from file1, file2)
column_pairs = [
    ('FLUX_AP06_G', 'A(G)_SFD'),
    ('FLUX_UNCERT_AP06_G', 'A(G)_SFD'),
    ('FLUX_AP06_R', 'A(R)_SFD'),
    ('FLUX_UNCERT_AP06_R', 'A(R)_SFD'), 
    ('FLUX_AP06_Z', 'A(Z)_SFD'),
    ('FLUX_UNCERT_AP06_Z', 'A(Z)_SFD'),
    ('FLUX_AP06_FUV', 'A(FUV)_SFD'),
    ('FLUX_UNCERT_AP06_FUV', 'A(FUV)_SFD'),
    ('FLUX_AP06_NUV', 'A(NUV)_SFD'),
    ('FLUX_UNCERT_AP06_NUV', 'A(NUV)_SFD'),
    ('FLUX_AP06_W1', 'A(W1)_SFD'),
    ('FLUX_UNCERT_AP06_W1', 'A(W1)_SFD'),
    ('FLUX_AP06_W2', 'A(W2)_SFD'),
    ('FLUX_UNCERT_AP06_W2', 'A(W2)_SFD'),
    ('FLUX_AP06_W3', 'A(W3)_SFD'),
    ('FLUX_UNCERT_AP06_W3', 'A(W3)_SFD'),
    ('FLUX_AP06_W4', 'A(W4)_SFD'),
    ('FLUX_UNCERT_AP06_W4', 'A(W4)_SFD')
    ]


output_file = tabledir + 'ephot_extinctionoutput.csv' #name of the output file

#print(df2_check.columns[df2_check.columns.duplicated()])

#print(pd.read_csv(file1).columns.tolist())
#print(pd.read_csv(file2).columns.tolist())

#multiply_columns_and_save is defined in CIGALEInputprep.py 
multiply_columns_and_save(file1, file2, column_pairs, output_file)


test = pd.read_csv(output_file)

#print(test.columns[test.columns.duplicated()])

# %%


#exporting the ephots table to only the columns we need
csv_input_file = tabledir + 'ephot_extinctionoutput.csv'
csv_output_file = tabledir + 'ephot_finalreal.csv'





#specify the column we want to output
columns_to_keep = [
    'VFID', 'GALAXY', 'RA_MOMENT', 'DEC_MOMENT',
    'FLUX_AP06_G',   'FLUX_UNCERT_AP06_G', 
    'FLUX_AP06_R',   'FLUX_UNCERT_AP06_R',
    'FLUX_AP06_Z',   'FLUX_UNCERT_AP06_Z',
    'FLUX_AP06_FUV', 'FLUX_UNCERT_AP06_FUV',
    'FLUX_AP06_NUV', 'FLUX_UNCERT_AP06_NUV',
    'FLUX_AP06_W1',  'FLUX_UNCERT_AP06_W1',
    'FLUX_AP06_W2',  'FLUX_UNCERT_AP06_W2',
    'FLUX_AP06_W3',  'FLUX_UNCERT_AP06_W3',
    'FLUX_AP06_W4',  'FLUX_UNCERT_AP06_W4'
    ]


#read the csv file
df = pd.read_csv(csv_input_file)

#another duplicate columns check --> there had been issues with the VFID column
print("Duplicate columns:")
print(df.columns[df.columns.duplicated()].tolist())

#another duplicate columns check --> make sure VFID exists only once
if 'VF_ID' in df.columns:
    df['VFID'] = df['VF_ID']

#remove duplicate columns (keeps first occurrence)
df = df.loc[:, ~df.columns.duplicated()]

#format VFIDs consistently --> this way we can use multiply_columns_and_save from CIGALEInputprep.py
df['VFID'] = (
    df['VFID']
    .astype(str)
    .str.strip()
    .str.replace('.0', '', regex=False)
    .str.replace('VFID', '', regex=False)
    .apply(lambda x: f"VFID{int(float(x)):04d}")
)

#convert nanomaggies to mJy:
#1 nanomaggy = 3.631e-3 mJy
factor = 3.631e-3
flux_columns = [col for col in df.columns if col.startswith('FLUX')]
df[flux_columns] = df[flux_columns] * factor



#iterate over each pair of FLUX and FLUX_UNCERT columns and apply the check
for flux_col in flux_columns:
    uncert_col = flux_col.replace('FLUX_', 'FLUX_UNCERT_')
    
    if uncert_col in df.columns:
        #replace the value in the FLUX_UNCERT column if it's smaller than 5% of the corresponding FLUX column
        #this is to keep the errors to a minimum of 5% to keep results somewhat accurate
        df[uncert_col] = df[[flux_col, uncert_col]].apply(lambda row: max(row[uncert_col], 0.05 * abs(row[flux_col])), axis=1)
        
        
 


#write the selected columns to a new CSV file
df.to_csv(csv_output_file, index=False)


# %%

#write out the flux data to text files for CIGALE

#define paths for north and south output files
north_path = cigaledir + '/vf_ephot_north/vf_ephot_north.txt'
south_path = cigaledir + '/vf_ephot_south/vf_ephot_south.txt'

#read the main CSV file
csv_input_file = tabledir + 'ephot_finalreal.csv'
flux_tab = pd.read_csv(csv_input_file)

#replace zeros with NaN
flux_tab.replace(0, np.nan, inplace=True)


south = flux_tab[flux_tab["DEC_MOMENT"] < 32]




good = flux_tab[flux_tab["FLUX_AP06_G"].notna()]


#same north/south flag as CIGALE-inputs-text-file.py
south = flux_tab[flux_tab["DEC_MOMENT"] < 32]
north = flux_tab[flux_tab["DEC_MOMENT"] > 32]



#create a new environment file with a redshift column
environment_input = tabledir + 'vf_v2_environment.csv'
environment_output = tabledir + 'vf_v2_environment_redshift.csv'

add_redshift_column(environment_input, environment_output)

#read the updated file
redshift_tab = pd.read_csv(environment_output)

#format VFIDs to match the photometry table
redshift_tab["VFID"] = redshift_tab["VFID"].apply(lambda x: f"VFID{int(x):04d}")

#merge flux table with redshift table
flux_tab = pd.merge(
    flux_tab,
    redshift_tab[['VFID', 'redshift']],
    on='VFID',
    how='left'
)



# photometry columns
phot_cols = [
    'FLUX_AP06_FUV',
    'FLUX_AP06_NUV',
    'FLUX_AP06_G',
    'FLUX_AP06_R',
    'FLUX_AP06_Z',
    'FLUX_AP06_W1',
    'FLUX_AP06_W2',
    'FLUX_AP06_W3',
    'FLUX_AP06_W4'
]


# remove galaxies with no photometry
flux_tab = flux_tab.dropna(subset=phot_cols, how='all')

# remove galaxies with missing redshift
flux_tab = flux_tab.dropna(subset=['redshift'])


#print("North =", (flux_tab["DEC_MOMENT"] > 32).sum())
#print("South =", (flux_tab["DEC_MOMENT"] < 32).sum())

#create north/south flags
north_flag = flux_tab['DEC_MOMENT'] > 32
south_flag = flux_tab['DEC_MOMENT'] < 32




#write north data
with open(north_path, 'w') as file:

    header = (
        "id redshift "
        "galex.FUV galex.FUV_err "
        "galex.NUV galex.NUV_err "
        "BASS-g BASS-g_err "
        "BASS-r BASS-r_err "
        "wise.W1 wise.W1_err "
        "wise.W2 wise.W2_err "
        "wise.W3 wise.W3_err "
        "wise.W4 wise.W4_err\n"
    )

    file.write(header)

    for _, n in flux_tab[north_flag].iterrows():

        s = (
            f"{n['VFID']} {n['redshift']} "
            f"{n['FLUX_AP06_FUV']} {n['FLUX_UNCERT_AP06_FUV']} "
            f"{n['FLUX_AP06_NUV']} {n['FLUX_UNCERT_AP06_NUV']} "
            f"{n['FLUX_AP06_G']} {n['FLUX_UNCERT_AP06_G']} "
            f"{n['FLUX_AP06_R']} {n['FLUX_UNCERT_AP06_R']} "
            f"{n['FLUX_AP06_W1']} {n['FLUX_UNCERT_AP06_W1']} "
            f"{n['FLUX_AP06_W2']} {n['FLUX_UNCERT_AP06_W2']} "
            f"{n['FLUX_AP06_W3']} {n['FLUX_UNCERT_AP06_W3']} "
            f"{n['FLUX_AP06_W4']} {n['FLUX_UNCERT_AP06_W4']}\n"
        )

        file.write(s)


#write south data
with open(south_path, 'w') as file:

    header = (
        "id redshift "
        "galex.FUV galex.FUV_err "
        "galex.NUV galex.NUV_err "
        "decamDR1-g decamDR1-g_err "
        "decamDR1-r decamDR1-r_err "
        "decamDR1-z decamDR1-z_err "
        "wise.W1 wise.W1_err "
        "wise.W2 wise.W2_err "
        "wise.W3 wise.W3_err "
        "wise.W4 wise.W4_err\n"
    )

    file.write(header)

    for _, n in flux_tab[south_flag].iterrows():

        s = (
            f"{n['VFID']} {n['redshift']} "
            f"{n['FLUX_AP06_FUV']} {n['FLUX_UNCERT_AP06_FUV']} "
            f"{n['FLUX_AP06_NUV']} {n['FLUX_UNCERT_AP06_NUV']} "
            f"{n['FLUX_AP06_G']} {n['FLUX_UNCERT_AP06_G']} "
            f"{n['FLUX_AP06_R']} {n['FLUX_UNCERT_AP06_R']} "
            f"{n['FLUX_AP06_Z']} {n['FLUX_UNCERT_AP06_Z']} "
            f"{n['FLUX_AP06_W1']} {n['FLUX_UNCERT_AP06_W1']} "
            f"{n['FLUX_AP06_W2']} {n['FLUX_UNCERT_AP06_W2']} "
            f"{n['FLUX_AP06_W3']} {n['FLUX_UNCERT_AP06_W3']} "
            f"{n['FLUX_AP06_W4']} {n['FLUX_UNCERT_AP06_W4']}\n"
        )

        file.write(s)




