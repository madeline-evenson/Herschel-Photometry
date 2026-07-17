#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 18:43:24 2026

@author: madeline.evenson
"""

#CIGALEInputprep.py

import pandas as pd
import numpy as np

from astropy.io import fits


# %%

#turn inverse variance into uncertainty
#fix from Tom's code --> need to create all uncertainty columns first, THEN drop the IVAR columns 
def inverse_variance_into_uncertainty(input_file, output_file):
    
    #load the csv file
    df = pd.read_csv(input_file)
    
    #find all the inverse-variance columns
    ivar_cols = [col for col in df.columns if col.startswith('FLUX_IVAR_AP06')]
    
    #crate uncertainty columns
    for col in ivar_cols:
        
        new_col_name = col.replace('FLUX_IVAR_AP06', 'FLUX_UNCERT_AP06')
        
        #avoid divide-by-zero warnings
        df[new_col_name] = np.where(df[col] > 0, 1.0/np.sqrt(df[col]), np.nan)
        
    #drop IVAR columns after all conversions are complete
    df.drop(columns=ivar_cols, inplace=True)
    
    #save to output_file
    df.to_csv(output_file, index=False)
        
    print(f'File successfully written to {output_file}')
        
        
# %%

def multiply_columns_and_save(file1, file2, column_pairs, output_file):

    import pandas as pd

    flux_df = pd.read_csv(file1)
    ext_df = pd.read_csv(file2)

    # Remove duplicate column names if they exist
    flux_df = flux_df.loc[:, ~flux_df.columns.duplicated()]
    ext_df = ext_df.loc[:, ~ext_df.columns.duplicated()]


    # Rename VF_ID -> VFID if needed ---> this step used so this function can be used in both
        #CIGALE-inputs-text-file.py and vf-ephots-text-file.py
    if 'VF_ID' in flux_df.columns and 'VFID' not in flux_df.columns:
        flux_df = flux_df.rename(columns={'VF_ID':'VFID'})

    if 'VF_ID' in ext_df.columns and 'VFID' not in ext_df.columns:
        ext_df = ext_df.rename(columns={'VF_ID':'VFID'})


    #standardize VFID formatting --> again so this function can be used in both .py files
    def format_vfid(series):

        return (
            series
            .astype(str)
            .str.strip()
            .str.replace('VFID', '', regex=False)
            .str.replace('.0', '', regex=False)
            .astype(int)
            .apply(lambda x: f"VFID{x:04d}")
        )


    flux_df['VFID'] = format_vfid(flux_df['VFID'])
    ext_df['VFID'] = format_vfid(ext_df['VFID'])

    #print statements to double check nothing weird happened after merge
    print("===== MERGE CHECK =====")
    print("Photometry rows:", len(flux_df))
    print("Extinction rows:", len(ext_df))
    print("Photometry VFID duplicates:", flux_df['VFID'].duplicated().sum())
    print("Extinction VFID duplicates:", ext_df['VFID'].duplicated().sum())


    #keep only extinction columns needed for correction --> mostly used in vf_ephots-text-file.py
    ext_keep = ['VFID']

    for _, ext_col in column_pairs:
        if ext_col not in ext_keep:
            ext_keep.append(ext_col)

    ext_df = ext_df[ext_keep]


    #merge extinction information
    result_df = flux_df.merge(
        ext_df,
        on='VFID',
        how='left',
        suffixes=('', '_ext')
    )


    print("Matched extinction values:")
    print(result_df['VFID'].isin(ext_df['VFID']).sum())


    #apply extinction correction
    #F_corrected = F * 10^(0.4*A_lambda)
    for flux_col, ext_col in column_pairs:

        if flux_col in result_df.columns and ext_col in result_df.columns:

            result_df[flux_col] = (
                result_df[flux_col] *
                10**(0.4 * result_df[ext_col])
            )


    #remove extinction columns after correction --> cleans up files because we don't need this info anymore
    drop_cols = [
        col for col in result_df.columns
        if col.startswith('A(')
    ]

    result_df = result_df.drop(columns=drop_cols)


    #final duplicate cleanup
    result_df = result_df.loc[:, ~result_df.columns.duplicated()]


    result_df.to_csv(output_file, index=False)

    print(f"File successfully written to {output_file}")
# %%

#function to convert a .fits file to a .csv file
def fits_to_csv(fits_file, csv_file):
    
    #open the .fits file
    with fits.open(fits_file) as hdul:
        #assuming the data is in the first extension (index 1) or primary HDU (index 0)
        data = hdul[1].data if len(hdul) > 1 else hdul[0].data
        
        #convert the FITS data to a pandas dataframe
        df = pd.DataFrame(data)
        
        #safe the dataframe to a csv file
        df.to_csv(csv_file, index=False)
        print(f"Data successfully written to {csv_file}")
        
        
        
# %%


#getting the redshift out of vcosmic
def add_redshift_column(input_file, output_file):
    
    #speed of light in km/s
    speed_of_light = 299792.458
    
    #load the csv file
    df = pd.read_csv(input_file)
    
    #check if the "Vcosmic" column exists
    if 'Vcosmic' not in df.columns:
        raise ValueError("Column 'Vcosmic' not found in the input file")
        
    #calculate the redshift and add it as a new column
    df['redshift'] = df['Vcosmic'] / speed_of_light
    
    #replace negative redshift values with NaN
    df['redshift'] = df['redshift'].apply(lambda x: np.nan if x < 0 else x)
    
    #round the redshift values to 3 decimal places
    df['redshift'] = df['redshift'].round(4)
    
    #save the modified DataFrame to a new CSV file
    df.to_csv(output_file, index=False)
    
    print(f'File successfully written to {output_file}')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    