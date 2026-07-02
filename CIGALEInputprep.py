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



#adding extinctino into the template table
def multiply_columns_and_save(file1, file2, column_pairs, output_file):
    
    #load the csv files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    
    #initialize a dataframe to store the results
    result_df = df1.copy()
    
    #multiply each pair of columns and store the result in the result_df 
    for col1, col2 in column_pairs:
        if col1 not in df1.columns:
            raise ValueError(f"Column '{col1}' not found in {file1}")
        if col2 not in df2.columns:
            raise ValueError(f"Column '{col2}' not found in {file2}")
            
            
        #perform the multiplication with extinction correction
        result_df[col1] = df1[col1] * (10 ** (df2[col2] / 2.5))
        
        #replace 0, infinity, or very large/small values with NaN
        result_df[col1] = result_df[col1].replace([0, np.inf, -np.inf], np.nan)
        
    #save the resulting dataframe to a new csv file
    result_df.to_csv(output_file, index=False)
    
    
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    