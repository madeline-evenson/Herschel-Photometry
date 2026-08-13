#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 15:34:59 2026

@author: madeline.evenson
"""

#CIGALE outputs

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy.stats import linregress
from astropy.table import Table

# %%

#function to load FITS data into a pandas DataFrame

def fits_to_dataframe(fits_file, columns=None):
    
    with fits.open(fits_file) as hdul:
        
        data = hdul[1].data
        df = pd.DataFrame(data.tolist(), columns=data.names)
        
        if columns:
            df = df[columns]
            
        return df
    
    
#function to calculate MH2 and errors
def calculate_mh2(row):
    
    if pd.isna(row['MHI']) or pd.isna(row['MH2_to_MHI']):
        return np.nan, np.nan, np.nan
    
    else:
        MH2 = row['MHI'] * row['MH2_to_MHI']
        
        #calculate upper and lower errors in log space for MH2
        error_term_up = np.sqrt((row['MHI_err_up'] / row['MHI'])**2 + (row['MH2_to_MHI_err_up'] / row['MH2_to_MHI'])**2)
        error_term_down = np.sqrt((row['MHI_err_down'] / row['MHI'])**2 + (row['MH2_to_MHI_err_down'] / row['MH2_to_MHI'])**2)

        MH2_err_up = np.log10(MH2 * (1 + error_term_up))
        MH2_err_down = np.log10(MH2 * (1 + error_term_down)) #ensure positive term in log

        return MH2, MH2_err_up, MH2_err_down        