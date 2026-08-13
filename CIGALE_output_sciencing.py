#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 15:29:54 2026

@author: madeline.evenson
"""


#Tom's notes: 'after getting the CIGALE outputs from Kim, I then do all of the analysis using this
#notebook. so this notebook wll mainly revolve around converting the dust mass from CIGALE to a total gas
#mass, as well as using that alongside other CIGALE outputs and data within the VFS to conduct 
#analysis and produce prots for the writeup

import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy.stats import linregress
from astropy.table import Table


#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are
cigaledir = homedir + '/cigale/'

n_ephot_dir = cigaledir + '/vf_ephot_north/out/'
s_ephot_dir = cigaledir + 'vf_ephot_south/out/'

n_herschel_dir = cigaledir + '/vf_north/out/'
s_herschel_dir = cigaledir + '/vf_south/out/'

from CIGALE_outputs import *

# %%

#opening each results.fits file and concatinating north/south results


n_ephot_file = n_ephot_dir + 'results.fits'
n_ephot_table = Table.read(n_ephot_file, format='fits')
n_ephot_df = n_ephot_table.to_pandas()

s_ephot_file = s_ephot_dir + 'results.fits'
s_ephot_table = Table.read(s_ephot_file, format='fits')
s_ephot_df = s_ephot_table.to_pandas()

n_herschel_file = n_herschel_dir + 'results.fits'
n_herschel_table = Table.read(n_herschel_file, format='fits')
n_herschel_df = n_herschel_table.to_pandas()

s_herschel_file = s_herschel_dir + 'results.fits'
s_herschel_table = Table.read(s_herschel_file, format='fits')
s_herschel_df = s_herschel_table.to_pandas()

kim_file = cigaledir + 'kim_results.txt'
kim_df = pd.read_csv(kim_file, sep=r'\s+')

ephot_df = pd.concat([n_ephot_df, s_ephot_df], ignore_index=True)
herschel_df = pd.concat([n_herschel_df, s_herschel_df], ignore_index=True)

# %%

#Tom's notes: 'the purpose of this next code is to convert the CIGALE dust mass into total gas mass
#using the ratio from Sandstrom. it first converst the dust masses from kg to solar masses, along with
#their uncertainties. then it converst all of these valeus into log space. then, to get the gas mass, 
#the program then calculates the gas mass in log space using the dust-to-gas ratio, alongside its
#uncertainties, and outputs it as a csv file

#link to sandstrom et al. (2013): https://ui.adsabs.harvard.edu/abs/2013ApJ...777....5S/abstract

# %%


#convert dust masses outputed by CIGALE into stellar mass, then calculate gas mass from dust to gas ratio

#solar mass in kg
solar_mass_kg = 1.98847e30

#dust to gas ratio parameters in log space
log_DGR = -1.86
DGR_err = 0.22 #DGR uncertainty in dex

# %%

##### ephot data #####

#converting bayes and best dust mass to msol
ephot_df['bayes.dust.mass.msol'] = ephot_df['bayes.dust.mass'] / solar_mass_kg
ephot_df['best.dust.mass.msol'] = ephot_df['best.dust.mass'] / solar_mass_kg

#convert uncertainties to solar masses  in log space
ephot_df['bayes.dust.mass.msol_err'] = ephot_df['bayes.dust.mass_err'] / (ephot_df['bayes.dust.mass'] * np.log(10)) #convert to log space

#convert bayes and best dust mass (msol) to log space
ephot_df['log.bayes.dust.mass'] = np.log10(ephot_df['bayes.dust.mass.msol'])
ephot_df['log.best.dust.mass'] = np.log10(ephot_df['best.dust.mass.msol'])

#calculate log gas masses and their uncertainties
ephot_df['log.bayes.gas.mass'] = ephot_df['log.bayes.dust.mass'] - log_DGR
ephot_df['log.bayes.gas.mass_err'] = np.sqrt(ephot_df['bayes.dust.mass.msol_err']**2 + DGR_err**2)

ephot_df['log.best.gas.mass'] = ephot_df['log.best.dust.mass'] - log_DGR
ephot_df['log.best.gas.mass_err'] = DGR_err #best dust mass assumed to have no additional error


#save to csv 
csv_file = tabledir + 'ephot-testoutput.csv' #output csv file

#create dataframe to store results
df = pd.DataFrame({
    'id': ephot_df['id'],
    'log(bayes.dust.mass)' : ephot_df['log.bayes.dust.mass'],
    'log(bayes.dust.mass_err)' : ephot_df['bayes.dust.mass.msol_err'],
    'log(best.dust.mass)' : ephot_df['log.best.dust.mass'],
    'log(bayes.gas.mass)' : ephot_df['log.bayes.gas.mass'],
    'log(bayes.gas.mass_err)' : ephot_df['log.bayes.gas.mass_err'],
    'log(best.gas.mass)' : ephot_df['log.best.gas.mass'],
    'log(best.gas.mass_err)' : ephot_df['log.best.gas.mass_err'],
    'bayes.sfh.sfr' : ephot_df['bayes.sfh.sfr'],
    'best.reduced_chi_square' : ephot_df['best.reduced_chi_square']})

#export to csv
df.to_csv(csv_file, index=False)
print(f'Ephot data successfully written to {csv_file}')

# %%

##### herschel data #####

#converting bayes and best dust mass to msol
herschel_df['bayes.dust.mass.msol'] = herschel_df['bayes.dust.mass'] / solar_mass_kg
herschel_df['best.dust.mass.msol'] = herschel_df['best.dust.mass'] / solar_mass_kg

#convert uncertainties to solar masses  in log space
herschel_df['bayes.dust.mass.msol_err'] = herschel_df['bayes.dust.mass_err'] / (herschel_df['bayes.dust.mass'] * np.log(10)) #convert to log space

#convert bayes and best dust mass (msol) to log space
herschel_df['log.bayes.dust.mass'] = np.log10(herschel_df['bayes.dust.mass.msol'])
herschel_df['log.best.dust.mass'] = np.log10(herschel_df['best.dust.mass.msol'])

#calculate log gas masses and their uncertainties
herschel_df['log.bayes.gas.mass'] = herschel_df['log.bayes.dust.mass'] - log_DGR
herschel_df['log.bayes.gas.mass_err'] = np.sqrt(herschel_df['bayes.dust.mass.msol_err']**2 + DGR_err**2)

herschel_df['log.best.gas.mass'] = herschel_df['log.best.dust.mass'] - log_DGR
herschel_df['log.best.gas.mass_err'] = DGR_err #best dust mass assumed to have no additional error


#save to csv 
csv_file = tabledir + 'herschel-testoutput.csv' #output csv file

#create dataframe to store results
df = pd.DataFrame({
    'id': herschel_df['id'],
    'log(bayes.dust.mass)' : herschel_df['log.bayes.dust.mass'],
    'log(bayes.dust.mass_err)' : herschel_df['bayes.dust.mass.msol_err'],
    'log(best.dust.mass)' : herschel_df['log.best.dust.mass'],
    'log(bayes.gas.mass)' : herschel_df['log.bayes.gas.mass'],
    'log(bayes.gas.mass_err)' : herschel_df['log.bayes.gas.mass_err'],
    'log(best.gas.mass)' : herschel_df['log.best.gas.mass'],
    'log(best.gas.mass_err)' : herschel_df['log.best.gas.mass_err'],
    'bayes.sfh.sfr' : herschel_df['bayes.sfh.sfr'],
    'best.reduced_chi_square' : herschel_df['best.reduced_chi_square']})

#export to csv
df.to_csv(csv_file, index=False)
print(f'Herschel data successfully written to {csv_file}')


# %%

##### kim data #####

#converting bayes and best dust mass to msol
kim_df['bayes.dust.mass.msol'] = kim_df['bayes.dust.mass'] / solar_mass_kg
kim_df['best.dust.mass.msol'] = kim_df['best.dust.mass'] / solar_mass_kg

#convert uncertainties to solar masses  in log space
kim_df['bayes.dust.mass.msol_err'] = kim_df['bayes.dust.mass_err'] / (kim_df['bayes.dust.mass'] * np.log(10)) #convert to log space

#convert bayes and best dust mass (msol) to log space
kim_df['log.bayes.dust.mass'] = np.log10(kim_df['bayes.dust.mass.msol'])
kim_df['log.best.dust.mass'] = np.log10(kim_df['best.dust.mass.msol'])

#calculate log gas masses and their uncertainties
kim_df['log.bayes.gas.mass'] = kim_df['log.bayes.dust.mass'] - log_DGR
kim_df['log.bayes.gas.mass_err'] = np.sqrt(kim_df['bayes.dust.mass.msol_err']**2 + DGR_err**2)

kim_df['log.best.gas.mass'] = kim_df['log.best.dust.mass'] - log_DGR
kim_df['log.best.gas.mass_err'] = DGR_err #best dust mass assumed to have no additional error


#save to csv 
csv_file = tabledir + 'kim-testoutput.csv' #output csv file

#create dataframe to store results
df = pd.DataFrame({
    'id': kim_df['id'],
    'log(bayes.dust.mass)' : kim_df['log.bayes.dust.mass'],
    'log(bayes.dust.mass_err)' : kim_df['bayes.dust.mass.msol_err'],
    'log(best.dust.mass)' : kim_df['log.best.dust.mass'],
    'log(bayes.gas.mass)' : kim_df['log.bayes.gas.mass'],
    'log(bayes.gas.mass_err)' : kim_df['log.bayes.gas.mass_err'],
    'log(best.gas.mass)' : kim_df['log.best.gas.mass'],
    'log(best.gas.mass_err)' : kim_df['log.best.gas.mass_err'],
    'bayes.sfh.sfr' : kim_df['bayes.sfh.sfr'],
    'best.reduced_chi_square' : kim_df['best.reduced_chi_square']})

#export to csv
df.to_csv(csv_file, index=False)
print(f'Kim data successfully written to {csv_file}')

# %%

#Tom's notes: 'once we have the CIGALE gas masses, we will also need something to compare it against
#this is where the VFS' CO measurements come in handy. this code takes in the HI measurements from 
#the HI reference table, along with its given HI to H2 ratio, converts it, and combines it to get
#the total gas mass

# %%

#calculate the H2 data from the HI reference table

#load the SV file
csv_file = tabledir + 'test.csv' #the table with HI measurements data
df = pd.read_csv(csv_file)

#apply the function to calculate the new MH2 columns
df[['MH2', 'MH2_err_up', 'MH2_err_down']] = df.apply(lambda row: pd.Series(calculate_mh2(row)), axis=1)

#convert all mass columns and errors to log10
mass_columns = ['MHI', 'MHI_err_up', 'MHI_err_down', 'MH2', 'MH2_err_up', 'MH2_err_down']

for col in mass_columns:
    df[col] = df[col].apply(lambda x: np.log10(x) if pd.notna(x) and x > 0 else np.nan)


df['Mgas'] = df.apply(
    lambda row: np.log10(10**row['MHI'] + 10**row['MH2'])
    if pd.notna(row['MHI']) and pd.notna(row['MH2']) else np.nan,
    axis=1)

#save the modified dataframe to a new CSV
output_csv = 'reference_gas_mass.csv'
df.to_csv(output_csv, index=False)

print(f"CSV file with MH2 and Mgas columns in log10 successfully written to {output_csv}")

# %%

#Tom's notes: once we have both of these tables, we can then perform analysis with it
#the first code, however, is a histogram of the dsut luminosity from CIGALE
#we did not convert the luminosity from the CIGALE's standard of Watts


def plot_dust_lum(df, title, output_file):
    dust_luminosity = df['best.dust.luminosity']
    
    #keep only positive, finite luminosities
    dust_luminosity = dust_luminosity[np.isfinite(dust_luminosity) & (dust_luminosity > 0)]
    
    #convert dust luminosity to log space
    log_dust_luminosity = np.log10(dust_luminosity)
    
    #plot histogram of log dust luminosity
    plt.figure(figsize=(8, 6))
    plt.hist(log_dust_luminosity, bins=30, color='violet', edgecolor='m')
    
    plt.xlabel('Log Dust Luminosity (Watts)')
    plt.ylabel('N')
    plt.title(title)
    
    output_path = plotdir + output_file
    plt.savefig(output_path, dpi=150)
    
    plt.show()
    
# %%

plot_dust_lum(ephot_df, 'Ephot Histogram of Log Dust Luminosity', output_file='Ephot-hist-log-dustlum.png')
plot_dust_lum(herschel_df, 'Herschel Histogram of Log Dust Luminosity', output_file='Herschel-hist-log-dustlum.png')
plot_dust_lum(kim_df, 'Kim Histogram of Log Dust Luminosity', output_file='Kim-hist-log-dustlum.png')

# %%


#function to plot Ha SFR vs the SFR from CIGALE to see whenever the two methods of calculating the SFR agree with each other

def Ha_vs_CIGALE_SFR(df, title, output_file):
    
    
    fits = tabledir + 'vf_v2_halpha.fits'
    
    #load both FITS tables as DataFrames

    df1 = df[['id', 'bayes.sfh.sfr', 'bayes.sfh.sfr_err']].copy()
    
    df2 = fits_to_dataframe(fits, columns=['VFID', 'GAL_LOG_SFR_HA', 'GAL_LOG_SFR_HA_ERR', 'GAL_LOG_SFR_HA_FLAG'])
    
    #ensure data types match for the merge (e.g. converting IDs to strings if needed)
    df1['id'] = df1['id'].astype(str)
    df2['VFID'] = df2['VFID'].astype(str)
    
    #cross-match the tables using the first table's 'id' and the second table's 'VFID'
    merged_df = pd.merge(df1, df2, left_on='id', right_on='VFID')
    
    #filter out rows where GAL_LOG_SFR_HA_FLAG is false
    filtered_df = merged_df[merged_df['GAL_LOG_SFR_HA_FLAG'] == True]
    
    #convert bayes.sfh.sfr to log space(handling non-positive values)
    filtered_df = filtered_df[filtered_df['bayes.sfh.sfr'] > 0] #remove non-positive values
    filtered_df['log_bayes_sfh_sfr'] = np.log10(filtered_df['bayes.sfh.sfr'])
    
    #check if there is valid data after filtering
    if not filtered_df.empty:
        #calculate the best-fit line using linear regression
        slope, intercept, _, _, _ = linregress(
            filtered_df['log_bayes_sfh_sfr'], filtered_df['GAL_LOG_SFR_HA'])
        
        #scatter plot for log_bayes_sfh_sfr vs GAL_LOG_SFR_HA
        plt.figure(figsize=(10, 6))
        plt.scatter(
            filtered_df['log_bayes_sfh_sfr'],
            filtered_df['GAL_LOG_SFR_HA'],
            color='blue', 
            label='Galaxy'
            )
        
        #plot the 1-to-1 line
        min_val = min(filtered_df['log_bayes_sfh_sfr'].min(), filtered_df['GAL_LOG_SFR_HA'].min())
        max_val = max(filtered_df['log_bayes_sfh_sfr'].max(), filtered_df['GAL_LOG_SFR_HA'].max())
        plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle=':', linewidth=1.5, label='1-to-1 Line')
        
        #set title and labels
        plt.title(title)
        plt.xlabel(r'Log Bayesian SFR (CIGALE) $(M_{\odot}yr^{-1})$')
        plt.ylabel(r'SFR Hα (VFS) $(M_{\odot}yr^{-1})$')
        plt.grid(True)
        plt.legend()
        
        #save and show the plot
        plt.tight_layout()
        
        output_path = plotdir + output_file
        plt.savefig(output_path, dpi=150)
        
        plt.show()

# %%

Ha_vs_CIGALE_SFR(ephot_df, 'Ephot SFR Hα vs SFR Cigale', output_file='Ephot-SFR-Ha-CIGALE.png')
Ha_vs_CIGALE_SFR(herschel_df, 'Herschel SFR Hα vs SFR Cigale', output_file='Herschel-SFR-Ha-CIGALE.png')
Ha_vs_CIGALE_SFR(kim_df, 'Kim SFR Hα vs SFR Cigale', output_file='Kim-SFR-Ha-CIGALE.png')


# %%

#this code takes the difference between CIGALE's outputted SFRs to see the affect Herschel data has on accurate CIGALE measurements

#load both dataframes
df1 = herschel_df[['id', 'bayes.sfh.sfr']].copy()
df2 = ephot_df[['id', 'bayes.sfh.sfr']].copy()


#ensure data types match for the merge (e.g. converting IDs to strings if needed)
df1['id'] = df1['id'].astype(str)
df2['id'] = df2['id'].astype(str)

#merge the tables by the 'id' column
merged_df = pd.merge(df1, df2, on='id', suffixes=('_herschel', '_ephot'))

#filter out rows where SFR values are non-positive to avoid log errors
merged_df = merged_df[(merged_df['bayes.sfh.sfr_herschel'] > 0) & (merged_df['bayes.sfh.sfr_ephot'] > 0)]

#convert SFR values to log10 space
merged_df['log_bayes_sfh_sfr_herschel'] = np.log10(merged_df['bayes.sfh.sfr_herschel'])
merged_df['log_bayes_sfh_sfr_diff'] = np.log10(merged_df['bayes.sfh.sfr_herschel']) - np.log10(merged_df['bayes.sfh.sfr_ephot'])

#calculate the median of the plotted SFR difference values
median_diff = merged_df['log_bayes_sfh_sfr_diff'].median()

#define plot limits for consistent scaling and plotting
x_min = merged_df['log_bayes_sfh_sfr_herschel'].min()
x_max = merged_df['log_bayes_sfh_sfr_herschel'].max()

#plot SFR with Herschel vs the difference in SFRs
plt.figure(figsize=(10, 6))
plt.scatter(
    merged_df['log_bayes_sfh_sfr_herschel'],
    merged_df['log_bayes_sfh_sfr_diff'],
    color='tab:blue',
    label='Galaxy'
    )

plt.axhline(0, color='black', linestyle='--', label='y = 0')
plt.axhline(median_diff, color='red', linestyle=':', label=f'Median = {median_diff:.2f}')

#set title and labels
plt.title('SFR Difference (With - Without Herschel)')
plt.xlabel(r'log SFR (with Herschel) $(M_{\odot}yr^{-1})$')
plt.ylabel(r'SFR Difference (log(with Herschel) - log(without Herschel))')
plt.xlim(x_min, x_max)
plt.grid(True)
plt.legend()

#save and show the plot
plt.tight_layout()

output_file = 'SFR-difference-with-and-without-Herschel.png'
output_path = plotdir + output_file

plt.savefig(output_path, dpi=150)
plt.show()



# %%

# this code plots the VFS gas mass calculated earlier
# against the CIGALE gas mass converted using the Sandstrom relation


# load CSV files
first_csv_file = tabledir + 'reference_gas_mass.csv'
df1 = pd.read_csv(first_csv_file)

second_csv_file = tabledir + 'herschel-testoutput.csv'
df2 = pd.read_csv(second_csv_file)


#clean galaxy ids before mergind 
#reference gas mass IDs are already normal strings
df1['VFID'] = df1['VFID'].astype(str).str.strip()

# CIGALE IDs are stored as byte strings, e.g.: b'VFID0018'
#convert them to normal strings, e.g.:  VFID0018

def clean_id(x):
    if isinstance(x, bytes):
        return x.decode('utf-8').strip()
    return str(x).strip().replace("b'", "").replace("'", "")

df2['id'] = df2['id'].apply(clean_id)


#check that the IDs now have the same format
print("Reference VFID examples:")
print(df1['VFID'].head(10).tolist())

print("\nCIGALE ID examples:")
print(df2['id'].head(10).tolist())


#check how many ids match
reference_ids = set(df1['VFID'])
cigale_ids = set(df2['id'])

matching_ids = reference_ids & cigale_ids

print("\nNumber of matching IDs:", len(matching_ids))
print("Example matching IDs:", list(matching_ids)[:10])



#merge the two datasets
merged_df = pd.merge(
    df1,
    df2,
    left_on='VFID',
    right_on='id'
)

print("\nReference gas mass rows:", len(df1))
print("Herschel rows:", len(df2))
print("Merged rows:", len(merged_df))

# Find the most extreme gas-mass values
print("\nLowest reference gas masses:")
print(
    merged_df[
        ['VFID', 'Mgas', 'log(bayes.gas.mass)']
    ].sort_values('Mgas').head(10)
)

print("\nLowest CIGALE gas masses:")
print(
    merged_df[
        ['VFID', 'Mgas', 'log(bayes.gas.mass)']
    ].sort_values('log(bayes.gas.mass)').head(10)
)


#stop if there are no matching galaxies
if len(merged_df) == 0:
    raise ValueError(
        "No galaxies matched between the reference gas mass "
        "table and the CIGALE Herschel table.")




#reference gas mass
Mgas = merged_df['Mgas'] #already in log space

#CIGALE Bayesian gas mass
bayes_gas_mass = merged_df['log(bayes.gas.mass)'] #already in log space


#check data
print("\nMgas:")
print(Mgas.describe())

print("\nBayesian gas mass:")
print(bayes_gas_mass.describe())



#after a few code revisions, 5 galaxies taken out of dust/gas mass fits 
#these give galaxies had zero-value best-fit dust masses, extremely small Bayesian dust masses, and huge chi-squared values
#indicates that these are poor/problematic CIGALE fits
bad_ids = [
    'VFID2143',
    'VFID4051',
    'VFID4311',
    'VFID4723',
    'VFID5625']


#remove nan and inf and -inf values
valid_mask = (
    np.isfinite(Mgas) &
    np.isfinite(bayes_gas_mass) &
    ~merged_df['VFID'].isin(bad_ids))

print(f"\nValid points for regression: {valid_mask.sum()} / {len(merged_df)}")


#make sure there are enough points for a regression
if valid_mask.sum() < 2:
    raise ValueError(
        "There are fewer than two valid galaxies available "
        "for the linear regression.")


#linear regression
slope_bayes, intercept_bayes, r_value, p_value, std_err = linregress(
    Mgas[valid_mask],
    bayes_gas_mass[valid_mask])

print("\nLinear regression results:")
print(f"Slope = {slope_bayes:.3f}")
print(f"Intercept = {intercept_bayes:.3f}")
print(f"R = {r_value:.3f}")
print(f"R^2 = {r_value**2:.3f}")
print(f"p-value = {p_value:.3e}")
print(f"Standard error = {std_err:.3f}")



##### plotting ######

plt.figure(figsize=(8, 6))


#scatter plot of valid galaxies
plt.plot(
    Mgas[valid_mask],
    bayes_gas_mass[valid_mask],
    'o',
    color='tab:blue',
    label='Galaxy')


#create x-values for the fit and the 1-to-1 lines
x_fit = np.linspace(
    Mgas[valid_mask].min(),
    Mgas[valid_mask].max(),
    100)


#plot best fit line
y_fit = slope_bayes * x_fit + intercept_bayes

plt.plot(x_fit, y_fit, color='tab:orange', linestyle='-',
    label=(f'Best Fit Line: y = {slope_bayes:.2f}x + {intercept_bayes:.2f}'))

#plot 1-to-1 line
plt.plot(x_fit, x_fit, color='gray', linestyle='--', label='1-to-1 Line')

plt.xlabel(r'log (VFS Gas Mass) $(M_{\odot})$')
plt.ylabel(r'log (CIGALE Gas Mass) $(M_{\odot})$')
plt.title('CIGALE Gas Mass vs VFS Gas Mass')

plt.legend()

plt.tight_layout()

output_file = 'CIGALE-gmass-vs-vfs-gmass.png'
output_path = plotdir + output_file

plt.savefig(output_path, dpi=150)
plt.show()



# %%

#Tom's notes: next we will calcualte the depletion time scales by using the respective gas masses and
#SFRs. for this current example we use the gas mass from the VFS and SFR from CIGALE and VFS

# calculate depletion times

csv_file = tabledir + 'herschel-testoutput.csv'
fits_file = tabledir + 'vf_v2_halpha.fits'
mgas_csv_file = tabledir + 'reference_gas_mass.csv'

#load cigale data
df1 = pd.read_csv(csv_file, usecols=['id', 'log(bayes.gas.mass)', 'bayes.sfh.sfr'])

#convert CIGALE gas mass from log space to normal space
df1['bayes.gas.mass'] = 10**df1['log(bayes.gas.mass)']


#load H-alpha data
df2 = fits_to_dataframe(fits_file,columns=[
    'VFID', 
    'GAL_LOG_SFR_HA', 
    'GAL_LOG_SFR_HA_ERR', 
    'GAL_LOG_SFR_HA_FLAG'])

#load VFS gas mass
df3 = pd.read_csv(
    mgas_csv_file,
    usecols=['VFID', 'Mgas'])


#clean ids so formatting is the same
df1['id'] = df1['id'].astype(str).str.strip()
df2['VFID'] = df2['VFID'].astype(str).str.strip()
df3['VFID'] = df3['VFID'].astype(str).str.strip()

df1['id'] = (
    df1['id']
    .astype(str)
    .str.strip()
    .str.replace("b'", "", regex=False)
    .str.replace("'", "", regex=False)
)

df2['VFID'] = (
    df2['VFID']
    .astype(str)
    .str.strip()
    .str.replace("b'", "", regex=False)
    .str.replace("'", "", regex=False)
)

df3['VFID'] = (
    df3['VFID']
    .astype(str)
    .str.strip()
    .str.replace("b'", "", regex=False)
    .str.replace("'", "", regex=False)
)



#merge dataframes
merged_df = pd.merge(df1, df2, left_on='id', right_on='VFID')

merged_df = pd.merge(merged_df, df3, on='VFID', how='left')

print("Merged rows:", len(merged_df))


#keep only valid h-alpha measurements
filtered_df = merged_df[merged_df['GAL_LOG_SFR_HA_FLAG'] == True].copy()

print("After H-alpha flag:", len(filtered_df))


#convert SFRs to normal space
filtered_df['SFR1'] = 10**filtered_df['GAL_LOG_SFR_HA']
filtered_df['SFR2'] = 10**filtered_df['bayes.sfh.sfr']


##### calculate depletion times #####

# H-alpha SFR
mask1 = (
    np.isfinite(filtered_df['Mgas']) &
    np.isfinite(filtered_df['SFR1']) &
    (filtered_df['SFR1'] > 0))

filtered_df.loc[mask1, 't_dep1(CO)'] = (
    10**filtered_df.loc[mask1, 'Mgas']
    / filtered_df.loc[mask1, 'SFR1'])


# CIGALE SFR
mask2 = (
    np.isfinite(filtered_df['Mgas']) &
    np.isfinite(filtered_df['SFR2']) &
    (filtered_df['SFR2'] > 0))

filtered_df.loc[mask2, 't_dep2(CIGALE)'] = (
    10**filtered_df.loc[mask2, 'Mgas']
    / filtered_df.loc[mask2, 'SFR2'])


#only keep galaxies with BOTH depletion times
plot_df = filtered_df[
    np.isfinite(filtered_df['t_dep1(CO)']) &
    np.isfinite(filtered_df['t_dep2(CIGALE)']) &
    (filtered_df['t_dep1(CO)'] > 0) &
    (filtered_df['t_dep2(CIGALE)'] > 0)
].copy()

print("Number of valid depletion-time points:", len(plot_df))

print(plot_df[['id', 't_dep1(CO)', 't_dep2(CIGALE)']].head())



# %%

plt.figure(figsize=(7, 7))

plt.scatter(plot_df['t_dep1(CO)'], plot_df['t_dep2(CIGALE)'], color='tab:blue',
    alpha=1, label='Data points')

min_val = min(plot_df['t_dep1(CO)'].min(), plot_df['t_dep2(CIGALE)'].min())

max_val = max(plot_df['t_dep1(CO)'].max(), plot_df['t_dep2(CIGALE)'].max())

#plot 1-to-1 line
plt.plot([min_val, max_val], [min_val, max_val],linestyle='--',
    color='red', label='1-to-1 Line')

plt.xlabel('t_dep (H-alpha) (years)', fontsize=12)
plt.ylabel('t_dep (CIGALE) (years)', fontsize=12)
plt.title('Depletion times (CO GAS)', fontsize=14)

plt.xscale('log')
plt.yscale('log')

plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig(plotdir + 't_dep1_vs_t_dep2_CO_gas.png', dpi=300)

plt.show()
















