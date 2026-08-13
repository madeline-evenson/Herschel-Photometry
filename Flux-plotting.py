#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 11:05:14 2026

@author: madeline.evenson
"""

# this notebook is to start plotting fluxes :)

# %%



#import necessary libraries

import matplotlib.pyplot as plt
import pandas as pd
import os
#import wget


mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are

#make directory if doesn't exist
os.makedirs(plotdir, exist_ok=True)


# %%

csv_file = tabledir + 'Photometrytesting2.csv'

df = pd.read_csv(csv_file)

#print(df.columns)


#super conservative cuts --> can maybe keep 2665/2667, 4177/4179
#overlap_VFID_list = [711, 712, 2665, 2667, 4177, 4179, 4961, 4962]

#df_clean = df.copy()

#removes rows where VF_ID column contains int value in overlap_VFID_list
#df_clean = df[~df['VF_ID'].isin(overlap_VFID_list)]

sma_labels = [f'SMA_AP0{i}' for i in range(1,9)]


# %%

#create and save all the 70 vs. 100 flux plots

for sma in sma_labels:

    fig, ax = plt.subplots(figsize=(8,6))

    x_axis = f'70Flux_{sma}'
    y_axis = f'100Flux_{sma}'

    ax.scatter(df[x_axis], df[y_axis], alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(f'{x_axis} (Jy)')
    ax.set_ylabel(f'{y_axis} (Jy)')
    ax.set_title(x_axis+' vs. '+ y_axis)


    savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

    print('Saving:', savepath)

    plt.savefig(savepath, dpi=150)

    plt.close()
    

#create and save all the 100 vs. 160 flux plots
for sma in sma_labels:

    fig, ax = plt.subplots(figsize=(8,6))

    x_axis = f'100Flux_{sma}'
    y_axis = f'160Flux_{sma}'

    ax.scatter(df[x_axis], df[y_axis], alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(f'{x_axis} (Jy)')
    ax.set_ylabel(f'{y_axis} (Jy)')
    ax.set_title(x_axis+' vs. '+ y_axis)


    savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

    print('Saving:', savepath)

    plt.savefig(savepath, dpi=150)

    plt.close()


#create and save all the 70 vs. 160 flux plots
for sma in sma_labels:

    fig, ax = plt.subplots(figsize=(8,6))

    x_axis = f'70Flux_{sma}'
    y_axis = f'160Flux_{sma}'

    ax.scatter(df[x_axis], df[y_axis], alpha=0.5)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel(f'{x_axis} (Jy)')
    ax.set_ylabel(f'{y_axis} (Jy)')
    ax.set_title(x_axis+' vs. '+ y_axis)


    savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

    print('Saving:', savepath)

    plt.savefig(savepath, dpi=150)

    plt.close()
    
 
    
# %%
 
#calculate total flux across all apertures
df['70_total_flux'] = 0
df['100_total_flux'] = 0
df['160_total_flux'] = 0

for sma in sma_labels:
    df['70_total_flux'] += df[f'70Flux_{sma}']
    df['100_total_flux'] += df[f'100Flux_{sma}']
    df['160_total_flux'] += df[f'160Flux_{sma}']



#plot 70 vs 100 total flux
fig, ax = plt.subplots(figsize=(8,6))

x_axis = '70_total_flux'
y_axis = '100_total_flux'

ax.scatter(df[x_axis], df[y_axis], alpha = 0.5)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel('70 Total Flux (Jy)')
ax.set_ylabel('100 Total Flux (Jy')
ax.set_title('70 vs. 100 Total Flux')

savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

print('Saving:', savepath)

plt.savefig(savepath, dpi=150)

plt.close()



#plot 100 vs 160 total flux
fig, ax = plt.subplots(figsize=(8,6))

x_axis = '100_total_flux'
y_axis = '160_total_flux'

ax.scatter(df[x_axis], df[y_axis], alpha = 0.5)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel('100 Total Flux (Jy)')
ax.set_ylabel('160 Total Flux (Jy')
ax.set_title('100 vs. 160 Total Flux')

savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

print('Saving:', savepath)

plt.savefig(savepath, dpi=150)

plt.close()



#plot 70 vs 160 total flux
fig, ax = plt.subplots(figsize=(8,6))

x_axis = '70_total_flux'
y_axis = '160_total_flux'

ax.scatter(df[x_axis], df[y_axis], alpha = 0.5)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel('70 Total Flux (Jy)')
ax.set_ylabel('160 Total Flux (Jy')
ax.set_title('70 vs. 160 Total Flux')

savepath = os.path.join(plotdir, f"{x_axis}-vs-{y_axis}.png")

print('Saving:', savepath)

plt.savefig(savepath, dpi=150)

plt.close()





