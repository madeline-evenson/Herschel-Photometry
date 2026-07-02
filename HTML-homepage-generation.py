#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 20:58:15 2026

@author: madeline.evenson
"""
# %%
#import the necessary packages (some might be unnecessary)

from astropy.table import Table
from matplotlib import pyplot as plt
import os
import numpy as np

mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

#set the home path


os.environ['HOME'] = '/Users/madeline.evenson/Research' #general path for all the coding space
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/Virgo/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are

from masking_funct import *
from photometry_funct import *



# %%

#this code block aims to generate the central home HTML page for the project
#located at https://herschelimages.ku.edu/p/Tom/output.html

    #username: herschel
    #password: xyvQHP,wHLoT

#website writes a standard HTML webpage, then iterates over the list of galaxies and generates their path to the specific pages,
#   the legacy survey postage stamps, as well as any of the 3 bands of Herschel for every galaxy
#there is also a B/G significant mask and R significant mask parameter that signifies whenever over 70% of a galaxy is masked or not

# %%

#writes the output html file


#open the HTML file in write mode
with open(os.path.join(homedir, "Virgo/HTML-building/output.html"), "w") as html:
    #write the HTML code line by line
    html.write('<html><body>\n')
    html.write('<title>Herschel images</title>\n')
    html.write('<style type="text/css">\n')
    html.write('table, td, th {padding: 5px; text-align: center; border: 2px solid black;}\n')
    html.write('p {display: inline-block;}\n')
    html.write('</style>\n')
    html.write('<table><tr><th>#</th><th>VFID</th><th>Name</th><th>Legacy Image</th><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th><th>RA</th><th>DEC</th><th>B/G significant mask</th><th>R significant mask</th></tr>')
    
    
    #read the CSV file containing the galaxy data
    galaxy = Table.read(tabledir + '/Photometrytesting2.csv')
    
    for i in range(len(galaxy)):
        n = str(i+1)
        galaxy_name_raw = str(galaxy['GALAXY'][i])
        galaxy_name = galaxy_name_raw.replace(" ", "")
        path = os.path.join(datadir, 'pipeline', galaxy_name_raw)
        VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"
        RA = f"{galaxy['RA_MOMENT'][i]:.3f}"
        DEC = f"{galaxy['DEC_MOMENT'][i]:.3f}"
        
        #compute the percentage of good pixels directly from the aperture columns
        tot_blue_green = 0
        tot_red = 0

        for j in range(1, 9):
            tot_blue_green += galaxy[f"70GoodPixels_SMA_AP0{j}"][i]
            tot_blue_green += galaxy[f"100GoodPixels_SMA_AP0{j}"][i]
            tot_red += galaxy[f"160GoodPixels_SMA_AP0{j}"][i]

        good_percentage_blue_green = tot_blue_green / 16   # 2 bands × 8 apertures
        good_percentage_red = tot_red / 8
        
        #check if values are less than 0.95 and set the mask accordingly
        bg_significant_mask = 'True' if good_percentage_blue_green < 0.95 else 'False'
        r_significant_mask = 'True' if good_percentage_red < 0.95 else 'False'
        
        if os.path.exists(path):
            html.write('<tr><td>'+n+'</td><td><a href="galaxy/html/'+VFID+'-'+galaxy_name+'.html">'+VFID+'</a></td><td>'+galaxy_name+'</td><td><a href="galaxy/png/'+VFID+'-'+galaxy_name+'-LS.jpg"><img src="galaxy/png/'+VFID+'-'+galaxy_name+'-LS.jpg" alt="'+VFID+'-'+galaxy_name+'-LS.jpg" height="auto" width="100%"></a></td><td><a href="galaxy/png/'+VFID+'-'+galaxy_name+'blue.png"><img src="galaxy/png/'+VFID+'-'+galaxy_name+'blue.png" alt="'+VFID+'-'+galaxy_name+'blue.jpg" height="auto" width="100%"></a></td><td><a href="galaxy/png/'+VFID+'-'+galaxy_name+'green.png"><img src="galaxy/png/'+VFID+'-'+galaxy_name+'green.png" alt="'+VFID+'-'+galaxy_name+'green.jpg" height="auto" width="100%"></a></td><td><a href="galaxy/png/'+VFID+'-'+galaxy_name+'red.png"><img src="galaxy/png/'+VFID+'-'+galaxy_name+'red.png" alt="'+VFID+'-'+galaxy_name+'red.jpg" height="auto" width="100%"></a></td><td>'+RA+'</td><td>'+DEC+'</td><td>'+bg_significant_mask+'</td><td>'+r_significant_mask+'</td></tr>\n')
        else:
            continue
        
    html.write('</table>\n')
    html.write('<br /><br />\n')
    html.write('</html></body>\n')
    
print("HTML code has been written to output.html")


# %%

#refer to the photometry notebook for generation of individual HTML files, as I want the images from there first
#there is a codeblock in there that both calculates the photometry AND generates a postage stamp of the galaxies









