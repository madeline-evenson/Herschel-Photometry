#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 21:07:02 2026

@author: madeline.evenson
"""


# %%



#import necessary libraries

import matplotlib.pyplot as plt


from astropy.table import Table
import os

import warnings
warnings.filterwarnings('ignore')


mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']


#set the home path
os.environ['HOME'] = '/Users/madeline.evenson/Research'
homedir = os.getenv('HOME')
tabledir = homedir+'/Virgo/tables/'
plotdir = homedir+'/Virgo/plots/'
htmldir = homedir+'/Virgo/HTML-building/galaxy/' #set to where the html resources should be (parent folder of the different html folders)
datadir = homedir+'/masking/' #set to where the completed mask fits are




# %%


#this is the code to write out individual HTML files for the galaxies. these can be accessed using the link in the HTML website generation notebook


#individual html writing
galaxy = Table.read(tabledir+'/Photometrytesting2.csv')

for i in range(len(galaxy)):
    #check bounds for i+1
    if i + 1 < len(galaxy):
        galaxy_name2 = str(galaxy['GALAXY'][i+1])
        VFID2 = f"VFID{int(galaxy['VF_ID'][i+1]):04d}"
    else:
        galaxy_name2 = '0'
        VFID2 = '0'
        
    #check bounds for i-1
    if i - 1 >= 0:
        galaxy_name3 = str(galaxy['GALAXY'][i-1])
        VFID3 = f"VFID{int(galaxy['VF_ID'][i-1]):04d}"
    else:
        galaxy_name3 = '0'
        VFID3 = '0'
        
    galaxy_name_raw = str(galaxy['GALAXY'][i])
    galaxy_name = galaxy_name_raw.replace(" ", "")
    path = os.path.join(datadir, 'pipeline', galaxy_name_raw)
    VFID = f"VFID{int(galaxy['VF_ID'][i]):04d}"
    RA = str(galaxy['RA_MOMENT'][i])
    DEC = str(galaxy['DEC_MOMENT'][i])
    
    filepath = os.path.join(htmldir, 'html', VFID + '-' + galaxy_name + '.html')
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    test_path = os.path.join(htmldir, 'png', VFID+'-'+galaxy_name+'blue.png')
    print("Checking:", test_path)
    
    #write the HTML code line by line
    with open(filepath, "w") as html:
        html.write('<html><body>\n')
        html.write('<h1>'+VFID+'</h1>')
        
        if VFID3 =='0':
            html.write('<a href="' + homedir + '/HTML-building/output.html">Home</a>')
            html.write('<br />')
        else:
            html.write('<a href="' + VFID3 + '-' + galaxy_name3 + '.html">Previous'+VFID3+'</a>')
            html.write('<br />')
            
        if VFID2 == '0':
            html.write('<a href="' + homedir + '/HTML-building/output.html">Home</a>')
            html.write('<br />')
        else:
            html.write('<a href="' + VFID2 + '-' + galaxy_name2 + '.html">Next' + VFID2 + '</a>')
            html.write('<br />')
            
        html.write('<h2>Galaxy Data:</h2>\n')
        html.write('<style type="text/css">\n')
        html.write('table, td, th {padding:5px; text-align:center; border:2px solid black;}\n')
        html.write('table {margin-top:0; margin-bottom:15px; border-collapse:collapse;}\n')
        html.write('h2 {margin-top:20px; margin-bottom:6px;}\n')
        html.write('p {display:inline-block;}\n')
        html.write('</style>\n')
        html.write('<table width="100%"><tr><th>VFID</th><th>Name</th><th>Legacy image</th><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th><th>RA</th><th>DEC</th></tr>')
        html.write('<tr><td>'+VFID+'</td><td>'+galaxy_name+'</td><td><a href="../png/'+VFID+'-'+galaxy_name+'-LS.jpg"><img src="../png/'+VFID+'-'+galaxy_name+'-LS.jpg" alt="No LS data.jpg" height="auto" width="100%"></a></td><td><a href="../png/'+VFID+'-'+galaxy_name+'blue.png"><img src="../png/'+VFID+'-'+galaxy_name+'blue.png" alt="Missing file 70microns.jpg" height="auto" width="100%"></a></td><td><a href="../png/'+VFID+'-'+galaxy_name+'green.png"><img src="../png/'+VFID+'-'+galaxy_name+'green.png" alt="Missing file 100microns.jpg" height="auto" width="100%"></a></td><td><a href="../png/'+VFID+'-'+galaxy_name+'red.png"><img src="../png/'+VFID+'-'+galaxy_name+'red.png" alt="Missing file 160microns.jpg" height="auto" width="100%"></a></td><td>'+RA+'</td><td>'+DEC+'</td>\n')
        html.write('</table>\n')
        html.write('<h2>AP06 Aperture on Masked Herschel Images:</h2>')
        html.write('<table width="100%"><tr><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th></tr>')
        html.write('<tr><td><a href="../png/'+VFID+'-'+galaxy_name+'blue.png"><img src="../png/'+VFID+'-'+galaxy_name+'blue.png" alt="Missing file 70microns.jpg" height="auto" width="100%"></a></td><td><a href="../png/'+VFID+'-'+galaxy_name+'green.png"><img src="../png/'+VFID+'-'+galaxy_name+'green.png" alt="Missing file 100microns.jpg" height="auto" width="100%"></a></td><td><a href="../png/'+VFID+'-'+galaxy_name+'red.png"><img src="../png/'+VFID+'-'+galaxy_name+'red.png" alt="Missing file 160microns.jpg" height="auto" width="100%"></a></td></tr> \n')
        html.write('</table>\n')
        html.write('<h2>All 8 Apertures on Legacy Survey Images:</h2>')
        html.write('<table width="100%"><tr><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th></tr>')
        html.write('<tr><td><a href="../AP/'+VFID+'-'+galaxy_name+'-70-AP.png"><img src="../AP/'+VFID+'-'+galaxy_name+'-70-AP.png" alt="Missing file 70microns.jpg" height="auto" width="100%"></a></td><td><a href="../AP/'+VFID+'-'+galaxy_name+'-100-AP.png"><img src="../AP/'+VFID+'-'+galaxy_name+'-100-AP.png" alt="Missing file 100microns.jpg" height="auto" width="100%"></a></td><td><a href="../AP/'+VFID+'-'+galaxy_name+'-160-AP.png"><img src="../AP/'+VFID+'-'+galaxy_name+'-160-AP.png" alt="Missing file 160microns.jpg" height="auto" width="100%"></a></td></tr> \n')
        html.write('</table>\n')
        html.write('<h2>All 8 Apertures on Masked Herschel Images:</h2>')
        html.write('<table width="100%"><tr><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th></tr>')
        html.write('<tr><td><a href="../mask/'+VFID+'-'+galaxy_name+'-70-AP.png"><img src="../mask/'+VFID+'-'+galaxy_name+'-70-AP.png" alt="Missing file 70microns.jpg" height="auto" width="100%"></a></td><td><a href="../mask/'+VFID+'-'+galaxy_name+'-100-AP.png"><img src="../mask/'+VFID+'-'+galaxy_name+'-100-AP.png" alt="Missing file 100microns.jpg" height="auto" width="100%"></a></td><td><a href="../mask/'+VFID+'-'+galaxy_name+'-160-AP.png"><img src="../mask/'+VFID+'-'+galaxy_name+'-160-AP.png" alt="Missing file 160microns.jpg" height="auto" width="100%"></a></td></tr> \n')
        html.write('</table>\n')
        html.write('<h2>Background Comparison:</h2>')
        html.write('<table width="100%"><tr><th>Comparison of Raw Background to Plane-Subtracted Background</th></tr>')
        html.write('<tr><td><a href="../bg-comparison-plots/'+VFID+'-'+galaxy_name+'_bg_comparison.png"><img src="../bg-comparison-plots/'+VFID+'-'+galaxy_name+'_bg_comparison.png" alt="No background data.png" height="auto" width="100%"></td></tr>')
        html.write('</table>\n')
        html.write('<h2>Background-Subtracted Flux Histograms:</h2>')
        html.write('<table width="100%"><tr><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th></tr>')
        html.write('<tr><td><a href="../psub-bg-flux-hists/'+galaxy_name+'_B_psub_bg_flux_hist.png"><img src="../psub-bg-flux-hists/'+galaxy_name+'_B_psub_bg_flux_hist.png" alt="Missing file 70microns.png" height="auto" width="100%"></a></td><td><a href="../psub-bg-flux-hists/'+galaxy_name+'_G_psub_bg_flux_hist.png"><img src="../psub-bg-flux-hists/'+galaxy_name+'_G_psub_bg_flux_hist.png" alt="Missing file 100microns.png" height="auto" width="100%"></a></td><td><a href="../psub-bg-flux-hists/'+galaxy_name+'_R_psub_bg_flux_hist.png"><img src="../psub-bg-flux-hists/'+galaxy_name+'_R_psub_bg_flux_hist.png" alt="Missing file 160microns.png" height="auto" width="100%"></a></td></tr> \n')
        html.write('</table>\n')
        html.write('<h2>Background-Subtracted Vertical Cut Profiles:</h2>')
        html.write('<table width="100%"><tr><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th></tr>')
        html.write('<tr><td><a href="../psub-bg-v-cut-profiles/'+galaxy_name+'_B_vertical_profiles.png"><img src="../psub-bg-v-cut-profiles/'+galaxy_name+'_B_vertical_profiles.png" alt="No good columns for '+galaxy_name+' B cuts" height="auto" width="100%"></a></td><td><a href="../psub-bg-v-cut-profiles/'+galaxy_name+'_G_vertical_profiles.png"><img src="../psub-bg-v-cut-profiles/'+galaxy_name+'_G_vertical_profiles.png" alt="No good columns for '+galaxy_name+' G cuts" height="auto" width="100%"></a></td><td><a href="../psub-bg-v-cut-profiles/'+galaxy_name+'_R_vertical_profiles.png"><img src="../psub-bg-v-cut-profiles/'+galaxy_name+'_R_vertical_profiles.png" alt="No good columns for '+galaxy_name+' R cuts" height="auto" width="100%"></a></td></tr> \n')
        html.write('</table>\n')
        html.write('<h2>Background-Subtracted Horizontal Cut Profiles:</h2>')
        html.write('<table width="100%"><tr><th>Herschel-UnimapBlue (70microns)</th><th>Herschel-UnimapGreen (100microns)</th><th>Herschel-UnimapRed (160microns)</th></tr>')
        html.write('<tr><td><a href="../psub-bg-h-cut-profiles/'+galaxy_name+'_B_horizontal_profiles.png"><img src="../psub-bg-h-cut-profiles/'+galaxy_name+'_B_horizontal_profiles.png" alt="No good rows for '+galaxy_name+' B cuts" height="auto" width="100%"></a></td><td><a href="../psub-bg-h-cut-profiles/'+galaxy_name+'_G_horizontal_profiles.png"><img src="../psub-bg-h-cut-profiles/'+galaxy_name+'_G_horizontal_profiles.png" alt="No good rows for '+galaxy_name+' G cuts" height="auto" width="100%"></a></td><td><a href="../psub-bg-h-cut-profiles/'+galaxy_name+'_R_horizontal_profiles.png"><img src="../psub-bg-h-cut-profiles/'+galaxy_name+'_R_horizontal_profiles.png" alt="No good rows for '+galaxy_name+' R cuts" height="auto" width="100%"></a></td></tr> \n')
        html.write('</table>\n')
        html.write('<h2>Photometry Results:</h2>')
        html.write('<table width="100%"><tr><th>Photometry profiles</th></tr>')
        html.write('<tr><td><a href="../photometry/'+VFID+'-'+galaxy_name+'.png"><img src="../profiles/'+VFID+'-'+galaxy_name+'_flux_sb_profile.png" alt="No photometry data.png" height="auto" width="100%"></td></tr>')
        html.write('</table>')
        html.write('<br /><br />\n')
        
        if VFID3 == '0':
            html.write('<a href="'+homedir+'/HTML-building/output.html">Home</a>')
            html.write('<br />')
        else:
            html.write('<a href="'+VFID3+'-'+galaxy_name3+'.html">Previous'+VFID3+'</a>')
            html.write('<br />')
        
        if VFID2 == '0':
            html.write('<a href="'+homedir+'/HTML-building/output.html">Home</a>')
            html.write('<br />')
        else:
            html.write('<a href="'+VFID2+'-'+galaxy_name2+'.html">Next'+VFID2+'</a>')
            html.write('<br />')
            
        html.write('</body></html>\n')
        html.close()
        
print("done")






