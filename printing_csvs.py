#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 01:19:55 2026

@author: madeline.evenson
"""

from astropy.table import Table
galaxy = Table.read('/Users/madeline.evenson/Research/Virgo/tables/Photometrytesting2.csv')
print(galaxy.colnames)