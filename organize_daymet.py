#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 19:22:37 2021

This code does two things:
1) Brings in monthly precip, tmin, and tmax gridded data to make raster map for 
   aridity index in Fig. 1
2) Gets site-specific AI values for data analysis of SEM results along the aridity 
   gradient

@author: ay394
"""

import os
import subprocess

import numpy as np
import pandas as pd

import gdal
import rasterio
import rasterio.plot
import geopandas as gpd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.legend_handler import HandlerPatch
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
from matplotlib import gridspec

from adjustText import adjust_text


wdir = "/Volumes/GoogleDrive/My Drive/Young_phenology_ET_analysis"
wdir_shell_str = "/Volumes/GoogleDrive/My\ Drive/Young_phenology_ET_analysis"

climdata_dir = "/data/raw_data/daymet"

os.chdir(wdir + "/data/ancillary_data/") 
site_info = pd.read_csv('pheno_flux_sites_to_use.csv')
site_info_geo = gpd.GeoDataFrame(site_info, 
                                       geometry=gpd.points_from_xy(site_info.lon, site_info.lat))

# First, make gridded maps of aridity index (AI). PET will be calculated using
# Thornwaite's equation (1948) 

os.chdir(wdir + climdata_dir)
fn = os.listdir()
