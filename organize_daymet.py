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
#import subprocess

import numpy as np
import pandas as pd

import gdal
import geopandas as gpd

wdir = "/Volumes/GoogleDrive/My Drive/Young_phenology_ET_analysis"
wdir_shell_str = "/Volumes/GoogleDrive/My\ Drive/Young_phenology_ET_analysis"

climdata_dir = "/data/raw_data/daymet"
climvars = ("prcp_monttl","tmin_monavg","tmax_monavg")
yrs = range(1980,2011)

os.chdir(wdir + "/data/ancillary_data/")
site_info = pd.read_csv('pheno_flux_sites_to_use.csv')
site_info_geo = gpd.GeoDataFrame(site_info, 
                                 geometry=gpd.points_from_xy(site_info.lon, site_info.lat))

# First, import NALCMS merged map for USA and Canada and then use this projection
# info to trim and reproject the 

os.chdir(wdir + "/data/raw_data/NALCMS/")
nalcms = gdal.Open("nalcms_usa_can_merge.tif",gdal.GA_ReadOnly)
projinfo = nalcms.GetProjection()
geoTransform = nalcms.GetGeoTransform()

# This will provide the x and y coordinate limits for the map
minx = geoTransform[0]
maxy = geoTransform[3]
maxx = minx + geoTransform[1] * nalcms.RasterXSize
miny = maxy + geoTransform[5] * nalcms.RasterYSize

# Second, make gridded maps of aridity index (AI = PET/P). PET will be calculated using
# Thornthwaite's equation (1948)
for y in yrs:
        
    os.chdir(wdir + climdata_dir + "/prcp")
    
    prcp_y = gdal.Open("daymet_v4_prcp_monttl_na_" + str(y) + ".tif",gdal.GA_ReadOnly)
    prcp_arr_y = prcp_y.ReadAsArray()
    prcp_arr_y[prcp_arr_y == -9999] = np.nan
    
    prcp_arr_annsum_y = np.sum(prcp_arr_y,axis=0)
    
    os.chdir(wdir + climdata_dir + "/prcp")
        
        fn_i = "daymet_v4_" + climvars[v] + "_na_" + str(y) + ".tif"
        rst_i = gdal.Open(fn_i,gdal.GA_ReadOnly)
        
        projinfo = rst_i.GetProjection()
        geoTransform = rst_i.GetGeoTransform()
        
        rst_arr_i = rst_i.ReadAsArray()
        rst_arr_i[rst_arr_i == -9999] = np.nan
                
        for m in range(0,12):
            
            rst_arr_im = rst_arr_i[m,:,:]

                
            
