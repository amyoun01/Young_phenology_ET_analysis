#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 19:22:37 2021

This code does two things:
1) Brings in monthly precip, tmin, and tmax gridded data to make raster map for 
   aridity index in Fig. 1
2) Gets site-specific AI values for data analysis of SEM results along the aridity 
   gradient. We do this by calculating PET using Thornthwaite's equation as presented
   here: https://upcommons.upc.edu/bitstream/handle/2117/89152/Appendix_10.pdf?sequence=3&isAllowed=y

@author: ay394
"""

import os
#import subprocess

import numpy as np
import pandas as pd

import gdal
import osr
import geopandas as gpd
from pyproj import Transformer

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
nalcms = gdal.Open("nalcms_usa_can_merge.tif")#,gdal.GA_ReadOnly)

# Second, make gridded maps of aridity index (AI = PET/P). PET will be calculated using
# Thornthwaite's equation (1948)

count = 0

for y in yrs[0]:
        
    # PRCP ------------------------
    
    os.chdir(wdir + climdata_dir + "/prcp")
    
    prcp_y = gdal.Open("daymet_v4_prcp_monttl_na_" + str(y) + ".tif",gdal.GA_ReadOnly)
    prcp_arr_y = prcp_y.ReadAsArray()
    prcp_arr_y[prcp_arr_y == -9999] = np.nan
    
    prcp_arr_annsum_y = np.sum(prcp_arr_y,axis=0)
    
    
    # TMAX ------------------------
    
    os.chdir(wdir + climdata_dir + "/tmax")
    
    tmax_y = gdal.Open("daymet_v4_tmax_monavg_na_" + str(y) + ".tif",gdal.GA_ReadOnly)
    tmax_arr_y = tmax_y.ReadAsArray()
    tmax_arr_y[tmax_arr_y == -9999] = np.nan
    
    # TMIN ------------------------

    os.chdir(wdir + climdata_dir + "/tmin")
    
    tmin_y = gdal.Open("daymet_v4_tmin_monavg_na_" + str(y) + ".tif",gdal.GA_ReadOnly)
    tmin_arr_y = tmin_y.ReadAsArray()
    tmin_arr_y[tmin_arr_y == -9999] = np.nan
    
    # Calculate monthly tavg
    
    tavg_arr_y = (tmax_arr_y + tmin_arr_y)/2
    del tmin_y, tmin_arr_y, tmax_y, tmax_arr_y
    
    # -------------------------------------------------------------------------
    # Calculate daylength for PET calculation. Only do this for first iteration 
    # of the loop (i.e. count == 0).
    # -------------------------------------------------------------------------
    
    if count == 0:
        
        projinfo = prcp_y.GetProjection()
        geoTransform = prcp_y.GetGeoTransform()
    
        # This will provide the x and y coordinate limits for the map
        minx = geoTransform[0]
        maxy = geoTransform[3]
        maxx = minx + geoTransform[1] * prcp_y.RasterXSize
        miny = maxy + geoTransform[5] * prcp_y.RasterYSize
        
        # Get Proj4 string for inverse projection using pyproj    
        srs = osr.SpatialReference()
        srs.ImportFromWkt(prcp_y.GetProjection())
        proj_str = srs.ExportToProj4()
        wgs84_str = "+proj=longlat +datum=WGS84 +no_defs"
        
        transformer = Transformer.from_proj(proj_str, wgs84_str)
        x_coords = np.arange(minx,maxx,1000)
        y_coords = np.arange(maxy,miny,-1000)
        
        coord_grid = np.meshgrid(x_coords,y_coords)
        
        x_coords_allpix = coord_grid[0].flatten()
        y_coords_allpix = coord_grid[1].flatten()
        
        latlong = transformer.transform(x_coords_allpix,y_coords_allpix)
        
        lat = latlong[1].reshape((prcp_y.RasterYSize,prcp_y.RasterXSize))
        lat = np.tile(lat,(12,1,1))
        
        tanLat = np.tan(lat/57.2957795); del lat
        tanDelta_values = np.array([-0.37012566, -0.23853358, -0.04679872,  0.16321764, 
                                     0.32930908,  0.40677729,  0.3747741,   0.239063, 
                                     0.04044485, -0.16905776, -0.33306377, -0.40743608])
        
        days_per_month = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
        d = np.zeros((12,prcp_y.RasterYSize,prcp_y.RasterXSize))
        
        for m in range(0,12):
            
            tanLat[m,:,:] = tanLat[m,:,:] * tanDelta_values[m]
            d[m,:,:] + days_per_month[m]        
        
        tanLat[tanLat < -1] = -1  
        tanLat[tanLat > 1] = 1
        
        omega = np.arccos(-tanLat)
        
        N = 24 * omega/np.pi
    
    # -------------------------------------------------------------------------
    # Calculate monthly PET
    # -------------------------------------------------------------------------
    
    T = tavg_arr_y; del tavg_arr_y
    T[T < 0] = 0
    i = (T/5.0)**1.514
    I = np.nansum(i,axis=0); del i
    I[I == 0.0] = np.nan
    
    alpha = 675 * 10**-9 * I**3 - 771 * 10**-7 * I**2 + 1792 * 10**-5 * I + 0.49239
    PET_nc = 16 * ((10*T)/np.tile(I,(12,1,1)))**np.tile(alpha,(12,1,1))
    PET = (N/12) * (d/30) * PET_nc; del PET_nc
    
    count = count + 1
    
    