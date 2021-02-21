#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 13:11:42 2021

@author: adam
"""
import os

import numpy as np
import pandas as pd

import osr, gdal
import rasterio
import pyproj

wdir = "/Volumes/GoogleDrive/My Drive/Young_phenology_ET_analysis"

wgs84_str = "+proj=longlat +datum=WGS84 +no_defs"

os.chdir(wdir + "/data/ancillary_data")
phenoflux_metadata = pd.read_csv("pheno_flux_sites_to_use.csv")

os.chdir(wdir + "/data/raw_data/daymet")

aridity_1980_2009 = rasterio.open("aridity_mean_1980_2009.tif")
aridity_1980_2009_arr = aridity_1980_2009.read()
aridity_1980_2009_arr = aridity_1980_2009_arr.squeeze()

crs = aridity_1980_2009.crs.to_proj4()

transformer = Transformer.from_proj(wgs84_str,crs)

x,y = transformer.transform(np.array(phenoflux_metadata.lon),np.array(phenoflux_metadata.lat))

row,col = aridity_1980_2009.index(x,y)

aridity_vals = aridity_1980_2009_arr[row,col]
