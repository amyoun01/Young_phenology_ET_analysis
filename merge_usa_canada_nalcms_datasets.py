#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:16:57 2021

@author: ay394
"""


import os
import subprocess

import gdal

wdir = "/Volumes/GoogleDrive/My Drive/Young_phenology_ET_analysis"
wdir_shell_str = "/Volumes/GoogleDrive/My\ Drive/Young_phenology_ET_analysis"

usa_nalcms_fn = wdir + "/data/raw_data/NALCMS/united_states_2015/USA_NALCMS_2015_LC_30m_LAEA_mmu5pix_/USA_NALCMS_2015_LC_30m_LAEA_mmu5pix_.tif"
can_nalcms_fn = wdir + "/data/raw_data/NALCMS/canada_2015/CAN_NALCMS_2015_LC_30m_LAEA_mmu5pix_/CAN_NALCMS_2015_LC_30m_LAEA_mmu5pix_.tif"

usa_1k_fn = wdir + "/data/raw_data/NALCMS/usa_1km.tif"
can_1k_fn = wdir + "/data/raw_data/NALCMS/can_1km.tif"

usa_nalcms = gdal.Open(usa_nalcms_fn,gdal.GA_ReadOnly)
can_nalcms = gdal.Open(can_nalcms_fn,gdal.GA_ReadOnly)

geoTransform = usa_nalcms.GetGeoTransform()
minx = geoTransform[0]
maxy = geoTransform[3]
maxx = minx + geoTransform[1] * usa_nalcms.RasterXSize
miny = maxy + geoTransform[5] * usa_nalcms.RasterYSize

can_1k_ds = gdal.Warp(can_1k_fn, \
                           can_nalcms, \
                           format="GTiff", \
                           xRes = 1000.0, \
                           yRes = 1000.0, \
                           outputBounds = (minx,miny,maxx,maxy), \
                           srcSRS = can_nalcms.GetProjection(), \
                           dstSRS = usa_nalcms.GetProjection(), \
                           resampleAlg="near")

can_1k_ds = None

usa_1k_ds = gdal.Warp(usa_1k_fn, \
                           usa_nalcms, \
                           format="GTiff", \
                           xRes = 1000.0, \
                           yRes = 1000.0, \
                           outputBounds = (minx,miny,maxx,maxy), \
                           srcSRS = can_nalcms.GetProjection(), \
                           dstSRS = usa_nalcms.GetProjection(), \
                           resampleAlg="near")

usa_1k_ds = None

os.chdir(wdir + "/data/raw_data/NALCMS")

nalcms_fn = "nalcms_usa_can_merge.tif"

rcalc_command = "gdal_calc.py -A usa_1km.tif " + " -B can_1km.tif " + " --outfile=" + nalcms_fn + " --NoDataValue=0.0 --quiet --calc=A+B"
subprocess.call(rcalc_command,shell=True)