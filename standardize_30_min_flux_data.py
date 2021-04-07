#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 08:26:30 2021

@author: ay394
"""

import os
import re

from datetime import datetime

import pandas as pd
import numpy as np

# Set Working Directory
wdir = "/Volumes/GoogleDrive/My Drive/Young_phenology_ET_analysis"

# Add path to load in custom written functions to workspace
# addpath([wdir,'/code/z_functions']);

os.chdir(wdir + "/data/ancillary_data")

# Load in metadata for each study site
phenocam_flux_metadata_table = pd.read_csv("pheno_flux_sites_to_use.csv")

sites = phenocam_flux_metadata_table.fluxsite
start_date = pd.to_datetime(phenocam_flux_metadata_table.start_date)
end_date = pd.to_datetime(phenocam_flux_metadata_table.end_date)

var_names = pd.read_csv("variables_to_import_for_fluxsites.csv")
vars_of_interest = var_names.columns[1:len(var_names.columns)]

AMF_file_names = sorted(os.listdir(wdir + "/data/raw_data/ameriflux/BASE"))

for i in range(0,len(sites)):
    
    os.chdir(wdir + "/data/raw_data/ameriflux/BASE")

    fluxdat = pd.read_csv(AMF_file_names[i],skiprows=[0,1])
    fluxdat[fluxdat == -9999] = np.nan
    
    datetime_start = pd.to_datetime(fluxdat.TIMESTAMP_START,format='%Y%m%d%H%M')
    datetime_end   = pd.to_datetime(fluxdat.TIMESTAMP_END,format='%Y%m%d%H%M')

    dt_id = (datetime_start > start_date[i]) & (datetime_end < end_date[i])
    
    fluxdat = fluxdat.loc[dt_id]
    
    raw_data_var_names = fluxdat.columns
    
    fluxdat_subset = []
    fluxdat_subset['datetime_start'] = np.nan datetime_start
    
    for k in range(0,len(vars_of_interest)):
        
        if vars_of_interest[k] == "precip":
            continue            
        
        column_id = k + 1
        
        varname_to_import = var_names.iloc[i,column_id]
        
        if varname_to_import == 'NA' :
            
            fluxdat_subset[varname_to_import] = np.nan
            
        else:
            
            fluxdat_subset[varname_to_import] = fluxdat.(varname_to_import);
 
            fluxdat_subset.(char(vars_of_interest(k))) = 
            
        end
        
    end
    
    fluxdat =  fluxdat_filtered; clear  fluxdat_filtered;
    
    sum_netrad_component_values = ...
        fluxdat.sw_in - fluxdat.sw_out + fluxdat.lw_in - fluxdat.lw_out;
    
    both_nan_idx = isnan(fluxdat.netrad) & (isnan(sum_netrad_component_values) == false);
    fluxdat.netrad(both_nan_idx) = sum_netrad_component_values(both_nan_idx);
            
    fluxdat.datetime_start = flux_datetime_start;
    fluxdat.datetime_end = flux_datetime_end;
    
    fluxdat = fluxdat(:,[end-1,end,1:end-2]);
    
    export_filename = ...
        sprintf('%s//results//2_filtered_flux_data//%s_%s.csv',...
                wdir,char(sites(i)),timestep);
    
    % Export data tables
    writetable(fluxdat,export_filename);
    
end
# End of script ----------------------------------------------------------------