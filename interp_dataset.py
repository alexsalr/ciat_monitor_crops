# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 15:00:06 2018

@author: ASALAZAR
"""

import sys
import os
import re
import xarray as xr
import numpy as np

sys.path.append('b_Temporal_Stack')
sys.path.append('c_Class_Models')

import regionstack
import interpolatets as itp


from dask.distributed import Client
#import eotempstack
#import xr_eotemp

def main_optical(dates):
    
    saldana = regionstack.regionStack('Saldana', attrs=['S2','LC08'])
    
    #saldana.harmonize_L8()
    
    boi = ['red','blue','green','nir','swir1','swir2']
    
    vars_loc = os.environ['WIN_SVR_DATA']+'Saldana/vars/'
    
    for band in boi:
        if not os.path.isfile(vars_loc+'opt_'+band+'.nc'):
            
            ## TO-DO add isel for relevant dates
            sentinel = saldana.S2[band]#.isel(time=slice(min(dates)-120?D-security_margin, max(dates)+16D+security_margin))
            
            landsat = saldana.LC08[band]#.isel(time=slice(min(dates)-120?D, max(dates)+16D))
            
            aligned = xr.align(sentinel, landsat, exclude={'time'})
            
            da = xr.concat([aligned[0], aligned[1]], dim='time')
            
            da.sortby('time').to_netcdf(vars_loc+'opt_'+band+'.nc')
    
    client = Client(n_workers=12)
    
    client.upload_file('c_Class_Models/interpolatets.py')
    
    files = list(filter(re.compile(r'^opt_.*').search, os.listdir(vars_loc)))
    files = list(map(lambda x: vars_loc+x, files))
    
    print('Reading concatenated Dataset')
    dataset = xr.open_mfdataset(files,chunks={'y':1000,'x':750,'time':-1})
    
    print('Concat Dataset\n')
    print(dataset)
    
    _location = os.environ['WIN_SVR_DATA']+'Saldana/features/'
    
    itp.interpolate_dataset(dataset,_location,boi,
                            date_of_analysis=dates, der=False)
    
    print('Dataset interpolation done!\n')
    ##
    
def der_optical(dates):

    client = Client(n_workers=12)
    
    boi = ['red','blue','green','nir','swir1','swir2']
    
    vars_loc = os.environ['WIN_SVR_DATA']+'Saldana/vars/'
    
    files = list(filter(re.compile(r'^opt_.*').search, os.listdir(vars_loc)))
    files = list(map(lambda x: vars_loc+x, files))
    
    dataset = xr.open_mfdataset(files,chunks={'y':1000,'x':750,'time':-1})
    
    print('Concat Dataset\n')
    print(dataset)
    
    _location = os.environ['WIN_SVR_DATA']+'Saldana/features/'
    
    itp.interpolate_dataset(dataset, _location, boi,
                            date_of_analysis=dates, der=True)
    
    print('Dataset interpolation done!\n')

def main_radar(dates):

    saldana = regionstack.regionStack('Saldana', attrs=['S1_ASCENDING','S1_DESCENDING'])
    
    boi = ['VV']
    
    vars_loc = os.environ['WIN_SVR_DATA']+'Saldana/vars/'
    
    for band in boi:
        
        if not os.path.isfile(vars_loc+'rad_'+band+'.nc'):
            asc = saldana.S1_ASCENDING[band]
            
            dsc = saldana.S1_DESCENDING[band]
            
            da = xr.concat([asc, dsc], dim='time')
            
            print('{} band was concatenated. Writing DataArray'.format(band))
            
            da.to_netcdf(vars_loc+'rad_'+band+'.nc')
    
    client = Client(n_workers=12)
    
    client.upload_file('b_Temporal_Stack/interpolatets.py')
    
    files = list(filter(re.compile(r'^rad_.*').search, os.listdir(vars_loc)))
    files = list(map(lambda x: vars_loc+x, files))
    
    print('Reading concatenated Dataset')
    dataset = xr.open_mfdataset(files,chunks={'y':1000,'x':750,'time':-1})
    
    print('Concat Dataset\n')
    print(dataset)
    
    _location = os.environ['WIN_SVR_DATA']+'Saldana/features/'
    
    itp.interpolate_dataset(dataset, _location, boi,
                            date_of_analysis=dates)
    
    itp.interpolate_dataset(dataset, _location, boi,
                            date_of_analysis=dates, der=True)
    

def main_radar_text(dates):
    
    
    
    saldana = regionstack.regionStack('Saldana', attrs=['S1_ASCENDING_GLCM','S1_DESCENDING_GLCM'])
    
    boi = ['VV_ASM','VV_Contrast','VV_Dissimilarity','VV_Energy','VV_Entropy',
           'VV_GLCMCorrelation','VV_GLCMMean','VV_GLCMVariance','VV_Homogeneity']
    
    vars_loc = os.environ['WIN_SVR_DATA']+'Saldana/vars/'
    
    for band in boi:
        
        if not os.path.isfile(vars_loc+'rad_'+band+'.nc'):
            asc = saldana.S1_ASCENDING_GLCM[band]
            
            dsc = saldana.S1_DESCENDING_GLCM[band]
            
            da = xr.concat([asc, dsc], dim='time')
            
            print('{} band was concatenated. Writing DataArray'.format(band))
            
            da.to_netcdf(vars_loc+'rad_'+band+'.nc')
    
    client = Client(n_workers=12)
    
    client.upload_file('b_Temporal_Stack/interpolatets.py')
    
    files = list(filter(re.compile(r'^rad_.*').search, os.listdir(vars_loc)))
    files = list(map(lambda x: vars_loc+x, files))
    
    print('Reading concatenated Dataset')
    dataset = xr.open_mfdataset(files,chunks={'y':1000,'x':750,'time':-1})
    
    print('Concat Dataset\n')
    print(dataset)
    
    _location = os.environ['WIN_SVR_DATA']+'Saldana/features/'
    
    itp.interpolate_dataset(dataset, _location, boi,
                            date_of_analysis=dates)
    
    itp.interpolate_dataset(dataset, _location, boi,
                            date_of_analysis=dates, der=True)
    
    #client.close()
    
if __name__ == '__main__':
    
    dates = [np.datetime64('2015-11-27'),np.datetime64('2015-12-21'),np.datetime64('2016-01-07')]
    
    #main_radar(dates)
    
    main_optical(dates)
    
    #der_optical(dates)
    
    #main_radar_text(dates)
    