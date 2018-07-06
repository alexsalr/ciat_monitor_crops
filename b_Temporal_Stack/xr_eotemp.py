# -*- coding: utf-8 -*-
"""
Created on Wed May 23 09:31:30 2018

@author: ASALAZAR
"""

import math
import xarray as xr
import pandas as pd
import numpy as np
import datetime as dt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

@xr.register_dataset_accessor('eotemp')
class EOTempDataset(object):
    
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        
    def plot_RGB(self, bands=['red', 'green', 'blue'], factor=0.0001, writefile = False, filename='rgb.pdf'):
        
        crs = ccrs.UTM('18N')
        av_dates = self._obj.coords['time'].data.tolist()
        
        # Set plot dimensions
        plt.figure(figsize=(15,3*math.ceil(len(av_dates)/5.0)))
        
        # Try to create the rgb xarray
        try:
            arrays = []
            for band in bands:
                arrays.append(self._obj[band])
            rgb = xr.concat(arrays, pd.Index(bands, name='band')).sortby('time')
            rgb.name = 'rgb'
        except:
            print('The EOTempArray does not contain rgb bands')
            return None
        
        # Plot every date
        for ix, date in enumerate(av_dates):
            ap = rgb.isel(time=ix)*factor
            ax = plt.subplot(math.ceil(len(av_dates)/5.0),5,1+ix, projection=crs)
            ap.plot.imshow(rgb='band', transform=crs)
        
        # Write plot
        if writefile:
            plt.savefig(filename, dpi=300)
        
        plt.show()
    
    def plot_raster(self, band, factor=1, writefile = False, filename='rgb.pdf', **kwargs):
    
        crs = ccrs.UTM('18N')
        av_dates = self._obj.coords['time'].data.tolist()
        
        # Set plot dimensions
        plt.figure(figsize=(15,3*math.ceil(len(av_dates)/5.0)))
        
        # Plot every date
        for ix, date in enumerate(av_dates):
            if factor is not 1:
                da = self._obj.sel(time=date)[band]*factor
            else:
                da = self._obj.sel(time=date)[band]
            #ap = rgb.isel(time=ix)*factor
            ax = plt.subplot(math.ceil(len(av_dates)/5.0),5,1+ix, projection=crs)
            #ap.plot.imshow(rgb='band', transform=crs)
            plt.imshow(da.data, **kwargs)
            
        # Write plot
        if writefile:
            plt.savefig(filename, dpi=300)
        
        plt.show()
        
#    def calcTempTrend(self, band, ndate=-1, tempwindowsize=1):
#        """
#        Calculates the temporal trend of a given band. Returns a xa DataArray with the calulated trend.
#        
#        @params
#            band (str): name of the band 
#            ndate (int): index of the date of analysis, default is -1, i.e. last available date
#            tempwindowsize (int): size of the temporal window in available dates distance
#                                  default value is 1, i.e. the previous available date
#        """
#        # Make temporal subset from last date and the tempwindowsize distance
#        try:
#            temp_subset = self._obj.sortby('time').isel(time=[ndate-tempwindowsize,ndate])
#        except:
#            raise IndexError('The specified time indices are not in the dataset')
#        try:
#            temp_subset = temp_subset[band]
#        except:
#            raise IndexError('Band {} is not a dataset variable'.format(band))
#        
#        # Extract dates
#        time = pd.DatetimeIndex(temp_subset.time.values).date
#        time = np.array([x.toordinal() for x in time])
#        
#        # Mask the subset when mask is available
#        try:
#            temp_subset = temp_subset.where(temp_subset.mask)
#        except:
#            pass
#        
#        # Extract band values
#        vals = temp_subset.values
#        # Reshape to an array with as many rows as dates and as many columns as there are pixels
#        vals_reshaped = vals.reshape(len(time), -1)
#        
#        # Declare array to store slope coefficient
#        slope = np.empty(vals_reshaped[0].shape)
#        slope[:] = np.nan
#        
#        # Make a mask of missing values (if any)in any of the dates
#        mask = np.any(np.isnan(vals_reshaped), axis=0)
#        # Use mask to get only valid pairs
#        masked_vals = vals_reshaped.T[~mask].T
#        
#        try:
#            # Do a first-degree polyfit with masked vals
#            polyfit = np.polyfit(time, masked_vals, 1)
#            # Put slope valid values in slope array
#            np.place(slope, ~mask, polyfit[0].tolist())
#        except np.linalg.LinAlgError:
#            print('Value pairs are all nAn. No good quality pixels are available for the specified dates')
#        
#        trends = slope.reshape(vals.shape[1], vals.shape[2])
#        
#        # Return xa dataarray
#        return trends#self._obj.copy().isel(time=ndate).assign(trend=(['y','x'],trends)).trend
        
#def determineTrendImages(regionstack, band, maxcloud=0.2):
#    
#    if band in ['NDVI', 'LSWI']:
#        s2_dates = np.intersect1d(regionstack.s2.time.values, regionstack.train.time.values)
#        l8_dates = np.intersect1d(regionstack.l8.time.values, regionstack.train.time.values)
#        int_dates = np.concatenate([s2_dates, l8_dates])
#        
#        quality = regionstack.train.mask.sel(time=int_dates).mean(dim=['x', 'y']).compute()
#        
#        return regionstack.train.where(quality>(1.0-maxcloud), drop=True)
#        
#    elif band in ['VV_ASC', 'VH_ASC']:
#        int_dates = np.intersect1d(regionstack.s1ASC.time.values, regionstack.train.time.values)
#        
#        return regionstack.train.sel(time=int_dates)
#        
#    elif band in ['VV_DSC', 'VH_DSC']:
#        int_dates = np.intersect1d(regionstack.s1DSC.time.values, regionstack.train.time.values)
#        
#        return regionstack.train.sel(time=int_dates)
#
#def calcAllTrends(regionstack, bands = ['NDVI','LSWI','VV_ASC','VV_DSC','VH_DSC'], maxcloud=0.1):
#    
#    for band in bands:
#        
#        valid = determineTrendImages(regionstack, band, maxcloud=maxcloud)
#        
#        first_date = np.empty(valid[band].isel(time=0).shape)
#        first_date[:] = np.nan
#        
#        c_arrays = [first_date]
#        
#        for idx, time in enumerate(valid.time[1:].values):
#            print(('Processing band {} for date {}'.format(band,time)))
#            
#            c_arrays.append(valid.eotemp.calcTempTrend(band,ndate=-idx-1))
#            
#        c = np.stack(c_arrays,axis=2)
#        
#        c_array = valid.assign(change=(['x','y','time'],c))
#        c_array = c_array.transpose('time','x','y').change
#        
#        regionstack.train[band+'_c'] = c_array
#        