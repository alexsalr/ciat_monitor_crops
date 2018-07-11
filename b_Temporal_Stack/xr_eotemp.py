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

import rasterio
import os

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
    
    
    def to_tiff(self, bands, filename, time, addbands=None, crs=None, transform=None):
        """
        Write tiff raster file with the selected variables as bands in
        xarray dataset.
        
        Args:
            bands ([str]): name of dataset variables to write as bands
            filename (str): file name to use to write the tif file
            time (np.datetime64 or int): used to select time in eots
            addbands ({str:xr.DataArray}): optional. If passed, all
                              variables in dataset are added as bands
            crs (str): e.g. +init=epsg:32618
            transform (tuple): gdal transform (c, a, b, f, d, e)
                e.g.(488900.0, 10.0, -0.0, 448460.0, -0.0, -10.0)
                |c| x-coord of the upper-left corner of the upper-left pixel
                |a| width of a pixel
                |b| row rotation (typically zero)
                |f| y-coord of the of the upper-left corner of the upper-left pixel
                |d| column rotation (typically zero)
                |e| height of a pixel (typically negative)
        """
        
        if type(time) is np.datetime64:
            temp_subset = self._obj.sel(time=time).drop(['time','mask'])[bands]
        elif type(time) is int:
            temp_subset = self._obj.isel(time=time).drop(['time','mask'])[bands]
        
        # Put bands in list
        bandlist = [temp_subset[band] for band in bands]
        bandnames = bands
        
        # Add additional bands to list
        if addbands is not None:
            bandlist = bandlist + list(addbands.values())
            bandnames = bandnames + list(addbands.keys())
        
        # Show the array that is being written
        for band in bandlist:
            print(band)
        
        # Load data to memory
        bandlist = [band.load() for band in bandlist]
        
        # Merge datasets
        xa = xr.merge([band.rename(bandnames[idx]) for idx, band in enumerate(bandlist)])
        
        # Reshape dataset to DataArray
        xa = xa.to_array(dim='bands')
        
        # Assign crs
        if crs is not None:
            xa.attrs['crs'] = crs
        else:
            try:
                xa.attrs['crs'] = temp_subset.attrs['crs']
            except KeyError:
                raise KeyError('crs not present in Dataset, please provide as argument')
        
        # Assign transform
        if transform is not None:
            xa.attrs['transform'] = transform
        else:
            try:
                xa.attrs['transform'] = temp_subset.attrs['transform']
            except KeyError:
                raise KeyError('transform not present in Dataset, please provide as argument')
        
        # Write to tif using rasterio
        xarray_to_rasterio(xa, filename)

def xarray_to_rasterio(xa, output_filename):
    """
    
    Converts the given xarray.DataArray object to a raster output file
    using rasterio.
    
    Args:
        xa (xr.DataArray) DataArray to convert
        output_filename (str): filename of output GeoTIFF file
    
    Notes:
        Converts the given xarray.DataArray to a GeoTIFF output file using rasterio.
        This function only supports 2D or 3D DataArrays, and GeoTIFF output.
        The input DataArray must have attributes (stored as xa.attrs) specifying
        geographic metadata, or the output will have _no_ geographic information.
        If the DataArray uses dask as the storage backend then this function will
        force a load of the raw data.
        
        Modified coordinates of np.array using the affine attribute, flipping coordinates
        when coordinates are decreasing. Requires transposed array (y, x).
        
    References:
        https://github.com/robintw/XArrayAndRasterio
    
    """
    
    # Forcibly compute the data, to ensure that all of the metadata is
    # the same as the actual data (ie. dtypes are the same etc)
    xa = xa.load()
    bands = xa.bands
    
    if len(xa.shape) == 2:
        count = 1
        height = xa.shape[0]
        width = xa.shape[1]
        band_indicies = 1
    else:
        count = xa.shape[0]
        height = xa.shape[1]
        width = xa.shape[2]
        band_indicies = np.arange(count) + 1

    processed_attrs = {}
    
    # Get values to write
    bands_to_write = xa.values
    
    try:
        val = xa.attrs['transform']
        processed_attrs['affine'] = rasterio.Affine.from_gdal(*val)
        # Test if axes need to be inversed
        if val[1] < 0:
            bands_to_write = np.flip(bands_to_write, axis=2)
        if val[5] < 0:
            bands_to_write = np.flip(bands_to_write, axis=1)
    except KeyError:
        raise KeyError('transform attribute is not present in dataset')
    try:
        val = xa.attrs['crs']
        processed_attrs['crs'] = rasterio.crs.CRS.from_string(val)
    except KeyError:
        raise KeyError('crs attribute is not present in dataset')
    
    with rasterio.open(output_filename, 'w',
                       driver='GTiff',
                       height=height, width=width,
                       dtype=str(xa.dtype), count=count,
                       **processed_attrs) as dst:
        dst.write(bands_to_write, band_indicies)
        


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