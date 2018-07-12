# -*- coding: utf-8 -*-
"""
A python module to interpolate time series of earth observation images in
xarray.Dataset objects.

Created on Fri Jun 29 13:59:21 2018

@author: ASALAZAR
"""

import os
import warnings

import xarray as xr
import numpy as np
import pandas as pd

from scipy import interpolate

def interpolate_dataset(dataset, location, bands=[], date_of_analysis=None):
    """
    
    Dataset must be chunked following the same: no chunks in time dimension,
    
    
    Args:
        dataset:
        
    
    """
    dof, times = calculate_time_periods(dataset, date_of_analysis=date_of_analysis)
    print('Date of analysis is {}, interpolating:\n{}'.format(dof, times))
    
    for band in bands:
        
        if len(dataset[band].chunks[dataset[band].dims.index('time')]) > 1:
            #make sure DataArray is not chunked in time dimension
            dataset[band] = dataset[band].chunk({'time':-1})
        
        # Subset and sort by increasing time
        xa_band = dataset[band].sortby('time').persist()
        
        try: #try to apply pixel quality mask
            xa_band = xa_band.where(xa_band.mask).persist()
            mask = True
        except AttributeError: #when mask is not available
            mask = False
            pass
        
        #Call ufunction to interpolate
        if mask:
            int_band = interpolate_band(xa_band, times)
        else:
            int_band = interpolate_band_vect(xa_band, times)
        
        # Parallelized writing now supported!
        print('Writing {} band to {}'.format(band, location))
        int_band.rename(band).to_netcdf(location+band+'.nc')
        print('Done!')
    
    return xr.open_mfdataset(list(map(lambda x: location+x, os.listdir(location))))
            #Reading fails when loading data netcdf IO Error
            #chunks={'ntime':1,'x':1000,'y':1000}, parallel=True)

def interpolate_band(dataarray, intdates):
    """
    Interpolate time series in a data array for the specified dates. Tries to
    parallelize using dask.
    
    Args:
        dataarray (xarray.DataArray): data array (time,y,x) with time series
        intdates (np.ndarray(np.datetime64)): array of dates to interpolate
        
    Returns:
        result (xarray.DataArray): data array with interpolated values stacked
                                    in new dimension itime (interpolated-time)
    """
    # Apply ufunc-- inputs xr.DataArray and dates for interpolation
    # returns data array with interpolated values for int_dates
    result = xr.apply_ufunc(ufunc_cubic_spline, dataarray,
                            input_core_dims=[['time']],
                            output_core_dims=[['time']],
                            exclude_dims=frozenset(['time']),
                            kwargs={'axis': -1,
                                    'orig_times': dataarray.time.values,
                                    'new_times': intdates},
                            dask='parallelized',
                            output_dtypes=[np.float32],
                            output_sizes={'time':intdates.shape[0]})
    result['time'] = ('time', intdates)
    return result

def interpolate_band_vect(dataarray, intdates):
    """"""
    result = xr.apply_ufunc(ufunc_cubic_spline_vect, dataarray,
                            input_core_dims=[['time']],
                            output_core_dims=[['time']],
                            exclude_dims=frozenset(['time']),
                            kwargs={'axis': -1,
                                    'orig_times': dataarray.time.values,
                                    'new_times': intdates},
                            dask='parallelized',
                            output_dtypes=[np.float32],
                            output_sizes={'time':intdates.shape[0]})
    result['time'] = ('time', intdates)
    return result

def ufunc_cubic_spline_vect(array, orig_times, new_times, axis):
    
    spl = interpolate.interp1d(orig_times,
                               array,
                               kind='cubic',
                               axis=axis,
                               assume_sorted=True)
    
    return spl(new_times)

def ufunc_cubic_spline(array, orig_times, new_times, axis):
    """
    Ufunc to fit a cubic spline on eo time series (y=f(x)) and interpolate
    values for the specified dates.
    
    Args:
        array (np.ndarray): 3-D array with shape (y,x,time) having the response
                                values (y), earth observation bands values
        axis (int): axis of the time dimension to apply the 
        orig_times: dates of the original earth observation images
        new_times: dates to interpolate
        
    Returns:
        interpolated (np.ndarray): 3-D array with interpolated values, stacked
                                    by interpolated dates in the third axis.
    """
    # Convert datetime objects to int
    time_base = min(pd.to_datetime(orig_times))
    
    int_orig_times = (pd.to_datetime(orig_times) - x_base).days.values
    
    int_new_times = (pd.to_datetime(new_times) - x_base).days.values
    
    # Fit cubic spline and interpolate dates
    interpolated = np.apply_along_axis(int_cubic_spline,
                                       axis,
                                       array,
                                       orig_times=int_orig_times,
                                       new_times=int_new_times)
    
    return interpolated

def int_cubic_spline(y, orig_times, new_times):
    """
    Cubic spline fitting and interpolation removing NaN values.
    
    Args:
        y (np.ndarray): 1D array with response to fit spline
        orig_times (np.ndarray): 1D array with predictor to fit spline
        new_times (np.ndarray): 1D array with predictor values to interpolate
    
    Returns:
        interpolated (np 1d-array): interpolated values
    """
    # Filter NaNs in response
    nans = np.isnan(y)
    
    # Try to fit cubic spline with filtered y values
    try:     
        spl = interpolate.CubicSpline(orig_times[~nans], y[~nans], extrapolate=False)
        
        interpolated = spl(new_times)
    
    except ValueError:
        e = sys.exc_info()
        print('{} {} {}'.format(e[0],e[1],e[2]))
        warnings.warn('CubicSpline could not be fitted for one or more pixels')
        ## When spline cannot be fitted(not enought data), return NaN
        interpolated = np.empty(new_times.shape[0])
        interpolated[:] = np.nan
        
    return interpolated

def calculate_time_periods(ds, time_delta=16, date_of_analysis=None):
    """
    Defines date of analysis and calculates times periods for remote sensing image dates
    interpolation.
    
    Args:
        ds (xr.DataArray or xr.Dataset): the dataset or data array with remote sensing images
        time_delta (int): period of time between dates for interpolation
        date_of_analysis(np.datetime64): default last date in dataset minus time_delta
    
    Returns:
        date_of_analysis (np.datetime64): the date to perform the analysis
        dates (np.ndarray, dtype:np.datetime64): 1D array with dates for interpolatation
    """
    # Extract time from Dataset/DataArray
    times = ds.sortby('time').time.values
    last_time = times[-1]
    first_time = times[0]
    
    # List to append dates for interpolation
    dfi = []
    
    # Define default date of analysis
    if date_of_analysis is None:
        date_of_analysis = last_time - np.timedelta64(time_delta, 'D')
    # Time post-analysis
    time = date_of_analysis + np.timedelta64(time_delta, 'D')
    # Append times to list
    while time > first_time:
        dfi.append(time)
        time -= np.timedelta64(time_delta, 'D')
    # Merge in a numpy array
    dates = np.array(dfi, dtype=np.datetime64)
    
    return date_of_analysis, dates
