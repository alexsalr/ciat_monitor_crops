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

from scipy import interpolate

def interpolate_dataset(dataset, location, bands=[], date_of_analysis=None):
    
    dds = dataset.sortby('time').chunk({'time':-1,'x':100,'y':100})#.persist()
    
    #output = []
    
    dof, times = calculate_time_periods(dds, date_of_analysis=date_of_analysis)
    print('Date of analysis is {}, interpolating:\n{}'.format(dof, times))
    
    for band in bands:    
        int_band = interpolate_band(dds[band].persist(), times)
        
        try: #try to apply pixel quality mask
            int_band = int_band.where(int_band.mask)
        except KeyError: #when mask is not available
            pass
        
        # Parallelized writting not supported, need to load in memory
        int_band.compute().to_netcdf(location+band+'.nc')
        print('Dataset for band {} was written in {}'.format(band, location))
    
    return xr.open_mfdataset(list(map(lambda x: location+x, os.listdir(location))),
                            chunks={'ntime':1,'x':1000,'y':1000}, parallel=True)

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
                            output_core_dims=[['itime']],
                            kwargs={'axis': -1,
                                    'orig_times': dataarray.time.values,
                                    'new_times': intdates},
                            dask='parallelized',
                            output_dtypes=[np.float32],
                            output_sizes={'itime':intdates.shape[0]})
    result['itime'] = ('itime', intdates)
    return result

def ufunc_cubic_spline(array, axis=-1, orig_times, new_times):
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
    # Fit cubic spline and interpolate dates
    interpolated = np.apply_along_axis(int_cubic_spline,
                                       axis,
                                       array,
                                       orig_times=orig_times,
                                       new_times=new_times)
    
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
        spl = interpolate.CubicSpline(orig_times.astype('d')[~nans], y[~nans])
        
        interpolated = spl(new_times.astype('d'))
    
    except ValueError:
        warnings.warn('CubicSpline could not be fitted for one or more pixels')
        ## When spline cannot be fitted(not enought data), return NaN
        interpolated = np.empty(new_times.shape[0])
        interpolated[:] = np.nan
        
    return interpolated

# =============================================================================
# def savitzky_golay(y, window_size, order, deriv=0, rate=1):
#     """http://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html"""
#     import numpy as np
#     from math import factorial
#     
#     try:
#         window_size = np.abs(np.int(window_size))
#         order = np.abs(np.int(order))
#     except ValueError:
#         raise ValueError("window_size and order have to be of type int")
#     if window_size % 2 != 1 or window_size < 1:
#         raise TypeError("window_size size must be a positive odd number")
#     if window_size < order + 2:
#         raise TypeError("window_size is too small for the polynomials order")
#     order_range = range(order+1)
#     half_window = (window_size -1) // 2
#     # precompute coefficients
#     b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
#     m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
#     # pad the signal at the extremes with values taken from the signal itself
#     firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
#     lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
#     y = np.concatenate((firstvals, y, lastvals))
#     return np.convolve( m[::-1], y, mode='valid')
# =============================================================================

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

