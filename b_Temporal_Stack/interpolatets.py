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

def interpolate_dataset(dataset, location, bands=[], date_of_analysis='default', time_delta=16, der=False):
    """
    Method to interpolate earth observation time series in a dataset.
    
    To be passed to the xr.apply_ufunc method, the dataset is chunked
    with no chunks along the time dimension.
    
    Args:
        dataset (xr.Dataset): xarray Dataset with earth observation time series
        location (str): location to write interpolated dataset
        bands ([str]): list of string names of bands in dataset
        date_of_analysis (np.datetime64 or 'default'): date of interest
        time_delta: time between interpolated dates, in day units
        der (bool): if the derivative is beign computed
    """
    dof, times = calculate_time_periods(dataset, date_of_analysis=date_of_analysis, time_delta=time_delta)
    print('Date of analysis is {}, interpolating:\n{}'.format(dof, times))
    
    for band in bands:
        
        if dataset[band].chunks is None:
            pass # test if DataArray is not a dask array
        #Rechunk (dask) DataArray
        elif len(dataset[band].chunks[dataset[band].dims.index('time')]) > 1:
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
        if mask or der:
            int_band = interpolate_band(xa_band, times, der=der)
        else:
            int_band = interpolate_band_vect(xa_band, times)
        
        #Set output file name
        if der:
            file_name = location+band+'_der.nc'
        else:
            file_name = location+band+'.nc'
        
        # Parallelized writing
        print('Writing {} band to {}'.format(band, file_name))
        int_band.rename(band).to_netcdf(file_name)
        print('Done!')
    
    #return #xr.open_mfdataset(list(map(lambda x: location+x, os.listdir(location))))
            #Reading fails when loading data netcdf IO Error
            #chunks={'ntime':1,'x':1000,'y':1000}, parallel=True)

def interpolate_band(dataarray, intdates, der):
    """
    Interpolate time series in a data array for the specified dates. Tries to
    parallelize using dask.
    
    Args:
        dataarray (xarray.DataArray): data array (time,y,x) with time series
        intdates (np.ndarray(np.datetime64)): array of dates to interpolate
        der (bool): if the first derivative is beign calculated
    Returns:
        result (xarray.DataArray): data array with interpolated values stacked
                                    in new dimension itime (interpolated-time)
    """
    # Convert datetime objects to int
    time_base = min(pd.to_datetime(dataarray.time.values))
    
    int_orig_times = (pd.to_datetime(dataarray.time.values) - time_base).days.values
    
    int_new_times = (pd.to_datetime(intdates) - time_base).days.values
    
    # Apply ufunc-- inputs xr.DataArray and dates for interpolation
    # returns data array with interpolated values for int_dates
    result = xr.apply_ufunc(ufunc_cubic_spline, dataarray,
                            input_core_dims=[['time']],
                            output_core_dims=[['time']],
                            exclude_dims=frozenset(['time']),
                            kwargs={'axis': -1,
                                    'orig_times': int_orig_times,
                                    'new_times': int_new_times,
                                    'der': der},
                            dask='parallelized',
                            output_dtypes=[np.float32],
                            output_sizes={'time':intdates.shape[0]})
    result['time'] = ('time', intdates)
    return result

def interpolate_band_vect(dataarray, intdates):
    """"""
    # Convert datetime objects to int
    time_base = min(pd.to_datetime(dataarray.time.values))
    
    int_orig_times = (pd.to_datetime(dataarray.time.values) - time_base).days.values
    
    int_new_times = (pd.to_datetime(intdates) - time_base).days.values
    
    result = xr.apply_ufunc(ufunc_cubic_spline_vect, dataarray,
                            input_core_dims=[['time']],
                            output_core_dims=[['time']],
                            exclude_dims=frozenset(['time']),
                            kwargs={'axis': -1,
                                    'orig_times': int_orig_times,
                                    'new_times': int_new_times},
                            dask='parallelized',
                            output_dtypes=[np.float32],
                            output_sizes={'time':intdates.shape[0]})
    result['time'] = ('time', intdates)
    return result

def ufunc_cubic_spline_vect(array, orig_times, new_times, axis):
    
    # Convert datetime objects to int
    #time_base = min(pd.to_datetime(orig_times))
    
    #int_orig_times = (pd.to_datetime(orig_times) - time_base).days.values
    
    #int_new_times = (pd.to_datetime(new_times) - time_base).days.values
    
    spl = interpolate.interp1d(orig_times,
                               array,
                               kind='cubic',
                               axis=axis,
                               assume_sorted=True)
    
    return spl(new_times)

def ufunc_cubic_spline(array, orig_times, new_times, axis, der):
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
                                       new_times=new_times,
                                       der=der)
    
    return interpolated

def int_cubic_spline(y, orig_times, new_times, der):
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
        
        if der:
            deriv = spl.derivative()
            interpolated = deriv(new_times)
        else:
            interpolated = spl(new_times)
    
    except ValueError:
        #e = sys.exc_info()
        #print('{} {} {}'.format(e[0],e[1],e[2]))
        warnings.warn('CubicSpline could not be fitted for one or more pixels')
        ## When spline cannot be fitted(not enought data), return NaN
        interpolated = np.empty(new_times.shape[0])
        interpolated[:] = np.nan
        
    return interpolated

def calculate_time_periods(ds, time_delta, date_of_analysis='default'):
    """
    Defines date of analysis and calculates times periods for remote sensing image dates
    interpolation.
    
    Args:
        ds (xr.DataArray or xr.Dataset): the dataset or data array with remote sensing images
        time_delta (int): period of time between dates for interpolation
        date_of_analysis([np.datetime64]): default last date in dataset minus time_delta
    
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
    if date_of_analysis == 'default':
        date_of_analysis = [last_time - np.timedelta64(time_delta, 'D')]
        
    for date_to_process in date_of_analysis:
        # Time post-analysis
        time = date_to_process + np.timedelta64(time_delta, 'D')
        # Append times to list
        while time > first_time:
            dfi.append(time)
            time -= np.timedelta64(time_delta, 'D')
    
    # Merge unique dates in a numpy array
    dates = np.array(list(set(dfi)), dtype=np.datetime64)
    
    return date_of_analysis, np.sort(dates, axis=0)
