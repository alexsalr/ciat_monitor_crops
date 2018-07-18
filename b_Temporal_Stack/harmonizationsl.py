# -*- coding: utf-8 -*-
"""
A python module to correct Landsat-8 earth observation datasets to match
Sentinel-2 images.

Created on Fri Jun 29 08:55:57 2018

@author: ASALAZAR
"""
import sys
import xarray as xr
import numpy as np

from skimage.feature import register_translation
from skimage.transform import SimilarityTransform
from skimage.transform import warp

def apply_bandpass_correction(lc08dataset):
    """
    Apply band pass adjustment of Landsat-8 images acording to linear regression
    coefficients calculated by NASAs Harmonized Landsat Sentinel-2 project. 
    See https://hls.gsfc.nasa.gov/algorithms/bandpass-adjustment/
    
    Please note that modifies the input dataset.
    
    Args:
        lc08dataset (xarray.Dataset): Landsat-8 dataset to correct
        offset ([float,float]): image offset (y,x) between S2 and LC08 images
    
    Returns:
        lc08dataset (xarray.Dataset): a reference to the corrected dataset
    """
    
    # Inversed linear regression parameters [B1, B0] for each band
    lr_harm = {'blue':[1.0/0.977, -0.00411/-0.977],
               'green':[1.0/1.005, -0.00093/-1.005],
               'red':[1.0/0.982, 0.00094/-0.982],
               'nir':[1.0/1.001, -0.00029/-1.001],
               'swir1':[1.0/1.001,-0.00015/-1.001],
               'swir2':[1.0/0.996,-0.00097/-0.996]}
    
    # Apply bandpass adjustment
    for band, coeff in lr_harm.items():
        try:
            #print('Correcting {} band'.format(band))
            lc08dataset[band] = (coeff[0]*lc08dataset[band]+coeff[1]).astype(np.uint16)
        except KeyError: #when band is not in dataset
            pass
    
    return lc08dataset

def apply_geo_correction(dataset, offset = [3.37, 1.18]):
    """
    Warps the dataset variables applying a transformation to correct the
    specified 2D-offset (y,x)
    
    Please note that modifies the input dataset and that each variable in the
    dataset is loaded to memory. Passes the y axis offset directy as the translation
    to apply as they are applied in an inverted y-axis raster within the xr.apply_ufunc
    method.
    
    Args:
        dataset (xarray.Dataset): dataset to correct
        offset ([float,float]): sub-pixel offset in y and x dims(in pixels)
    
    Returns:
        dataset (xarray.Dataset): a reference to the original (now corrected)
                                    dataset.
    """
    # Correct pixel offset translation parameters (tx, ty)
    tform = SimilarityTransform(translation=(-offset[1],offset[0])) #inverted y-axis rasters
    # Read variables names in dataset, exclusing dimensions (y,x,time)
    vars_in_ds = [band for band in list(dataset.variables.keys()) if band not in ['time','y','x','mask']]
    
    # Iterate variables
    for band in vars_in_ds:
        try:
            #print('Warping {} band to y={},x={}'.format(band, -offset[0], -offset[1])) #Cannot be chunked in x,y as data is lost in chunk limits
            dataset[band] = xr.apply_ufunc(__warp_clip__, dataset[band].load(),#dask='parallelized',
                   input_core_dims=[['time']], output_core_dims=[['time']],#output_dtypes=[dataset[band].dtype],
                   kwargs={'inverse_map':tform,'dtype_':dataset[band].dtype,'order':0,'preserve_range':True})
        except:
            e = sys.exc_info()
            print('Warping {} band failed with error {}, {}, {}'.format(band, e[0], e[1], e[2]))
    
    return dataset
    
def __warp_clip__(da,dtype_,**kwargs):
    """
    Apply skimage.transform warp to numpy ndarray.
    Third dimension (band) is in this case the temporal dimension.
    
    Args:
        da (np.ndarray): 2D or 3D array with independent
                    rasters stacked in the 3rd dimension
        dtype_ (np.dtype): Type to cast the resulting array
        **kwargs: optional arguments to pass to skimage.transform.warp
        
    Returns:
        da (np.ndarray): Warped image, clipped and casted to specified numpy type
    """
    return np.clip(warp(da, **kwargs), a_min=0, a_max=65535).astype(dtype_)

def calculate_offset(sentinel2_image, landsat8_image):
    """
    Calculates offset . Recieves single band images as xarray.DataArray objects
    with only 'y' and 'x' dimensions and coordinates.
    
    Args:
        sentinel2_image (xarray.DataArray): reference image (y,x)
        landsat8_image (xarray.DataArray): secondary image with offset (y,x)
    
    Returns:
        shift ([float,float]): list with (y,x) subpixel offset
    """
    
    ## TO-DO: test if images have the same dimension
    
    image = np.nan_to_num(sentinel2_image)
    offset_image = np.nan_to_num(landsat8_image)
    
    # 1-pixel precison (sub pixel correction requires order-1+ warp with interpolation
    shift, error, diffphase = register_translation(image, offset_image)
    
    print("Detected pixel offset (y, x): {}".format(shift))
    
    # check correction
    tform = SimilarityTransform(translation=(-shift[1],shift[0])) #inverted y-axis rasters
    warped = xr.apply_ufunc(__warp_clip__,landsat8_image.load(),
                   kwargs={'inverse_map':tform,
                           'dtype_':landsat8_image.dtype,
                           'order':0,
                           'preserve_range':True})
    
    corr_shift, corr_error, corr_diffphase = register_translation(image,
                                                                  warped,
                                                                  100)
    
    print("Sub-pixel offset after correction in reference image (y, x): {}".format(corr_shift))
    
    return shift
    