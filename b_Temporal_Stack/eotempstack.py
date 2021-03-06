# -*- coding: utf-8 -*-
"""
Module for constructing a multi-temporal earth observation dataset from pre-processed Sentinel-1, Sentinel-2, Landsat-7/8 raster files

@author: asalazar
"""

import os
import re
import sys
import rasterio
import dask
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from collections import defaultdict
from datetime import datetime as dt
from calendar import monthrange
from rasterio.merge import merge
#from rasterio.rio import stack

class eoTempStack(object):
    """
    Class to construct and analyze a multi temporal stack of earth observation images. Supports Sentinel-1, Sentinel-2,
    Landsat 7 and Landsat 8 pre-processed products. All images in the same object have the same coordinate system, spatial
    extent and resolution. Operates on a dictionary of band names for each of the supported products. Required pre-processing:
        - For Sentinel-1... 
        - For Sentinel-2... L2A
        - For Landsat(7/8)...
    """
    def __init__(self, sourcedir, outdir, prodtype):
        """
        Class constructor.
        
        @sourcedir (str): location of products as described in @prodtype
        @outdir (str): directory to store resulting data products (raster stack)
        @prodtype (str): accepts 'S1', 'S2', 'LE07', 'L8' where
                            'S1': set of .img files with the same extent and polarization, product from SNAP collocation Op
                            'S2': set of directories resulting from sen2cor processing Sentinel-2 level 2A product
        """
        ## Set parameters as object attributes
        self.source_directory = check_dir(sourcedir)
        self.prod_type = prodtype
        self.out_directory = outdir
        self.bands_loc = {} #Declare empty dictionary. Set in setBandsLoc
        self.bands_temporal_range = {} #Declare empty dictionary. Set in setTempData
        
        # Set the minimum required variables to construct the object. Implementation varies by product subclass.
        self.setBandsLoc()
        self.setTempData()
        
        print(('{} object initialized from {}'.format(self.prod_type, self.source_directory)))
        
    def getBandsLoc(self, band=None):
        if band is not None:
            return self.bands_loc[band]
        else:
            return self.bands_loc
    
    def getTempData(self, band=None):
        if band is not None:
            return self.bands_temporal_range[band]
        else:
            return self.bands_temporal_range
    
    def getSourceDir(self):
        return self.source_directory
    
    def getOutDir(self):
        return self.out_directory
    
    def getProdType(self):
        return self.prod_type
    
    def getBand(self, band, tempid=None, date=None):
        """Operates on all bands in object unless an tempid or date is provided.
        Returns a band or a list of bands"""
        bandsloc = self.getBandsLoc(band)
        if tempid is not None:
            with rasterio.open(bandsloc[tempid]) as src:
                readband = src.read(1)
        elif date is not None:
            with rasterio.open(bandsloc[self.getBandIndex(date)]) as src:
                readband = src.read(1)
        else:
            readband = []
            for loc in bandsloc:
                with rasterio.open(loc) as src:
                    readband.append(src.read(1))
        return readband
    
    def getBandIndex(self, date):
        """Return the tempid in of the date in the  
        define if date is a string or a date type
        """
        for id, layer in enumerate(next(iter(self.getTempData().values()))):
            if date == layer:
                return id
    
    def getBandXarray(self, band, time_range):
        """Returns an xarray.Array object from all dates in a given band in a given time_range tuple (start_date, end_date)"""
        start_date = time_range[0]
        end_date = time_range[1]
        # Get the times within the time_range
        time = xr.Variable('time', pd.DatetimeIndex([pd.Timestamp(f) for f in self.getTempData(band) if start_date <= f <= end_date]))
        # Use numpy to filter band locations and return in band
        band_locations = np.array(self.getBandsLoc(band))[np.array(list([start_date <= f <= end_date for f in self.getTempData(band)]))].tolist()
        # Read all bands using open_rasterio xarray method
        arlist = [xr.open_rasterio(f, chunks={'x':5000, 'y':5000}) for f in band_locations]
        # Concatenate all dates
        da = xr.concat(arlist, dim=time)
        return da
    
    def mosaicBand(self, band):
        """
        Checks if dates are stored more than once in the temporal_data for one
        band. If it does, read and mosaic all duplicated dates. Then update
        band location and dates data.
        
        Restrictions: if any date is duplicated, only mosaiced data will remain
        in object dictionaries.
        
        Args:
            band (str): name of band (keys in object loc/temp_range dicts)
        """
        ls = list_indices(self.getTempData(band))
        
        # By default updating dictionaries is disabled
        update_locs = False
        
        for date, ilist in ls.items():
            
            # If more than one image is on the same date
            if len(ilist)>1:
                # If any date is duplicated, enable updating dictionaries
                update_locs = True
                
                mosaicloc = self.getSourceDir()+self.getProdType()+'_'+date.strftime('%Y%m%d')+'_'+band+'_mosaic.tif'
                
                tiles_loc = list(map(lambda x: self.getBandsLoc(band)[x], ilist))
                
                if os.path.isfile(mosaicloc):
                    pass
                else:
                    try:
                        print('Mosaicing {} band for {}'.format(band, date.strftime('%Y-%m-%d')))
                        merge_tiles(tiles_loc, mosaicloc)
                    except:
                        print('{} band mosaic for {} failed, will be excluded from {} eotempstack'.format(band, date, self.prod_type))
        
        if update_locs:
            self.__updateBandsMosaic__(band)
        
    def __updateBandsMosaic__(self, band):
        """
        Update mosaiced band information
        """
        
        fl = os.listdir(self.getSourceDir())
        ffilter = r'^'+self.getProdType()+'.*'+band+'.*mosaic.tif$'
        
        file_names = list(filter(re.compile(ffilter).search, fl))
        
        self.bands_loc[band] = list(map(lambda x: self.getSourceDir()+x, file_names))
        self.bands_temporal_range[band] = list(map(lambda x: dt.strptime(x[3:11],
                                 '%Y%m%d').date(), file_names))
        
        
##Helper methods
def merge_tiles(tilesloc, mosaicloc, **args):
    """
    Mosaic images using rasterio merge
    
    Args:
        tilesloc([str]): list of string paths of tile files
        mosaicloc(str): path to store mosaic image file
        **args: arguments to pass to rasterio.merge
    """
    
    # Retrieve the meta data from the first in list
    try:
        kwargs = rasterio.open(tilesloc[0]).meta
        # Merge images in list
        arr, out_trans = merge(list(map(lambda x: rasterio.open(x), tilesloc)), **args)
    except RasterioIOError:
        warning.warn('One or more of the images could not be read.')
        raise
    # Update meta data
    kwargs.update(width=arr.shape[2],height=arr.shape[1],transform=out_trans.to_gdal(),affine=out_trans)
    # Write mosaic file
    with rasterio.open(mosaicloc, 'w', **kwargs) as dst:
        for band in range(int(kwargs['count'])):
            dst.write_band(band+1, arr[band,:,:])

def list_indices(seq):
    """
    Make a dictionary of values of a sequence with a list of their indices in the original sequence
    
    Args:
        seq (list): a list of elements
    
    Returns:
        counter (defaultdic): dictionary with elements as keys and indices
                                as values
    """
    counter = defaultdict(list)
    for i,item in enumerate(seq):
        counter[item].append(i)
    return counter
    
class S1TempStack(eoTempStack):
    bands_of_interest = ['VV', 'VH']
    
    def __init__(self, sourcedir, outdir, orbit, prod_type = 'S1'):
        # Call eo_tempstack initialization method
        self.orbit = orbit
        super(S1TempStack, self).__init__(sourcedir, outdir, prod_type)
    
    def setBandsLoc(self):
        ## Get directory names of pre-processed S1 products to read
        proddirs = list(filter(re.compile(r'^S1_SPF_'+self.orbit+'.*data$').search, os.listdir(self.source_directory)))
        
        ## Declare dictionary to store file location by polarizations
        prodlist = {}
        for prod_dir in proddirs:
            # List and filter img files for Sigma0 bands, store full path to file
            image_files = list(filter(re.compile(r'Sigma0.*img$').search, os.listdir(self.source_directory+prod_dir)))
            prodlist[prod_dir.split('_')[3]] = list([self.getSourceDir()+prod_dir+'/'+x for x in image_files])
        # Store the dictionary of files location as instance variable
        self.bands_loc = prodlist
        
    def setTempData(self, key=None, tempdata=None):
        ## Declare dictionary to store dates by polarization
        temp_range = {}
        for key, value in self.getBandsLoc().items():
            temp_range[key] = list([dt.strptime(x.split('/')[-1].split('_')[3],'%d%b%Y').date() for x in value])
            # Store dictionary as instance variable
        self.bands_temporal_range = temp_range
    
    #DEPRECATED
    #def getXarray(self):
    #    xarrays = []
    #    for band in self.bands_of_interest:
    #        x = self.getBandXarray(band).isel(band=0)
    #        xarrays.append(x)
    #    ## Modify this to match standard_band_dict
    #    band = xr.Variable('band', pd.Index(self.bands_of_interest))
    #    da = xr.concat(xarrays, band)
    #    return da
    
    def createXDataset(self):
        try:
            # Get unique date values in dict
            date_values = list(set([date for bandlist in list(self.getTempData().values()) for date in bandlist]))
            
            # Calculate ranges of the available data. We consider monthly ranges, from first to last day of each month
            time_ranges = list([(dt.strptime('01'+x.strftime('%m%Y'), '%d%m%Y').date(),
                            dt.strptime(str(monthrange(x.year, x.month)[1])+x.strftime('%m%Y'), '%d%m%Y').date()) for x in date_values])
            
            # Iterate over unique date ranges (i.e. per month)            
            for time_range in list(set(time_ranges)):
                # Extract the month in format %Y%m
                month = time_range[0].strftime('%Y%m')
                # Put arrays for each band in list
                xarrays = []
                # Retrieve all bands
                for band in self.bands_of_interest:
                    try:
                        x = self.getBandXarray(band, time_range).isel(band=0)
                        x.name = band
                        xarrays.append(x)
                    except:
                        print(('{} band for S1 {} {} could not be processed'.format(band, self.orbit, month)))
                # Merge all bands in a single array
                xa = xr.merge(xarrays)
                # Drops unused bands
                xa = xa.drop('band')
                # Filename
                file_name = self.out_directory+self.prod_type+'_'+self.orbit+'_'+month+'.nc'
                # Write array as netcdf
                if not os.path.isfile(file_name):
                    xa.to_netcdf(file_name)
                else:
                    warnings.warn('The dataset was not written. The file {} already exists.'.format(file_name))
                
        except:
            e = sys.exc_info()
            print('Stacking {} failed: {} {} {}'.format(self.prod_type, e[0], e[1], e[2]))

class S1TextureTempStack(S1TempStack):
    bands_of_interest = ['VV_ASM',
                         'VV_Contrast',
                         'VV_Dissimilarity',
                         'VV_Energy',
                         'VV_Entropy',
                         'VV_GLCMCorrelation',
                         'VV_GLCMMean',
                         'VV_GLCMVariance',
                         'VV_Homogeneity',
                         'VH_ASM',
                         'VH_Contrast',
                         'VH_Dissimilarity',
                         'VH_Energy',
                         'VH_Entropy',
                         'VH_GLCMCorrelation',
                         'VH_GLCMMean',
                         'VH_GLCMVariance',
                         'VH_Homogeneity',]
    
    # 'GLCM_S1_DESCENDING_20180410_'.
    def __init__(self, sourcedir, outdir, orbit):
        # Call eo_tempstack initialization method
        self.orbit = orbit
        super(S1TextureTempStack, self).__init__(sourcedir, outdir, orbit, 'GLCM_S1')
        
    def setBandsLoc(self):
    
        ## Get directory names of pre-processed S1 products to read
        prodlist = list(filter(re.compile(r'^GLCM_S1_'+self.orbit+'.*data$').
                               search, os.listdir(self.source_directory)))
        
        ## Get names of files to stack in raster
        prodloclist = {}
        for band in self.bands_of_interest:
            prodloclist[band] = list([self.getSourceDir()+x+'/Sigma0_'+band+
                       '_S.img' for x in prodlist])
        
        self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        temp_range = {}
        
        for key, value in self.getBandsLoc().items():
            temp_range[key] = list([dt.strptime(x.split('/')[-2].split('_')[3],
                      '%Y%m%d').date() for x in value])
        
        self.bands_temporal_range = temp_range
        
class opticalTempStack(eoTempStack):
    ## Init method to include indices calculation
    def __init__(self, sourcedir, outdir, prodtype):
        # Call eo_tempstack initialization method
        super(opticalTempStack, self).__init__(sourcedir, outdir, prodtype)
        
        # Check if dates are duplicated, and mosaic
        for band in self.bands_of_interest: self.mosaicBand(self.standard_band_dict[band])
        
        # Calculate indices at initialization TODO determine if is appropiate now
        for band in self.calculated_bands: self.calcIndex(band)
        
    def calcQualityPixels(self):
        """Calculates cloud cover and shade by date from original Sentinel-2 files"""
        qualitypixels = []
        for date in self.getTempData('qa_cloud'):
            cloud_mask = self.getBand('qa_cloud', date=date) > self.cloud_quality_limit['qa_cloud']
            shadow_mask = self.getBand('qa_class', date=date).isin(self.cloud_quality_limit['qa_class'])
            rmask = np.logical_or(cloud_mask, shadow_mask)
            qualitypixels.append(1 - float(np.sum(rmask)) / float(rmask.size))
        return qualitypixels
        
    def calcIndex(self, index):
        """index(str): accepts NDVI and LSWI"""
        #calc_band = []
        locations = []
        
        ## Iterate to read bands and calculate (first check if file exists)
        for id, date in enumerate(self.getTempData('nir')):
            # define location of index raster file
            bandLocation = self.source_directory+self.prod_type+'_'+index+'_'+self.bands_temporal_range['nir'][id].strftime('%Y%m%d')+'_.tif'
            # check if raster file already exists
            if os.path.isfile(bandLocation):
                locations.append(bandLocation)
            # if not, calculate index and put in raster file
            else:
                nir = self.getBand('nir', date=date)
                # Allow division by zero
                np.seterr(divide='ignore', invalid='ignore')
                # Calculate index
                if index == 'NDVI': #(nir-red)/(nir+red)
                    red = self.getBand('red', date=date)
                    calc_band = (nir.astype(float) - red.astype(float)) / (nir + red)
                if index == 'LSWI': #(nir-swir1)/(nir+swir1)
                    swir1 = self.getBand('swir1', date=date)
                    calc_band = (nir.astype(float) - swir1.astype(float)) / (nir + swir1)
                locations.append(bandLocation)
                with rasterio.open(self.getBandsLoc('nir')[id]) as src:
                    kwargs = src.meta
                kwargs.update(dtype=rasterio.float32)
                with rasterio.open(bandLocation, 'w', **kwargs) as dst:
                    dst.write_band(1, calc_band.astype(rasterio.float32))
        
        # Update object variables
        self.setBandsLoc(index, locations)
        self.setTempData(index, self.getTempData('nir'))
    
    def getMaskedBand(self, key, tempid=None, date=None):
        if tempid is not None or date is not None:
            cloud_mask = self.getBand('qa_cloud', tempid=tempid, date=date) > self.cloud_quality_limit['qa_cloud']
            shadow_mask = self.getBand('qa_class', tempid=tempid, date=date).isin(self.cloud_quality_limit['qa_class'])
            rmask = np.logical_or(cloud_mask, shadow_mask)
            return np.ma.masked_array(self.getBand(key, tempid=tempid, date=date), mask=rmask)
        else:
            raise NotImplementedError('Not Implemented. Please specify date or tempid.')
    
    def getMaskArray(self, tempid=None, date=None):
        if tempid is not None or date is not None:
            cloud_mask = self.getBand('qa_cloud', tempid=tempid, date=date) > self.cloud_quality_limit['qa_cloud']
            shadow_mask = self.getBand('qa_class', tempid=tempid, date=date).isin(self.cloud_quality_limit['qa_class'])
            rmask = np.logical_or(cloud_mask, shadow_mask)
            return rmask
        else:
            raise NotImplementedError('Not Implemented. Please specify date or tempid.')
    
    def getMaskXarray(self, time_range):
        Q1 = self.getBandXarray('qa_cloud', time_range) > self.cloud_quality_limit['qa_cloud']
        Q2 = self.getBandXarray('qa_class', time_range).isin(self.cloud_quality_limit['qa_class'])
        xrmask = xr.ufuncs.logical_not(xr.ufuncs.logical_or(Q1, Q2))
        return xrmask.isel(band=0).drop('band')
    
    #DEPRECATED
    #def getXarray(self):
    #    xarrays = []
    #    bandnames = []
    #    for band in self.getBandsLoc().keys():
    #        bandnames.append(band)
    #        if self.prod_type == 'S2' and band in ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']:
    #            x = self.getBandXarray(band).isel(band=0).drop('wavelength')
    #        else:
    #            x = self.getBandXarray(band).isel(band=0)
    #        xarrays.append(x)
    #    band = xr.Variable('band', pd.Index(bandnames))
    #    xa = xr.concat(xarrays, band)
    #    return xa
    
    def createXDataset(self):
        """
        Reads images from all bands in objects dictionary into xarray.DataArray objects using rasterio and concatenates all available dates
        in fixed periods (implemented, monthly) to store as netcdf files of earth observation time series.
        """
        #try:
        # Calculate ranges of the available data. We consider monthly ranges, from first to last day of each month
        time_ranges = list([(dt.strptime('01'+x.strftime('%m%Y'), '%d%m%Y').date(),
                        dt.strptime(str(monthrange(x.year, x.month)[1])+x.strftime('%m%Y'), '%d%m%Y').date()) for x in next(iter(self.getTempData().values()))])
        
        # Iterate over unique date ranges (i.e. per month)            
        for time_range in list(set(time_ranges)):
            # Extract the month in format %Y%m
            month = time_range[0].strftime('%Y%m')
            # Put arrays for each band in list
            xarrays = []
            # Retrieve all bands
            for band in list(self.getBandsLoc().keys()):
                if self.prod_type == 'S2' and band in ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']:
                    x = self.getBandXarray(band, time_range).isel(band=0).astype(np.uint16, copy=False)
                    try: # Try to remove unused dimension for Sentinel-2 .img files (fails when mosaicing)
                        x = x.drop('wavelength')
                    except:
                        pass
                elif self.prod_type in ['LE07', 'LC8'] and band in ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']:
                    x = self.getBandXarray(band, time_range).isel(band=0).astype(np.int16, copy=False)
                else:
                    x = self.getBandXarray(band, time_range).isel(band=0)
                x.name = band
                xarrays.append(x)
            # Merge all bands in a single array
            xa = xr.merge(xarrays)
            # Include quality pixel band
            xa.coords['mask'] = (('time', 'y', 'x'), self.getMaskXarray(time_range))
            # Drops unused bands when present
            try:
                xa = xa.drop('band')
            except:
                pass
            try:
                xa = xa.drop(['qa_class', 'qa_cloud'])
            except:
                pass
            file_name = self.out_directory+self.prod_type+'_'+month+'.nc'
            # Write array as netcdf
            if not os.path.isfile(file_name):
                xa.to_netcdf(file_name)
            else:
                warnings.warn('The dataset was not written. The file {} already exists.'.format(file_name))
        #except:
        #    e = sys.exc_info()
        #    print(('Stacking {} failed: {} {} {}'.format(self.prod_type, e[0], e[1], e[2])))

class S2TempStack(opticalTempStack):
    cloud_quality_limit = {'qa_cloud':9,'qa_class':[3]}     #quality_cloud_confidence>9, quality_scene_classification == 3
    bands_of_interest = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'quality_cloud_confidence','quality_scene_classification']
    
    standard_band_dict = {'B2':'blue', 'B3':'green', 'B4':'red', 'B8':'nir', 'B11':'swir1', 'B12':'swir2', 
                          'quality_cloud_confidence':'qa_cloud','quality_scene_classification':'qa_class'}
    calculated_bands = ['NDVI', 'LSWI'] 
    
    def __init__(self, sourcedir, outdir):
        super(S2TempStack, self).__init__(sourcedir, outdir, 'S2')
        
    def setBandsLoc(self, key=None, bandloc=None):
        if key is not None and bandloc is not None:
            self.bands_loc[key] = bandloc
        else:
            ## Get names of files to stack in raster
            prodlist = list(filter(re.compile(r'^S2.*[MSIL2A].*data$').search, os.listdir(self.getSourceDir())))
            
            prodloclist = {}
            for band in self.bands_of_interest:
                prodloclist[self.standard_band_dict[band]] = list([self.getSourceDir()+x+'/'+band+'.img' for x in prodlist])
                
            self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        if key is not None and tempdata is not None:
            self.bands_temporal_range[key] = tempdata
        else:
            temp_range = {}
            for key, value in self.getBandsLoc().items():
                try:
                    temp_range[key] = list([dt.strptime(x.split('/')[-2][11:19], 
                                                                                '%Y%m%d').date() for x in value])
                except:
                    try:
                        temp_range[key] = list([dt.strptime(x.split('/')[-2][47:55], 
                                                                                '%Y%m%d').date() for x in value])
                    except ValueError:
                        raise Exception('The S2 L2A product does not follow naming conventions (dates %Y%m%d at 11:19 or 47:55).')
            self.bands_temporal_range = temp_range
        
    
class L8TempStack(opticalTempStack):
    cloud_quality_limit = {'qa_cloud':323, 'qa_class':[16, 80, 144, 208]}     #pixel_qa>323, sr_aerosol==16/80/144/208
    bands_of_interest = ['ref_srband2', 'ref_srband3', 'ref_srband4', 'ref_srband5', 'ref_srband6', 'ref_srband7',
                         'ref_pixelqa', 'ref_sraerosol']
    standard_band_dict = {'ref_srband2':'blue', 'ref_srband3':'green', 'ref_srband4':'red', 'ref_srband5':'nir',
                          'ref_srband6':'swir1', 'ref_srband7':'swir2', 'ref_pixelqa':'qa_cloud',
                          'ref_sraerosol':'qa_class'}
    calculated_bands = ['NDVI', 'LSWI']
    
    def __init__(self, sourcedir, outdir):
        super(L8TempStack, self).__init__(sourcedir, outdir, 'LC08')
        
    def setBandsLoc(self, key=None, bandloc=None):
        if key is not None and bandloc is not None:
            self.bands_loc[key] = bandloc
        else:
            ## Get names of files to stack in raster
            prodlist = list(filter(re.compile(r'^LC08.*[0-9]$').search, os.listdir(self.getSourceDir())))
            prodloclist = {}
            for band in self.bands_of_interest:
                prodloclist[self.standard_band_dict[band]] = list([self.getSourceDir()+x+'/'+band+'.tif' for x in prodlist])
            self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        if key is not None and tempdata is not None:
            self.bands_temporal_range[key] = tempdata
        else:
            temp_range = {}
            for key, value in self.getBandsLoc().items():
                temp_range[key] = list([dt.strptime(x.split('/')[-2][10:18], 
                                                                                '%Y%m%d').date() for x in value])
            self.bands_temporal_range = temp_range


class L7TempStack(opticalTempStack):
    cloud_quality_limit = {'qa_cloud':67, 'qa_class':[4, 12, 20, 36, 52]} #pixel_qa>67 cloud_qa ==4,12,20,36,52
    bands_of_interest = ['ref_srband1', 'ref_srband2', 'ref_srband3', 'ref_srband4', 'ref_srband5', 'ref_srband7',
                         'ref_pixelqa', 'ref_cloudqa']
    standard_band_dict = {'ref_srband1':'blue', 'ref_srband2':'green', 'ref_srband3':'red', 'ref_srband4':'nir',
                          'ref_srband5':'swir1', 'ref_srband7':'swir2', 'ref_pixelqa':'qa_cloud',
                          'ref_cloudqa':'qa_class'}
    calculated_bands = ['NDVI', 'LSWI']
    
    ## Init method to include indices calculation
    def __init__(self, sourcedir, outdir):
        super(L7TempStack, self).__init__(sourcedir, outdir, 'LE07')
        
    def setBandsLoc(self, key=None, bandloc=None):
        if key is not None and bandloc is not None:
            self.bands_loc[key] = bandloc
        else:
            ## Get names of files to stack in raster
            prodlist = list(filter(re.compile(r'^LE07.*[0-9]$').search, os.listdir(self.getSourceDir())))
            prodloclist = {}
            for band in self.bands_of_interest:
                prodloclist[self.standard_band_dict[band]] = list([self.getSourceDir()+x+'/'+band+'.tif' for x in prodlist])
            self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        if key is not None and tempdata is not None:
            self.bands_temporal_range[key] = tempdata
        else:
            temp_range = {}
            for key, value in self.getBandsLoc().items():
                temp_range[key] = list([dt.strptime(x.split('/')[-2][10:18], 
                                                                                '%Y%m%d').date() for x in value])
            self.bands_temporal_range = temp_range
        
def check_dir(directory):
    """ Helper method. Check if directory ends in slash, if not append one"""
    if not directory[-1] == '/':
        return directory+'/'
    else:
        return directory
    