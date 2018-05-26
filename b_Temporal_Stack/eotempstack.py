# -*- coding: utf-8 -*-
"""
Module for constructing a multi-temporal earth observation dataset from pre-processed Sentinel-1, Sentinel-2, Landsat-7/8 raster files

@author: asalazar
"""

import os
import re
import datetime
import sys
import rasterio, dask
import numpy as np
import pandas as pd
import xarray as xr

from calendar import monthrange

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
        ## TODO assert that sourcedir finished in slash 
        self.source_directory = check_dir(sourcedir)
        self.prod_type = prodtype
        self.out_directory = outdir
        self.bands_loc = {} #Declare dictionary. Set in setBandsLoc
        self.bands_temporal_range = {} #Declare dictionary. Set in setTempData
        #self.orbit = None
        
        # Set the minimum required variables to construct the object. Implementation varies by product (subclass).
        self.setBandsLoc()
        self.setTempData()
        
        print('{} object initialized from {}'.format(self.prod_type, self.source_directory))
        
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
        for id, layer in enumerate(self.getTempData().itervalues().next()):
            if date == layer:
                return id
    
    def getBandXarray(self, band, time_range):
        """Returns an xarray object from all dates in a given band in a given timerange"""
        start_date = time_range[0]
        end_date = time_range[1]
        # Get the times within the time_range
        time = xr.Variable('time', pd.DatetimeIndex([pd.Timestamp(f) for f in self.getTempData(band) if start_date <= f <= end_date]))
        # Use numpy to filter band locations and return in band
        band_locations = np.array(self.getBandsLoc(band))[np.array(list(map(lambda f: start_date <= f <= end_date, self.getTempData(band))))].tolist()
        # Read all bands using open_rasterio xarray method
        arlist = [xr.open_rasterio(f) for f in band_locations]
        # Concatenate all dates
        da = xr.concat(arlist, dim=time)
        return da
    
class S1TempStack(eoTempStack):
    bands_of_interest = ['VV', 'VH']
    
    def __init__(self, sourcedir, outdir, orbit):
        # Call eo_tempstack initialization method
        self.orbit = orbit
        super(S1TempStack, self).__init__(sourcedir, outdir, 'S1')
    
    def setBandsLoc(self):
        ## Get directory names of pre-processed S1 products to read
        proddirs = filter(re.compile(r'S1_'+self.orbit+'.*data$').search, os.listdir(self.source_directory))
        
        ## Declare dictionary to store file location by polarizations
        prodlist = {}
        for prod_dir in proddirs:
            # List and filter img files for Sigma0 bands, store full path to file
            image_files = filter(re.compile(r'Sigma0.*img$').search, os.listdir(self.source_directory+prod_dir))
            prodlist[prod_dir.split('_')[2]] = list(map(lambda x: self.getSourceDir()+prod_dir+'/'+x, image_files))
        # Store the dictionary of files location as instance variable
        self.bands_loc = prodlist
        
    def setTempData(self, key=None, tempdata=None):
        ## Declare dictionary to store dates by polarization
        temp_range = {}
        for key, value in self.getBandsLoc().iteritems():
            temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-1].split('_')[3],'%d%b%Y').date(), value))
            # Store dictionary as instance variable
        self.bands_temporal_range = temp_range
    
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
        #try:
            # Calculate ranges of the available data. We consider monthly ranges, from first to last day of each month
            # print(self.getTempData())
            
            time_ranges = list(map(lambda x: (datetime.datetime.strptime('01'+x.strftime('%m%Y'), '%d%m%Y').date(),
                            datetime.datetime.strptime(str(monthrange(x.year, x.month)[1])+x.strftime('%m%Y'), '%d%m%Y').date()), self.getTempData().itervalues().next()))
            
            # Iterate over unique date ranges (i.e. per month)            
            for time_range in set(list(time_ranges)):
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
                        print('{} band for S1 {} {} could not be processed'.format(band, self.orbit, month))
                # Merge all bands in a single array
                xa = xr.merge(xarrays)
                # Drops unused bands
                xa = xa.drop('band')
                # Write array as netcdf
                xa.to_netcdf(self.out_directory+self.prod_type+'_'+self.orbit+'_'+month+'.nc')
                
        #except:
        #    e = sys.exc_info()
        #    print('Stacking {} failed: {} {} {}'.format(self.prod_type, e[0], e[1], e[2]))

class opticalTempStack(eoTempStack):
    ## Init method to include indices calculation
    def __init__(self, sourcedir, outdir, prodtype):
        # Call eo_tempstack initialization method
        super(opticalTempStack, self).__init__(sourcedir, outdir, prodtype)
        # Calculate indices at initialization
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
        # TODO make possible to update when new eo images are obtained
        # Check if dataset already exists
        #if os.path.isfile(self.out_directory+self.prod_type+'.nc'):
        #    ds = xr.open_dataset(self.out_directory+self.prod_type+'.nc', chunks={'time':1})
        # Otherwise create
        try:
            # Calculate ranges of the available data. We consider monthly ranges, from first to last day of each month
            time_ranges = list(map(lambda x: (datetime.datetime.strptime('01'+x.strftime('%m%Y'), '%d%m%Y').date(),
                            datetime.datetime.strptime(str(monthrange(x.year, x.month)[1])+x.strftime('%m%Y'), '%d%m%Y').date()), self.getTempData().itervalues().next()))
            
            # Iterate over unique date ranges (i.e. per month)            
            for time_range in set(list(time_ranges)):
                # Extract the month in format %Y%m
                month = time_range[0].strftime('%Y%m')
                # Put arrays for each band in list
                xarrays = []
                # Retrieve all bands
                for band in self.getBandsLoc().keys():
                    if self.prod_type == 'S2' and band in ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']:
                        x = self.getBandXarray(band, time_range).isel(band=0).drop('wavelength').astype(np.uint16, copy=False)
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
                # Drops unused bands
                xa = xa.drop('band').drop(['qa_class', 'qa_cloud'])
                # Write array as netcdf
                xa.to_netcdf(self.out_directory+self.prod_type+'_'+month+'.nc')
                
        except:
            e = sys.exc_info()
            print('Stacking {} failed: {} {} {}'.format(self.prod_type, e[0], e[1], e[2]))

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
            prodlist = filter(re.compile(r'^S2.*data$').search, os.listdir(self.getSourceDir()))
            prodloclist = {}
            for band in self.bands_of_interest:
                prodloclist[self.standard_band_dict[band]] = list(map(lambda x: self.getSourceDir()+x+'/'+band+'.img', prodlist))
            self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        if key is not None and tempdata is not None:
            self.bands_temporal_range[key] = tempdata
        else:
            temp_range = {}
            for key, value in self.getBandsLoc().iteritems():
                try:
                    temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-2][11:19], 
                                                                                '%Y%m%d').date(),value))
                except:
                    try:
                        temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-2][47:55], 
                                                                                '%Y%m%d').date(),value))
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
            prodlist = filter(re.compile(r'^LC08.*[0-9]$').search, os.listdir(self.getSourceDir()))
            prodloclist = {}
            for band in self.bands_of_interest:
                prodloclist[self.standard_band_dict[band]] = list(map(lambda x: self.getSourceDir()+x+'/'+band+'.tif',
                                                                 prodlist))
            self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        if key is not None and tempdata is not None:
            self.bands_temporal_range[key] = tempdata
        else:
            temp_range = {}
            for key, value in self.getBandsLoc().iteritems():
                temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-2][10:18], 
                                                                                '%Y%m%d').date(),value))
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
            prodlist = filter(re.compile(r'^LE07.*[0-9]$').search, os.listdir(self.getSourceDir()))
            prodloclist = {}
            for band in self.bands_of_interest:
                prodloclist[self.standard_band_dict[band]] = list(map(lambda x: self.getSourceDir()+x+'/'+band+'.tif',
                                                                 prodlist))
            self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        if key is not None and tempdata is not None:
            self.bands_temporal_range[key] = tempdata
        else:
            temp_range = {}
            for key, value in self.getBandsLoc().iteritems():
                temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-2][10:18], 
                                                                                '%Y%m%d').date(),value))
            self.bands_temporal_range = temp_range
        
def check_dir(directory):
    """
    Check if directory ends in slash, if not append one
    """
    if not directory[-1] == '/':
        return directory+'/'
    else:
        return directory