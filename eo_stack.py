# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 15:32:31 2018

@author: asalazar
"""

import rasterio, os, re, datetime
import numpy as np
import pandas as pd
import xarray as xr

#from rasterio.rio import stack

class eoTempStack:
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
        @prodtype (str): accepts 'S1', 'S2', 'L7', 'L8' where
                            'S1': set of .img files with the same extent and polarization, product from SNAP collocation Op
                            'S2': set of directories resulting from sen2cor processing Sentinel-2 level 2A product
        """        
        ##toDO assert that sourcedir finished in slash 
        self.source_directory = sourcedir
        self.prod_type = prodtype
        self.out_directory = outdir
        self.bands_loc = {} #Declare dictionary. Set in setBandsLoc
        self.bands_temporal_range = {} #Declare dictionary. Set in setTempData
        #self.tif_stack_location = {} #Declare dictionary. Set in setStackLoc
        
        # Set the minimum required variables to construct the object. Implementation varies by product (subclass).
        self.setBandsLoc()
        #self.setStackLoc()
        self.setTempData()
		
		#buildAllStack implements temporal stacking of each band as GeoTIFF. Currently not used.
        #self.buildAllStack()
		   
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
				readband = src.read()
		elif date is not None:
			with rasterio.open(bandsloc[self.getBandIndex(date)]) as src:
				readband = src.read()
		else:
			readband = []
			for loc in bandsloc:
				with rasterio.open(loc) as src:
					readband.append(src.read())
        return readband
    
    def getBandIndex(self, date):
        """Return the tempid in of the date in the  
        define if date is a string or a date type
        """
        for id, layer in enumerate(self.getTempData().itervalues().next()):
            if date == layer:
                return id
    
    def getBandXarray(self, band):
		"""Returns an xarray object from all dates in a given band"""
		time = xr.Variable('time', pd.DatetimeIndex([pd.Timestamp(f) for f in self.getTempData(band)]))
		arlist = [xr.open_rasterio(f) for f in self.getBandsLoc(band)]
		da = xr.concat(arlist, dim=time)
		return da
    
class S1TempStack(eoTempStack):
	bands_of_interest = ['VV', 'VH']
	
    def setBandsLoc(self):
        ## Get directory names of pre-processed S1 GRD products to process
        proddirs = filter(re.compile(r'S1.*data$').search, os.listdir(self.source_directory))
        
        ## Declare dictionary to store file location by polarizations
        prodlist = {}
        for prod_dir in proddirs:
            # List and filter img files for Sigma0 bands, store full path to file
            image_files = filter(re.compile(r'Sigma0.*img$').search, os.listdir(self.source_directory+prod_dir))
            prodlist[prod_dir.split('_')[1]] = list(map(lambda x: self.getSourceDir()+prod_dir+'/'+x, image_files))
        # Store the dictionary of files location as instance variable
        self.bands_loc = prodlist
        
    def setTempData(self, key=None, tempdata=None):
        #if key is not None and tempdata is not None:
        #    self.bands_temporal_range[key] = tempdata
        #else:
        ## Declare dictionary to store dates by polarization
        temp_range = {}
        for key, value in self.getBandsLoc().iteritems():
			temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-1].split('_')[3],'%d%b%Y').date(), value))
            # Store dictionary as instance variable
        self.bands_temporal_range = temp_range
	
	def getXarray(self):
		xarrays = []
		for band in self.bands_of_interest:
			x = self.getBandXarray(band).isel(band=0)
			xarrays.append(x)
		## Modify this to match standard_band_dict
		band = xr.Variable('band', pd.Index(self.bands_of_interest))
		da = xr.concat(xarrays, band)
		return da
	
class S2TempStack(eoTempStack):
    cloud_quality_limit = {'quality_cloud_confidence':9,'quality_scene_classification':3} 	#quality_cloud_confidence>9,
																							#quality_scene_classification == 3
    bands_of_interest = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'quality_cloud_confidence','quality_scene_classification']
    #TODO implement change names of bands to standard
	#standard_band_dict = {'B2':'blue', 'B3':'green', 'B4':'red', 'B8':'nir', 'B11':'swir1', 'B12':'swir2', 
	#'quality_cloud_confidence':'CloudQA1','quality_scene_classification':'CloudQA2'}
    calculated_bands = ['NDVI', 'LSWI'] 
    
	## TODO extend init method to include indices calculation
	
    def setBandsLoc(self, key=None, bandloc=None):
		if key is not None and bandloc is not None:
            self.bands_loc[key] = bandloc
		else:
			## Get names of files to stack in raster
			prodlist = filter(re.compile(r'^S2.*data$').search, os.listdir(self.getSourceDir()))
			prodloclist = {}
			for band in self.bands_of_interest:
				prodloclist[band] = list(map(lambda x: self.getSourceDir()+x+'/'+band+'.img', prodlist))
			self.bands_loc = prodloclist
        
    def setTempData(self, key=None, tempdata=None):
        if key is not None and tempdata is not None:
            self.bands_temporal_range[key] = tempdata
        else:
            temp_range = {}
            for key, value in self.getBandsLoc().iteritems():
                temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-2][11:19], 
                                                                                '%Y%m%d').date(),value))
            self.bands_temporal_range = temp_range
        
    def calcCloudCover(self):
        cloudcover = []
        for date in self.getTempData('quality_cloud_confidence'):
            quality_band = self.getBand('quality_cloud_confidence', date=date)
            cloudpixelsmask = quality_band > self.cloud_quality_limit
            cloudcover.append(float(np.sum(cloudpixelsmask)) / float(cloudpixelsmask.size))
        return cloudcover
        
    def calcIndex(self, index):
		"""index(str): accepts NDVI and LSWI"""
        calc_band = []
		locations = []
        
		for date in self.getTempData('B8'):
			nir = self.getBand('B8', date=date)
			# Allow division by zero
			np.seterr(divide='ignore', invalid='ignore')
			# Calculate index
			if index == 'NDVI': #(nir-red)/(nir+red)
				red = self.getBand('B4', date=date)
				calc_band.append((nir.astype(float) - red.astype(float)) / (nir + red))
			if index == 'LSWI': #(nir-swir1)/(nir+swir1)
				swir1 = self.getBand('B11', date=date)
                calc_band.append((nir.astype(float) - swir1.astype(float)) / (nir + swir1))
		
        # Define spatial characteristics of output object ## CHECK WHAT ARE THEY
		for id, date in enumerate(self.getBandsLoc('B8')):
			bandLocation = self.out_directory + index + '_' + id + '_.tif'
			locations.append(bandLocation)
			with rasterio.open(date) as src:
				kwargs = src.meta
			kwargs.update(dtype=rasterio.float32)
			with rasterio.open(bandLocation, 'w', **kwargs) as dst:
				dst.write_band(1, calc_band[id].astype(rasterio.float32))
		
        # Update object variables
        self.setStackLoc(index, locations)
        self.setTempData(index, self.getTempData('B8'))
    
    def getMaskedBand(self, key, tempid=None, date=None):
        if tempid is not None or date is not None:
            rmask = self.getBand('quality_cloud_confidence', tempid=tempid, date=date) > self.cloud_quality_limit
            return np.ma.masked_array(self.getBand(key, tempid=tempid, date=date), mask=rmask)
        else:
            raise NotImplementedError('Not Implemented. Please specify date or tempid.')
            
    def getXarray(self):
		xarrays = []
		for id, band in enumerate(self.bands_of_interest):
			if id<6:
				x = self.getBandXarray(band).isel(band=0).drop('wavelength')
			else:
				x = self.getBandXarray(band).isel(band=0)
			xarrays.append(x)
		## Modify this to match standard_band_dict
		band = xr.Variable('band', pd.Index(self.bands_of_interest))
		da = xr.concat(xarrays, band)
		return da