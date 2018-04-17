# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 15:32:31 2018

@author: asalazar
"""

import rasterio, os, re, datetime
import numpy as np
#from rasterio.rio import stack

class eoTempStack:
    """
    Class to construct and analyze a multi temporal stack of earth observation images. Supports Sentinel-1, Sentinel-2,
    Landsat 7 and Landsat 8 pre-processed products. All images in the same object have the same coordinate system, spatial
    extent and resolution. Required pre-processing:
        - For Sentinel-1...
        - For Sentinel-2...
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
        self.bands_orig_files = {} #Declare dictionary. Set in setOrigBandsLoc
        self.bands_temporal_range = {} #Declare dictionary. Set in setTempData
        self.stack_location = {} #Declare dictionary. Set in setStackLoc
        
        print('{} empty object initialized from {}'.format(self.prod_type, self.source_directory))
        
        # Set the minimum required variables to construct the object. Implementation varies by product (subclass).
        self.setOrigBandsLoc()
        self.setStackLoc()
        self.setTempData()
        self.buildAllStack()
        
    def getOrigBandsLoc(self, key=None):
        if key is not None:
            return self.bands_orig_files[key]
        else:
            return self.bands_orig_files
    
    def getTempData(self, key=None):
        if key is not None:
            return self.bands_temporal_range[key]
        else:
            return self.bands_temporal_range
    
    def buildAllStack(self):
        # Makes sure stack location is set to save the results
        for key in self.getStackLoc():
            self.buildStack(key)
    
    def buildStack(self, key):
        prodlist = self.getOrigBandsLoc(key)
        with rasterio.open(prodlist[0]) as src0:
            meta = src0.meta
            # update meta to reflect the number of layers
            meta.update(count = len(prodlist))
            # read each layer and write it to stack
            with rasterio.open(self.getStackLoc(key), 'w', **meta) as dst:
                for id, layer in enumerate(prodlist):
                    with rasterio.open(layer) as src1:
                        # band numbering in rasterio goes from 1
                        dst.write_band(id+1, src1.read(1))
                        # update metadata with date of image and polarization
                        dst.update_tags(id+1, image_date=self.getTempData(key)[id])
                        dst.update_tags(id+1, band_name=key)
    
    def getSourceDir(self):
        return self.source_directory
    
    def setStackLoc(self):
        # Declare a dictionary to store stack location by key
        stackLoc = {}
        for key in self.getOrigBandsLoc():
            stackLoc[key] = self.out_directory + '/' + key + '.tif'
        # Store dictionary as instance variable
        self.stack_location = stackLoc
    
    def getStackLoc(self, key=None):
        if key is not None:
            return self.stack_location[key]
        else:
            return self.stack_location
    
    def getBand(self, key, index=None, date=None):
        with rasterio.open(self.getStackLoc(key)) as src:
            if date is not None:
                band = src.read(self.getBandIndex(date)+1)
            elif index is not None:
                band = src.read(index+1)
            else:
                band = src.read()
        return band
    
    def getTags(self, key, index=None, date=None):
        with rasterio.open(self.stack_location) as src:
            if date is not None:
                tags = src.tags(self.getBandIndex(date)+1)
            elif index is not None:
                tags = src.tags(index+1)
            else:
                tags = src.tags()
        return tags
    
    def getBandIndex(self, date):
        """Return the index in of the date in the  
        define if date is a string or a date type
        """
        for id, layer in enumerate(self.getTempData().itervalues().next()):
            if date == layer:
                return id
        
class S1TempStack(eoTempStack):
                
    def setOrigBandsLoc(self):
        ## Get directory names of pre-processed S1 GRD products to process
        proddirs = filter(re.compile(r'S1.*data$').search, os.listdir(self.source_directory))
        
        ## Declare dictionary to store file location by polarizations
        prodlist = {}
        for prod_dir in proddirs:
            # List and filter img files for Sigma0 bands, store full path to file
            image_files = filter(re.compile(r'Sigma0.*img$').search, os.listdir(self.source_directory+prod_dir))
            prodlist[prod_dir.split('_')[1]] = list(map(lambda x: self.getSourceDir()+prod_dir+'/'+x, image_files))
        # Store the dictionary of files location as instance variable
        self.bands_orig_files = prodlist
        
    def setTempData(self):
        ## Declare dictionary to store dates by polarization
        temp_range = {}
        for key, value in self.getOrigBandsLoc().iteritems():
            temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-1].split('_')[3], '%d%b%Y').date(),
                                                 value))
        # Store dictionary as instance variable
        self.bands_temporal_range = temp_range

        
class S2TempStack(eoTempStack):
    cloud_quality_limit = 9
    bands_of_interest = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'quality_cloud_confidence',
                         'quality_scene_classification']
    calculated_bands = ['NDVI', 'NDWI']
    
    def setOrigBandsLoc(self):
        ## Get names of files to stack in raster
        prodlist = filter(re.compile(r'^S2.*data$').search, os.listdir(self.getSourceDir()))
        
        prodloclist = {}
        
        for band in self.bands_of_interest:
            prodloclist[band] = list(map(lambda x: self.getSourceDir()+x+'/'+band+'.img', prodlist))
        
        self.bands_orig_files = prodloclist
            
    def setTempData(self):
        temp_range = {}
        for key, value in self.getOrigBandsLoc().iteritems():
            temp_range[key] = list(map(lambda x: datetime.datetime.strptime(x.split('/')[-2][11:19], '%Y%m%d').date(),
                                                 value))
        self.bands_temporal_range = temp_range
        
    def calcCloudQuality(self):
        """"Recieves a rasterio read product.. what is the name?"""
        cloudcover = []
        for date in self.getTempData('quality_cloud_confidence'):
            quality_band = self.getBand('quality_cloud_confidence')
            cloudpixelsmask = quality_band > self.cloud_quality_limit
            cloudcover.append(1 - float(np.sum(cloudpixelsmask)) / float(cloudpixelsmask.size))
        self.cloud_cover = cloudcover
        
    def calcBand(self, band):
        # We handle the connections with "with"
        #with rasterio.open(bands[0]) as src:
        #    b3 = src.read(1)
        #with rasterio.open(bands[1]) as src:
        #    b4 = src.read(1)
        calc_band = []
        
        if band == 'NDVI':
            for date in self.getTempData('B3'):
                b3 = self.getBand('B3', date=date)
                b4 = self.getBand('B4', date=date)
                # Allow division by zero
                np.seterr(divide='ignore', invalid='ignore')
                # Calculate NDVI
                calc_band.append((b4.astype(float) - b3.astype(float)) / (b4 + b3))
            # Define spatial characteristics of output object (basically they are analog to the input)
            with rasterio.open(self.getBand('B3', index=0)) as src:
                kwargs = src.meta
            # Update kwargs (change in data type)
            kwargs.update(dtype=rasterio.float32)
            # Write raster stack of results
            with rasterio.open(self.out_directory + '/' + band + '.tif', 'w', **kwargs) as dst:
                for id, layer in enumerate(calc_band):
                        # band numbering in rasterio goes from 1
                        dst.write_band(id+1, layer.astype(rasterio.float32))
                        # update metadata with date of image and band name
                        dst.update_tags(id+1, image_date=self.getTempData('B3')[id])
                        dst.update_tags(id+1, band_name=band)
            
        