# -*- coding: utf-8 -*-
"""
Created on Wed May 23 09:31:30 2018

@author: ASALAZAR
"""

import os
import re
import shutil
import rasterio
import dask
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd

from eotempstack import *
from affine import Affine
from rasterio import features

class regionStack(object):
    """
    Class to interface the earth observation images stacking into multidimensional labelled arrays.
    Operates on the directory structure created by the pre-processing methods, calling the eoTempArray
    constructors.
    
    At instantiation tries to process any available preprocessed images for the specifies region and
    references all available netCDF files as dask arrays for each of the supported products (S1,S2,L07,L08).
    
    @params:
        regname (str): regions in ciat_monitor_crops project e.g. Saldana, Ibague, Casanare, Huila.
        dataserver (str): as environmental variable in system. default 'WIN_SVR_DATA' @ dapadfs
        
    """
    
    def __init__(self, regname, dataserver = 'WIN_SVR_DATA'):
        self.region_name = regname
        self.data_directory = os.environ[dataserver]+self.region_name+'/'
    
        #Loads all available products as dask arrays
        self.s1ASC = self.__setDataset('S1',orbit='ASCENDING')
        
        self.s1DSC = self.__setDataset('S1',orbit='DESCENDING')
        
        self.s2 = self.__setDataset('S2')
        
        self.l7 = self.__setDataset('LE07')
        
        self.l8 = self.__setDataset('LC08')
        
    def __setDataset(self, prodtype, orbit=None):
        
        # Try to process new products, if available
        try:
            self.__createDataset(prodtype, orbit)
        except:
            print('Creation of new {} stacks failed, reading existing stacks'.format(prodtype))
        
        # Try to read existing stacks
        try:
            datadir = self.data_directory+'stack/'
            if orbit is None:
                files = filter(re.compile(r'^'+prodtype+'.*').search, os.listdir(datadir))
            else:
                files = filter(re.compile(r'^'+prodtype+'_'+orbit+'.*').search, os.listdir(datadir))
            print('Reading {} {} stack files'.format(len(files), prodtype))
            ds = xr.open_mfdataset(list(map(lambda x: datadir+x, files)),
                  chunks={'time':1})
            return ds.sortby('time')
        except:
            print('No stacks available for {}'.format(prodtype))
            return None
    
    def __createDataset(self, prodtype, orbit=None):
        
        sourcedir = self.data_directory+'pre/'
        outdir = self.data_directory+'stack/'
        stackeddir = self.data_directory+prodtype+'/stacked/'
        
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            print('New directory {} was created'.format(outdir))
        
        if prodtype == 'S1':
            S1TempStack(sourcedir, outdir, orbit = orbit).createXDataset()
        elif prodtype == 'S2':
            S2TempStack(sourcedir, outdir).createXDataset()
        elif prodtype == 'LE07':
            L7TempStack(sourcedir, outdir).createXDataset()
        elif prodtype == 'LC08':
            L8TempStack(sourcedir, outdir).createXDataset()
        
        # Relocate files to avoid reprocessing
        if not os.path.exists(stackeddir):
            os.makedirs(stackeddir)
            print('New directory {} was created'.format(stackeddir))
        
        if prodtype == 'S1':
            files = filter(re.compile(r'^'+prodtype+'_'+orbit+'.*').search, os.listdir(sourcedir))
            for f in files:
                shutil.move(sourcedir+f, stackeddir)
        else:
            files = filter(re.compile(r'^'+prodtype+'.*').search, os.listdir(sourcedir))
            for f in files:
                shutil.move(sourcedir+f, stackeddir)
    
    def regionTrainingClasses(self, shapefile=None, testset = 0.0):
        """
        @params
            testset (float): fraction of shapes to be reserved for testing
        """
        # Definition of location to write files
        out_dir = getattr(self, 'data_directory')+'training/'
        
        # Create directory if does not exist
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        
        if shapefile is not None:
            # Take any Sentinel-2 index calculation as reference file
            reference_file = out_dir[:-9]+'S2/stacked/'+filter(re.compile(r'.*tif$').search,os.listdir(out_dir[:-9]+'S2/stacked/'))[0]
            
            fields_shp = gpd.read_file(shapefile).to_crs({'init': 'epsg:32618'})
            
            # Check if some shapes will be reserved for testing 
            if testset > 0.0:
                train=fields_shp.sample(frac=(1-testset),random_state=200)
                test=fields_shp.drop(train.index)
                
                try:
                    phc_dates = rasterizeClassdf(train, reference_file, out_dir, nametag = '')
                    train_dates = rasterizeClassdf(train, reference_file, out_dir, nametag = 'train')
                    test_dates = rasterizeClassdf(test, reference_file, out_dir, nametag = 'test')
                except:
                    print('The shapefile could not be rasterized')
                    
            else:
                try:
                    phc_dates = rasterizeClassdf(train, reference_file, out_dir, nametag = '')
                except:
                    print('The shapefile could not be rasterized')
            
            #phc_dates = []
        
        products_to_read = ['train']#, 'test']
        
        # Try to create if does not exist
        for eo_ds in products_to_read:
            
            netcdf_filename = out_dir+eo_ds+'.nc'
            
            if os.path.isfile(netcdf_filename):
                pass 
            else:
                #try:
                print('Attempting {} class dataset creation'.format(eo_ds))
                self.__mergeTrainingClasses(out_dir, phc_dates).to_netcdf(netcdf_filename)
                #except:
                #    print('{} class dataset creation failed'.format(eo_ds))
        
        # Try to read
        for eo_ds in products_to_read:
            
            netcdf_filename = out_dir+eo_ds+'.nc'
            
            try:
                self.__readTrainingDataset(netcdf_filename, eo_ds)
                
            except:
                print('{} class dataset could not be read'.format(eo_ds))
    
    def __readTrainingDataset(self, netcdf_filename, eods):
        tds = xr.open_mfdataset(netcdf_filename,chunks={'time':1})
        
        try:
            tds['mask'] = tds.mask.astype(bool)
        except:
            pass
        
        print('Adding dataset with classes as {} attribute of regionStack'.format(eods))
        setattr(self, eods, tds)
        
    def __mergeTrainingClasses(self, outdir, classdates):
        
        s1a = getattr(self, 's1ASC')
        s1d = getattr(self, 's1DSC')
        s2 = getattr(self, 's2')
        l7 = getattr(self, 'l7')
        l8 = getattr(self, 'l8')
        
        opt_bands = ['NDVI','LSWI']
        
        optical_dataset = xr.merge([s2[opt_bands],l7[opt_bands],l8[opt_bands]])
        
        # Extract all dates with any eo image available
        #time = np.array(pd.DatetimeIndex(optical_dataset.time.values).date)
        #time = list(set(list(map(lambda x: x.strftime('%Y%m%d'), time.tolist()))))
        
        # Read all dates of class data
        time = xr.Variable('time', pd.DatetimeIndex([pd.Timestamp(f) for f in classdates]))
        arlist = [xr.open_rasterio(f) for f in list(map(lambda x: outdir+'class_'+x+'.tif', classdates))]
        phc_data = xr.concat(arlist, dim=time).isel(band=0).drop('band')
        
        # Extract the intersect dates
        intersect_dates = np.intersect1d(phc_data.time.values, optical_dataset.time.values)
        
        final_dates = chooseBestQuality(optical_dataset, intersect_dates)
        
        radar_dataset = self.assignAproxRadar(optical_dataset.sel(time=final_dates))
        
        sat_dataset = optical_dataset.sel(time=final_dates).merge(radar_dataset).transpose('time','x','y')
        
        # Temporal subset of class data
        class_data = phc_data.sortby('time').sel(time=final_dates).transpose('time','x','y')
        
        # Convert class to uint8 dtype
        class_data = class_data.astype(np.uint8, copy=False)
        # Merge to dataset
        class_data.name = 'class'
        merged_ds = xr.merge([sat_dataset, class_data])
        
        # Make a mask for pixels for training and pixels for testing
        #self.        
        
        ## look for train rasters
        #try:
        #train_arlist = [xr.open_rasterio(f) for f in list(map(lambda x: outdir+'class_'+x+'train.tif', classdates))]
        #train_data = xr.concat(train_arlist, dim=time).isel(band=0).drop('band').values
        merged_ds.coords['train'] = (('time', 'x', 'y'), make_train_test_mask(outdir, final_dates, 'train'))
        #except:
        #    print('No train tags for pixels available')
        
        #try:
        #test_arlist = [xr.open_rasterio(f) for f in list(map(lambda x: outdir+'class_'+x+'test.tif', classdates))]
        #test_data = xr.concat(train_arlist, dim=time).isel(band=0).drop('band').values
        merged_ds.coords['test'] = (('time', 'x', 'y'), make_train_test_mask(outdir, final_dates, 'test'))
        #except:
        #    print('No test tags for pixels available')
        
        return merged_ds.chunk({'time':1})
        
    def assignAproxRadar(self, optical_dataset):
        """"""
        orbits = []
        for orbit in ['ASC', 'DSC']:
            renamed = rename_variables(getattr(self, 's1'+orbit), orbit)
            
            candidates = []
            
            for idx, time in enumerate(optical_dataset.time.values):
                try:
                    asc = renamed.sel(time=np.datetime64(time), method='nearest', 
                                      tolerance=np.timedelta64(10, 'D'))
                    asc['time'] = time
                    candidates.append(asc)
                
                except KeyError:
                    print('Closest available image to {} for {} is more than 10 days away'.format(time, orbit))
                    
            assig = xr.concat(candidates, dim='time')
            
            orbits.append(assig)
        
        return orbits[0].merge(orbits[1])
        

def make_train_test_mask(outdir, dates, classtag):
    """
    @params
        outdir (str): location of the rasterized images
        dates ([np.datetime64]): list of dates of dates of the raster
        classtag (str): tag of the rasterized image name
    """
    classdates = list(map(lambda t: str(t)[:4]+str(t)[5:7]+str(t)[8:10], dates))
    arlist = [xr.open_rasterio(f) for f in list(map(lambda x: outdir+'class_'+x+classtag+'.tif', classdates))]
    time = xr.Variable('time', pd.DatetimeIndex([pd.Timestamp(f) for f in dates]))#pd.DatetimeIndex(dates)
    data_array = xr.concat(arlist, dim=time).isel(band=0).drop('band').transpose('time','x','y').values
    #train_values = train_data.values
    mask = np.zeros(data_array.shape)
    mask[data_array>0.0]=1
    
    return mask.astype(np.int16)
    

def chooseBestQuality(opticalds, subset_dates=None, maxcloudcover = 0.2, mintempdistance = 5):
    """
    Selects the best dates complying  in an optical dataset
    
    @params
    opticalds (xarray.Dataset): dataset with coordinates time,x,y
    subset_dates ([np.datetime64]): list of dates to subset the dataset
    maxcloudcover (float): proportion of acceptable cloud covered pixels
    mintempdistance (int): min temporal distance in days for between images
    """
    if subset_dates is None:
        subset_dates = opticalds.time.values
    
    quality = opticalds.mask.sel(time=subset_dates).mean(dim=['x', 'y']).compute()
    
    # Get only the images with sufficient quality
    q_images = quality.where(quality>(1.0-maxcloudcover), drop=True).sortby('time')
    
    # Calculate distance in time when less than 5 days
    temp_distance = q_images.time.diff('time')<np.timedelta64(mintempdistance, 'D')
    # List to append resulting times
    times = []
    
    for idx, value in enumerate(q_images.time.values):
        #print '{}, {}'.format(idx,value)
        if idx>0 and temp_distance.isel(time=idx-1):
            # Get the max value
            window = q_images.isel(time=slice(-1+idx,1+idx))
            # Append the date of the max value
            times.append(window.where(window==window.max(), drop=True).time.values[0])
        else:
            times.append(value)
    
    return list(set(times))

def rename_variables(dataset,suffix,sep='_',ommit_vars=['time','x','y',]):
    variables = dataset.variables.keys()
    for var in ommit_vars:
        variables.remove(var)
    r_dict = {}
    map(lambda x: r_dict.update({x:x+sep+suffix}), variables)
    return dataset.rename(r_dict)
    
def rasterizeClassdf (geopandasdf, referencefile, location, nametag = ''):
    """ Rasterizes a pandas dataframe, writing the raster files as tif. Returns a list of dates
    @ params
    geopandasdf (geopandas): dataframe with classes as attributes and 
                  dates as column names in format 'X%Y&m&d'
    referencefile (str): location of reference file to rasterize
    location (str): location to write rasterized tif files
    nametag (str): nametag to add to the tif file name
    """
    
    ## Rasterize shapefile
    rst = rasterio.open(referencefile)
    meta = rst.meta.copy()
    meta.update(compress='lzw', dtype=rasterio.float64)
    
    phc_dates = []
    
    for column in geopandasdf:
        if geopandasdf[column].name[0] == 'X':
        
            # Extract the date
            date = geopandasdf[column].name[1:]
            
            out_file = location+'class_'+date+nametag+'.tif'
            
            phc_dates.append(date)
            
            if not os.path.isfile(out_file):
                with rasterio.open(out_file, 'w', **meta) as out:
                    out_arr = out.read(1)
                    shapes = ((geom,value) for geom, value in zip(geopandasdf.geometry, geopandasdf[column].astype('float64')))
                    burned = rasterio.features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
                    out.write_band(1, burned.astype(rasterio.float64))
    
    # Return a list of dates
    return phc_dates
    
    