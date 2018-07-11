# -*- coding: utf-8 -*-
"""
Created on Wed May 23 09:31:30 2018

@author: ASALAZAR
"""

import os
import re
import sys
import shutil
import rasterio
import dask

import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import eotempstack as eots

from datetime import datetime as dt
from calendar import monthrange

#from affine import Affine

class regionStack(object):
    """
    Class to interface the earth observation images stacking into multidimen-
    sional labelled arrays. Operates on the directory structure created by the
    pre-processing methods, calling the eoTempArray constructors.        
    """
    
    def __init__(self, regname, dataserver = 'WIN_SVR_DATA', attrs=['S1_ASCENDING','S1_DESCENDING','S1_ASCENDING_GLCM','S1_DESCENDING_GLCM','S2','LE07','LC08'], ignore_harm=False):
        """
        At instantiation tries to process any available preprocessed images for
        the specified region and references all available netCDF files as dask
        arrays for each of the supported products (S1,S2,L07,L08).
    
        @params:
            regname (str): regions in ciat_monitor_crops project e.g. Saldana, 
                                                    Ibague, Casanare, Huila.
            dataserver (str): directory where data is stored following pre-processing conventions
                                default environmental variable in system. default
                                os.environ['WIN_SVR_DATA'] @ dapadfs
            attrs ([str]): earth observation images to set as attributes in object. defaults to all supported
                          Supports: 'S1_ASCENDING' - Sentinel-1 ascending orbit intensity in dB units
                                    'S1_DESCENDING' - Sentinel-1 descending orbit intensity in dB units
                                    'S1_ASCENDING_GLCM' - GLCM texture stats for Sentinel-1 in asc orbit
                                    'S1_DESCENDING_GLCM' - GLCM texture stats for Sentinel-1 in asc orbit
                                    'S2' - Sentinel-2 Level-2 BOA reflectance
                                    'LE07' - Landsat-7 Level-2
                                    'LC08' - Landsat-8 Level-2
        
        """
        self.region_name = regname
        self.data_directory = os.environ[dataserver]+self.region_name+'/'
    
        #Loads all available products as dask arrays
        #Iterate attributes to set
        for product in attrs:
            if product[:2] == 'S1' and product[-4:] == 'GLCM':
                setattr(self, product, self.__setDataset(product[-4:]+'_'+product[:2],orbit=product[3:-5]))
            elif product[:2] == 'S1':
                setattr(self, product, self.__setDataset(product[:2] ,orbit=product[3:]))
            elif product[:4] == 'LC08':
                if not ignore_harm:
                    setattr(self, product, self.__setDataset('h'+product)) #try to read harmonized products, if exist
                try:
                    if getattr(self, product) is not None: #if harmonized dataset was loaded, notify user
                        print('Harmonized Landsat-8 dataset was read')
                    else:
                        raise AttributeError('Harmonized dataset was not properly read')
                except AttributeError: # then, read the original dataset
                    setattr(self, product, self.__setDataset(product))
            else:
                setattr(self, product, self.__setDataset(product))
    
    def __setDataset(self, prodtype, orbit=None):
        """
        Wrapper to access earth observation temporal stacks. Calls method to
        create new stacks of pre-processed data in pre/ directory, then reads
        existing stack files.
        
        params:
            prodtype (str): S1, S2, LE07, LC08
            orbit (str): Sentinel-1 orbit, ASCENDING or DESCENDING
        """
        # Try to process new products, if available
        #try:
        self.__createDataset(prodtype, orbit)
        #except:
        #    print(('Creation of new {} stacks failed, reading existing stacks'
        #           .format(prodtype)))
        
        # Try to read existing stacks
        try:
            # Default directory of stacks
            datadir = self.data_directory+'stack/'
        
            # Lookup files to read and open Dataset with dask
            if orbit is None:
                files = list(filter(re.compile(r'^'+prodtype+'.*').search,
                                    os.listdir(datadir)))
            else:
                files = list(filter(re.compile(r'^'+prodtype+'_'+orbit+'.*').
                                    search, os.listdir(datadir)))
            print('Reading {} {} stack files'.format(len(files), prodtype))
            ds = xr.open_mfdataset(list([datadir+x for x in files]),
                  chunks={'time':1,'x':1000,'y':1000}, parallel=True)
            return ds.sortby('time')
        except:
            e = sys.exc_info()
            if orbit is not None:
                pm = prodtype+orbit+' orbit'
            else:
                pm = prodtype
            print('Reading {} stack files failed: {} {} {}'.format(pm, e[0], e[1], e[2]))
            
            return None
    
    def __createDataset(self, prodtype, orbit=None):
        """
        Class to create eoTempStack objects for the specified product. Processes
        and relocated files in pre/ directory to a stacked/ directory.
        
        params:
            prodtype (str): 
            orbit (str): 
        """
        
        sourcedir = self.data_directory+'pre/'
        outdir = self.data_directory+'stack/'
        stackeddir = self.data_directory+prodtype+'/stacked/'
        
        #By default, do not move files (to be moved only if create dataset succeeds
        move=False
        
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            print(('New directory {} was created'.format(outdir)))
        
        if prodtype == 'S1':
            eots.S1TempStack(sourcedir, outdir, orbit = orbit).createXDataset()
        elif prodtype == 'S2':
            eots.S2TempStack(sourcedir, outdir).createXDataset()
        elif prodtype == 'LE07':
            eots.L7TempStack(sourcedir, outdir).createXDataset()
        elif prodtype == 'LC08':
            eots.L8TempStack(sourcedir, outdir).createXDataset()
        elif prodtype == 'GLCM_S1':
            eots.S1TextureTempStack(sourcedir, outdir, orbit = orbit).createXDataset()
        elif prodtype == 'hLC08':
            pass
        else:
            raise NotImplementedError('Chosen product type is not implemented.')
        
        # Relocate files to avoid reprocessing
        move=True
        
        # create directory
        if not os.path.exists(stackeddir):
            os.makedirs(stackeddir)
            print(('New directory {} was created'.format(stackeddir)))
        # get files names
        if prodtype == 'S1':
            files = list(filter(re.compile(r'^'+'S1_SPF'+'_'+orbit+'.*').
                                search, os.listdir(sourcedir)))
        elif prodtype == 'GLCM_S1':
            files = list(filter(re.compile(r'^'+prodtype+'_'+orbit+'.*').
                                search, os.listdir(sourcedir)))
        else:
            files = list(filter(re.compile(r'^'+prodtype+'.*').search,
                                os.listdir(sourcedir)))
        # move files
        if move is True:
            for f in files:
                shutil.move(sourcedir+f, stackeddir)
    
    ##########################################################################
    ##########################################################################
    #################### Sentinel-Landsat Harmonization ######################
    ##########################################################################
    ##########################################################################
    
    def harmonize_L8(self, band='red', rewrite=False):
        """
        Implement harmonization with Sentinel-2
        
        Args:
            s2ds (xr.DataArray): Sentinel-2 data array
            band (str): Band to use to perform the correction
        """
        
        import harmonizationsl as h
        
        #Get the best dates (for both images)
        best_s2_image = choose1BestQuality(getattr(self, 'S2'))[band].drop('time').drop('mask')
        best_l8_image = choose1BestQuality(getattr(self, 'LC08'))[band].drop('time').drop('mask')
        
        #Merge the datasets
        merged_dataset = xr.Dataset({'s2': best_s2_image, 'l8': best_l8_image})
        
        #Calculate offset
        offset = h.calculate_offset(merged_dataset.s2, merged_dataset.l8)
        
        #Get the dataset to correct
        lc08dataset = getattr(self, 'LC08')
        
        #Calculate time ranges
        time_ranges = list(map(lambda x: getMonthRange(x), lc08dataset.time.values))
        
        for drange in list(set(time_ranges)):
            
            month = drange[0].strftime('%Y%m')
            
            file_name = self.data_directory+'stack/hLC08_'+month+'.nc'
            
            if not os.path.isfile(file_name) or rewrite:
                
                print('Correcting {} dataset'.format(month))
                                
                #Apply bandpass correction
                cds = h.apply_bandpass_correction(lc08dataset.sel(time=slice(drange[0],drange[1])))
                
                #Apply geo correction
                cds = h.apply_geo_correction(cds, offset)
                
                #Write files
                cds.chunk({'time':1}).to_netcdf(file_name)
        
        setattr(self, 'LC08', xr.open_mfdataset(list([self.data_directory+'stack/hLC08_'+x[0].strftime('%Y%m')+'.nc' for x in time_ranges]),
                  chunks={'time':1,'x':1000,'y':1000}, parallel=True))
    
    ##########################################################################
    ##########################################################################
    ####################### Training dataset building ########################
    ##########################################################################
    ##########################################################################
    
    def regionTrainingClasses(self, shapefile=None, testset = 0.0):
        """
        Constructs a model training dataset based on a shapefile.
        
        @params
            shapefile (str): location of shapefile with classes. Shapefile
                        contains polygons with an attribute table of dates
                        as columns with format X%Y%m%d as column name and
                        dummy ints representing the classes for each date.
            testset (float): fraction of polygons to be reserved for testing.
        """
        # Definition of location to write files
        out_dir = getattr(self, 'data_directory')+'training/'
        
        # Create directory if does not exist
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
        
        if shapefile is not None:
            # Take any Sentinel-2 index calculation as reference file
            rfd = out_dir[:-9]+'S2/stacked/'
            
            reference_file = rfd+list(filter(re.compile(r'.*tif$').search,
                                             os.listdir(rfd)))[0]
            
            fields_shp = gpd.read_file(shapefile).to_crs({'init': 'epsg:32618'})
            
            # Check if some shapes will be reserved for testing 
            if testset > 0.0:
                train=fields_shp.sample(frac=(1-testset),random_state=200)
                test=fields_shp.drop(train.index)
            else:
                train=fields_shp
            
            # Rasterize train and test polygons
            #try:
            class_dates = rasterizeClassdf(train, reference_file,
                                             out_dir, classtag = 'train')
            if testset > 0.0:
                rasterizeClassdf(test, reference_file,out_dir,classtag = 'test')
            
            #except:
            #    print('The shapefile classes not be rasterized')
            
        # Define raster files to be read
        if testset > 0.0:
            raster_to_read = ['train','test']
        else:
            raster_to_read = ['train']
        
        if not os.path.isfile(out_dir+'class.nc'):
            print('Attempting {} class dataset creation'.format(raster_to_read))
            self.__mergeTrainingClasses(out_dir,class_dates,raster_to_read)
            
        #### WRITE FILES #### TODO
        
        #delayed_obj = mergedds.to_netcdf(out_dir+'optical_dataset.nc')
        
        #from dask.diagnostics import ProgressBar
        #with ProgressBar():
        #    results = delayed_obj.compute()
        
        #writeParquetDataset(mergedds, out_dir, raster_to_read)
        
        # Try to read
        self.__readTrainingDataset(out_dir)
    
    def calcTempVariation(self, ds_type, bands):
        """Calculate slope of the specified bands between two dates"""
        
        # Get the reference to the dataset
        dataset = getattr(self, ds_type).sortby('time')
        # Reconstruct the name of the file using methods in eotempstack class
        
        # Empty list to store new band names
        new_bands = []
        # Append new bands to dataset
        for band in bands:
            # Make first date of nan
            first_date = np.empty(dataset[band].sortby('time').isel(time=0).shape)
            first_date[:] = np.nan
            # List to append values by date
            c_arrays = [first_date]
            # Process all dates for a band
            for idx, time in enumerate(dataset.time[1:].values):
                print('Processing band {} for date {}'.format(band,time))
                # Calculate temporal trend
                c_arrays.append(calcTempTrend(dataset,band,ndate=-idx-1))
            # Stack all dates in np array    
            c = np.stack(c_arrays,axis=2)
            # Merge with xr dataset
            dataset[band+'_c'] = (['y','x','time'],c)
            # Append to new bands list
            new_bands.append(band+'_c')
        
        # Calculate monthly time ranges for each date
        time_ranges = list(map(lambda x: getMonthRange(x), dataset.time.values))
        # Iterate unique time ranges
        for time_range in set(list(time_ranges)):
            # Format month
            month = time_range[0].strftime('%Y%m')
            # Select bands and time range
            xa = dataset[new_bands].transpose('time','y','x').sel(time=slice(time_range[0], time_range[1]))
            # Write to file
            xa.to_netcdf(self.data_directory+'stack/'+ds_type+'_'+month+'c.nc')
            
    
    def __readTrainingDataset(self, outdir):
        
        try:
            cds = xr.open_mfdataset(outdir+'class.nc',chunks={'time':1}, parallel=True)
            
            print('Adding dataset with classes as class_ds attribute of regionStack')
            
            setattr(self, 'class_ds', cds)
        except:
            print('Class netcdf dataset could not be read')
        
    

    
    def __mergeTrainingClasses(self, outdir, classdates, raster_to_read,
                               maxcloudcover = 0.2):
        """
        Merge optical dates with max cloud cover and assigns closest radar dates
        """
        
        #dataset = getattr(self, ds_type)
        
        # Read all dates with class data
        #time = xr.Variable('time', pd.DatetimeIndex([pd.Timestamp(f) for f in
        #                                             classdates]))
        
        # Extract the intersect dates
        #intersect_dates = np.intersect1d(time.values,
        #                                 dataset.time.values)
        
        #if ds_type in ['S2','LE07','LC08']:
        #    intersect_dates = chooseBestQuality(dataset, intersect_dates, maxcloudcover=maxcloudcover)
        
        d_merge = []
        
        for raster in raster_to_read:
            d_merge.append(readClassArray(outdir,classdates,raster))
        
        dataset = xr.merge(d_merge)
        
        #time_ranges = list(map(lambda x: getMonthRange(x), intersect_dates))
        #print(time_ranges)
        #print(dataset.time.values)
        
        #for time_range in set(list(time_ranges)):
        #    
        #    month = time_range[0].strftime('%Y%m')
        #    
        xa = dataset[raster_to_read].transpose('time','y','x').sortby('time')#.sel(time=slice(time_range[0], time_range[1]))
        #
        xa.to_netcdf(outdir+'class.nc')
        
        #return intersect_dates
        
    def assignAproxRadar(self, optical_dataset):
        """Assign the closest S1 dates using the dates in provided dataset"""
        orbits = []
        for orbit in ['ASCENDING', 'DESCENDING']:
            renamed = rename_variables(getattr(self, 'S1_'+orbit), orbit)
            
            candidates = []
            
            for idx, time in enumerate(optical_dataset.time.values):
                try:
                    asc = renamed.sel(time=np.datetime64(time), method='nearest', 
                                      tolerance=np.timedelta64(10, 'D'))
                    asc['time'] = time
                    candidates.append(asc)
                
                except KeyError:
                    print(('Closest available image to {} for {} is more \
                           than 10 days away'.format(time, orbit)))
                    
            assig = xr.concat(candidates, dim='time')
            
            orbits.append(assig)
        
        return orbits[0].merge(orbits[1])
    
    def buildTrainingDataFrame(self, datadict, classsets=['train','test']):
        """Build a training dataframe for the data and variables specified in datadict"""
        
        try:
            classds = getattr(self, 'class_ds')
        except AttributeError:
            print('No training dataset has been assigned')
        
        intersect_dates = np.array([], dtype=np.datetime64)
        for key in datadict:
            intersect_dates = np.append(intersect_dates, getattr(self, key).time.values)
        #intersect_dates = np.array(list(set(intersect_dates)))
        intersect_dates = np.intersect1d(classds.time.values, intersect_dates)
        
        dataset = []        
        for key, bands in datadict.items():
            dataset.append(getattr(self, key)[bands].sel(time=np.intersect1d(getattr(self, key).time.values, intersect_dates)))
                
        dataset = xr.merge(dataset+[classds.sel(time=intersect_dates)]).chunk({'time':1,'y':1000,'x':1000})
        
        writeParquetDataset(dataset, self.data_directory, classsets)    
    
def writeParquetDataset(dataset, datadir, classsets):
    """Write a parquet file of the dataset using a dask array to subset only
    values with class data.
    
    params:
        dataset (xarray.Dataset): dataset to write as parquet
        datadir (str): parent folder to write class sets
        classets ([str]): names of colums in dataset with class values
                          with valid values are >0 e.g. ['train','test']
    """
    
    ds_location = datadir+'parquet/'
    
    for classn in classsets:
        # Iterate individual dates
        for index, date in np.ndenumerate(dataset.time.values):
            data_pd = dataset.sel(time=date).to_dask_dataframe()
            
            # Filter the values with class data
            data_values = data_pd[data_pd[classn]>0]
            file_location = ds_location+classn+str(index[0])+'/'
            
            # Write to parquet
            print('Writing {} parquet dataset to {} . {} of {}...'.format(classn, file_location, index[0]+1, dataset.time.values.size))
            dask.dataframe.to_parquet(data_values, file_location)
        

def calcTempTrend(dataset, band, ndate=-1, tempwindowsize=1):
    """
    Calculates the temporal trend of a given band. Returns a xa DataArray with the calulated trend.
    
    @params
        band (str): name of the band 
        ndate (int): index of the date of analysis, default is -1, i.e. last available date
        tempwindowsize (int): size of the temporal window in available dates distance
                              default value is 1, i.e. the previous available date
    """
    # Make temporal subset from last date and the tempwindowsize distance
    try:
        temp_subset = dataset.sortby('time').isel(time=[ndate-tempwindowsize,ndate])
    except:
        raise IndexError('The specified time indices are not in the dataset')
    try:
        temp_subset = temp_subset[band]
    except:
        raise IndexError('Band {} is not a dataset variable'.format(band))
    
    # Extract dates
    time = pd.DatetimeIndex(temp_subset.time.values).date
    time = np.array([x.toordinal() for x in time])
    
    # Mask the subset when mask is available
    try:
        temp_subset = temp_subset.where(temp_subset.mask)
    except:
        pass
    
    # Extract band values
    vals = temp_subset.values
    # Reshape to an array with as many rows as dates and as many columns as there are pixels
    vals_reshaped = vals.reshape(len(time), -1)
    
    # Declare array to store slope coefficient
    slope = np.empty(vals_reshaped[0].shape)
    slope[:] = np.nan
    
    # Make a mask of missing values (if any)in any of the dates
    mask = np.any(np.isnan(vals_reshaped), axis=0)
    # Use mask to get only valid pairs
    masked_vals = vals_reshaped.T[~mask].T
    
    try:
        # Do a first-degree polyfit with masked vals
        polyfit = np.polyfit(time, masked_vals, 1)
        # Put slope valid values in slope array
        np.place(slope, ~mask, polyfit[0].tolist())
    except np.linalg.LinAlgError:
        print('Value pairs are all nan. No good quality pixels are available for the specified dates')
    
    trends = slope.reshape(vals.shape[1], vals.shape[2])
    
    # Return np array with slope values
    return trends

def readClassArray(outdir, dates, classtag):
    """
    @params
        outdir (str): location of the rasterized images
        dates ([str]): list of dates of dates of the raster in format %Y%m%d
        classtag (str): tag of the rasterized image name
    """
    classdates = dates #list([str(t)[:4]+str(t)[5:7]+str(t)[8:10] for t in dates]) ##  this was for np.datetime64
    arlist = [xr.open_rasterio(f) for f in list([outdir+'class_'+x+classtag+
              '.tif' for x in classdates])]
    time = xr.Variable('time', pd.DatetimeIndex([pd.Timestamp(f) for f in dates]))
    data_array = xr.concat(arlist, dim=time).isel(band=0).drop('band').transpose('time','y','x')#.values
    
    data_array.name = classtag
    
    return data_array.astype(np.uint8, copy=False)#mask.astype(np.int16)
    

def choose1BestQuality(opticalds):
    """
    Selects the date with the least clouds in an optical dataset
    
    Parameters:
        opticalds (xr.Dataset): multi-temp dataset of optical rs images
    """
    try:
        quality = opticalds.mask.mean(dim=['x', 'y']).compute()
    except AttributeError:
        print('Only optical datasets with bool pixel quality mask are supported')
        raise
    
    return opticalds.isel(time=np.argmax(quality.values))
    

def chooseBestQuality(opticalds, subset_dates=None, maxcloudcover = 0.2, mintempdistance = 5):
    """
    Selects the best dates complying  in an optical dataset
    
    params:
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
    
    # Bool df, True when distance in time when less than mintempdistance
    temp_distance = q_images.time.diff('time')<np.timedelta64(mintempdistance, 'D')
    # List to append resulting times
    times = []
    
    # Iterate dates to select best date within the mintempdistance period
    for idx, value in enumerate(q_images.time.values):
        # If the time distance was calculated (from second index) and temp_distance true
        if idx>0 and temp_distance.isel(time=idx-1):
            # Get the max value
            window = q_images.isel(time=slice(-1+idx,1+idx))
            # Append the date of the max value
            times.append(window.where(window==window.max(), drop=True).time.values[0])
        else:
            times.append(value)
    
    return list(set(times))

def rename_variables(dataset,suffix,sep='_',ommit_vars=['time','x','y',]):
    """
    Method to rename variables of a dataset, ommiting the dimensions
    
    params:
        dataset (xarray.Dataset): dataset to rename
        suffix (str): suffix to append to the original name
        sep (str): string separator before suffix, defs to underscore
        
    returns:
        renamed dataset (xarray.Dataset)
    """
    # Get the original variable names
    variables = list(dataset.variables.keys())
    # Remove the dimension names    
    for var in list(dataset.dims.keys()):
        variables.remove(var)
    # Dictionary to map original and new names
    r_dict = {}
    # Populate dictionary with variable names to change
    list(map(lambda x: r_dict.update({x:x+sep+suffix}), variables))
    # Return renamed dataset
    return dataset.rename(r_dict)
    
def rasterizeClassdf (geopandasdf, referencefile, location, classtag = ''):
    """
    Rasterizes a pandas dataframe, writing the raster files as tif. Returns a
    list of dates.
    
    params:
        geopandasdf (geopandas): dataframe with classes as attributes and dates
                                 as column names ('X%Y%m%d')
        referencefile (str): location of reference file to rasterize
        location (str): location to write rasterized tif files
        classtag (str): name tag to add to the tif file name
    """
    
    ## Rasterize shapefile
    rst = rasterio.open(referencefile)
    meta = rst.meta.copy()
    meta.update(compress='lzw', dtype=rasterio.float64)
    
    class_dates = []
    
    for column in geopandasdf:
        if geopandasdf[column].name[0] == 'X':
        
            # Extract the date
            date = geopandasdf[column].name[1:]
            
            out_file = location+'class_'+date+classtag+'.tif'
            
            class_dates.append(date)
            
            if not os.path.isfile(out_file):
                with rasterio.open(out_file, 'w', **meta) as out:
                    out_arr = out.read(1)
                    shapes = ((geom,value) for geom, value in zip(geopandasdf.geometry,
                              geopandasdf[column].astype('float64')))
                    burned = rasterio.features.rasterize(shapes=shapes,
                                                         fill=0,
                                                         out=out_arr,
                                                         transform=out.transform)
                    out.write_band(1, burned.astype(rasterio.float64))
    
    # Return a list of dates
    return class_dates
    
def getMonthRange(date):
    """"""
    x = pd.to_datetime(str(date))
    first_date = dt.strptime('01'+x.strftime('%m%Y'), '%d%m%Y').date()
    last_date = dt.strptime(str(monthrange(x.year, x.month)[1])+x.strftime('%m%Y'), '%d%m%Y').date()
    
    return (first_date, last_date)
    