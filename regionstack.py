# -*- coding: utf-8 -*-
"""
Created on Wed May 23 09:31:30 2018

@author: ASALAZAR
"""

from eo_stack import *
import os, re, shutil
import xarray as xa

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
        self.s1_ASC = self.__setDataset__('S1',orbit='ASCENDING')
        
        self.s1_DSC = self.__setDataset__('S1',orbit='DESCENDING')
        
        self.s2 = self.__setDataset__('S2')
        
        self.l7 = self.__setDataset__('LE07')
        
        self.l8 = self.__setDataset__('LC08')
        
    def __setDataset__(self, prodtype, orbit=None):
        
        # Try to process new products, if available
        try:
            self.__createDataset__(prodtype, orbit)
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
            ds = xa.open_mfdataset(list(map(lambda x: datadir+x, files)),
                  chunks={'time':1})
            return ds.sortby('time')
        except:
            print('No stacks available for {}'.format(prodtype))
            return None
    
    def __createDataset__(self, prodtype, orbit=None):
        
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
        