# -*- coding: utf-8 -*-
"""
Created on Fri May 25 15:59:36 2018

@author: ASALAZAR
"""

import h2o
import xarray as xr
import pandas as pd

def prepare_dataset(dataframe):
    
    data_pd = dataframe.to_dataframe()
    
    # Leave only data with class tags
    data_pd = data_pd.loc[data_pd['class_mask'] == True]
    
    # Remove not valid data (cloudy pixels)
    try:
        data_pd = data_pd.loc[data_pd['mask'] == True]
    except:
        pass #Data does not have mask dimension e.g. Sentinel-1
        
    
    # Load data to h2o frame
    data = h2o.H2OFrame(data_pd)
    
    