# -*- coding: utf-8 -*-
"""
Module to extract dataset variables for polygons in a shapefile

Created on Tue Jul 17 14:55:02 2018

@author: ASALAZAR
"""

import os

import geopandas as gpd
import dask.dataframe as dd
import pandas as pd
import numpy as np
import xarray as xr

from shapely import geometry

def prepare_dataset(inloc, outloc, date_class, shapefile, crs='epsg:32618', seed=33):
    """
    Extract xarray dataset variables for locations in the shapefile
    
    Args:
        inloc (str): path of the location of datasets to process
        outloc (str): path of location to write resulting dataframe as parquet
        crs (str): crs corresponding to dataset coordinates, default 'epsg:32618'
        date_class ([(str, np.datetime64)]): list of tuples with the first
                        element indicating the column name with the phenology obs
                        data in the shapefile and the second the corresponding
                        date. e.g ('X20151222', np.datetime64('2015-12-22'))
        shapefile (str): path of location of the shapefile with class data
                        regions in shapefile must have a lot identifier,
                        coded as (IDLote)
        seed (int): seed for the generation of reproducible datasets, default 33
        
    """
    
    #Read shapefile
    fields_shp = gpd.read_file(shapefile).to_crs({'init': crs})
    
    #Indicate if used as train or test
    # based on dictionary for each dic[date][polygon_index]:test/train
    ## eg {np.datetime('2015-12-22'):'2WKSD7':'train'} 
    class_dict = assign_polygons_to_class(date_class, fields_shp, seed)
    
    dss = os.listdir(inloc)
    #Open the datasets
    dss = list(map(lambda ds: xr.open_dataset(inloc+ds), dss))
    
    date_dataframes = []
    
    for date in date_class:
        date_dataframes.append(create_date_dataframe(dss,
                                                     fields_shp,
                                                     date[1],
                                                     date[0],
                                                     class_dict[date[1]]))
    
    data = dd.concat(date_dataframes, axis=0, interleave_partitions=True)
    
    #Write results with valid polygon class data outloc
    dd.to_parquet(data[data['vclass']>0], outloc)

def assign_polygons_to_class(passed_list, fields_shp, seed):
    """
    Performs stratiefied dataset split based on number of phenology observations in the shapefile
    
    passed_list ([(str, np.datetime64)]): list of tuples with the first
                        element indicating the column name with the class
                        data in the shapefile and the second the corresponding
                        date. e.g ('X20151222', np.datetime64('2015-12-22'))
        fields_shp (str): path of location of the shapefile with class data
    """
    
    classes = np.empty(0)
    polygons = np.empty(0).astype(str)
    
    #Get the classes and polygons to be used
    for date in passed_list:
        classes = np.append(classes, fields_shp[date[0]])
        polygons = np.append(polygons, fields_shp['IDLote'])
    
    #Use this to know the date for assigning at the end
    npolygons = fields_shp.shape[0]
    
    #Use method to map the number of polygons for each class and their indices
    seqidx = list_indices(classes)
    
    
    #Randomly shuffle indices in a new dictionary
    unsorted = {}
    
    for k,val in seqidx.items():
        if not np.isnan(k):
            unsorted[k] = np.random.RandomState(seed=seed).permutation(np.array(val)).tolist()
    
    #Split the train and test datasets in around 70% train and 30% test
    train = []
    test = []
    
    for k,val in unsorted.items():
        split_loc = int(round(len(val)*0.7,0))
        train = train + val[:split_loc]
        test = test + val[split_loc:]
        
    # Declare final dictionary    
    final_dict = {}
    for element in passed_list:
        final_dict[element[1]] = {}
    
    # Populate final dictionary using split indices
    for idx in train:
        date = passed_list[idx // npolygons][1]
        final_dict[date][polygons[idx]] = 'train'
    
    for idx in test:
        date = passed_list[idx // npolygons][1]
        final_dict[date][polygons[idx]] = 'test'
        
    return final_dict

def list_indices(seq):
    """
    Make a dictionary of values of a sequence with a list of their indices in the original sequence
    
    Args:
        seq (list): a list of elements
    
    Returns:
        counter (defaultdic): dictionary with elements as keys and indices
                                as values
    """
    from collections import defaultdict
    counter = defaultdict(list)
    for i,item in enumerate(seq):
        counter[item].append(i)
    return counter

def create_date_dataframe(dss, fields_shp, doa, class_col_name, class_use_dict):
    
    #Put intermediate data frames in a list
    field_dataframes = []
    
    #Iterate geo-dataframe
    for idx, shape in fields_shp.iterrows():
        
        # Indicate column of field ID
        IDLote = shape['IDLote']
        # Indicate column of phenology class value
        class_value  = shape[class_col_name]
        # Indicate the geometry column
        polygon = shape['geometry']
        
        if not np.isnan(class_value):
            df = get_field_dataset(dss, polygon, class_value, doa, IDLote)
            
            df['IDLote'] = IDLote
            
            df['time'] = doa
            
            try:
                df['tt'] = class_use_dict[IDLote]
            except KeyError:
                df['tt'] = 'nd'
                
            field_dataframes.append(df)
    
    #Concatenate dask dataframes for all fields 
    data = dd.concat(field_dataframes, axis=0, interleave_partitions=True)
    
    return data
    

def get_field_dataset(datasets, polygon, class_value, doa, pol_id):
    """
    Open xarray Datasets and subset data using a polygon.
    
    Args:
        datasets([str]): xr.Dataset objects for the region of the polygons
        polygon (shapely.geometry.polygon.Polygon): field polygon
        class_value (int): field class value
        doa (np.datetime64): date of analysis
        pol_id (str): polygon ID
    Returns:
        df (dask.DataFrame): a dask dataframe of all variables in datasets
    """
    
    #Select to polygon bounds
    bounds = polygon.bounds
    
    #Subset the datasets using polygon bounds
    var_data_list = list(map(lambda ds: ds.sel(x=slice(bounds[0],bounds[2]),
                                               y=slice(bounds[3],bounds[1])),
                             datasets))
    
    # Reshape and rename variables
    var_data_list = list(map(lambda ds: reshape_variables(ds, doa=doa), var_data_list))
    
    # Merge all dataset variables files
    final_ds = xr.merge(var_data_list)
    
    # Assign vegetation class
    assign_class_tag(final_ds, polygon, class_value)
    
    #
    df = final_ds.to_dask_dataframe()
    
    return df

def reshape_variables(dataset, doa, var='all'):
    """
    Reshape and rename dataset variables with names relative to a date of interest
    
    Args:
        dataset(xr.Dataset): dataset with bariables to be reshaped
        doa(np.datetime64): date to be used to calculate names
        var([str] or 'all'): list with string names of variables to be included
                        if not specified all variables in dataset are considered
    
    Returns:
        result(xr.Dataset): all variables renamed to Var_#Days relative to doa
                            a letter 'n' is used to replace negative signs for
                            dates before doa
    """
    
    list_for_ints = []
    
    if var == 'all':
        var = list(dataset.data_vars.keys())
    elif type(var) is not list:
        raise ValueError('The var parameter must be a list')
        
    for variable in var:
        interm = dataset[variable]
        
        try:
            relative_time = (pd.to_datetime(interm.time.values) - doa).days.values
            
            # Add 2 lines to only keep data for dates divisible by 16
            #interm = interm. TO-DO expose the the selected interval, to pass as parameter
            nr = np.where(relative_time%16 == 0)
            
            interm = interm.isel(time=nr[0])
            
            interm['time'] = list(map(lambda x: variable+'_'+str(x).replace('-','n'),relative_time[nr[0]]))
            
            list_for_ints.append(interm.to_dataset(dim='time'))
        
        except AttributeError:
            list_for_ints.append(interm)
    
    result = xr.merge(list_for_ints)
    
    return result

def assign_class_tag(dataset, polygon, class_value):
    """
    Assign new dataset variable in points intersected by a polygon
    to the specified class_value, otherwihise assigns zero.
    
    Args:
        dataset(xr.Dataset): dataset to be modified
        polygon (shapely.geometry.polygon.Polygon): field polygon
        class_value (uint8): value to be assigned to polygon
    
    """
    
    dataset['vclass'] = (['y','x'],
                         np.zeros((dataset.y.shape[0],
                                   dataset.x.shape[0])).astype(np.uint8))
    
    for x in dataset.x.values:
        for y in dataset.y.values:
            if geometry.point.Point(x,y).within(polygon):
                dataset['vclass'].loc[dict(y=y,x=x)] = class_value