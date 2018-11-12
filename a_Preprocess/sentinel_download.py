# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:13:18 2018

@author: ASALAZAR
"""

def download_sentinel(platform, prod_type, scihub_user, scihub_pass, start_date, end_date, region=None, filename=None, down_dir=None):
    from sentinelsat.sentinel import SentinelAPI
    import os
    # change the working directory to the location of files
    if down_dir!=None:
        os.chdir(down_dir)
    print(region)
    # connect to the API
    api = SentinelAPI(scihub_user, scihub_pass, 'https://scihub.copernicus.eu/dhus')
    
    # search by polygon, time, and Hub query keywords
    if region is not None and filename is not None:
      products = api.query(region, date = (start_date, end_date), filename = filename, producttype = prod_type, platformname = platform)
    elif region is not None:
      products = api.query(region, date = (start_date, end_date), producttype = prod_type, platformname = platform)
    elif filename is not None:
      products = api.query(date = (start_date, end_date), filename = filename, producttype = prod_type, platformname = platform)
    else:
      products = api.query(date = (start_date, end_date), producttype = prod_type, platformname = platform)
    
    # download all results from the search
    print("Files will be downloaded to {}".format(os.getcwd()))
    api.download_all(products)
    