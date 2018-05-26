import os, re, parmap, sys

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 14:56:16 2018

Functions with pre-processing workflow for Sentinel-1 data for CIAT
crop_monitoring_project using ESA SNAP api tools.

@author: ASALAZAR
sen2cor functions based on:     
Rodrigo Almeida, Maya Situnayake, Kees Baake, Mortimer Werther, Timon Weitkamp, Arnan Araza. Wageningen University
Academic Consultancy Project for EagleSensing. Remote Sensing and GIS Integration course, Period 6, 2016-2017.
https://github.com/rodrigoalmeida94/ACT_EagleSensing.
"""

## TODO: remove hard coding of sen2cor local installation
def sen2cor_L2A (res, prod): 
    """
    Function to call sen2cor L2A_Process for obtaining L2A products
    
    res (str/num): resolution, accepts 10, 20, 60, or all 
    prod (str): location of S1 L1C product
    """
    # Hard-code location of sen2cor installation
    os.chdir(os.environ["SEN2COR_BIN"])
    # Coerce resolution to string
    res = str(res)
    # Execute L2A_Process with resolution parameter when specified
    if res == 'all':
        os.system("bash L2A_Process" + " " + prod)
    else:
        os.system("bash L2A_Process --resolution=" + res + " " + str(prod))

def sen2cor_L2A_batch (res, L1Cdir):
    """
    Function to batch processing of S1-L1C files in a directory to S1-L2A
    
    res (str/num): resolution, accepts 10, 20, 60, or all
    L1Cdir (str): location of S1 L1C products
    """
    # Put S1 L1C directory names in list
    L1C_files = filter(re.compile(r'^S2.*L1C.*SAFE$').search, os.listdir(L1Cdir))
    print("{} L1C files found in directory".format(str(len(L1C_files))))
    l1cList = []
    for L1C_file in L1C_files: # Iterate over directory names
        # Check if the file exists
        checker = r'S2._MSIL2A_' + L1C_file[11:]
        checker_list = filter(re.compile(checker).search, os.listdir(L1Cdir))
        if len(checker_list) == 0:
            # Call sen2cor function for individual product
            print("{} is set for processing".format(L1C_file))
            #sen2cor_L2A(res, L1Cdir+L1C_file)
            l1cList.append((res, L1Cdir+L1C_file))
        else:
            print("{} was already processed, removing from list".format(L1C_file))
    
    parmap.starmap(sen2cor_L2A, l1cList, pm_chunksize=12)
    
## Pre-processing of Sentinel-2 L2A products

def pre_process_s2(data_dir, out_dir, area_of_int):
    """
    
    data_dir (list): list of location of directory of L2A products (.SAFE dir)
    area_of_int (geoJSON): extent of region of interest
    
    """
    import snappy
    from sentinelsat.sentinel import SentinelAPI, read_geojson, geojson_to_wkt
    from snappy import ProductIO
    from snappy import HashMap
    from snappy import GPF
    GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
    HashMap = snappy.jpy.get_type('java.util.HashMap')
    
    #/home/azalazar/data/
    os.chdir(data_dir)
    
    # Check location for saving results
    out_direc = out_dir
    if not out_direc.endswith('/'):
        out_direc += '/'
    if not os.path.exists(out_direc):
        os.makedirs(out_direc)
        print("New directory {} was created".format(out_direc))
    
    # Get a list of S2 L2A product directory names
    prdlist = filter(re.compile(r'^S2.*L2A.*SAFE$').search, os.listdir(data_dir))
    
    # Filter products that have already been processed
    pre_files = filter(re.compile(r'^S2.*L2A.*data$').search, os.listdir(out_dir))
    checker = list(map(lambda x: x[0:-5] + r'.SAFE', pre_files))
    prdlist = [i for i in prdlist if i not in checker]
    
    # Create a dictionary to read Sentinel-2 L2A products
    product = {}
    for element in prdlist:
        product[element[:-5]] = {}
        #print(element)
    
    for key, value in product.iteritems():
        try:
            # Read the product
            #print('Reading {}'.format(key+'.SAFE/MTD_MSIL2A.xml'))
            reader = filter(re.compile(r'MTD_.*xml$').search, os.listdir(data_dir+key+'.SAFE/'))
            print('Reading {}'.format(key+'.SAFE/'+reader[0]))
            value['GRD'] = ProductIO.readProduct(data_dir+key+'.SAFE/'+reader[0])
        
            # Resample all bands to 10m resolution
            resample_subset = HashMap()
            resample_subset.put('targetResolution', 10)
            print('Resampling {}'.format(key))
            value['res10'] = GPF.createProduct('Resample', resample_subset, value['GRD'])
        
            # Subset to area of interest
            param_subset = HashMap()
            param_subset.put('geoRegion', area_of_int)
            param_subset.put('outputImageScaleInDb', False)
            param_subset.put('bandNames', 'B2,B3,B4,B8,B11,B12,quality_cloud_confidence,quality_scene_classification')
            print('Subsetting {}'.format(key))
            value['sub'] = GPF.createProduct("Subset", param_subset, value['res10'])
            
            # Write product
            print('Writing {} subset resampled to 10m'.format(key))
            ProductIO.writeProduct(value['sub'], out_dir+key, 'BEAM-DIMAP')
        
            # Dispose all the intermediate products
            value['GRD'].dispose()
            value['res10'].dispose()
            value['sub'].dispose()
            
        except:
            e = sys.exc_info()
            print("{} could not be processed: {} {} {}".format(key, e[0], e[1], e[2]))
            with open(data_dir+"L2A_corrupt.txt", "a") as efile:
                efile.write(key+'\n')

            
            