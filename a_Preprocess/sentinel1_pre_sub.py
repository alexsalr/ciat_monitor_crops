# -*- coding: utf-8 -*-
"""
Created on Wed May 16 10:10:44 2018

@author: ASALAZAR
"""
import re
import sys
import datetime
import json
from snappy import ProductIO
from snappy import GPF
from snappy import jpy
from snappy import HashMap as hashp

ConcisePM = jpy.get_type('com.bc.ceres.core.PrintWriterConciseProgressMonitor')
System = jpy.get_type('java.lang.System')

from sentinel1_pre_utils import getBandNames
from sentinel1_pre_utils import write_product

def process_date (prod_list, out_dir, orbit, area_of_int, ref_raster, polarizations):
    """products is a list of products of the same orbit (ASCENDING/DESCENDING) and date (STATE_VECTOR_TIME[0:11]). Returns product subset for a given polarizaton"""
    if len(prod_list) > 1:
        # 0 Slice assembly
        param_sliceass = hashp()
        param_sliceass.put('selectedPolarisations','VH,VV')
        product = GPF.createProduct("SliceAssembly", param_sliceass, prod_list)
    else:
        product = prod_list[0]
    
    # 1 Apply orbit file
    param = hashp()
    product = GPF.createProduct("Apply-Orbit-File", param, product)
    
    # 2 Radiometric calibration
    param = hashp()
    param.put('outputSigmaBand', True)
    param.put('sourceBands', getBandNames(product, 'Intensity_'))
    param.put('selectedPolarisations', 'VH,VV')
    param.put('outputImageScaleInDb', False)
    
    product = GPF.createProduct("Calibration", param, product)
    
    # 3 Terrain correction
    param = hashp()
    param.put('demResamplingMethod', 'NEAREST_NEIGHBOUR')
    param.put('imgResamplingMethod', 'NEAREST_NEIGHBOUR')
    param.put('applyRadiometricNormalization', True)
    param.put('demName', 'SRTM 3Sec')
    param.put('pixelSpacingInMeter', 10.0)
    param.put('sourceBands', getBandNames(product, 'Sigma0_'))
    param.put('mapProjection', 'WGS84(DD)')
    
    product = GPF.createProduct("Terrain-Correction", param, product)
    
    # 4 Subset
    if area_of_int is not None:
        param = hashp()
        param.put('geoRegion', area_of_int)
        param.put('outputImageScaleInDb', False)
        param.put('sourceBandNames', getBandNames(product, 'Sigma0_'))
        
        product = GPF.createProduct("Subset", param, product)
    
    out_name = out_dir+'S1_'+orbit+'_'+product.getName()[24:32]+'_'
    
    write_product(product, out_name, ConcisePM(System.out))
    

def main(data_dir, out_dir, orbit, area_of_int, ref_raster, polarizations, write_int, bkey, batch):
    
    # Read products
    for key, value in batch.iteritems():
        # Read the product
        value['S1GRD'] = ProductIO.readProduct(data_dir+key+'.SAFE/manifest.safe')
        print('Reading '+key)
    
    # Make a list of dates of the products
    prod_dates = list(map(lambda x: str(batch[x]['S1GRD'].getMetadataRoot().getElement('Abstracted_Metadata').getAttribute('STATE_VECTOR_TIME').getData())[0:11], batch.keys()))
    
    print('Processing {} dates: {}.'.format(len(set(prod_dates)), list(set(prod_dates))))
    
    # Declare a(nother) dictionary to put the sliced products by date (then use this dictionary to perform th rest of operations)
    dates = {}
    # Declare lists for intermediate products in the two possible polarizations
    inter_prods = []
    
    # Iterate unique dates
    for date in list(set(prod_dates)):
        
        # Put all S1 products in dates dictionary
        dates[date] = []
        for idx, prod_date in enumerate(prod_dates):
            if prod_date == date:
                dates[date].append(batch[batch.keys()[idx]]['S1GRD'])
        
        # Process each date
        process_date(dates[date], out_dir, orbit, area_of_int, ref_raster, polarizations)

if __name__ == '__main__':
    
    param_str = sys.argv[1]
    # Decode json string into dictionnary
    param = json.loads(param_str)
    
    # Retrieve parameters
    data_dir = param['data_dir']
    out_dir= param['out_dir']
    orbit = param['orbit']
    area_of_int= param['area_of_int']
    ref_raster = param['ref_raster']
    polarizations = param['polarizations']
    write_int = param['write_int']
    bkey = param['bkey']
    batch = param['batch']
    
    # Pre-process batch
    main(data_dir, out_dir, orbit, area_of_int, ref_raster, polarizations, write_int, bkey, batch)
    