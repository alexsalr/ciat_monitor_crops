# -*- coding: utf-8 -*-
"""
Created on Wed June 15 10:10:44 2018

@author: ASALAZAR
"""
import os, re, sys, datetime, json
from snappy import ProductIO, GPF, jpy #, HashMap, GPF, jpy
from snappy import HashMap as hashp

ConcisePM = jpy.get_type('com.bc.ceres.core.PrintWriterConciseProgressMonitor')
System = jpy.get_type('java.lang.System')

from sentinel1_pre_utils import getBandNames
from sentinel1_pre_utils import write_product
from sentinel1_pre_utils import collocateToRef

def GLCM_Textures(product):
    
    params = hashp()
    
    params.put('sourceBandNames', getBandNames(product, "Sigma0_"))
    params.put('windowSizeStr', '5x5')
    params.put('quantizerStr', 'Probabilistic Quantizer')
    params.put('quantizationLevelsStr', '16')
    params.put('displacement','4' )
    params.put('outputContrast','true')
    params.put('outputDissimilarity','true')
    params.put('outputHomogeneity','true')
    params.put('outputASM','true')
    params.put('outputEnergy','true')
    params.put('outputMean','true')
    params.put('outputVariance','true')
    params.put('outputCorrelation','true')
    
    return GPF.createProduct("GLCM", params, product)

def main(data_dir, out_dir, orbit, ref_raster, product_name):
    
    # Check if output exists
    if not os.path.isfile(out_dir+'GLCM_'+product_name):
        # Read int product
        product = ProductIO.readProduct(out_dir+product_name)
    
        # Write prods
        write_product(collocateToRef(GLCM_Textures(product), ref_raster), out_dir+'GLCM_'+product_name, ConcisePM(System.out))
    else:
        pass

if __name__ == '__main__':
    
    param_str = sys.argv[1]
    # Decode json string into dictionnary
    param = json.loads(param_str)
    
    # Retrieve parameters
    data_dir = param['data_dir']
    out_dir= param['out_dir']
    orbit = param['orbit']
    #area_of_int= param['area_of_int']
    ref_raster = param['ref_raster']
    #polarizations = param['polarizations']
    #write_int = param['write_int']
    #bkey = param['bkey']
    #batch = param['batch']
    product_name = param['product_name']
    
    # Pre-process product
    main(data_dir, out_dir, orbit, ref_raster, product_name)
    