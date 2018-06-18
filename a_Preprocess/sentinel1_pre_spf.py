# -*- coding: utf-8 -*-
"""
Created on Wed June 15 10:10:44 2018

@author: ASALAZAR
"""
import re, sys, datetime, json
from snappy import ProductIO, GPF, jpy #, HashMap, GPF, jpy
from snappy import HashMap as hashp

ConcisePM = jpy.get_type('com.bc.ceres.core.PrintWriterConciseProgressMonitor')
System = jpy.get_type('java.lang.System')

from sentinel1_pre_utils import getBandNames
from sentinel1_pre_utils import write_product
from sentinel1_pre_utils import collocateToRef

def stacking(product_set):
    """
    Takes a list of SNAP products and returns a stacked product with all the bands named with
    the products acquisition dates.
    Args:
        product_set: a list of products to be stacked
        ref_raster (str): location of a raster product with same extent and target resolution.
    Output: returns an individual product with the bands of the other products 
    """
    # check if products contain any bands, discard when not
    prod_set = [product for product in product_set if not product.getNumBands() == 0]
    
    # define the stack parameters
    params = hashp()
    params.put('resamplingType', 'NEAREST_NEIGHBOUR')
    params.put('initialOffsetMethod', 'Product Geolocation')
    params.put('extent', 'Master')
    
    # create the stack
    print("Creating stack of {} products...".format(str(len(prod_set))))
    create_stack = GPF.createProduct('CreateStack', params, prod_set)
    return create_stack

def mtspeckle_sigma0 (stacked_prod, pol):
    """
    Applies the a multi-temporal speckle filter to the a corregistered calibrated product stack. Takes the product bands
    which name starts with 'Sigma0'.
    
    Args:
        stacked_prod (): product with all the bands to be used for the multi-temporal speckle filter operation
        pol (str): polarization to apply the speckle filter (VV or VH)
    Output:
    """
    param_specklefilter = hashp()
    pol_bands = getBandNames(stacked_prod, "Sigma0_"+pol)
    if pol_bands is None:
        raise IndexError('No bands of {} polarization available'.format(pol))
    else:
        param_specklefilter.put('sourceBandNames', pol_bands)
        param_specklefilter.put('filter', 'Lee Sigma')
        sf_product = GPF.createProduct("Multi-Temporal-Speckle-Filter", param_specklefilter, stacked_prod)
        return sf_product

def Sigma0_todB (product):
    """
    Transforms the product bands to a logaritmic scale in dB (10*log10[band]).
    
    Args:
        product: product with Sigma0 bands in linear units
    Output:
    """
    param_logdB = hashp()
    param_logdB.put('sourceBandNames', getBandNames(product, 'Sigma0_'))
    db_product = GPF.createProduct("LinearToFromdB", param_logdB, product)
    return db_product

def main(data_dir, out_dir, orbit, area_of_int, ref_raster, polarizations, write_int, bkey, batch):
    
    inter_prods = []
    
    for key, value in batch.iteritems(): # Read products
        inter_prods.append(ProductIO.readProduct(out_dir+key+'_.dim'))
        print('Reading '+key)
    
    stack = stacking(inter_prods) # Stack all products
    
    for pol in polarizations: # for each polarization
        output_name = out_dir + 'S1_SPF_' + orbit + '_' + pol + '_P' + datetime.datetime.now().strftime("%Y%m%d") + '_' + str(bkey)
        # Multitemporal speckle filter, log transform, write
        write_product(Sigma0_todB(collocateToRef(mtspeckle_sigma0(stack, pol),ref_raster)), output_name, ConcisePM(System.out))
    
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
    