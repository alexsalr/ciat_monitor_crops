# -*- coding: utf-8 -*-
"""
Created on Wed May 16 10:10:44 2018

@author: ASALAZAR
"""
import re, sys, datetime, json
from snappy import ProductIO, GPF #, HashMap, GPF, jpy
from snappy import HashMap as hashp

def getBandNames (product, sfilter = ''):
    """
    Produces a string to use in the sourceBandNames parameter specification of SNAP operators.
    Args:
        product (SNAP product): the product to get band names of
        sfilter (string): regular expression to filter the name of the bands
    Output:
        returns a string with comma-separated band names
    """
    band_names = product.getBandNames()
    if sfilter != '':
        band_names = filter(re.compile(r''+sfilter).search, band_names)
    if len(band_names) > 0:
        band_names = ','.join(band_names)
    else:
        band_names = None
    return band_names

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
    params = hashp()#HashMap()
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
    param_specklefilter = hashp()#HashMap()
    param_specklefilter.put('sourceBandNames', getBandNames(stacked_prod, "Sigma0_"+pol))
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
    param_logdB = hashp()#HashMap()
    param_logdB.put('sourceBandNames', getBandNames(product, 'Sigma0_'))
    db_product = GPF.createProduct("LinearToFromdB", param_logdB, product)
    return db_product

def write_product (product, out_name):
    """
    Writes a GDF product in BEAM-DIMAP format (.dim). Prints informative text with product name and
    names of the bands.
    
    Args:
        product (): product to be written
        out_name (str): name/location of the output file
    """
    print('Writing {}, with bands: {}.'.format(out_name, getBandNames(product)))
    ProductIO.writeProduct(product, out_name, 'BEAM-DIMAP')#, pm = createProgressMonitor())

def collocateToRef(product, ref_raster):
        if ref_raster is not None:
            cparams = hashp()#HashMap()
            sourceProducts = hashp()#HashMap()
            sourceProducts.put("master", ProductIO.readProduct(ref_raster))
            sourceProducts.put("slave", product)
            return GPF.createProduct('Collocate', cparams, sourceProducts)
        else:
            print('No reference raster was provided')
            return None

def process_date (prod_list, area_of_int):
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
        
    return product

def main(data_dir, out_dir, area_of_int, ref_raster, polarizations, write_int, bkey, batch):
    
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
        inter_prods.append(dict([('product', process_date(dates[date], area_of_int)), ('outname', out_dir+'S1_'+date)])) 
    
    # Write products
    if write_int == True:
        map(lambda x: write_product(x['product'], x['outname']), inter_prods)
    
    # stack, apply multi-temporal speckle filter and logaritmic transform
    # map(lambda x: ProductIO.writeProduct(Sigma0_todB(mtspeckle_sigma0(stacking(polprods), x)))
    
    stack = stacking(list(map(lambda x: x['product'], inter_prods)))
    
    ## Make stack of polarizations, apply mt speckle filter, log transform and write
    for pol in polarizations:
        output_name = out_dir + 'S1_' + pol + '_dB_P' + datetime.datetime.now().strftime("%Y%m%d") + '_' + str(bkey)
        # stack, apply multi-temporal speckle filter and logaritmic transform
        write_product(Sigma0_todB(collocateToRef(mtspeckle_sigma0(stack, pol),ref_raster)), output_name)

if __name__ == '__main__':
    
    param_str = sys.argv[1]
    # Decode json string into dictionnary
    param = json.loads(param_str)
    
    # Retrieve parameters
    data_dir = param['data_dir']
    out_dir= param['out_dir']
    area_of_int= param['area_of_int']
    ref_raster = param['ref_raster']
    polarizations = param['polarizations']
    write_int = param['write_int']
    bkey = param['bkey']
    batch = param['batch']
    
    # Pre-process batch
    main(data_dir, out_dir, area_of_int, ref_raster, polarizations, write_int, bkey, batch)
    