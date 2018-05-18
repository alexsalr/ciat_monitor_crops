# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 14:56:16 2018

Functions with pre-processing workflow for Sentinel-1 data for CIAT
crop_monitoring_project using ESA SNAP api tools.

@author: ASALAZAR
"""

# Import packages
import os, shutil, re, sys, subprocess, json
from math import ceil

#from snappy import ProductIO, GPF #, HashMap, GPF, jpy
#from snappy import HashMap as hashp

#GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
#HashMap = jpy.get_type('java.util.HashMap')
#WKTReader = jpy.get_type('com.vividsolutions.jts.io.WKTReader')
#now = datetime.datetime.now()

def pre_process_s1_by_orbit(data_dir, out_dir, area_of_int=None, ref_raster=None, polarizations=['VV','VH'], write_int=False):
    #import os, shutil, re
    from snappy import ProductIO
    #polarizations_asc = ['VV','VH']
    #polarizations_desc = ['VV','VH']
    
    all_products = filter(re.compile(r'^S1.....GRD.*SAFE$').search, os.listdir(data_dir))
    
    orbits = list(map(lambda x: str(ProductIO.readProduct(data_dir+x+'/manifest.safe').getMetadataRoot().getElement('Abstracted_Metadata').getAttribute('PASS').getData()), all_products))
    
    ## TODO avoid reprocessing of already processed products, possible approach, reading dates of results?, other is moving to an already processed dir
    
    # Move the files to new directory
    for idx, product in enumerate(all_products):
        #TODO check and create directory
        shutil.move(data_dir+product,check_dir(data_dir+orbits[idx]+'/'))
        shutil.move(data_dir+product[:-4]+'zip',check_dir(data_dir+orbits[idx]+'/'))
    
    # Process individually each directory
    for orbit in ['ASCENDING', 'DESCENDING']:#list(set(orbits)):
        #if orbit == 'ASCENDING':
        #    polar_pro = polarizations_asc
        #elif orbit == 'DESCENDING':
        #    polar_pro = polarizations_desc
        # Pass new data directory
        data_dir_orbit = data_dir + orbit + '/'
        out_dir_orbit = check_dir(out_dir + orbit + '/')
        try:
            #Call pre-process function
            pre_process_s1(data_dir_orbit, out_dir_orbit, area_of_int=area_of_int, ref_raster=ref_raster, polarizations=polarizations, write_int=write_int)
            
            # Move the processed files to avoid reprocessing
            map(lambda x: shutil.move(data_dir_orbit+x, check_dir(data_dir_orbit+'processed/')), os.listdir(data_dir_orbit))
        except:
            e = sys.exc_info()
            print('Sentinel-1 with {} orbit could not be processed: {} {} {}'.format(orbit, e[0], e[1], e[2]))

def check_dir(direc):
    #import os
    if not os.path.exists(direc):
        os.makedirs(direc)
        print 'New directory {} was created'.format(direc)
    return direc

def pre_process_s1(data_dir, out_dir, area_of_int=None, ref_raster=None, polarizations=['VV','VH'], write_int=False):
    """
    Args:
        data_dir (str): The location of Sentinel-1 unzipped products (.SAFE dir)
        our_dir (str): Directory for saving the results
        area_of_int (geojson):  A geoJSON file specifying coordinates of a
                                regions within the Sentinel-1 images (WGS84)
        ref_raster (str): location of a raster product with same extent and target resolution
        polarizations ([str]): strings of polarizations to consider in processing chain
    """
    print('Pre-processing Sentinel-1 images...')
    batches = make_batches(sorted(filter(re.compile(r'^S1.....GRD.*SAFE$').search, os.listdir(data_dir))))
    
    print('Processing {} batches from {}. Results will be saved in {}'.format(str(len(batches.keys())), data_dir, out_dir))
    
    for bkey, batch in batches.iteritems():
        params_dict = dict(data_dir = data_dir, out_dir=out_dir, area_of_int=area_of_int, ref_raster = ref_raster, polarizations = polarizations, write_int = write_int, bkey = bkey, batch = batch)
        batch_json = json.dumps(params_dict)
        
        subprocess.call(['python', 'sentinel1_pre_sub.py', batch_json])
        
#    for bkey, batch in batches.iteritems():
#        
#        # List for storing location if intermediate results
#        results = {}
#        for pol in polarizations:
#            results[pol] = []
#        
#        for key, value in batch.iteritems():
#            # Read the product
#            value['GRD'] = ProductIO.readProduct(data_dir+key+'.SAFE/manifest.safe')
#            print('Reading '+key)
#            
#            # Apply orbit
#            param_orbit = hashp()#HashMap()
#            value['orbit'] = GPF.createProduct("Apply-Orbit-File", param_orbit, value['GRD'])
#            
#            # The following operations are specific for each polarization
#            for pol in polarizations:
#                # Radiometric calibration
#                param_calibration = hashp()#HashMap()
#                param_calibration.put('outputSigmaBand', True)
#                param_calibration.put('sourceBands', getBandNames(value['orbit'], 'Intensity_'+pol))
#                param_calibration.put('selectedPolarisations', pol)
#                param_calibration.put('outputImageScaleInDb', False)
#                value['calibration_'+pol] = GPF.createProduct("Calibration", param_calibration, value['orbit'])
#                
#                # Terrain correction
#                param_terraincor = hashp()#HashMap()
#                param_terraincor.put('demResamplingMethod', 'NEAREST_NEIGHBOUR')
#                param_terraincor.put('imgResamplingMethod', 'NEAREST_NEIGHBOUR')
#                param_terraincor.put('applyRadiometricNormalization', True)
#                param_terraincor.put('demName', 'SRTM 3Sec')
#                param_terraincor.put('pixelSpacingInMeter', 10.0)
#                param_terraincor.put('sourceBands', getBandNames(value['calibration_'+pol], 'Sigma0_'+pol))
#                param_terraincor.put('mapProjection', 'WGS84(DD)')
#                value['terraincor_'+pol] = GPF.createProduct("Terrain-Correction", param_terraincor, value['calibration_'+pol])
#                
#                # Subset to area of interest
#                if area_of_int is not None:
#                    param_subset = hashp()#HashMap()
#                    param_subset.put('geoRegion', area_of_int)
#                    param_subset.put('outputImageScaleInDb', False)
#                    param_subset.put('sourceBandNames', getBandNames(value['terraincor_'+pol], 'Sigma0_'+pol))
#                    value['subset_'+pol] = GPF.createProduct("Subset", param_subset, value['terraincor_'+pol])
#                
#                # define the name of the intermediate output
#                output_name = out_dir + pol + "_" + key
#                
#                # store immediate results
#                if write_int == True:
#                    results[pol].append(output_name+'.dim')
#                    # Write the results to files
#                    if area_of_int is not None:
#                        write_product(value['subset_'+pol], output_name)
#                    else:
#                        write_product(value['terraincor_'+pol], output_name)
#                    # dispose all the intermediate products
#                    value['calibration_'+pol].dispose()
#                    value['terraincor_'+pol].dispose()
#                    value['subset_'+pol].dispose()
#                else:
#                    results[pol].append(value['subset_'+pol])
#                    
#                #dispose all the intermediate products
#                if write_int == True:
#                    value['GRD'].dispose()
#                    value['orbit'].dispose()
#                    
#        ## Make stack of polarizations, apply mt speckle filter, log trasnform and write
#        for pol in polarizations:
#            # retrieve int products
#            if write_int == True:
#                polprods = []
#                for result in results[pol]:
#                    polprods.append(ProductIO.readProduct(result))
#            else:
#                polprods = results[pol]
#                
#            # stack, apply multi-temporal speckle filter and logaritmic transform
#            stack = Sigma0_todB(mtspeckle_sigma0(stacking(polprods), pol))
#            
#            if ref_raster is not None:
#                cparams = hashp()#HashMap()
#                sourceProducts = hashp()#HashMap()
#                sourceProducts.put("master", ProductIO.readProduct(ref_raster))
#                sourceProducts.put("slave", stack)
#                stack = GPF.createProduct('Collocate', cparams, sourceProducts)
#            
#            # define the name of the output
#            output_name = out_dir + 'S1_' + pol + '_dB_P' + now.strftime("%Y%m%d") + '_' + str(bkey)
#            
#            # write results
#            write_product(stack, output_name)
#    
#        results = None
#            
def make_batches(prdlist):
     # Calculate number of batches to allow no more than 10 products
     max_prods = 8
     #prdlist = filter(re.compile(r'^S1.....GRD.*SAFE$').search, os.listdir(data_dir))
     nbatch = ceil(len(prdlist)/float(max_prods))
     
     # Get the indices of products to save
     batchindex = {}
     
     # Create a dictionary to read Sentinel-1 L1 GRD products
     batches = {}
     
     counter = len(prdlist)
     for batch in range(int(nbatch)):
         if counter >max_prods:
             batchindex[batch] = range(int(batch*max_prods), max_prods+max_prods*batch)
         else:
             #print(range(int(batch*10), counter+10*batch))
             batchindex[batch] = range(int(batch*max_prods), counter+max_prods*batch)
         counter = counter - max_prods
         
     for key, value in batchindex.iteritems():
         batches[key] = {}
         for index in value:
             batches[key][prdlist[index][:-5]] = {}
             
     return batches


# def getBandNames (product, sfilter = ''):
#     """
#     Produces a string to use in the sourceBandNames parameter specification of SNAP operators.
#     Args:
#         product (SNAP product): the product to get band names of
#         sfilter (string): regular expression to filter the name of the bands
#     Output:
#         returns a string with comma-separated band names
#     """
#     band_names = product.getBandNames()
#     if sfilter != '':
#         band_names = filter(re.compile(r''+sfilter).search, band_names)
#     if len(band_names) > 0:
#         band_names = ','.join(band_names)
#     else:
#         band_names = None
#     return band_names
# 
# def stacking(product_set):
#     """
#     Takes a list of SNAP products and returns a stacked product with all the bands named with
#     the products acquisition dates.
#     Args:
#         product_set: a list of products to be stacked
#         ref_raster (str): location of a raster product with same extent and target resolution.
#     Output: returns an individual product with the bands of the other products 
#     """
#     # check if products contain any bands, discard when not
#     prod_set = [product for product in product_set if not product.getNumBands() == 0]
#             
#     # define the stack parameters
#     params = hashp()#HashMap()
#     params.put('resamplingType', 'NEAREST_NEIGHBOUR')
#     params.put('initialOffsetMethod', 'Product Geolocation')
#     params.put('extent', 'Master')
#     
#     # create the stack
#     print("Creating stack of {} products...".format(str(len(prod_set))))
#     create_stack = GPF.createProduct('CreateStack', params, prod_set)
#     return create_stack
# 
# def mtspeckle_sigma0 (stacked_prod, pol):
#     """
#     Applies the a multi-temporal speckle filter to the a corregistered calibrated product stack. Takes the product bands
#     which name starts with 'Sigma0'.
#     
#     Args:
#         stacked_prod (): product with all the bands to be used for the multi-temporal speckle filter operation
#         pol (str): polarization to apply the speckle filter (VV or VH)
#     Output:
#     """
#     param_specklefilter = hashp()#HashMap()
#     param_specklefilter.put('sourceBandNames', getBandNames(stacked_prod, "Sigma0_"+pol))
#     param_specklefilter.put('filter', 'Lee Sigma')
#     sf_product = GPF.createProduct("Multi-Temporal-Speckle-Filter", param_specklefilter, stacked_prod)
#     return sf_product
# 
# def Sigma0_todB (product):
#     """
#     Transforms the product bands to a logaritmic scale in dB (10*log10[band]).
#     
#     Args:
#         product: product with Sigma0 bands in linear units
#     Output:
#     """
#     param_logdB = hashp()#HashMap()
#     param_logdB.put('sourceBandNames', getBandNames(product))
#     db_product = GPF.createProduct("LinearToFromdB", param_logdB, product)
#     return db_product
# 
# def write_product (product, out_name):
#     """
#     Writes a GDF product in BEAM-DIMAP format (.dim). Prints informative text with product name and
#     names of the bands.
#     
#     Args:
#         product (): product to be written
#         out_name (str): name/location of the output file
#     """
#     print('Writing {}, with bands: {}.'.format(out_name, getBandNames(product)))
#     ProductIO.writeProduct(product, out_name, 'BEAM-DIMAP')#, pm = createProgressMonitor())
#
#def createProgressMonitor():
#    PWPM = jpy.get_type('com.bc.ceres.core.PrintWriterProgressMonitor')
#    JavaSystem = jpy.get_type('java.lang.System')
#    monitor = PWPM(JavaSystem.out)
#    return monitor