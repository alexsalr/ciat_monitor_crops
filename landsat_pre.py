import snappy, os, re
from sentinelsat.sentinel import read_geojson, geojson_to_wkt
from snappy import ProductIO, HashMap, GPF, jpy

GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
HashMap = snappy.jpy.get_type('java.util.HashMap')
WKTReader = snappy.jpy.get_type('com.vividsolutions.jts.io.WKTReader')

def pre_process_landsat(data_dir, out_dir, area_of_int, ref_raster):

    # Check location for saving results
    out_direc = out_dir
    if not out_direc.endswith('/'):
        out_direc += '/'
    if not os.path.exists(out_direc):
        os.makedirs(out_direc)
        print("New directory {} was created".format(out_direc))
    
    # Get a list of S2 L2A product directory names
    landsat_products = filter(re.compile(r'^L.*[0-9]$').search, os.listdir(data_dir))
    
    for key in landsat_products:
        
        ldir = filter(re.compile('tif$').search, os.listdir(data_dir+key))
        eo_files = sorted(ldir, key=lambda x: x.split('_')[-1])
        order = map(lambda x: x.split('_')[-1], eo_files)
        
        bands = list(map(lambda x: ProductIO.readProduct(data_dir+key+'/'+x), eo_files))
        
        # Stack
        prod_set = [product for product in bands if not product.getNumBands() == 0]
        stack_params = HashMap()
        stack_params.put('resamplingType', 'NEAREST_NEIGHBOUR')
        stack_params.put('initialOffsetMethod', 'Product Geolocation')
        stack_params.put('extent', 'Master')
        stack = GPF.createProduct('CreateStack', stack_params, prod_set)
        
        # Collocate
        cparams = HashMap()
        sourceProducts = HashMap()
        sourceProducts.put("master", ProductIO.readProduct(ref_raster))
        sourceProducts.put("slave", stack)
        collocated = GPF.createProduct('Collocate', cparams, sourceProducts)
        
        # Resample all bands to 10m resolution
        resample_subset = HashMap()
        resample_subset.put('targetResolution', 10)
        print('Resampling {}'.format(key))
        resampled = GPF.createProduct('Resample', resample_subset, collocated)
        
        # Subset to area of interest
        param_subset = HashMap()
        param_subset.put('geoRegion', geojson_to_wkt(read_geojson(area_of_int)))
        param_subset.put('outputImageScaleInDb', False)
        #param_subset.put('bandNames', 'B2,B3,B4,B8,B11,B12,quality_cloud_confidence,quality_scene_classification')
        print('Subsetting {}'.format(key))
        subset = GPF.createProduct("Subset", param_subset, resampled)
        
        # Write product
        print('Writing {} with bands: {}'.format(key, ', '.join(order)))
        ProductIO.writeProduct(subset, out_dir+key, 'BEAM-DIMAP')
        