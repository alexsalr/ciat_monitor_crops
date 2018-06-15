import re
from snappy import HashMap as hashp
from snappy import ProductIO
from snappy import GPF

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

def collocateToRef(product, ref_raster):
        if ref_raster is not None:
            cparams = hashp()
            sourceProducts = hashp()
            sourceProducts.put("master", ProductIO.readProduct(ref_raster))
            sourceProducts.put("slave", product)
            return GPF.createProduct('Collocate', cparams, sourceProducts)
        else:
            print('No reference raster was provided')
            return None

def write_product (product, out_name, pm):
    """
    Writes a GDF product in BEAM-DIMAP format (.dim). Prints informative text with product name and
    names of the bands.
    
    Args:
        product (): product to be written
        out_name (str): name/location of the output file
    """
    print('Writing {}, with bands: {}.'.format(out_name, getBandNames(product)))
    ProductIO.writeProduct(product, out_name, 'BEAM-DIMAP', pm)
