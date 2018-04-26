import os, re, parmap


def pre_landsat_batch(data_dir, ref_raster):

    # Get a list of S2 L2A product directory names
    landsat_products = filter(re.compile(r'^L.*[0-9]$').search, os.listdir(data_dir))
    
    for key in landsat_products:
        
        ldir = filter(re.compile('tif$').search, os.listdir(data_dir+key))
        
        bands = list(map(lambda x: (data_dir+key+'/'+x, data_dir+key+'/ref_'+x.split('_')[-2]+x.split('_')[-1], ref_raster), ldir))
                
        parmap.starmap(pre_process_landsat, bands)
        
def pre_process_landsat(orig_tif, dest_tif, ref_raster):
    
    os.system("rio warp {} {} --like {}".format(orig_tif, dest_tif, ref_raster))
    