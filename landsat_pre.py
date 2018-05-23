import os, re, parmap


def pre_landsat_batch(data_dir, out_dir, ref_raster):

    # Get a list of S2 L2A product directory names
    landsat_products = filter(re.compile(r'^L.*[0-9]$').search, os.listdir(data_dir))
    
    for key in landsat_products:
        
        ldir = filter(re.compile('tif$').search, os.listdir(data_dir+key))
        
        # Create directories in dest location
        if not os.path.exists(out_dir+key):
            os.makedirs(out_dir+key) 
        
        # Put bands in list
        bands = list(map(lambda x: (data_dir+key+'/'+x, out_dir+key+'/ref_'+x.split('_')[-2]+x.split('_')[-1], ref_raster), ldir))
        
        parmap.starmap(pre_process_landsat, bands)
        
def pre_process_landsat(orig_tif, dest_tif, ref_raster):
    
    os.system("rio warp {} {} --like {}".format(orig_tif, dest_tif, ref_raster))
    