import os, re, csv, shutil, zipfile, tarfile, parmap
from sentinel1_pre import *
from sentinel2_pre import *
from landsat_pre import *

def pre_process_region(region, prods, download=False, start_date=None, end_date=None, tile=None, ntry=1, data_server = 'LOCAL_DATA', polarizations=['VV','VH']):
    """
    region (str): Huila, Casanare, Saldana, Ibague, Valledupar
    prods ([str]): S1, S2, Landsat
    download (boolean): If True tries to download Sentinel data for the specified time period
    start_date (str): Date to use for data download, format YYYYmmdd. TODO support processing of dates windows
    end_date (str): Date to use for data download, format YYYYmmdd.
    tile (str): 
    ntry (int): number of times to try processing Sentinel-2 data, as sometimes sen2cor fails
    """
    
    area_of_int = read_aoi(region, data_server)
    # Try to read or build reference rasters, depends on S2-data availability
    try:
        ref_raster_dim = read_ref_raster(region, data_server)
        ref_raster_img = read_ref_raster(region, data_server)[:-3]+'data/B1.img'
    except:
        e = sys.exc_info()
        print("Reference raster for {} could not be generated: {} {} {}".format(region, e[0], e[1], e[2]))
        
    # Read 
    for prod in prods:
        
        data_dir = set_data_dir(region, prod, data_server)
        out_dir = set_out_dir(region, data_server)
        
        if prod is 'S1':
            
            if download is True:
                from sentinel_dndl import download_sentinel
                download_sentinel('Sentinel-1', 'GRD', 'asalazarr', 'tila8sude', start_date, end_date, region=area_of_int, down_dir=data_dir)                    
            
            #uncompress_files(data_dir)
            
            # Try again to get the reference raster
            try:
                ref_raster_dim = read_ref_raster(region, data_server)
            except:
                e = sys.exc_info()
                print("Reference raster for {} could not be generated: {} {} {}. Processing complete Sentinel-1 tile.".format(region, e[0], e[1], e[2]))
                pass
            
            pre_process_s1_by_orbit(data_dir, out_dir, area_of_int, ref_raster_dim, polarizations=polarizations)
            
        elif prod is 'S2':
        
            if download is True:
                from sentinel_dndl import download_sentinel
                download_sentinel('Sentinel-2', 'S2MSI1C', 'asalazarr', 'tila8sude', start_date, end_date, region=area_of_int, down_dir=data_dir, filename=tile)
            
            uncompress_files(data_dir)
            
            # Process data ntry times
            trycounts = 0
            while trycounts < ntry:
                # Remove corrupted L2A files if exist
                if os.path.isfile(data_dir+'L2A_corrupt.txt'):
                    with open(data_dir+'L2A_corrupt.txt', 'r') as f:
                        for line in f:
                            # Remove L2A products
                            try:
                                print('Removing directory {}'.format(data_dir+line[:-2]+'.SAFE/'))
                                shutil.rmtree(data_dir+line[:-1]+'.SAFE/')
                            except:
                                e = sys.exc_info()
                                print("{} product can not be removed: {} {} {}.".format(region, e[0], e[1], e[2]))
                    # Remove L2A-corrupted text file
                    try:
                        os.remove(data_dir+'L2A_corrupt.txt')
                    except:
                        pass
                # Update counter
                trycounts += 1
                # Process sen2cor (processes all files that do not have a L2A product)
                sen2cor_L2A_batch('all', data_dir)
                # Check and subset
                pre_process_s2(data_dir, out_dir, area_of_int)
            
        elif prod in ['LE07','LC08']:
            # Uncompress files
            uncompress_files(data_dir)
            
            # Try again TODO make code not to break if not read
            ref_raster_img = read_ref_raster(region, data_server)[:-3]+'data/B1.img'
            
            # Crop and resample landsat images
            pre_landsat_batch(data_dir, out_dir, ref_raster_img)
            
            
def pre_process_s1_by_orbit(data_dir, out_dir, area_of_int=None, ref_raster=None, polarizations=['VV','VH'], write_int=False):
    
    from snappy import ProductIO
    
    all_products = filter(re.compile(r'^S1.....GRD.*SAFE$').search, os.listdir(data_dir))
    
    orbits = list(map(lambda x: str(ProductIO.readProduct(data_dir+x[:-4]+'zip').getMetadataRoot().getElement('Abstracted_Metadata').getAttribute('PASS').getData()), all_products))
    
    # Move the files to new directory
    for idx, product in enumerate(all_products):
        #shutil.move(data_dir+product,check_dir(data_dir+orbits[idx]+'/'))
        shutil.move(data_dir+product[:-4]+'zip',check_dir(data_dir+orbits[idx]+'/'))
    
    # Process individually each directory
    for orbit in ['ASCENDING', 'DESCENDING']:
        data_dir_orbit = data_dir + orbit + '/'
        out_dir_orbit = check_dir(out_dir)# + orbit + '/')
        try:
            #Unzip files
            uncompress_files(data_dir_orbit)
            
            #Call pre-process function
            pre_process_s1(data_dir_orbit, out_dir_orbit, orbit=orbit, area_of_int=area_of_int, ref_raster=ref_raster, polarizations=polarizations, write_int=write_int)
            
            # Move the processed files to avoid reprocessing. TODO avoid moving when pre_process_s1 fails
            map(lambda x: shutil.move(data_dir_orbit+x, check_dir(data_dir_orbit+'processed/')), os.listdir(data_dir_orbit))
            
        except:
            e = sys.exc_info()
            print('Sentinel-1 with {} orbit could not be processed: {} {} {}'.format(orbit, e[0], e[1], e[2]))

def check_dir(direc):
    if not os.path.exists(direc):
        os.makedirs(direc)
        print 'New directory {} was created'.format(direc)
    return direc        

def set_data_dir(region, prod, data_server):
    data_dir = os.environ[data_server]+region+'/'+prod+'/'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print('New directory {} was created'.format(data_dir))
    return data_dir

def set_out_dir(region, data_server):
    out_direc = os.environ[data_server]+region+'/pre/'
    #if not out_dir.endswith('/'):
    #    out_direc += '/'
    if not os.path.exists(out_direc):
        os.makedirs(out_direc)
        print("New directory {} was created".format(out_direc))
    return out_direc

def read_aoi(region, data_server):
    #os.chdir('~/data/spatial_ref/')
    with open(os.environ[data_server]+'spatial_ref/regions.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                area_wkt = row[region]
            except KeyError as err:
                print('The area of interest is not in spatial_ref regions file. {}'.format(err.args))
            break
    return area_wkt

def read_ref_raster(region, data_server):
    """
    Checks if reference raster for region exists, if not, creates it using B1 of a Sentinel-2 L1C product.
    """
        
    # Requires snappy, currently imported in sentinel1_pre/sentinel2_pre
    loc_raster = os.environ[data_server]+'spatial_ref/'+region+'.dim'
    if os.path.isfile(loc_raster):
        return loc_raster
    else:
        # Import required packages
        import snappy
        from snappy import ProductIO, HashMap, GPF, jpy
        GPF.getDefaultInstance().getOperatorSpiRegistry().loadOperatorSpis()
        HashMap = snappy.jpy.get_type('java.util.HashMap')
        WKTReader = snappy.jpy.get_type('com.vividsolutions.jts.io.WKTReader')
        #Construct ref_raster
        
        # Get any s2 product available
        prdlist = filter(re.compile(r'^S2.*L1C.*SAFE$').search, os.listdir(set_data_dir(region, 'S2', data_server)))
        first_product = prdlist[0]
        
        # Read reference product
        ref_product = ProductIO.readProduct(set_data_dir(region, 'S2', data_server)+first_product)
        
        # Resample all bands to 10m resolution
        resample_subset = HashMap()
        resample_subset.put('targetResolution', 10)
        resampled = GPF.createProduct('Resample', resample_subset, ref_product)
        
        # Subset to area of interest
        param_subset = HashMap()
        param_subset.put('geoRegion', read_aoi(region, data_server))
        param_subset.put('outputImageScaleInDb', False)
        param_subset.put('bandNames', 'B1')
        subset = GPF.createProduct("Subset", param_subset, resampled)
        
        # Write file
        ProductIO.writeProduct(subset, loc_raster, 'BEAM-DIMAP')
        
        return loc_raster

# TODO enable pattern as parameter for uncompressing specific files
def uncompress_files(eo_dir, unzip_dir = None):
    """
    Unzips every zipfile in the path, and stores in directory with zipfile name+.SAFE
    Args:
        eo_dir (string): string of directory where zipfiles are located
        unzip_dir (string): directory where files are to be unzipped, default relative path
                            uz_data in working directory
    """
    
    
    # List all zip files in directory
    eo_zip_files = filter(re.compile('zip$').search, os.listdir(eo_dir))
    
    # List all tar files files in directory
    eo_tar_files = filter(re.compile('tar.gz$').search, os.listdir(eo_dir))
    
    # Check if a data folder exist
    if unzip_dir is None:
        unzip_direc = eo_dir
    else:
        unzip_direc = unzip_dir
    
    if not os.path.exists(unzip_direc):
        os.makedirs(unzip_direc)
        print('New directory {} was created'.format(unzip_direc))
    
    # Make sure uncompress path ends with slash
    if unzip_direc[-1] != '/':
        unzip_direc = unzip_direc + '/'
    
    # Put parameter sets in tuples
    eo_zip_files = list(map(lambda x: (unzip_direc, eo_dir, x), eo_zip_files))
    eo_tar_files = list(map(lambda x: (unzip_direc, eo_dir, x), eo_tar_files))
    
    parmap.starmap(unzip_eo, eo_zip_files)
    parmap.starmap(untar_eo, eo_tar_files)

def unzip_eo(unzip_direc, eo_dir, im_id):
    ## Unzip only if a folder with the same name does not exist
    if not os.path.exists(unzip_direc+im_id[:-3]+'SAFE'):
        print('Unzipping ' + im_id)
        zip_ref = zipfile.ZipFile(eo_dir+im_id, 'r')
        zip_ref.extractall(unzip_direc)
        zip_ref.close()
    else:
        print(im_id[:-4] + ' was already uncompressed')

def untar_eo(unzip_direc, eo_dir, im_id):
    if not os.path.exists(unzip_direc+im_id[:-7]):
        print('Uncompressing ' + im_id)
        tar = tarfile.open(eo_dir+im_id, 'r')
        tar.extractall(unzip_direc+im_id[:-7])
        tar.close()
    else:
        print(im_id[:-7] + ' was already uncompressed')
