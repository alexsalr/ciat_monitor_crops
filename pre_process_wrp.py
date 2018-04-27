import os, re, csv
from sentinel1_pre import *
from sentinel2_pre import *
from landsat_pre import *

def pre_process_region(region, prods, download=False, start_date=None, end_date=None, tile=None):
    """
    region (str): Huila, 
    region 
    """
    
    area_of_int = read_aoi(region)
    ref_raster_dim = '/home/azalazar/data/spatial_ref/'+region+'.dim'
    ref_raster_img = '/home/azalazar/data/spatial_ref/'+region+'.data/ref.img'
    
    # Read 
    for prod in prods:
        
        data_dir = set_data_dir(region, prod)
        out_dir = set_out_dir(region)
        
        if prod is 'S1':
            
            if download is True:
                from sentinel_dndl import download_sentinel
                download_sentinel('Sentinel-1', 'GRD', 'asalazarr', 'tila8sude', start_date, end_date, region=area_of_int, down_dir=data_dir)                    
            
            uncompress_files(data_dir)
            
            pre_process_s1(data_dir, out_dir, area_of_int, ref_raster_dim, polarizations=['VV','VH'])
            
        elif prod is 'S2':
        
            if download is True:
                from sentinel_dndl import download_sentinel
                download_sentinel('Sentinel-2', 'S2MSI1C', 'asalazarr', 'tila8sude', start_date, end_date, region=area_of_int, down_dir=data_dir, filename=tile)
            
            uncompress_files(data_dir)
            sen2cor_L2A_batch('all', data_dir)
            pre_process_s2(data_dir, out_dir, area_of_int, ref_raster, polarizations=['VV','VH'])
            
        elif prod is 'Landsat':
            uncompress_files(data_dir)
            pre_landsat_batch(data_dir, ref_raster_img)
            pre_process
            
            

def set_data_dir(region, prod):
    data_dir = '/home/azalazar/data/'+region+'/'+prod+'/'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print ' {} dir was created'.format(data_dir)
    return data_dir

def set_out_dir(region):
    out_direc = '/home/azalazar/data/'+region+'/pre/'
    #if not out_dir.endswith('/'):
    #    out_direc += '/'
    if not os.path.exists(out_direc):
        os.makedirs(out_direc)
        print("New directory {} was created".format(out_direc))
    return out_direc

def read_aoi(area):
    #os.chdir('~/data/spatial_ref/')
    with open('/home/azalazar/data/spatial_ref/regions.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                area_wkt = row[area]
            except KeyError as err:
                print('The area of interest is not in spatial_ref regions file. {}'.format(err.args))
            break
    return area_wkt

def uncompress_files(eo_dir, unzip_dir = None):
    """
    Unzips every zipfile in the path, and stores in directory with zipfile name+.SAFE
    Args:
        eo_dir (string): string of directory where zipfiles are located
        unzip_dir (string): directory where files are to be unzipped, default relative path
                            uz_data in working directory
    """
    import zipfile, tarfile
    
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
        os.makedirs(unzip_dir)
        print unzip_dir + ' folder' + ' was created'
    
    # Make sure uncompress path ends with slash
    if unzip_direc[-1] != '/':
        unzip_dir = unzip_direc + '/'
    
    ## Loop over list of zip files
    for im_id in eo_zip_files:
        ## Unzip only if a folder with the same name does not exist
        if not os.path.exists(unzip_direc+im_id[:-3]+'SAFE'):
            print('Unzipping ' + im_id)
            zip_ref = zipfile.ZipFile(eo_dir+im_id, 'r')
            zip_ref.extractall(unzip_direc)
            zip_ref.close()
        else:
            print(im_id[:-4] + ' was already uncompressed')
    
    ## Loop over list of tar files
    for im_id in eo_tar_files:
        ## Unztar only if a folder with the same name does not exist
        if not os.path.exists(unzip_direc+im_id[:-7]):
            print('Uncompressing ' + im_id)
            tar = tarfile.open(eo_dir+im_id, 'r')
            tar.extractall(unzip_direc+im_id[:-7])
            tar.close()
        else:
            print(im_id[:-7] + ' was already uncompressed')