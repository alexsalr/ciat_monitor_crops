import tarfile, os, re

def untarfiles(eo_dir, unzip_dir = None):
    """
    Unzips every zipfile in the path, and stores in directory with zipfile name+.SAFE
    Args:
        eo_dir (string): string of directory where zipfiles are located
        unzip_dir (string): directory where files are to be unzipped, default relative path
                            uz_data in working directory
    """
    # List all zip files in directory
    eo_files = filter(re.compile('tar.gz$').search, os.listdir(eo_dir))
    
    # Check if a data folder exist
    if unzip_dir is None:
        unzip_direc = eo_dir
    else:
        unzip_direc = unzip_dir
    
    if not os.path.exists(unzip_direc):
        os.makedirs(unzip_dir)
        print unzip_dir + ' folder' + ' was created'
    
    # Make sure unzip direction ends with slash
    if unzip_direc[-1] != '/':
        unzip_dir = unzip_direc + '/'
    
    ## Loop over list of zip files
    for im_id in eo_files:
        ## Unzip only if a folder with the same name does not exist
        if not os.path.exists(unzip_direc+im_id[:-7]):
            print('Uncompressing ' + im_id)
            tar = tarfile.open(eo_dir+im_id, 'r')
            tar.extractall(unzip_direc+im_id[:-7])
            tar.close()
        else:
            print(im_id[:-7] + ' was already uncompressed')
    