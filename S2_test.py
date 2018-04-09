
from sentinel2_pre import *

if __name__ == '__main__':
    
    data_dir = '/home/azalazar/data/'
    out_dir = '/home/azalazar/data/pre/'
    area_of_int = '/home/azalazar/data/spatial_ref/ibagueextent.geojson'
        
    print('Pre-processing Sentinel-2 L1C data...')
    # Check and unzip files if needed
    
    unzipfiles(data_dir)
    
    # Sen2Cor processing
    sen2cor_L2A_batch('all', data_dir)
    
    #
    pre_process_s2(data_dir, out_dir, area_of_int)