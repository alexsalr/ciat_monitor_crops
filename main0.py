import sys ; sys.path.append('a_Preprocess')
from pre_process_wrp import *
#from sentinel_download import download_sentinel

if __name__ == '__main__':
    
    #pre_process_region('Huila', ['S1','Landsat'])
    
    #pre_process_region('Ibague', ['S2'], download=True, start_date='20180401', end_date='20180430', tile='*_T18NVK_*')
    
    #pre_process_region('Huila', ['S1'], data_server = 'WIN_SVR_DATA', download=False)#True, start_date='20180401', end_date='20180430', data_server = 'WIN_SVR_DATA')
    
    pre_process_region('Saldana', ['LE07'], download=False, data_server = 'WIN_SVR_DATA', start_date='20150401', end_date='20160228')#, tile='*_T18NZL_*')
    #pre_process_region('Casanare', ['S2'], download=True, start_date='20180101', end_date='20180430', tile='*_T18NYL_*')
    
    #pre_process_region('Casanare', ['S1'], download=True, start_date='20180101', end_date='20180430')
    
    #download_sentinel('Sentinel-2', 'S2MSI1C', 'asalazarr', 'tila8sude', '20151210', '20151212', region=read_aoi('Saldana'), down_dir='/home/azalazar/data/Saldana/S2_try2/')
    #uncompress_files('/home/azalazar/data/Saldana/S2_try2/')
    #sen2cor_L2A_batch('all', '/home/azalazar/data/Saldana/S2_try2/')
    