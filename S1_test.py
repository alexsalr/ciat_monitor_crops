from pre_process_wrp import *

if __name__ == '__main__':
    
    pre_process_region('Huila', ['S1','Landsat'])
    
    pre_process_region('Huila', ['S2'], download=True, start_date='20180311', end_date='20180331', tile='*_T18NVJ_*')
    
    

    
