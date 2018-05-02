from pre_process_wrp import *

if __name__ == '__main__':
    
    #pre_process_region('Huila', ['S1','Landsat'])
    
    pre_process_region('Huila', ['S2'], download=False)#, start_date='20151201', end_date='20160115')#, tile='*_T18NVJ_*')
    
    

    
