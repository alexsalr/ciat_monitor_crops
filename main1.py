# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 13:11:55 2018

@author: ASALAZAR
"""

import sys ; sys.path.append('a_Preprocess')
from pre_process_wrp import *

#
#Commented lines were used for moving files
#
#import os, re, zipfile, parmap
#
#def zipprod(proddir, outname):
#    print('Zipping {}'.format(proddir))
#    zipf = zipfile.ZipFile(outname, 'w', zipfile.ZIP_DEFLATED)
#    zipdir(proddir, zipf)
#    zipf.close()
#
#def zipdir(path, ziph):
#    # ziph is zipfile handle
#    for root, dirs, files in os.walk(path):
#        for file in files:
#            ziph.write(os.path.join(root, file))

if __name__ == '__main__':
    
    pre_process_region('Saldana', ['S2'], download=True,  start_date='20180501', end_date='20180727', data_server = 'SVR_DATA_2018')
    
#    cdir = '/mnt/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/Imagenes_Satelitales/Sentinel_1/Saldana_GRDH/0.Raw/'
#    
#    os.chdir(cdir)
#    
#    #files = filter(re.compile(r'^S1.....GRD.*SAFE$').search, os.listdir(cdir))
#    
#    files = ['S1A_IW_GRDH_1SDV_20151127T231327.SAFE', 'S1A_IW_GRDH_1SSV_20160309T104258.SAFE', 'S1A_IW_GRDH_1SSV_20160402T104259.SAFE', 'S1A_IW_GRDH_1SSV_20160419T231329.SAFE']
#    
#    params = list(map(lambda x: (x, os.environ['WIN_SVR_DATA']+'Saldana/S1/'+x[:-4]+'zip'), files))
#    
#    parmap.starmap(zipprod, params)