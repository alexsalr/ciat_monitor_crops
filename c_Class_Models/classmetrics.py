# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 13:17:25 2018

@author: ASALAZAR
"""


def calculate_metrics(model_performance):
    """
    Calculate remote sensing classification metrics
    
    Args:
        model_performance (h2o...): h2o model performance obj to extract confusion matrix
    
    """
    
    #Extract confusion matrix to pd dataframe
    cm = model_performance.confusion_matrix().as_data_frame()
    
    metrics = {}
    
    metrics['confmatrix'] = cm #dataframe
    
    metrics['upaccuracy'] = None #dataframe
    
    metrics['kappa_val'] = None #float64
    
    return metrics
    