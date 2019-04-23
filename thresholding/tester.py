# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:50:57 2019

@author: Xue_Lab
"""
#I still need to update this, but I just wanted to get it initially moved over
#to its own file

import signal_processing.thresholding.threshold as th
import Deep_Seize.core.EEG as eeg
import numpy as np

def test_stat_estimator():
    #loads the data into an EEG object
    path='C:\\Users\\Xue_Lab\\Desktop\\edf\\20150105RPSP_100155.edf'
    stat=th.fbr
    shortedf=eeg.EEG(path,ch_nums=[0])
    #calculates the estimated statistic
    est_stat=th.stat_estimator(shortedf,stat,1250,0,300)
    #splits the data into a matrix with each chunk being a row
    chunked_data=shortedf.load().reshape(int(len(shortedf)/1250),1250)
    #calculates the mean of the rms values of all of the chunks
    calc_stat=np.mean(stat(np.transpose(chunked_data)))
    #I had to round the answers to make sure they went out to the same number
    #of digits
    assert round(est_stat[0],5) == round(calc_stat,5)
    '''
    
    all of the coasts work with eeg object correctly except for when it is
    chunked like this
    
            est_stat     calc stat
    1-false 55.24927,   12.8054
    2-false 3075.84224, 12.8054
    3-false 21,         12.2975343
    4-true 12.8054,     12.8054
    5-false 22.59973,   3014.3782
    fbr-false 1.10738,  1.14301
    
    edfdata=shortedf.load()
    edfdata
    1-12.8054,12.8054
    2-12.8054,12.8064
    3-12.8054,12.8064
    4-12.8054,12.8064
    5-12.8054.3014.3782
    fbr-1.10738,1.14301
    
    
    chunking
    eeg object and coast 1,2,or3
    
    chunking
    loaded data from eeg
    
    '''
def random_data_test(): 
    #this also works with coast, fbr does not work with this data right now
    exarray=np.random.random(50) 
    stat=th.coast5
    est_stat=th.stat_estimator(exarray,stat,5,0,50)
    chunked_data=exarray.reshape(10,5)
    calc_stat=np.mean(stat(np.transpose(chunked_data)))
    assert round(est_stat,5)==round(calc_stat,5)
    
    '''
    this works for all five coasts
    '''