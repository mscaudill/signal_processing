# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:50:57 2019

@author: Xue_Lab
"""
#I still need to update this, but I just wanted to get it initially moved over
#to its own file

def test_stat_estimator():
    #loads the data into an EEG object
    path='C:\\Users\\Xue_Lab\\Desktop\\edf\\20150105RPSP_100155.edf'
    shortedf=EEG(path,ch_nums=[0])
    #calculates the estimated statistic
    est_stat=stat_estimator(shortedf,rms,1250,0,50)
    #splits the data into a matrix with each chunk being a row
    chunked_data=shortedf.load().reshape(int(len(shortedf)/1250),1250)
    #calculates the mean of the rms values of all of the chunks
    calc_stat=np.mean(rms(np.transpose(chunked_data)))
    #I had to round the answers to make sure they went out to the same number
    #of digits
    assert round(est_stat[0],5) == round(calc_stat,5)

def random_data_test(): 
    #this also works with coast, fbr does not work with this data right now
    stat=rms
    exarray=np.random.random(50)
    est_stat=stat_estimator(exarray,stat,5,0,50)
    chunked_data=exarray.reshape(10,5)
    calc_stat=np.mean(stat(np.transpose(chunked_data)))
    assert round(est_stat,5)==round(calc_stat,5)