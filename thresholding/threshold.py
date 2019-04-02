import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
import pytest

from core.EEG import EEG


def fbr(data,band1=100,band2=200,band3=300):
    """
    returns the frequency band ratio for the specified bins
    
    Args:
        data(array_like):       array whose data will be anaylzed
        band1(int):             starting point for the first bin
        band2(int):             ending place for the first bin and 
                                starting for the second bin
    returns the frequency band ratio
    """
    N=len(data)
    datafft=fft(data)
    datafft2=2.0/N*np.abs(datafft[0:N//2])
    return (sum(datafft2[band1:band2])/sum(datafft2[band2:band3]))


def rms(samples,axis=0):
    """Returns the root mean square statistic for samples.
    Args:
        samples (array_like):     array containing numbers whose rms is
                                  desired  
        axis (int):               axis along which rms is computed
                                  (Default=0)
    Returns: RMS (float)
    """
    result=np.mean(np.power(samples,2),axis=axis)
    return np.sqrt(result)

def coast(sample):
    """Returns the coast statistic for a sample.
    Args:
        data (array_like):        array containing numbers whose coast is desired
    Returns: coast (float)
    """
    suma=0
    lgth=len(sample)
    for segment,old_segment in zip(sample[1:],sample):
        suma+=abs(segment-old_segment)/lgth
    return suma

def ampcorr(data1,data2,data3):
    """Returns the amplitude correlation statistic for three segments of data.
    ***this will not work in the stat_estimator funtion***
    Args:
        data1(array_like:        array containing the first chunk of data
        data2(array_like):       array containing the second chunk of data
        data3(array_like):       array containing the third chunk of data
    Returns: amplitude correlation (float)
    """
    return min(max(data1),max(max(data2),max(data3)))-max(min(data1),min(min(data2),min(data3)))


def threshold(data, statistics, thresholds, stat_estimator, window_size=None,
              num_cores=None):
    """Returns windows in which data statistic exceeds a threshold.

    Args:
        data (array_like):                  array containing numbers whose 
                                            statistics are desired.
        statistics (func or seq of funcs):  callables returning statistic(s)
        thresholds (seq floats) 
        stat_estimator (func):              callable returning initial stats
        window_size (int):                  window_size in samples over
                                            which statistic will be determined
        num_cores (int):                    number of processing cores used
        

    Returns: boolean array where data statistics

    """
    pass

def window_maker(data,window_size,group=1):
    """
    takes a large set of data and breaks it into windows and randomizes the order 
    Args:
        data(array_like):            data to iterate over
        window_size(int):            size the of the windows the data will be
                                     broken into
        group(int):                  number of windows returned
    
    yields windows in a random order
    """
    window_count=len(data)//window_size#-group+1
    print(window_count)
    windows=np.random.choice(range(0,window_count),window_count,replace=False)
    #print(windows)
    for win_idx in windows:
         yield data[win_idx*window_size:(win_idx+1)*window_size]
        

def stat_estimator(data,statistic,window_size,percent_tolerance, threshold_win
                   ):
    """Caluclates an estimated statistic by randomly sampling chunks 
        until the difference the the mean for the current block of windows is
        within the tolereance percentage for the first statistic
        **works with eeg object if it only has one channel
    Args:
        data(array_like):               array
        statistics(funtion):            list of statistics to calculate
        window_size(int):               the data will be split into chunks 
                                        this size
        percent_tolerance(int or float):tolerance value for the mean using
                                        percentage, stops the calculations
        threshold_win(int):             how large of a sample to look at when
                                        calculating the threshold
    
    Prints the difference and the number of samples used 
        *this is only temporary while I work on it.
    Returns: returns estimated statistics in an array (float)
    """
    oldchunk=0
    old_stat=0
    samples_taken=0
    #loop going through the chunks in the data in a random order
    for window in window_maker(data,window_size):
        
        #calculates the running mean for each statistic
        mean_stat=(old_stat*(samples_taken)+statistic(window))/(samples_taken+1)
        old_stat=mean_stat
        if(samples_taken%threshold_win==0 and samples_taken>threshold_win):
            compare=abs(oldchunk-mean_stat)
            #compares the new and old mean to the threshold value
            if (np.all(compare<percent_tolerance/100*mean_stat)):
                #because I am still testing it I am having it print out these values
                print("samples_taken",samples_taken,"difference",compare)
                break
            oldchunk=mean_stat
        samples_taken+=1
    return mean_stat


'''
def stabilize(data,window_size,stat,samples_taken):
    """
    Works with array or eeg object with one channel, makes a graph of the 
    running mean for a statistic of randomly selected windows, used in order 
    to see how long it takes for one particular statistic to stabilize
    
    args:
        data(array or eeg object):     data used for ananlysis
        window_size(int):              how large the chunks will be to
                                       calculate a statistic from
        stat(funtion):                 statistic used for analysis
        samples_taken:                 number of random sample to be taken 
                                        before stop
    
    Returns: graph x-axis is number of samples, x-axis is running mean 
    of the statistic
    """
    old_stat=0
    run_means=[]

    for num_samples,sample_idx in enumerate(np.random.choice(
                range(2,int(len(data)/window_size)),int(len(data)/window_size-2),replace=False)):
        if(type(data)==EEG):
            sample=generator(data,sample_idx,window_size).load()
        else:
            sample=generator(data,sample_idx,window_size)
        mean_stat=(old_stat*num_samples+stat(sample))/(num_samples+1)
        old_stat=mean_stat
        run_means.append(mean_stat)
        if num_samples>samples_taken:
            break
    plt.plot(run_means)
'''
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

if __name__ == '__main__':
    np.random.seed(0)
    samples = np.random.random((100,4))
    #print(rms(samples))
    data=np.random.random(900000)*500
    stat_estimator(data,rms,100,.5,20)




'''
I have the window maker working properly with the stat_estimator. However,
it does not yet work with the challenge that amplitude correlation creates.

While working on this, I realized that in order to test the stat estimator, we 
cannot compare the estimated statistic to the statistic for the entire data.
For example if we were to take find the maximum value for all of the chunks and
find the mean of those, that would not be the same as the maximum value for the 
entire data. The same principle applies to our attributes.

The way that I tested it instead was by using the reshape funtion to break the
data into chunks and find the mean of the rms of all of the chunks that way. I
still need to test this for coast and fbr.

I created frequency band ratio. The problem with this is that when you use the
funtion you need to be able to specify the bands. How would I be able to make
the stat_estimator do that?

I have the same question with amplutede correlation. Once I have the window
maker working properly, how would I tell the stat estimator to use it that way
when it recieves amplitude correlation.


'''








