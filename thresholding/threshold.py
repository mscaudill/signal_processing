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
        band3(int):             ending place for the second bin
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

def coast2(data):
    """Returns the coast statistic for a sample.
    Args:
        data (array_like):        array containing numbers whose coast is desired
    Returns: coast (float)
    """
    diff_sum=0
    prior=data[0]
    for sample in data[1:]:
        diff_sum+=abs(sample-prior)
        prior=sample
    return diff_sum/len(data)

def coast3(data,window_size=100):
    window_generator=window_maker(data,window_size)
    num=0
    window_sums=[]
    window_firsts=[]
    window_lasts=[]
    for window in window_generator:
       # window[0]
       # window[-1]
        window_sums.append(np.sum(abs(np.diff(window,axis=0)),axis=0))
        window_firsts.append(window.load()[0])
        window_lasts.append(window.load()[-1])
        num+=1
    diff_sum=np.sum(window_sums,axis=0)
    for first,last in zip(window_firsts[1:],window_lasts):
        diff_sum+=abs(first-last)
    return diff_sum/len(data)

def coast4(data):
    return np.sum(abs(np.diff(data,axis=0)),axis=0)/len(data)

#currently all four coasts work properly with varying degrees of effeciency
#coast 3 currently only works for edf objects
#in order to make coast 3 work only for non edf objects remove the .load()

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
        thresholds (seq floats): 
        stat_estimator (func):              callable returning initial stats
        window_size (int):                  window_size in samples over
                                            which statistic will be determined
        num_cores (int):                    number of processing cores used
        

    Returns: boolean array where data statistics

    """
        

def window_maker(data,window_size,windows_yielded=1,random=False):
    """
    takes a large set of data and breaks it into windows and yields them in 
    order or randomizes the order 
    
    Args:
        data(array_like):            data to iterate over
        window_size(int):            size the of the windows the data will be
                                     broken into
        windows_yielded(int):        number of windows returned
        random(boolean):             if false, yields the windows in order, if 
                                     true randomizes the order
    
    yields windows in a random order
    """
    start=windows_yielded-1
    window_count=len(data)//window_size#-group+1
    if(random==True):
        windows=np.random.choice(range(start,window_count),window_count,replace=False)
    else:
        windows=range(start,window_count)
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
    for window in window_maker(data,window_size,random=True):
        
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



def stabilize(data,window_size,stat,total_samples,random=True):
    """
    Works with array or eeg object with one channel, makes a graph of the 
    running mean for a statistic of randomly selected windows, used in order 
    to see how long it takes for one particular statistic to stabilize
    
    args:
        data(array or eeg object):     data used for ananlysis
        window_size(int):              how large the chunks will be to
                                       calculate a statistic from
        stat(funtion):                 statistic used for analysis
        total_samples:                 number of random sample to be taken 
                                        before stop
        random(boolean):               determines whether sampling is random
    
    Returns: graph x-axis is number of samples, x-axis is running mean 
    of the statistic
    """
    old_stat=0
    run_means=[]
    samples_taken=0
    for window in window_maker(data,window_size,random=random):
        mean_stat=(old_stat*samples_taken+stat(window))/(samples_taken+1)
        old_stat=mean_stat
        run_means.append(mean_stat)
        if samples_taken>total_samples:
            break
        samples_taken+=1
    plt.plot(run_means)

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


if __name__ == '__main__':
    np.random.seed(0)
    samples = np.random.random((100,4))
    #print(rms(samples))
    data=np.random.random(900000)*500
    stat_estimator(data,rms,100,.5,20)




'''
I have the window maker working properly with the stat_estimator. However,
it does not yet work with the challenge that amplitude correlation creates. I
also added an option to yeld the windows in order so that it is more versatile

While working on this, I realized that in order to test the stat estimator, we 
cannot compare the estimated statistic to the statistic for the entire data.
For example if we were to take find the maximum value for all of the chunks and
find the mean of those, that would not be the same as the maximum value for the 
entire data. The same principle applies to our attributes.

The way that I tested it instead was by using the reshape funtion to break the
data into chunks and find the mean of the rms of all of the chunks that way.
When I used a random array to test it, the rms and coast worked.
When I tested it using a short edf file, the rms worked but the coast and the 
fbr did not. I am not sure why.

I created frequency band ratio. The problem with this is that when you use the
funtion you need to be able to specify the bands. How would I be able to make
the stat_estimator do that?

I have the same question with amplutede correlation. Once I have the window
maker working properly, how would I tell the stat estimator to use it that way
when it recieves amplitude correlation.

'''
class FBR():
    
    def __init__(self,band1,band2,band3):
        self.band1=band1
        self.band2=band2
        self.band3=band3
    
    def cfbr(self,data):
        N=len(data)
        datafft=fft(data)
        datafft2=2.0/N*np.abs(datafft[0:N//2])
        return (sum(datafft2[self.band1:self.band2])/sum(datafft2[self.band2:self.band3]))


#FBR(100,200,300).cfbr(shortedf)
#Out[27]: array([0.85931694])
        
