import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft


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

def coast1(sample):
    """Returns the coast statistic for a sample.This is the original one I
    created
    Args:
        data (array_like):        array containing numbers whose coast is desired
    Returns: coast (float)
    """
    suma=0
    for segment,old_segment in zip(sample[1:],sample):
        suma+=abs(segment-old_segment)
    return suma/len(sample)

def coast2(data):
    """Returns the coast statistic for a sample.This is faster than the first
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

def coast4(data):
    """Returns the coast statistic for a sample. This is the np.diff method
    This is the fastest method. It probably would not work with an EEG file 
    that is too big
    Args:
        data (array_like):        array containing numbers whose coast is desired
    Returns: coast (float)
    """
    return np.sum(abs(np.diff(data,axis=0)),axis=0)/len(data)

def coast3(data,window_size=50):
    """Returns the coast statistic for a sample. This is the window method
    Args:
        data (array_like):        array containing numbers whose coast is desired
        window_size (int):        size of windows CAUTION: if you set this to be
                                  larger than the sample or the sample does not
                                  divide evenly into the window size, 
                                  this will not work
    Returns: coast (float)
    """
    window_generator=window_maker(data,window_size)
    window_sums=[]
    window_firsts=[]
    window_lasts=[]
    for window in window_generator:
        window_sums.append(np.sum(abs(np.diff(window,axis=0)),axis=0))
        window_firsts.append(window[0])
        window_lasts.append(window[window_size-1])
    diff_sum=np.sum(window_sums,axis=0)
    for first,last in zip(window_firsts[1:],window_lasts):
        diff_sum+=abs(first-last)
    return diff_sum/len(data)



def coast5(data,window_size=50):
    """Returns the coast statistic for a sample. This is another variation of
    the window method. It is slightly slower but I think it will use less
    memory and the code is cleaner
    Args:
        data (array_like):        array containing numbers whose coast is desired
        window_size(int):         size of windows to break the data into
    Returns: coast (float)
    """
    #this makes sure that the window size is not larger than the sample and
    #sets the window size to the length of the data if it is
    if(window_size>len(data)):
        window_size=len(data)
    window_generator=window_maker(data,window_size)
    last=data[0]
    diff_sum=0
    for window in window_generator:
        diff_sum+=np.sum(abs(np.diff(window,axis=0)),axis=0) + abs(window[0]-last)
        last = window[window_size-1]
    return diff_sum/len(data)


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




if __name__ == '__main__':
    np.random.seed(0)
    samples = np.random.random((100,4))
    #print(rms(samples))
    data=np.random.random(900000)*500
    stat_estimator(data,rms,100,.5,20)




'''
One of the things that I realized with coast is that I have not tested it with
multiple channels. Coast 4 works well with everything because it only uses numpy
to perform calculations. The others do not work like numpy and we need to fix
them so they do work with multiple channels

Explanation of the graphs. Both show that as you increase the window_size the
time descreases exponentially. X axis is window_size. Y axis is time
The first one is each number on the x axis is really 50 times itself.

I have added more options for coast. Each one has its own description. 
I need to compare the memory usage for each.

this is the order from fastest to slowest with an array of data
coast4,coast5,coast3,coast2,coast1

this is the order from fastest to slowest with a small edf file
coast2,coast4,coast1,coast3,coast5

I discovered a problem with EEG object. I was not able to reference it backwards.
When I tried using window[-1] in order to retreive the last value in the window,
it raised an error. I am not sure if you had a reason for not allowing this. If
you do not it could be something to eventually fix. It was not a big deal here
because I was able to use window[window_length-1] here instead.
'''




    