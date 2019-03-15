import numpy as np
import matplotlib.pyplot as plt

def rms(samples,axis=0):
    """Returns the root mean square statistic for samples.
    Args:
        samples (array_like):     array containing numbers whose rms is
                                  desired  
        axis (int):               axis along which rms is computed
                                  (Default=0)
    Returns: RMS (float)
    """
    result=np.mean(samples**2,axis=axis)
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

def ampcorr(data,place):
    """Returns the amplitude correlation statistic for three segments of data.
    ***this will not work in the stat_estimator funtion***
    Args:
        data (array_like):        matrix containing numbers whose 
                                  coast is desired
        place (int):              place allows you to select
                                  an array from the matrix
    Returns: coast (float)
    """
    if(place>=2):
        #uses the array and the two arrays prior in the matrix
        #in order to do the calculation
        return min(max(data[place]),max(max(data[place-1]),max(data[place-2])))
        -max(min(data[place]),min(min(data[place-1]),min(data[place-2])))


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


def generator(data,place,window_size):
    return data[place*window_size:(place+1)*window_size]

def stat_estimator(data,statistics,window_size,percent_tolerance, threshold_win
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
    #reshapes the data into an array and cuts off the uneven segment at the end
  #  data=data[0:len(data)-len(data)%window_size].reshape(
  #          int(len(data)/window_size),window_size)
    oldchunk=0
    mean_stat=old_stat=np.zeros(len(statistics))
    #loop going through the chunks in the data in a random order
    for num_samples,sample_idx in enumerate(np.random.choice(
            range(2,int(len(data)/window_size)),int(len(data)/window_size-2),replace=False)):
        #loads the data from the random sample using generator
        if(type(data)==EEG):
            sample=generator(data,sample_idx,window_size).load()
        else:
            sample=generator(data,sample_idx,window_size)
        #calculates the running mean for each statistic
        for idx,stat in enumerate(statistics):
            mean_stat[idx]=(old_stat[idx]*num_samples+stat(sample))/(num_samples+1)
            old_stat[idx]=mean_stat[idx]
        if(num_samples%threshold_win==0 and num_samples>threshold_win):
            compare=abs(oldchunk-mean_stat[0])
            #compares the new and old mean to the threshold value
            if (compare<percent_tolerance/100*mean_stat[0]):
                #because I am still testing it I am having it print out these values
                print("samples_taken",num_samples,"difference",compare)
                break
            oldchunk=mean_stat[0]
    return mean_stat

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



if __name__ == '__main__':
    np.random.seed(0)
    samples = np.random.random((100,4))
    print(rms(samples))
    data=np.random.random(900000)*500
    stat_estimator(data,[rms,coast],100,.5,20)  
    
