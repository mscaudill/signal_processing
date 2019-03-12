import numpy as np


def rms(data,place=-1,axis=0):
    """Returns the root mean square statistic for samples.
    Args:
        samples (array_like):     array containing numbers whose rms is
                                  desired
        place (int):              if given matrix, and desired 
                                  place allows you to select
                                  an array from the matrix  
        axis (int):               axis along which rms is computed
                                  (Default=0)
    Returns: RMS (float)
    """
    if place<0:
        samples=data
    else:
        samples=data[place]
    result=np.mean(samples**2,axis=axis)
    return np.sqrt(result)

def coast(data,place=-1):
    """Returns the coast statistic for a sample.
    Args:
        data (array_like):        array or matrix containing numbers whose 
                                  coast is desired
        place (int):              if given matrix, and desired 
                                  place allows you to select
                                  an array from the matrix
    Returns: coast (float)
    """
    if place<0:
        sample=data
    else:
        sample=data[place]
    suma=0
    lgth=len(sample)
    for segment,old_segment in zip(sample[1:],sample):
        suma+=abs(segment-old_segment)/lgth
    return suma

def ampcorr(data,place):
    """Returns the amplitude correlation statistic for three segments of data.
    Args:
        data (array_like):        array or matrix containing numbers whose 
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
        initer (func):                      callable returning initial stats
        window_size (int):                  window_size in samples over
                                            which statistic will be determined
        num_cores (int):                    number of processing cores used
        

    Returns: boolean array where data statistics

    """
    pass


def stat_estimator(data,statistics,window_size,threshold_mean, threshold_win
                   ):
    """Caluclates an estimated statistic by randomly sampling chunks 
        until reaching a threshold range for the mean of that statistic
        my next step will be to make this work with multiple statistics
    Args:
        data(array):                    array
        statistics(funtion):            statistics to calculate
        window_size(int):                    the data will be split into chunks 
                                        this size
        threshold_mean(int or float):   threshold value for the mean, 
                                        stops the calculations
        threshold_win(int):             how large of a sample to look at when
                                        calculating the threshold
    
    Prints the difference and the number of samples used 
        *this is only temporary while I work on it.
    Returns: returns estimated statistic (float)
    """
    #reshapes the data into an array and cuts off the uneven segment at the end
    data=data[0:len(data)-len(data)%window_size].reshape(
            int(len(data)/window_size),window_size)
    old=0
    chunktotal=0
    oldchunk=0
    #loop going through the chunks in the data in a random order
    for i,h in enumerate(np.random.choice(range(2,len(data)),len(data)-2,replace=False)):
        #calculates the running mean
        new=(old*i+statistics(data,h))/(i+1)
        old=new
        chunktotal+=h
        if(i%threshold_win==0 and i>threshold_win):
            compare=abs(oldchunk-new)
            #compares the new and old mean to the threshold value
            if (compare<threshold_mean):
                #because I am still testing it I am having it print out these values
                print("samples_taken:",i,"difference",compare)
                break
            oldchunk=new
    return new



if __name__ == '__main__':
    np.random.seed(0)
    samples = np.random.random((100,4))
    print(rms(samples))
    data=np.random.random(900000)*500
    stat_estimator(data,rms,100,.5,20)
    

#the reason I added the place is so that the the statistics would also work
#with the amplitude correlation needing three segments of data

