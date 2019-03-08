import numpy as np

def rms(samples, axis=0):
    """Returns the root mean square statistic for samples.

    Args:
        samples (array_like):     array containing numbers whose rms is
                                  desired
        axis (int):               axis along which rms is computed
                                  (Default=0)

    Returns: RMS (float)
    """

    result  = np.mean(samples**2, axis=axis)
    return np.sqrt(result)

def threshold(data, statistics, thresholds, initer, window_size=None,
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


if __name__ == '__main__':
    np.random.seed(0)
    samples = np.random.random((100,4))
    print(rms(samples))

