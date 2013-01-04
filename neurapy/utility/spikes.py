"""Some methods for analysing spike trains."""
import pylab

def poisson_train(rate, duration):
  """Return time stamps from a simulated poisson process at given rate and over given duration. The granularity is 1ms
  Inputs:
    rate - In Hz
    duration - in s
  Output:
    Array of time stamps in s
  """
  p = rate/1000.0 #Probablity is rate (s-1) * bin size (s)
  N = int(duration*1000) #One bin per ms
  r = pylab.random(N)
  idx = pylab.find(r < p)
  return idx/1000.

def correlated_poisson_train(rate, duration, r):
  """Return time stamps from two simulated poisson processes at given rates and over given duration.

  The granularity is 1ms
  Inputs:
    rate - In Hz
    duration - in s
  Output:
    Array of time stamps in s
  """



def window_spike_train(timestamps, start_time=0, zero_times=0, end_time=None, window_len=1, subwindow_len=None):
  """Break up a spike train into windows and subwindows marching outwards from zero_time.

  * = zero
  | = window boundary
  - = subwindow

  --|----|----|*|----|----|----|----|--


  Inputs:
    Note: all time units should be in the same units as the timestamp time units

    timestamps - timestamps of the spikes
    start_time - time rel to zero_time we end our windows (needs to be <= 0).
                 If zero, means no pre windows.
                 If None, means prewindows stretch to begining of data
                 The start_time is extended to include an integer number of windows
    zero_times  - reference time. Can be a zx1 array, in which case will give us an array of windows. If scalar
                 will only give one set of windows.
    end_time   - time rel to zero_time we end our windows (needs to >= 0)
                 If zero, means no post-windows
                 If None, means post-windows stretch to end of data
                 The end_time is extended to include an integer number of windows
    window_len - Length of each window.
    subwindow_len - There can be multiple subwindows within a window
                    If subwindow_len is none, subwindow_len is made the same as window_len
  Output:
    window_edges - (w + 1 ) array of times, where w is the number of windows
    windows      - (z x w x 2) array of times.
                    z - number of epochs (size of zero_time array)
                    w - number of windows
                    0 - start idx of timestamps
                    1 - end idx of timestamps
                    If zero_time is a scalar, windows is 2D
    sub_windows  - (z x w x s x 2) array of times.
                    s - number of sub windows
                    If zero_time is a scalar, sub_windows is 3D
  """
  from math import ceil
  def march(timestamps, zero_time, start_idx, window_edges, subwindow_len, n_windows, n_sub_windows):
    """Repeatedly finding indexes from the same array that fit within the current window would be criminally
    inefficient. We make things a little less time consuming by keeping a running index.
    """
    window_edges_abs = window_edges + zero_time

    windows = pylab.ones((n_windows,2)) * -1
    subwindows = pylab.ones((n_windows,n_sub_windows,2)) * -1

    while timestamps[start_idx] <= window_edges_abs[0]:
      start_idx += 1
      if start_idx == timestamps.size:
        break

    end_idx = start_idx
    for widx in xrange(n_windows):
      if start_idx == timestamps.size:
        break
      window_start = window_edges_abs[widx]
      window_end = window_edges_abs[widx+1]
      twst = start_idx
      for swidx in xrange(n_sub_windows):
        if end_idx == timestamps.size:
          break
        subwindow_end = window_start + subwindow_len*(swidx+1)
        subwindow_end = min(window_end, subwindow_end) #This is needed so we synchronise the end of subwindows and windows
        while timestamps[end_idx] <= subwindow_end:
          end_idx += 1
          if end_idx == timestamps.size:
            break
        subwindows[widx, swidx, :] = [twst, end_idx]
        twst = end_idx

      windows[widx,:] = [start_idx,end_idx]
      start_idx = end_idx

    return window_edges, windows, subwindows, end_idx

  if subwindow_len is None:
    subwindow_len = window_len

  if pylab.isscalar(zero_times):
    zero_times = [zero_times]

  if start_time is None:
    start_time = timestamps[0] - max(zero_times)
  if end_time is None:
    end_time = timestamps[-1] - min(zero_times)

  n_pre_windows = int(ceil(-start_time / float(window_len)))
  n_post_windows = int(ceil(end_time / float(window_len)))
  n_windows = n_pre_windows + n_post_windows
  n_sub_windows = int(ceil(window_len / float(subwindow_len)))
  window_edges = (pylab.arange(n_windows+1) - n_pre_windows) * window_len

  all_windows = []
  all_subwindows = []
  start_idx = 0
  for this_zero_time in zero_times:
    window_edges, windows, subwindows, start_idx = march(timestamps, this_zero_time, start_idx, window_edges, subwindow_len, n_windows, n_sub_windows)
    all_windows.append(windows)
    all_subwindows.append(subwindows)

  all_windows = pylab.array(all_windows)
  all_subwindows = pylab.array(all_subwindows)

  return window_edges, all_windows, all_subwindows


def isi_histogram(timestamps, start_time=0, zero_times=0, end_time=None, window_len=1, range=.2, nbins=11):
  """Given the time stamps compute the isi histogram with a jumping window.
  Inputs:
    timestamps - the spike timestamps
    start_time - time rel to zero_time we end our windows (needs to be <= 0).
                 If zero, means no pre windows.
                 If None, means prewindows stretch to begining of data
                 The start_time is extended to include an integer number of windows
    zero_times  - reference time. Can be a zx1 array, in which case will give us an array of windows. If scalar
                 will only give one set of windows.
    end_time   - time rel to zero_time we end our windows (needs to >= 0)
                 If zero, means no post-windows
                 If None, means post-windows stretch to end of data
                 The end_time is extended to include an integer number of windows
    window_len - length of window to look at isi (in same units as time stamps)
    range - maximum isi
    nbins - number of bins in the histogram
  Outputs:
    t - time vector
    be - bin edges
    isihist - the histogram matrix. Normalized for each time slice

  (Can be plotted by doing pylab.pcolor(t, be, isihist.T,vmin=0,vmax=1, cmap=pylab.cm.gray_r))
  """
  window_edges, windows, subwindows = window_spike_train(timestamps, start_time, zero_times, end_time, window_len=window_len)
  isi = pylab.diff(timestamps)

  if windows.shape[1]:
    windows[:,-1,1] -= 1 #we have one less isi sample than timestamps

  bins = pylab.linspace(0, range, nbins+1) #We are doing bin edges

  isihist = pylab.zeros((windows.shape[1], nbins))
  for m in xrange(windows.shape[0]):
    for n in xrange(windows.shape[1]):
      hist,be = pylab.histogram(isi[windows[m,n,0]:windows[m,n,1]], bins, density=True)
      isihist[n,:] += hist

  t = (window_edges[1:] + window_edges[:-1])/2
  return t, be, isihist/windows.shape[0]

def spikecv(timestamps, start_time=0, zero_times=0, end_time=None, window_len=.1):
  """Given the time stamps compute the coefficient of variation with a jumping window.
  Returns cv and rate as an array.
  Inputs:
    timestamps - the spike timestamps
    window_len - length of window to look at isi (in same units as time stamps)

  Outputs:
    t  - time of the center of the window
    cv
    rate - in inverse units of timestamp
  """

  window_edges, windows, subwindows = window_spike_train(timestamps, start_time, zero_times, end_time, window_len=window_len)
  isi = pylab.diff(timestamps)
  if windows.shape[1]:
    windows[:,-1,1] -= 1 #we have one less isi sample than timestamps

  t = pylab.zeros(windows.shape[1])
  cv = pylab.zeros(windows.shape[1])
  rate = pylab.zeros(windows.shape[1])

  for n in xrange(windows.shape[1]):
    collected_isi = pylab.array([])
    for m in xrange(windows.shape[0]):
      #CV computation
      collected_isi = pylab.concatenate((collected_isi, isi[windows[m,n,0]:windows[m,n,1]]))

    mean = collected_isi.mean()
    std = collected_isi.std()
    cv[n] = std/mean
    rate[n] = 1./mean
    t[n] = window_len * (n + .5)

  return t, cv, rate

def spikefano(timestamps, start_time=0, zero_times=0, end_time=None, window_len=.1, subwindow_len=None):
  """Given the time stamps compute the fano factor with a jumping window.
  Inputs:
    timestamps - the spike timestamps
    window_len - length of window to look at ff (same units as timestamps). One window gets us one ff estimate
                 The fano factor is the LS fit of fano_windows (variance,mean) points
    subwindow_len - length of one spike count computation window

  Outputs:
    t   - time array
    ff  - fano factors
  """
  window_edges, windows, subwindows = window_spike_train(timestamps, start_time, zero_times, end_time, window_len=window_len, subwindow_len=subwindow_len)
  t = pylab.zeros(windows.shape[1])
  ff = pylab.zeros(windows.shape[1])
  for n in xrange(windows.shape[1]):
    spk_count = pylab.zeros(subwindows.shape[1]*subwindows.shape[2])
    for m in xrange(subwindows.shape[0]):#Epochs
      for l in xrange(subwindows.shape[2]):#Subwindows
        #FF computation
        sbw0 = subwindows[m,n,l,0]
        sbw1 = subwindows[m,n,l,1]
        spk_count[m*subwindows.shape[2]+l] = sbw1 - sbw0

    mean = spk_count.mean()
    std = spk_count.std()
    ff[n] = std**2/mean
    t[n] = window_len * (n+.5)

  return t, ff

def spike_triggered_histogram(tsA, tsB, window_len, range, nbins):
  """Given two spike trains compute their spike triggered spike histogram.
  Inputs:
    tsA, tsB - the two spike train timestamps. tsA is the reference train
    window_len - length of window over which to compute the histogram
    range - the time range over which to compute the histogram
    nbins - how many bins to make the histogram
  Outputs:
    t - time vector
    be - bin edges
    spkhist - the histogram matrix. Normalized for each time slice

  (Can be plotted by doing pylab.pcolor(t, be, spkhist.T,vmin=0,vmax=1, cmap=pylab.cm.gray_r))
  """


  #t = pylab.arange(300)*60
  #be = pylab.arange(-.1,.1,.01)
  #return t, be, pylab.rand(t.size-1, be.size-1)

  bins = pylab.linspace(range[0], range[1], nbins+1) #We are doing bin edges
  epochs = window_spike_train(tsA, window_len)

  idxBlo = 0 #Just a device to help us speed up operations
  idxBhi = 0
  t = pylab.zeros(len(epochs))
  spkhist = pylab.zeros((len(epochs), nbins))
  for n, epoch in enumerate(epochs):
    #Cycle through timestamps in this window
    thishist = pylab.zeros(nbins)
    these_ts = tsA[epoch[0]:epoch[1]]
    for j in xrange(these_ts.size):
      this_ts = these_ts[j]
      while (tsB[idxBlo] - this_ts > range[0]) & (idxBlo > 0):
        idxBlo -= 1
      while (tsB[idxBlo] - this_ts < range[0]) & (idxBlo < tsB.size - 1):
        idxBlo += 1
      while (tsB[idxBhi] - this_ts > range[1]) & (idxBhi > 0):
        idxBhi -= 1
      while (tsB[idxBhi] - this_ts < range[1]) & (idxBhi < tsB.size - 1):
        idxBhi += 1
      h,be = pylab.histogram(tsB[idxBlo:idxBhi+1] - this_ts, bins)
      thishist += h
    thishist /= thishist.max()
    spkhist[n,:] = thishist
    t[n] = window_len * (n + .5)

  return t, be, spkhist

def spikecount_correlation(tsA, tsB, start_time=0, zero_times=0, end_time=None, window_len=.1, subwindow_len=None):
  """Given the time stamps compute the spike count correlation with a jumping window.
  Inputs:
    tsA, tsB - the two spike train timestamps.
    start_time - time rel to zero_time we end our windows (needs to be <= 0).
                 If zero, means no pre windows.
                 If None, means prewindows stretch to begining of data
                 The start_time is extended to include an integer number of windows
    zero_times  - reference time. Can be a zx1 array, in which case will give us an array of windows. If scalar
                 will only give one set of windows.
    end_time   - time rel to zero_time we end our windows (needs to >= 0)
                 If zero, means no post-windows
                 If None, means post-windows stretch to end of data
                 The end_time is extended to include an integer number of windows
    window_len - length of window to look at isi (in same units as time stamps)
    subwindow_len - length of one spike count computation window

  Outputs:
    t  - time array
    r  -  the correlation vector
  """
  def correl_single_window(sbwA, sbwB):
    """sbwA and sbwB are obtained as subwindowsA[:,n,:,:], with n marching through all the windows.
    Indexes - 0 -> epochs
              1 -> subwindows
              2 -> start/stop

    """
    N0 = sbwA.shape[0]
    N1 = sbwA.shape[1]
    spk_countA = pylab.zeros((N0,N1))
    spk_countB = pylab.zeros((N0,N1))
    for n in xrange(N0): #Epochs
      for m in xrange(N1): #subwindows
        spk_countA[n,m] = sbwA[n,m,1] - sbwA[n,m,0]
        spk_countB[n,m] = sbwB[n,m,1] - sbwB[n,m,0]



    spk_countA = spk_countA.flatten() - spk_countA.mean()
    spk_countB = spk_countB.flatten() - spk_countB.mean()
    R = pylab.corrcoef(spk_countA, spk_countB)
    if R.size:
      r = R[0,1]
    else:
      r = 0
    return r

  window_edgesA, windowsA, subwindowsA =\
    window_spike_train(tsA, start_time=start_time, zero_times=zero_times, end_time=end_time, window_len=window_len, subwindow_len=subwindow_len)
  window_edgesB, windowsB, subwindowsB =\
    window_spike_train(tsB, start_time=start_time, zero_times=zero_times, end_time=end_time, window_len=window_len, subwindow_len=subwindow_len)

  #The only time we will have different numbers of windows is when zero_times is scalar, end_time is None and the two
  #trains are of different lengths
  wmin = min(windowsA.shape[1],windowsB.shape[1])
  #print wmin
  r = pylab.zeros(wmin)
  for n in xrange(wmin):
    #computation of r
    sbwA = subwindowsA[:,n,:,:]
    sbwB = subwindowsB[:,n,:,:]
    r[n] = correl_single_window(sbwA, sbwB)

  t = (window_edgesA[1:] + window_edgesA[:-1])/2
  return t, r



if __name__ == "__main__":
  zeros = pylab.arange(2, 40, 4)
  ts = pylab.array([])
  #Test window_spike_train, repeated several times
  #1 second silence
  #1 second 100Hz fixed
  #1 second silence
  #1 second 100Hz poisson
  for n in xrange(zeros.size):
    tsf = pylab.arange(0,1,.01) + zeros[n] - 1 #100 points spaced at .01 with 500ms of silence at the start
    tsp = poisson_train(rate=100, duration=1) + zeros[n] + 1
    ts = pylab.concatenate((ts,tsf,tsp))

  for start_time,zero_time,end_time in zip([0,-2],[0, zeros],[None,2]):
    pylab.figure()
    t, be, isihist = isi_histogram(ts, start_time=start_time, zero_times=zero_time, end_time=end_time, window_len=.1, range=.02, nbins=21)
    pylab.subplot(4,1,1)
    pylab.pcolor(t, be, isihist.T,vmin=0,vmax=1, cmap=pylab.cm.gray_r)
    pylab.ylabel('isi')
    t, cv, rate = spikecv(ts, start_time=start_time, zero_times=zero_time, end_time=end_time, window_len=.1)
    pylab.subplot(4,1,2)
    pylab.plot(t,cv)
    pylab.ylabel('cv')
    pylab.setp(pylab.gca(), 'xlim', [0, t.max()], 'ylim', [0,2])
    pylab.subplot(4,1,3)
    pylab.plot(t,rate)
    pylab.ylabel('Hz')
    pylab.setp(pylab.gca(), 'xlim', [0, t.max()], 'ylim', [0,200])
    t, ff = spikefano(ts, start_time=start_time, zero_times=zero_time, end_time=end_time, window_len=.1, subwindow_len=.005)
    pylab.subplot(4,1,4)
    pylab.plot(t,ff)
    pylab.ylabel('FF (s)')
    pylab.setp(pylab.gca(), 'xlim', [0, t.max()], 'ylim', [0,2])


