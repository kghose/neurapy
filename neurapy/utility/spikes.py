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

def window_spike_train(timestamps, window_len, subwindow_len=None):
  """Break up a spike train into epochs and windows (within epochs).
  Inputs:
    timestamps - timestamps of the spikes
    window_len - length of the windows. Must be in same units as that of the timestamps
    subwindow_len - There can be mulitple subwindows within a window
                 * The subwindow_len must be less than the window_len
                 * If the last subwindow goes past the end of the epoch, it will be discarded
                 If subwindow_len is none, subwindow_len is made the same as window_len
  Output:
    windows - list of tuples (start_idx, end_idx)
    sub_windows - list of list of tuples (start_idx, end_idx). Each inner list constitutes an epoch
  """
  if subwindow_len is None:
    subwindow_len = window_len

  start_idx = 0
  end_idx = 0
  window_end = window_len
  subwindow_end = subwindow_len
  windows = []
  subwindows = []
  while start_idx < timestamps.size:
    this_subwindow = []
    twst = start_idx
    while timestamps[end_idx] < window_end:
      subwindow_end = min(window_end, subwindow_end)
      while timestamps[end_idx] < subwindow_end:
        end_idx += 1
        if end_idx == timestamps.size:
          break

      this_subwindow.append([twst, end_idx])
      twst = end_idx
      subwindow_end += subwindow_len

      if end_idx == timestamps.size:
        break

    windows.append([start_idx,end_idx])
    subwindows.append(this_subwindow)
    start_idx = end_idx
    subwindow_end = window_end + subwindow_len
    window_end += window_len

  return windows, subwindows


def spikecv(timestamps, window_len):
  """Given the time stamps compute the coefficient of variation with a sliding window.
  Returns cv and rate as an array.
  Inputs:
    timestamps - the spike timestamps
    window_len - length of window to look at isi (in same units as time stamps)

  Outputs:
    t  - time of the center of the window
    cv
    rate - in inverse units of timestamp
  """

  windows, subwindows = window_spike_train(timestamps, window_len=window_len)
  isi = pylab.diff(timestamps)
  if len(windows):
    windows[-1][1] -= 1 #we have one less isi sample than timestamps

  t = pylab.zeros(len(windows))
  cv = pylab.zeros(len(windows))
  rate = pylab.zeros(len(windows))
  for n, window in enumerate(windows):
    #CV computation
    mean = isi[window[0]:window[1]].mean()
    std = isi[window[0]:window[1]].std()
    cv[n] = std/mean
    rate[n] = 1./mean
    t[n] = window_len * (n + .5)

  return t, cv, rate

def spikefano(timestamps, window_len, subwindow_len):
  """Given the time stamps compute the coefficient of variation with a sliding window.
  Returns cv and rate as an array.
  Inputs:
    timestamps - the spike timestamps
    window_len - length of window to look at ff (same units as timestamps). One window gets us one ff estimate
                 The fano factor is the LS fit of fano_windows (variance,mean) points
    spike_count_windows - How many spike count windows we fit inside window_len (which is an epoch window)

  Outputs:
    t   - time array
    ff  - fano factors
  """
  windows, subwindows = window_spike_train(timestamps, window_len=window_len, subwindow_len=subwindow_len)
  t = pylab.zeros(len(windows))
  ff = pylab.zeros(len(windows))
  for n in xrange(len(windows)):
    #FF computation
    sbw = subwindows[n]
    spk_count = pylab.zeros(len(sbw))
    for m in xrange(len(sbw)):
      spk_count[m] = sbw[m][1] - sbw[m][0]

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

if __name__ == "__main__":
  #Test window_spike_train
  ts = pylab.arange(0,1,.01) #100 points spaced at .01
  windows1, subwindows1 = window_spike_train(ts, window_len=.5, subwindow_len=None)
  windows2, subwindows2 = window_spike_train(ts, window_len=.5, subwindow_len=.45)
  t, cv, rate = spikecv(ts, window_len=.1)


  tsp = poisson_train(rate=100, duration=1)
  epochs3, windows3 = window_spike_train(tsp, window_len=.5, subwindow_len=.45)
  tp, cvp, ratep = spikecv(tsp, window_len=.2)

  tsff = pylab.concatenate((ts,tsp+1))
  tff, ff = spikefano(tsff, window_len=.2, subwindow_len=.01)