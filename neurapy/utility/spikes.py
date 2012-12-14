"""Some methods for analysing spike trains."""
import pylab

def poisson_train(rate, duration):
  """Return time stamps from a poisson process at given rate and over given duration.
  Inputs:
    rate - In Hz
    duration - in ms
  Output:
    Array of time stamps in ms
  """
  p = rate/1000.0 #Probablity is rate (s-1) * bin size (s)
  N = int(duration) #One bin per ms
  r = pylab.random(N)
  idx = pylab.find(r < p)
  return idx

def window_spike_train(timestamps, window_len):
  """Break up a spike train into windows (epochs).
  Inputs:
    timestamps - timestamps of the spikes
    window_len - length of the window. Must be in same units as that of the timestamps
  Output:
    epochs - list of tuples (start_idx, end_idx)
  """
  start_idx = 0
  end_idx = 0
  window_end = window_len
  epochs = []
  fragments = []
  while start_idx < timestamps.size:
    while timestamps[end_idx] < window_end:
      end_idx += 1
      if end_idx == timestamps.size:
        break
    epochs.append([start_idx,end_idx])
    start_idx = end_idx
    window_end += window_len

  return epochs


def spikecv(timestamps, window_len):
  """Given the time stamps compute the coefficient of variation with a sliding window.
  Returns cv and rate as an array.
  Inputs:
    timestamps - the spike timestamps
    window_len - length of window to look at isi in same units as timestamps

  Outputs:
    t  - time of the center of the window
    cv
    rate - in inverse units of timestamp
  """

  isi = pylab.diff(timestamps)

  start_idx = 0
  end_idx = 0
  window_end = window_len
  epochs = []
  while start_idx < isi.size:
    while timestamps[end_idx] < window_end:
      end_idx += 1
      if end_idx == isi.size:
        break
    epochs.append([start_idx,end_idx])
    start_idx = end_idx
    window_end += window_len

  t = pylab.zeros(len(epochs))
  cv = pylab.zeros(len(epochs))
  rate = pylab.zeros(len(epochs))
  for n, epoch in enumerate(epochs):
    #CV computation
    mean = isi[epoch[0]:epoch[1]].mean()
    std = isi[epoch[0]:epoch[1]].std()
    cv[n] = std/mean
    rate[n] = 1./mean
    t[n] = window_len * (n + .5)

  return t, cv, rate

def spikefano_singletrace(timestamps, window_len, spike_count_windows=100):
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

  fano_window_len = window_len / float(spike_count_windows)

  start_idx = 0
  end_idx = 0
  window_end = window_len
  epochs = []
  while start_idx < timestamps.size:
    while timestamps[end_idx] < window_end:
      end_idx += 1
      if end_idx == timestamps.size:
        break
    epochs.append([start_idx,end_idx])
    start_idx = end_idx
    window_end += window_len

  t = pylab.zeros(len(epochs))
  ff = pylab.zeros(len(epochs))
  for n, epoch in enumerate(epochs):
    #FF computation
    these_ts = timestamps[epoch[0]:epoch[1]]

    start_idx = 0
    end_idx = 0
    window_end = n * window_len + fano_window_len
    spk_count = pylab.zeros(spike_count_windows)
    spk_windows=0
    while start_idx < these_ts.size:
      while these_ts[end_idx] < window_end:
        end_idx += 1
        if end_idx == these_ts.size:
          break
      spk_count[spk_windows] = end_idx - start_idx
      start_idx = end_idx
      window_end +=  fano_window_len
      spk_windows += 1

    mean = spk_count[:spk_windows].mean()
    std = spk_count[:spk_windows].std()
    ff[n] = std**2/mean
    t[n] = window_len * (n+.5)

  return t, ff