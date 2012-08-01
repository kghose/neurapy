"""Some useful functions for data analysis. Requires pylab (matplotlib)
and scipy and some other more common python modules."""

import pylab
import scipy.signal as ss #For filtering
import datetime #for dates for the impedance file
from neurapy.cerebus import nsx, nev

#filfilt implementation from http://www.scipy.org/Cookbook/FiltFilt
def lfilter_zi(b,a):
    #compute the zi state from the filter parameters. see [Gust96].

    #Based on:
    # [Gust96] Fredrik Gustafsson, Determining the initial states in forward-backward 
    # filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
    # Volume 44, Issue 4

    n=max(len(a),len(b))

    zin = (  pylab.eye(n-1) - pylab.hstack( (-a[1:n,pylab.newaxis],
                                 pylab.vstack((pylab.eye(n-2), pylab.zeros(n-2))))))

    zid=  b[1:n] - a[1:n]*b[0]

    zi_matrix=pylab.linalg.inv(zin)*(pylab.matrix(zid).transpose())
    zi_return=[]

    #convert the result into a regular array (not a matrix)
    for i in range(len(zi_matrix)):
      zi_return.append(float(zi_matrix[i][0]))

    return pylab.array(zi_return)

def filtfilt(b,a,x):
    #For now only accepting 1d arrays
    ntaps=max(len(a),len(b))
    edge=ntaps*3

    if x.ndim != 1:
        raise ValueError, "Filiflit is only accepting 1 dimension arrays."

    #x must be bigger than edge
    if x.size < edge:
        raise ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b)."

    if len(a) < ntaps:
        a=pylab.r_[a,zeros(len(b)-len(a))]

    if len(b) < ntaps:
        b=pylab.r_[b,zeros(len(a)-len(b))]

    zi=lfilter_zi(b,a)

    #Grow the signal to have edges for stabilizing 
    #the filter with inverted replicas of the signal
    s=pylab.r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
    #in the case of one go we only need one of the extrems 
    # both are needed for filtfilt

    (y,zf)=ss.lfilter(b,a,s,-1,zi*s[0])

    (y,zf)=ss.lfilter(b,a,pylab.flipud(y),-1,zi*y[-1])

    return pylab.flipud(y[edge-1:-edge+1])

#Band choices from http://www.scholarpedia.org/article/Spike_sorting#Step_i.29_Filtering
def filter_waveforms(waveform, 
                     Fstop_lo = 800, Fpass_lo = 1000,
                     Fpass_hi = 3000, Fstop_hi = 3500, Fs = 30000.):
  """
  waveform - m x n array. m waveforms each of n samples
  Fstop - stop band for high pass filter
  Fpass - pass band for high pass filter
  Fs - sampling frequency of spike waveform."""
  ws = [2*Fstop_lo/Fs, 2*Fstop_hi/Fs]#2* because ws is in terms of nyquist freq which is .5*Fs
  wp = [2*Fpass_lo/Fs, 2*Fpass_hi/Fs]
  b,a = ss.iirdesign(wp, ws, gpass=1, gstop=10)
  for n in range(waveform.shape[0]):
    waveform[n,:] = filtfilt(b,a,waveform[n,:])#ss.lfilter(b,a, waveform[n,:])
  
  return waveform

def discard_artifacts(waveform, threshold_mv):
  idx = pylab.find(waveform.max(axis=1) < threshold_mv)
  waveform = waveform[idx,:]
  return waveform

def align_spike_peaks(waveform, threshold_mv):
  """Align each spike by its peak and scrunch down length apropriately"""
  #waveform -= pylab.matrix(waveform[:,:10].mean(axis=1)).T*pylab.matrix(pylab.ones((1,waveform.shape[1])))
  idx = pylab.find(waveform.max(axis=1) < threshold_mv)
  waveform = waveform[idx,:]
  peak_idx = waveform.argmin(axis=1)#a row vector of the peak indices for each spike
  idx = pylab.find((peak_idx > 5) & (peak_idx < 15))
  wv = pylab.zeros((idx.size,30))
  for n in range(idx.size):
    wv[n,:] = waveform[idx[n],peak_idx[idx[n]]-6:peak_idx[idx[n]]+24]
  
  return wv

def inspect_lfp(nsx_fname, channel = 1):
  """Plot the whole lfp for the given file."""
  f_nsx = open(nsx_fname)
  nsx_basic_header = nsx.read_basic_header(f_nsx)
  
  this_lfp = nsx.read_channel(f_nsx, nsx_basic_header,
                              channel, 
                              tstart_ms = 0,
                              tdur_ms = -1)
  Fs = float(nsx_basic_header['Fs Hz'])
  t_ms = 1000*pylab.arange(this_lfp.size)/Fs
  pylab.plot(t_ms, this_lfp)

def get_spikes_in_window(cerebus_times_ms = None,
               t1_ms = -50,
               t2_ms = 150,
               nev_basic_header = None,
               nev_extended_header = None,
               frag_dir = None,
               channel = 0,
               unit = 0):
  """A reading function that reads in data from the cerebus spike files. 
  (Replaces get_spikes)
  
  Inputs:
    cerebus_times_ms : an array of cerebus times indicating the start of the 
                       stimuli
    t1_ms : start time relative to cerebus_times_ms
    t2_ms : end time relative to cerebus_times_ms
    
    nev_basic_header 
    nev_extended_header
    frag_dir - directory where fragmented nev file is (absolute)
    channel - the channel we worry about
    unit - the unit  
                 
  Outputs:
    spike_time_ms : list (length same as cerebus_times) of spike times, given
                    relative to the cerebus_time 
  """
  
  # Read in the whole spike train - warning, if you read in the spike 
  # waveforms as well, you may run outta memory
  data, dummy = nev.read_frag_unit(
                 frag_dir, nev_basic_header, nev_extended_header,
                 channel = channel, 
                 unit = unit,
                 tstart_ms = 0.0,
                 tdur_ms = -1,
                 load_waveform = False,
                 buffer_increment_size = 1000)    
  all_spike_time_ms = data['spike time ms']
    
  spike_time_ms = [None] * cerebus_times_ms.size
  for n in range(cerebus_times_ms.size):
    tstart_ms = cerebus_times_ms[n] + t1_ms #all times in ms
    tstop_ms = cerebus_times_ms[n] + t2_ms
    idx = pylab.find((all_spike_time_ms >= tstart_ms) & (all_spike_time_ms <= tstop_ms))
    if idx.size > 0:
      rel_spike_time_ms = all_spike_time_ms[idx] - cerebus_times_ms[n]
      spike_time_ms[n] = rel_spike_time_ms
  
  return spike_time_ms

def spike_psth(spike_time_ms, t1_ms = -50., t2_ms = 250., bin_ms = 1):
  """."""
  N_trials = len(spike_time_ms)
  t2_ms = pylab.ceil((t2_ms - t1_ms) / bin_ms)*bin_ms + t1_ms
  N_bins = (t2_ms - t1_ms) / bin_ms
  
  spike_count_by_trial = pylab.zeros((N_trials,N_bins),dtype=float)
  if N_trials > 0:
    all_spikes_ms = pylab.array([],dtype=float)
    for trial in range(len(spike_time_ms)):
      if spike_time_ms[trial] is None:
        continue
      idx = pylab.find((spike_time_ms[trial] >= t1_ms) & 
                       (spike_time_ms[trial] <= t2_ms))
      spike_count_by_trial[trial,:], bin_edges = \
        pylab.histogram(spike_time_ms[trial][idx], bins = N_bins, 
                        range = (t1_ms, t2_ms))
      
    spike_rate = 1000*spike_count_by_trial.mean(axis=0)/bin_ms
  else:
    spike_rate = pylab.nan

  dummy, bin_edges = \
    pylab.histogram(None, bins = N_bins, range = (t1_ms, t2_ms))
  bin_center_ms = (bin_edges[1:] + bin_edges[:-1])/2.0

  return spike_rate, spike_count_by_trial, bin_center_ms

def get_lfp_in_window(cerebus_times_ms = None,
              t1_ms = -50,
              t2_ms = 150,
              p_thresh = +400,
              n_thresh = -400,
              hp_hz = None,
              lp_hz = None,
              nsx_basic_header = None,
              f_nsx = None,
              channel = 0):
  """A reading function that reads in data from the cerebus ns3 files to grab the
  lfp. 
  
  Inputs:
    cerebus_times_ms : an array of cerebus times indicating the start of the 
                      stimuli
    t1_ms : start time relative to cerebus_times_ms
    t2_ms : end time relative to cerebus_times_ms 

    p_thresh : If waveform exceeds this positive threshold, discard
    n_thresh : If waveform exceeds this negative threshold, discard
    
    hp_hz : high pass filter characteristic (no filtering if none)
    lp_hz : low pass filter characteristic (no filtering if none)
    
    nsx_basic_header :
    f_nsx : file handle
    channel = 0
                               
  Outputs:
    lfps : a list of arrays the same length as cerebus_times. Traces that have
           been thrown are marked by Nones
    mean_lfp : the mean lfp
    
  """
    
  trace_count = 0
  tdur_ms = t2_ms - t1_ms
  N_lfp = nsx.length_of_lfp(nsx_basic_header, tdur_ms)
  
  lfps = []
  mean_lfp = pylab.zeros(N_lfp, dtype=float)
  for n in range(cerebus_times_ms.size):
    tstart_ms = cerebus_times_ms[n] + t1_ms #all times in ms
    this_lfp = nsx.read_channel(f_nsx, nsx_basic_header,
                                channel, 
                                tstart_ms = tstart_ms,
                                tdur_ms = tdur_ms)
    if (this_lfp.max() < p_thresh) and (this_lfp.min() > n_thresh):
    #if threshold == None or pylab.absolute(this_lfp.max()) < threshold:
      lfps.append(this_lfp)
      mean_lfp += this_lfp
      trace_count += 1
    else:
      lfps.append(None)
      
  mean_lfp /= float(trace_count)#cerebus_times_ms.size
  Fs = float(nsx_basic_header['Fs Hz'])
  t_ms = 1000*pylab.arange(N_lfp)/Fs + t1_ms
  return lfps, mean_lfp, t_ms

  
# To remove --------------------------------------------------------------------

def get_spikes(cerebus_times_ms = None,
               t1_ms = -50,
               t2_ms = 150,
               fun_args = None):
  """A reading function that reads in data from the cerebus spike files. 
  
  Inputs:
    cerebus_times_ms : an array of cerebus times indicating the start of the 
                       stimuli
    t1_ms : start time relative to cerebus_times_ms
    t2_ms : end time relative to cerebus_times_ms
    fun_args : dictionary:
          'nev_basic_header', 
          'nev_extended_header', 
          'frag_dir', - directory where fragmented nev file is (absolute)
          'channel', - the channel we worry about
          'unit' the unit  
                 
  Outputs:
    neural_response : a 7D list with ultimate cells containing a list, array 
                         or any other data structure
    comment : a string that should explain what the data is  
  """
  nev_basic_header = fun_args['nev basic header']
  nev_extended_header = fun_args['nev extended header']
  frag_dir = fun_args['frag dir']
  channel = fun_args['channel']
  unit = fun_args['unit']
  
  # Read in the whole spike train - warning, if you read in the spike 
  # waveforms as well, you may run outta memory
  data, dummy = nev.read_frag_unit(
                 frag_dir, nev_basic_header, nev_extended_header,
                 channel = channel, 
                 unit = unit,
                 tstart_ms = 0.0,
                 tdur_ms = -1,
                 load_waveform = False,
                 buffer_increment_size = 1000)    
  all_spike_time_ms = data['spike time ms']
    
  spike_counts = pylab.zeros((cerebus_times_ms.size),dtype=float)
  #mean_latency_ms = pylab.ones((cerebus_times_ms.size),dtype=float) * 1000
  spike_time_ms = [None] * cerebus_times_ms.size
  for n in range(cerebus_times_ms.size):
    tstart_ms = cerebus_times_ms[n] + t1_ms #all times in ms
    tstop_ms = cerebus_times_ms[n] + t2_ms
    idx = pylab.find((all_spike_time_ms >= tstart_ms) & (all_spike_time_ms <= tstop_ms))
    spike_counts[n] = idx.size
    if idx.size > 0:
      rel_spike_time_ms = all_spike_time_ms[idx] - cerebus_times_ms[n]
      #mean_latency_ms[n] = rel_spike_time_ms.mean()
      spike_time_ms[n] = rel_spike_time_ms
  
  data = {'trials': cerebus_times_ms.size,
          'spike counts': spike_counts, 
          #'mean latency ms': mean_latency_ms, 
          'spike times ms': spike_time_ms}
  return data

def get_lfp(cerebus_times_ms = None,
              t1_ms = -50,
              t2_ms = 150,
              fun_args = None):
  """A reading function that reads in data from the cerebus ns3 files to grab the
  lfp. 
  
  Inputs:
    cerebus_times_ms : an array of cerebus times indicating the start of the 
                      stimuli
    t1_ms : start time relative to cerebus_times_ms
    t2_ms : end time relative to cerebus_times_ms 
    fun_args : dict:
          'basic nsx header', 
          'fnsx', - nsx file handle 
          'channel', - channel we are interested in
          'pos threshold', - 
          'neg threshold' - discard traces with values outside this range 
          (if not give, defaults to +/1 400 uV) 
                 
  Outputs:
    neural_response : a 7D list with ultimate cells containing a list, array 
                         or any other data structure
    comment : a string that should explain what the data is  
  """
  basic_header = fun_args['basic nsx header']
  f_nsx = fun_args['fnsx']
  channel = fun_args['channel']
  p_thresh = fun_args.get('pos threshold', +400)
  n_thresh = fun_args.get('neg threshold', -400)
    
  trace_count = 0
  tdur_ms = t2_ms - t1_ms
  N_lfp = nsx.length_of_lfp(basic_header, tdur_ms)
  
  mean_lfp = pylab.zeros(N_lfp, dtype=float)
  for n in range(cerebus_times_ms.size):
    tstart_ms = cerebus_times_ms[n] + t1_ms #all times in ms
    this_lfp = nsx.read_channel(f_nsx, basic_header,
                                channel, 
                                tstart_ms = tstart_ms,
                                tdur_ms = tdur_ms)
    if (this_lfp.max() < p_thresh) and (this_lfp.min() > n_thresh):
    #if threshold == None or pylab.absolute(this_lfp.max()) < threshold:
      mean_lfp += this_lfp
      trace_count += 1
  
  mean_lfp /= float(trace_count)#cerebus_times_ms.size
  
  data = {'total trials':cerebus_times_ms.size, 
          'trials': trace_count, #The number of valid trials
          'mean lfp': mean_lfp,
          'Fs Hz': float(basic_header['Fs Hz'])}
    
  return data

def old_spike_psth(data, t1_ms = -250., t2_ms = 0., bin_ms = 10):
  """Uses data format returned by get_spikes"""
  spike_time_ms = data['spike times ms']
  N_trials = data['trials']
  t2_ms = pylab.ceil((t2_ms - t1_ms) / bin_ms)*bin_ms + t1_ms
  N_bins = (t2_ms - t1_ms) / bin_ms
  
  if N_trials > 0:
    all_spikes_ms = pylab.array([],dtype=float)
    for trial in range(len(spike_time_ms)):
      if spike_time_ms[trial] is None:
        continue
      idx = pylab.find((spike_time_ms[trial] >= t1_ms) & 
                       (spike_time_ms[trial] <= t2_ms))
      all_spikes_ms = \
        pylab.concatenate((all_spikes_ms, spike_time_ms[trial][idx]))
    spike_n_bin, bin_edges = \
      pylab.histogram(all_spikes_ms, bins = N_bins, 
                      range = (t1_ms, t2_ms), new = True)

    spikes_per_trial_in_bin = spike_n_bin/float(N_trials) 
    spike_rate = 1000*spikes_per_trial_in_bin/bin_ms
  else:
    spike_rate = pylab.nan
  
  bin_center_ms = (bin_edges[1:] + bin_edges[:-1])/2.0

  return spike_rate, bin_center_ms

def old_get_lfp(cerebus_times_ms = None,
              t1_ms = -50,
              t2_ms = 150,
              fun_args = None):
  """A reading function that reads in data from the cerebus ns3 files to grab the
  lfp. 
  
  Inputs:
    cerebus_times_ms : an array of cerebus times indicating the start of the 
                      stimuli
    t1_ms : start time relative to cerebus_times_ms
    t2_ms : end time relative to cerebus_times_ms 
    fun_args : dict:
          'basic nsx header', 
          'fnsx', - nsx file handle 
          'channel', - channel we are interested in
          'pos threshold', - 
          'neg threshold' - discard traces with values outside this range 
          (if not give, defaults to +/1 400 uV) 
                 
  Outputs:
    neural_response : a 7D list with ultimate cells containing a list, array 
                         or any other data structure
    comment : a string that should explain what the data is  
  """
  basic_header = fun_args['basic nsx header']
  f_nsx = fun_args['fnsx']
  channel = fun_args['channel']
  p_thresh = fun_args.get('pos threshold', +400)
  n_thresh = fun_args.get('neg threshold', -400)
    
  trace_count = 0
  tdur_ms = t2_ms - t1_ms
  N_lfp = nsx.length_of_lfp(basic_header, tdur_ms)
  
  mean_lfp = pylab.zeros(N_lfp, dtype=float)
  for n in range(cerebus_times_ms.size):
    tstart_ms = cerebus_times_ms[n] + t1_ms #all times in ms
    this_lfp = nsx.read_channel(f_nsx, basic_header,
                                channel, 
                                tstart_ms = tstart_ms,
                                tdur_ms = tdur_ms)
    if (this_lfp.max() < p_thresh) and (this_lfp.min() > n_thresh):
    #if threshold == None or pylab.absolute(this_lfp.max()) < threshold:
      mean_lfp += this_lfp
      trace_count += 1
  
  mean_lfp /= float(trace_count)#cerebus_times_ms.size
  
  data = {'total trials':cerebus_times_ms.size, 
          'trials': trace_count, #The number of valid trials
          'mean lfp': mean_lfp,
          'Fs Hz': float(basic_header['Fs Hz'])}
    
  return data

def impedance(fname = '20090430impedance.txt', electrode_count = 96):
  """
  This method can read the impedance file saved by the cerebus impedance checker
  in text format. It expects the file to be in the format
  
  * AutoImpedance measurements
  *  8\10\2010     15:28:58
  *
  * m = 5.40
  * expected = 4460.00
  * b = 34.00
  *
  * Chan    RMS/Impedance
  *******  **************
   elec78     424 kOhm
   elec88     170 kOhm
 
  Given a filename and a list of electrodes to read, this method will return 
  a dictionary with the following key value pairs:
    'date': date time object
    'impedances': list of floats of values in kOhm electrode numbering starts
                  from 0
  """
  
  file = open(fname)
  lines = file.readlines()

  #Second line is the date and time in this format
  month = int(lines[1][2:4])
  day = int(lines[1][5:7])
  year = int(lines[1][8:13])
  hour = int(lines[1][17:19])
  minute = int(lines[1][20:22])
  second = int(lines[1][23:25])
  
  impedances = [0]*electrode_count
  
  #From the 10th line on each line is a cluster of three strings
  for n in range(9, len(lines)):
    values = lines[n].split()
    this_trode = int(values[0][4:])
    if this_trode <= electrode_count:
      impedances[this_trode-1] = int(values[1])
    
  data = {'date': datetime.datetime(year, month, day, hour, minute, second),
          'impedances': impedances}  
    
  return data