"""Functions to read the flotilla of files produced by the Neuralynx system."""

from struct import unpack as upk, pack as pk, calcsize as csize
import logging
logger = logging.getLogger(__name__)

def read_header(fin):
  """Standard 16 kB header."""
  return fin.read(16*1024).strip('\00')

def read_csc(fin):
  """Read a continuous record file.
  We assume that the sampling frequency is fixed throughout the trace. The sampling frequency is computed from the
  average time period between packets.
  We assume that no samples are dropped during the recording (within a segment)
  """
  hdr = read_header(fin)
  fmt = '=QIII512h' #= prevents padding and alignment nonsense 512 based on the Neuralynx we have
  sz = csize(fmt)

  real_Fs_list = [] #Calculated Fs

  segments = [] #We might have several 'runs' if we start and stop the recording during the same file
  segment_t = [] #Start time of the segments
  current_segment = None

  last_packet_t_us = -1
  predicted_next_packet_t_us = -1
  while fin:
    dain = fin.read(sz)
    if len(dain) < sz:
      break
    daup = upk(fmt, dain)
    current_packet_t_us = daup[0] #Time stamp of the packet
    if last_packet_t_us > 0:
      real_Fs_list.append(n_valid_samps * 1e6/(current_packet_t_us - last_packet_t_us))

    nominal_fs = daup[2] #The sampling freq the hardware claims
    n_valid_samps = daup[3] #Number of valid samples in the packet
    if current_packet_t_us - predicted_next_packet_t_us > 1e6/nominal_fs: #Start a new packet if the gap is larger than we expect
      last_packet_t_us = -1 #We have a discontinuity in the packets. We should not use this as part of the Fs computation
      if current_segment is not None:#close out the old segment
        segments.append(current_segment)
      current_segment = daup[4:4+n_valid_samps]
      segment_t.append(current_packet_t_us)
    else:
      current_segment += daup[4:4+n_valid_samps]
      last_packet_t_us = current_packet_t_us
    predicted_next_packet_t_us = current_packet_t_us + 1e6*n_valid_samps/nominal_fs


  if current_segment is not None:#close out the old packet
    segments.append(current_segment)

  return {'header': hdr, 'segments': segments, 't': segment_t, 'Fs': sum(real_Fs_list)/len(real_Fs_list)}

def read_nev(fin):
  """Read an event file."""
  system_id = []
  time_stamp = []
  event_id = []
  ttl_code = []
  extra_bits = []
  event_string = []

  hdr = read_header(fin)
  fmt = '=hhhQhHhhh8i128s'
  sz = csize(fmt)
  while fin:
    dain = fin.read(sz)
    if len(dain) < sz:
      break
      #data.append(upk(fmt, dain))
    daup = upk(fmt, dain)
    dum1, npkt_id, dum2, qwTimeStamp, evid, nttl, dum3, dum4, dum5, = daup[:9]
    dnExtra = daup[9:17]
    evstr = daup[17].replace('\00','').strip()

    system_id.append(npkt_id)
    time_stamp.append(qwTimeStamp)
    event_id.append(evid)
    ttl_code.append(nttl)
    extra_bits.append(dnExtra)
    event_string.append(evstr)

  return {'header': hdr, 'system id': system_id, 'time stamp': time_stamp, 'event id': event_id, 'event code': ttl_code,
          'extra bits': extra_bits, 'event string': event_string}

def read_nse(fin, only_timestamps=True):
  """Read single electrode spike record.
  Inputs:
    fin - file handle
    only_timestamps - if true, only load the timestamps, ignoring the waveform and feature data

  Output: Dictionary with fields


  Notes:
    0. What is spike acquizition entity number? Ask neuralynx
    1. Removing LOAD_ATTR overhead by defining time_stamp_append = [time_stamp[n].append for n in xrange(100)] and using
       time_stamp_append[dwCellNumber](qwTimeStamp) in the loop does not seem to improve performance.
       It reduces readability so I did not use it.
    2. Disabling garbage collection did not help
    3. In general, dictionary lookups slow things down
    4. using numpy arrays, with preallocation is slower than the dumb, straightforward python list appending
  """

  hdr = read_header(fin)
  fmt = '=QII8I32h'
  sz = csize(fmt)

  max_units = 0
  time_stamp = [[] for n in xrange(100)] #if we have more than 100 units on the wire we have other problems
  saen = [[] for n in xrange(100)]
  features = [[] for n in xrange(100)]
  waveform = [[] for n in xrange(100)]

  while fin:
    dain = fin.read(sz)
    if len(dain) < sz:
      break
    daup = upk(fmt, dain)
    qwTimeStamp, dwScNumber, dwCellNumber = daup[:3]
    time_stamp[dwCellNumber].append(qwTimeStamp)
    saen[dwCellNumber].append(dwScNumber)
    if not only_timestamps:
      dnParams = daup[3:12]
      snData = daup[12:]
      features[dwCellNumber].append(dnParams)
      waveform[dwCellNumber].append(snData)

    if dwCellNumber > max_units:
      max_units = dwCellNumber #Keeps track of maximum sorted units

  dict = {'header': hdr, 'time stamp': time_stamp[:max_units+1], 'spike acquisition entity': saen[:max_units+1]}
  if not only_timestamps:
    dict['features'] = features[:max_units+1]
    dict['waveform'] = waveform[:max_units+1]
  return dict

def write_nse(fname, time_stamps, remarks=''):
  """Write out the given time stamps into a nse file."""
  with open(fname,'wb') as fout:
    fmt = '=QII8I32h'
    dwScNumber = 1
    dwCellNumber = 1
    #dnParams = [0]*8
    #snData = [0]*32
    garbage = [0]*40

    header = remarks.ljust(16*1024,'\x00')
    fout.write(header)
    for ts in time_stamps:
      fout.write(pk(fmt, ts, dwScNumber, dwCellNumber, *garbage))

def test_read_nse(fname='/Users/kghose/Research/2013/Projects/Workingmemory/Data/NeuraLynx/2012-10-31_12-46-18/SE2.nse'):
  import time
  fin = open(fname)
  t0 = time.time()
  data = read_nse(fin)
  t1 = time.time()
  print t1 - t0
  return data
