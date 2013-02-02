"""Functions to read the flotilla of files produced by the Neuralynx system."""

import pylab
from struct import unpack as upk, pack as pk, calcsize as csize
import logging
logger = logging.getLogger(__name__)

def read_header(fin):
  """Standard 16 kB header."""
  return fin.read(16*1024).strip('\00')

def read_csc(fin):
  """Read a continuous record file.
  Input:
    fin - file handle
  Ouput:
    Dictionary with fields
      'header' - the file header
      'packets' - the actual packets as read. This is a new pylab dtype with fields:
        'timestamp' - timestamp (us)
        'chan' - channel
        'Fs' - the sampling frequency
        'Ns' - the number of valid samples in the packet
        'samp' - the samples in the packet.
          e.g. x['packets']['samp'] will return a 2D array, packets long and 512 wide (since each packet carries 512 wave points)
          similarly x['packets']['timestamp'] will return an array packets long
      'Fs': the average frequency computed from the timestamps (can differ from the nominal frequency the device reports)
      'trace': the concatenated data from all the packets
      't0': the timestamp of the first packet.
  NOTE: while 'packets' returns the exact packets read, 'Fs' and 'trace' assume that the record has no gaps and that the
  sampling frequency has not changed during the recording
  """
  hdr = read_header(fin)
  csc_packet = pylab.dtype([
    ('timestamp', 'Q'),
    ('chan', 'I'),
    ('Fs', 'I'),
    ('Ns', 'I'),
    ('samp', '512h')
  ])

  data = pylab.fromfile(fin, dtype=csc_packet, count=-1)
  avgFs = (data['Ns'][:-1]/(pylab.diff(data['timestamp'])*1e-6)).mean() #This will be incorrect if the trace has pauses
  trace = data['samp'].ravel()
  return {'header': hdr, 'packets': data, 'Fs': avgFs, 'trace': trace, 't0': data['timestamp'][0]}

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

  logger.info('{:d} events'.format(len(time_stamp)))
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
