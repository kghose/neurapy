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
          e.g. x['packets']['samp'] will return a 2D array, number of packets long and 512 wide (since each packet carries 512 wave points)
          similarly x['packets']['timestamp'] will return an array number of packets long
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


def read_nev(fin, parse_event_string=False):
  """Read an event file.
  Input:
    fin - file handle
    parse_event_string - If set to true then parse the eventstrings nicely. This takes extra time. Default is False
  Ouput:
    Dictionary with fields
      'header' - the file header
      'packets' - the events. This is a new pylab dtype with fields corresponding to the event packets.
        'nstx'
        'npkt_id'
        'npkt_data_size'
        'timestamp' - timestamp (us)
        'eventid'
        'nttl' - value of the TTL port
        'ncrc'
        'ndummy1'
        'ndummy2'
        'dnExtra'
        'eventstring' - The alphanumeric string NeuraLynx attaches to this event

      'eventstring' - Only is parse_event_string is set to True. This is a nicely formatted eventstring
  """
  hdr = read_header(fin)
  nev_packet = pylab.dtype([
    ('nstx', 'h'),
    ('npkt_id', 'h'),
    ('npkt_data_size', 'h'),
    ('timestamp', 'Q'),
    ('eventid', 'h'),
    ('nttl', 'H'),
    ('ncrc', 'h'),
    ('ndummy1', 'h'),
    ('ndummy2', 'h'),
    ('dnExtra', '8i'),
    ('eventstring', '128c')
  ])
  data = pylab.fromfile(fin, dtype=nev_packet, count=-1)
  logger.info('{:d} events'.format(data['timestamp'].size))
  if parse_event_string:
    logging.info('Packaging the event strings. This makes things slower.')
    # Makes things slow. Often this field is not needed
    evstring = [None]*data['timestamp'].size
    for n in xrange(data['timestamp'].size):
      str = ''.join(data['eventstring'][n])
      evstring[n] = str.replace('\00','').strip()
    return {'header': hdr, 'packets': data, 'eventstring': evstring}
  else:
    return {'header': hdr, 'packets': data}

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
