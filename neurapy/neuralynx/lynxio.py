"""Functions to read the flotilla of files produced by the Neuralynx system."""

from struct import unpack as upk, pack as pk, calcsize as csize
import logging, pylab
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

def read_nse(fin):
  """Read single electrode spike record.
  Inputs:
    fin - file handle
    only_timestamps - if true, only load the timestamps, ignoring the waveform and feature data

  Output: Dictionary with fields
    'header'  - header info
    'packets' - pylab data structure with the spike data

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
  nse_packet = pylab.dtype([
    ('timestamp', 'Q'),
    ('saen', 'I'),
    ('cellno', 'I'),
    ('Features', '8I'),
    ('waveform', '32h')
  ])
  data = pylab.fromfile(fin, dtype=nse_packet, count=-1)
  return {'header': hdr, 'packets': data}


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


def extract_nrd_ec(fname, ftsname, fttlname, fchanname, channel_list, channels=64, max_pkts=-1, buffer_size=10000, error_bugout=1000000000):
  """Read and write out selected raw traces from the .nrd file with error checking.
  Inputs:
    fname - name of nrd file
    ftsname - name under which timestamp vector will be saved
    fttlname - name under which the events will be saved
    fchanname - a list of file names for the
    channel_list - Which AD channels to convert.
    channels - total channels in the system
    max_pkts - total packets to read. If set to -1 then read all packets
    buffer_size   - how many chunks to read at a time.
    error_bugout - If the sum of stx, crc and timestamp errors exceed this value quit reading the file
  Outputs:
    Data are written to file

  e.g.
  ----------------------------------------------------------------------------------------------------------------------
  from neurapy.neuralynx import lynxio
  import logging
  logging.basicConfig(level=logging.DEBUG)

  channels = 64
  fname = '/Users/kghose/Research/2013/Projects/Workingmemory/Data/NeuraLynx/2013-01-25_14-53-04/DigitalLynxRawDataFile.nrd'
  channel_list = [0,1,2]

  ftsname = 'timestamps.raw'
  fttlname = 'ttl.raw'
  fchanname = ['chan_{:000d}.raw'.format(ch) for ch in channel_list]
  lynxio.extract_nrd(fname, ftsname, fttlname, fchanname, channel_list, channels, max_pkts=1000)
  ----------------------------------------------------------------------------------------------------------------------

  Data are written as a pure stream of binary data and can be easily and efficiently read using the numpy read function.
  For convenience, a function that reads the timestamps, events and channels (read_extracted_data) is included in the library.

  """
  def seek_packet(f):
    """Skip forward until we find the STX magic number."""
    #Read in 32bit increments until the magic number is found
    start = f.tell()
    pkt = f.read(4)
    while len(pkt) == 4:
      if pkt == '\x00\x08\x00\x00': #Magic number 2048 0x0800
        f.seek(-4,1) #Realign
        break
      pkt = f.read(4)
    stop = f.tell()
    return stop - start

  logger.info('Notice: you are using the slow version of the extractor. All error checks are done')

  #nrd packet format
  nrd_packet = pylab.dtype([
    ('stx', 'i'),
    ('pkt_id', 'i'),
    ('pkt_data_size', 'i'),
    ('timestamp high', 'I'), #Neuralynx timestamp is ... in its own 32 bit world
    ('timestamp low', 'I'),
    ('status', 'i'),
    ('ttl', 'I'),
    ('extra', '10i'),
    ('data', '{:d}i'.format(channels)),
    ('crc', 'i')
  ])
  packet_size = nrd_packet.itemsize

  pkt_cnt = 0
  garbage_bytes = 0
  stx_err_cnt = 0
  pkt_id_err_cnt = 0
  pkt_size_err_cnt = 0
  pkt_ts_err_cnt = 0
  pkt_crc_err_cnt = 0

  if max_pkts != -1: #An insidious bug was killed here!
    if buffer_size > max_pkts:
      buffer_size = max_pkts

  #The files we will write to.
  fts = open(ftsname,'wb')
  fttl = open(fttlname,'wb')
  fchan = [open(fcn,'wb') for fcn in fchanname]

  last_ts = 0L
  with open(fname,'rb') as f:
    hdr = read_header(f)
    logger.info('File header: {:s}'.format(hdr))

    garbage_bytes += seek_packet(f)
    these_packets = pylab.fromfile(f, dtype=nrd_packet, count=buffer_size)
    while these_packets.size > 0:
      all_packets_good = True
      packets_read = these_packets.size

      idx = pylab.find(these_packets['stx'] != 2048)
      if idx.size > 0:
        stx_err_cnt += 1
        all_packets_good = False
        max_good_packets = idx[0]
        these_packets = these_packets[:max_good_packets]

      if these_packets.size > 0:
        idx = pylab.find(these_packets['pkt_id'] != 1)
        if idx.size > 0:
          pkt_id_err_cnt += 1
          all_packets_good = False
          max_good_packets = idx[0]
          these_packets = these_packets[:max_good_packets]

      if these_packets.size > 0:
        idx = pylab.find(these_packets['pkt_data_size'] != 10 + channels)
        if idx.size > 0:
          pkt_size_err_cnt += 1
          all_packets_good = False
          max_good_packets = idx[0]
          these_packets = these_packets[:max_good_packets]

      if these_packets.size > 0:
        #crc computation
        field32 = pylab.vstack([these_packets[k].T for k in nrd_packet.fields.keys()]).astype('I')
        crc = pylab.zeros(these_packets.size,dtype='I')
        for idx in xrange(field32.shape[0]):
          crc ^= field32[idx,:]
        idx = pylab.find(crc != 0)
        if idx.size > 0:
          pkt_crc_err_cnt += 1
          all_packets_good = False
          max_good_packets = idx[0]
          these_packets = these_packets[:max_good_packets]

      if these_packets.size > 0:
        ts = pylab.array((these_packets['timestamp high']<<32) | (these_packets['timestamp low']), dtype='uint64')
        bad_idx = -1
        if last_ts > ts[0]:#Time stamps out of order at buffer boundary
          bad_idx = 0
        else:
          idx = pylab.find(ts[:-1] > ts[1:])
          if idx.size > 0:
            bad_idx = idx[0] + 1
        if bad_idx > -1:
          logger.info('Out of order timestamp {:d}'.format(ts[bad_idx]))
          pkt_ts_err_cnt += 1
          all_packets_good = False
          max_good_packets = bad_idx
          these_packets = these_packets[:max_good_packets]
          ts = ts[:max_good_packets]

      if these_packets.size > 0:
        last_ts = ts[-1] #Ready for the next read
        ts.tofile(fts)
        these_packets['ttl'].tofile(fttl)
        for idx,ch in enumerate(channel_list):
          these_packets['data'][:,ch].tofile(fchan[idx])

      pkt_cnt += these_packets.size
      if max_pkts != -1:
        if pkt_cnt >= max_pkts: #NOTE: This may give us upto buffer_size -1 more packets than we want.
          break

      if not all_packets_good:
        f.seek((these_packets.size-packets_read)*packet_size+4,1) #Rewind all the way except 32 bits
        garbage_bytes += seek_packet(f)

      if pkt_ts_err_cnt + pkt_crc_err_cnt + stx_err_cnt > error_bugout:
        logger.warning('Too many errors, bugging out')
        break

      these_packets = pylab.fromfile(f, dtype=nrd_packet, count=buffer_size)

  fts.close()
  fttl.close()
  [fch.close() for fch in fchan]

  logger.info('Extracted {:d} packets'.format(pkt_cnt))
  logger.info('{:d} garbage words'.format(garbage_bytes))
  logger.info('{:d} packets had bad stx'.format(stx_err_cnt))
  logger.info('{:d} packets had bad pkt id'.format(pkt_id_err_cnt))
  logger.info('{:d} packets had bad crc'.format(pkt_crc_err_cnt))
  logger.info('{:d} packets had out of order timestamps'.format(pkt_ts_err_cnt))




def extract_nrd_fast(fname, ftsname, fttlname, fchanname, channel_list, channels=64, max_pkts=-1, buffer_size=10000):
  """Read and write out selected raw traces from the .nrd file.
  Inputs:
    fname - name of nrd file
    ftsname - name under which timestamp vector will be saved
    fttlname - name under which the events will be saved
    fchanname - a list of file names for the
    channel_list - Which AD channels to convert.
    channels - total channels in the system
    max_pkts - total packets to read. If set to -1 then read all packets
    buffer_size   - how many chunks to read at a time.
  Outputs:
    Data are written to file

  e.g.
  ----------------------------------------------------------------------------------------------------------------------
  from neurapy.neuralynx import lynxio
  import logging
  logging.basicConfig(level=logging.DEBUG)

  channels = 64
  fname = '/Users/kghose/Research/2013/Projects/Workingmemory/Data/NeuraLynx/2013-01-25_14-53-04/DigitalLynxRawDataFile.nrd'
  channel_list = [0,1,2]

  ftsname = 'timestamps.raw'
  fttlname = 'ttl.raw'
  fchanname = ['chan_{:000d}.raw'.format(ch) for ch in channel_list]
  lynxio.extract_nrd(fname, ftsname, fttlname, fchanname, channel_list, channels, max_pkts=1000)
  ----------------------------------------------------------------------------------------------------------------------

  Data are written as a pure stream of binary data and can be easily and efficiently read using the numpy read function.
  For convenience, a function that reads the timestamps, events and channels (read_extracted_data) is included in the library.

  In my experience STX, CRC, timestamp errors and garbage bytes between packets are extremely rare in a properly working system. This function eschews any kind of checks on the data read and just converts the packets. If you suspect that your data has dropped packets, crc or other issues you should try the regular version of this function. You can note if you have packet errors from your Cheetah software.

  Personally, I recommend using the _ec version of the code. It runs fast enough.
  """
  logger.info('Notice: you are using the fast version of the extractor. No error checks are done')

#nrd packet format
  nrd_packet = pylab.dtype([
    ('stx', 'i'),
    ('pkt_id', 'i'),
    ('pkt_data_size', 'i'),
    ('timestamp high', 'I'), #Neuralynx timestamp is ... in its own 32 bit world
    ('timestamp low', 'I'),
    ('status', 'i'),
    ('ttl', 'I'),
    ('extra', '10i'),
    ('data', '{:d}i'.format(channels)),
    ('crc', 'i')
  ])
  #packet_size = nrd_packet.itemsize

  pkt_cnt = 0
  if max_pkts != -1: #An insidious bug was killed here!
    if buffer_size > max_pkts:
      buffer_size = max_pkts

  #The files we will write to. fixme: test for properly opened?
  fts = open(ftsname,'wb')
  fttl = open(fttlname,'wb')
  fchan = [open(fcn,'wb') for fcn in fchanname]

  with open(fname,'rb') as f:
    hdr = read_header(f)
    logger.info('File header: {:s}'.format(hdr))

    #Read in 32bit increments until the magic number is found
    pkt = f.read(4)
    while len(pkt) == 4:
      if pkt == '\x00\x08\x00\x00': #Magic number 2048 0x0800
        f.seek(-4,1) #Realign
        break
      pkt = f.read(4)

    these_packets = pylab.fromfile(f, dtype=nrd_packet, count=buffer_size)
    while these_packets.size > 0:
      ts = pylab.array((these_packets['timestamp high']<<32) | (these_packets['timestamp low']), dtype='uint64')
      ts.tofile(fts)
      these_packets['ttl'].tofile(fttl)
      for idx,ch in enumerate(channel_list):
        these_packets['data'][:,ch].tofile(fchan[idx])

      pkt_cnt += these_packets.size
      if max_pkts != -1:
        if pkt_cnt >= max_pkts: #NOTE: This may give us upto buffer_size -1 more packets than we want.
         break
      these_packets = pylab.fromfile(f, dtype=nrd_packet, count=buffer_size)

  fts.close()
  fttl.close()
  [fch.close() for fch in fchan]

  logger.info('Extracted {:d} packets'.format(pkt_cnt))



def read_extracted_data(fname, type='addata'):
  """Reads data file extracted by extract_nrd.
  Inputs:
    fname - name of the file we want to read.
    type  - type of the data. Has to be one of 'ts','ttl' or 'addata'
      'ts' - time stamps which are uint64 and give values in microseconds
      'ttl' - the parallel port input which is uint32
      'addata' - the continuous A/D channel data which is int32
  Output:
    data - pylab array of appropriate type
  """
  if type == 'ts':
    fmt = 'Q'
  elif type == 'ttl':
    fmt = 'I'
  elif type == 'addata':
    fmt = 'i'
  else:
    logger.error('Unrecognized data type {:s}'.format(type))
    return None

  fin = open(fname,'rb')
  data = pylab.fromfile(fin, dtype=fmt, count=-1)
  fin.close()

  return data
