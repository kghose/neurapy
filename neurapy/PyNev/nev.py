"""

Kaushik Ghose (kaushik_ghose@med.harvard.edu)
March 2009

Contains routines that allow us to read in NEV files as saved by the
I2S/Cyberkinetics/Blackrock/Whatelse Cerebus system

Path:
Under
/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/
Add a file called nev.pth with just one line
/Users/kghose/Research/Analysis/Python/Nev
or whatever your path to PyNev is
(From http://www.python.org/doc/2.5.2/inst/search-path.html#SECTION000410000000000000000)

Some other resources for nev files
http://nev2lkit.sourceforge.net/  - can't read 2.1 format files
http://neuroshare.sourceforge.net/index.shtml

"""

import logging
logger = logging.getLogger(__name__)

import struct
#For unpacking binary data

import numpy
#I wanted to make this module independent of 'non-standard' python modules, but
#decided that numpy.array would save space and time

import os
#Because we need to create directories if needed

import datetime
# because we set the date and time for the impedance measurement

# Read headers -----------------------------------------------------------------

def read_basic_header(f):
  """Given a freshly opened file handle, read us the basic nev header"""

  f.seek(0, 2)#skip to end
  file_length_in_bytes = f.tell()
  f.seek(0)#skip back to start
  basic_header = {}
  file_type_id = basic_header['file type id'] = f.read(8)
  logger.debug('File id %s' %(file_type_id))
  if file_type_id != 'NEURALEV':
    logger.warning('Not a neural events file : %s' %(file_type_id))
  file_spec = basic_header['file spec'] = f.read(2)
  logger.debug('File spec %d.%d' %(ord(file_spec[0]), ord(file_spec[1]))) 
  additional_flags = f.read(2)
  spike_waveform_is_16bit = basic_header['spike waveform is 16bit'] \
                = bool(ord(additional_flags[0]))
  if spike_waveform_is_16bit:
    logger.debug('Spike waveform is 16 bit')
  else:
    logger.debug('Spike waveform is mixed - look at NEUEVWAV')
  bytes_in_headers, = struct.unpack('I',f.read(4))
  basic_header['bytes in headers'] = bytes_in_headers
  logger.debug('Bytes in headers %d' %(bytes_in_headers))
  bytes_in_data_packets, = struct.unpack('I',f.read(4))
  basic_header['bytes in data packets'] = bytes_in_data_packets
  logger.debug('Bytes in data packets %d' %(bytes_in_data_packets))
  if (bytes_in_data_packets < 12) or \
     (bytes_in_data_packets > 256) or \
     (bytes_in_data_packets % 4 > 0):
    logger.warning('Bytes in data packets not to spec')
  time_stamp_resolution_hz, = struct.unpack('I',f.read(4))
  basic_header['time stamp resolution Hz'] = time_stamp_resolution_hz
  logger.debug('Time stamp resolution %d Hz' %(time_stamp_resolution_hz))
  neural_sample_resolution_hz, = struct.unpack('I',f.read(4))
  basic_header['neural sample resolution Hz'] = neural_sample_resolution_hz
  logger.debug('Neural sample resolution %d Hz' %(neural_sample_resolution_hz))
  time_origin = struct.unpack('8H', f.read(16))
  basic_header['time origin'] = time_origin
  logger.debug('File date %d-%02d-%02d, %02d:%02d:%02d' %(time_origin[0],time_origin[1],time_origin[3],time_origin[4],time_origin[5],time_origin[6]))
  creator = basic_header['creator'] = f.read(32)
  logger.debug('Creator: %s' %(creator))
  comment = basic_header['comment'] = f.read(256)
  logger.debug('Comment: %s' %(comment))
  n_extended_headers, = struct.unpack('I',f.read(4))
  basic_header['number of extended headers'] = n_extended_headers
  logger.debug('Number of extended headers %d' %(n_extended_headers))

  #This is for us, not stored in the actual file
  basic_header['file size'] = file_length_in_bytes
  total_packets = (basic_header['file size'] - basic_header['bytes in headers'])/bytes_in_data_packets
  basic_header['total packets'] = total_packets
  
  return basic_header

def read_extended_header(f, basic_header):
  """Given a file handle and the basic_header, read the extended header. File
  should be spun forward past the basic header"""

  n_extended_headers = basic_header['number of extended headers']
     
  # Extended header
  extended_header = {
    'comments': [],
    'neural event waveform': {}, #NEUEVWAV
    'neural event label': {}, #NEUEVLBL
    'neural event filter': {}, #NEUEVFLT
    'digital label': {}, #DIGLABEL
    'NSAS expt information channels': {}, #NSASEXEV
    'other packets': []
  }

  for nhdr in range(n_extended_headers):
    packet_id = f.read(8)
    payload = f.read(24)
    if packet_id == 'ARRAYNME':
      extended_header['electrode array'] = payload
    elif packet_id == 'ECOMMENT':
      extended_header['comments'].append(payload)
    elif packet_id == 'CCOMMENT':
      extended_header['comments'][-1] += payload
    elif packet_id == 'MAPFILE':
      extended_header['map file'] = payload
    elif packet_id == 'NEUEVWAV':
      parse_neuevwav(extended_header['neural event waveform'], 
                     packet_id, payload)
    else:
      extended_header['other packets'].append([packet_id, payload])

  return extended_header

def parse_neuevwav(neural_event_waveform_dict, packet_id, payload):
  """Given that we know this is a NEUEVWAV packet, parse it and insert it into
  the dict"""
  electrode_info_dict = {}
  offset = 0
  electrode_id, = struct.unpack_from('H', payload)
  offset += 2
  electrode_info_dict['physical connector'], = struct.unpack_from('B', payload, offset)
  offset += 1
  electrode_info_dict['connector pin'], = struct.unpack_from('B', payload, offset)
  offset += 1
  electrode_info_dict['nV per LSB'], = struct.unpack_from('H', payload, offset)
  offset += 2
  electrode_info_dict['energy threshold'], = struct.unpack_from('H', payload, offset)
  offset += 2
  electrode_info_dict['high threshold uV'], = struct.unpack_from('h', payload, offset)
  offset += 2
  electrode_info_dict['low threshold uV'], = struct.unpack_from('h', payload, offset)
  offset += 2
  electrode_info_dict['number of sorted units'], = struct.unpack_from('B', payload, offset)
  offset += 1
  electrode_info_dict['bytes per waveform sample'], = struct.unpack_from('B', payload, offset)
  offset += 1
  
  neural_event_waveform_dict[electrode_id] = electrode_info_dict

# Misc utility functions -------------------------------------------------------

def rewind(f, basic_header):
  """Position file pointer at start of packet data"""
  bytes_in_headers = basic_header['bytes in headers']
  f.seek(bytes_in_headers, 0) #we are now positioned at the start of the data packets


def seek(f, basic_header, extended_header, N = -20):
  """Skip forward(backword) N packets from the current position"""
  bytes_in_data_packets = basic_header['bytes in data packets']
  bytes_to_seek = N * bytes_in_data_packets
  current_pos = f.tell()
  if current_pos + bytes_to_seek < basic_header['bytes in headers']:
    logger.debug('seek_packets: put pointer at start of data packets')
    rewind(f, basic_header)
    return
  if current_pos + bytes_to_seek > basic_header['file size']:
    logger.debug('seek_packets: seek will put pointer beyond eof, not doing')
    return
      
  f.seek(bytes_to_seek, 1)
  
  
def skim_packets(f, basic_header, extended_header, N = 1000, packet_id = None):
  """Read in N sequential packets and return their time stamps and packet ids.
  If packet_id is not None, then pick specific packets with given id"""
  
  Fs = float(basic_header['time stamp resolution Hz'])
  Ts = 1.0/Fs
  bytes_in_data_packets = basic_header['bytes in data packets']
  time_stamps = numpy.zeros(N, dtype='float32')
  packet_ids = numpy.zeros(N, dtype='uint16')
  
  eof = False
  premature_eof = False
  counter = 0
  while not eof and counter < N:
    header_buffer = f.read(6)
    if len(header_buffer) < 6:
      eof = True
      if len(header_buffer) > 0:
        #This means we got cut off in an odd manner
        logger.warning('skim_packets : premature end of file : rewinding')
        premature_eof = True
        rewind(f, basic_header)
    else:
      pi, = struct.unpack('H', header_buffer[4:])
      if pi == packet_id or packet_id is None:
        timestamp, = struct.unpack('I', header_buffer[0:4])
        time_stamps[counter] = timestamp * Ts
        packet_ids[counter] = pi
        counter += 1
     
      f.seek(bytes_in_data_packets - 6, 1) #skip appropriate number of bytes
  
  return time_stamps[:counter], packet_ids[:counter]

def read_next_marker(f, basic_header, extended_header):
  """From the current place in the file find the next marker (non neural) packet
  and return its payload (the value of the digital input uint16) along with its
  time."""

  Fs = float(basic_header['time stamp resolution Hz'])
  bytes_in_data_packets = basic_header['bytes in data packets']
  
  timestamp = None
  code = None
  eof = False
  premature_eof = False
  found = False
  while not eof and not found:
    header_buffer = f.read(6)
    if len(header_buffer) < 6:
      eof = True
      if len(header_buffer) > 0:
        #This means we got cut off in an odd manner
        logger.warning('read_next_marker : premature end of file : rewinding')
        premature_eof = True
        rewind(f, basic_header)
    else:
      pi, = struct.unpack('H', header_buffer[4:])
      if not pi:
        timestamp, = struct.unpack('I', header_buffer[0:4])
        t = timestamp/Fs
        buffer = f.read(bytes_in_data_packets - 6)
        code, = struct.unpack('H', buffer[2:4])
        found = True
      else:  
        f.seek(bytes_in_data_packets - 6, 1) #skip appropriate number of bytes
  
  return t, code

import resource
#to set open file limits

# Code to dump data enmasse in preperation for further analysis ----------------
def fragment(f, basic_header, extended_header,
             frag_dir = 'myspikes/',
             channel_list = numpy.arange(1,97),
             ignore_spike_sorting = True):
  """Given a list of electrodes this will start from the beginning of the file
  and simply dump the spike times and the waveform for each electrode in a 
  separate file. 
  This automatically includes the non-neural events.
  Electrode numbering follows Cerebrus conventions i.e. starting from 1
  
  We might want to rewrite this in C and wrap a python module round it:
  http://starship.python.net/crew/mwh/toext/less-trivial.html
  
  Inputs:
  f - pointer to nev file
  basic_header
  extended_header - both read from the nev file
  channel_list - all the required channels
  ignore_spike_sorting - if true, ignore any online sorted units and dump 
                         everything to unit 0
  """

  if not os.path.exists(frag_dir):
    os.makedirs(frag_dir)
  
  #Open file
  fnonneural = open(frag_dir + '/nonneural.bin', 'wb') 
  
  #Up the limit on open files if needed
  fsoftlimit = resource.getrlimit(resource.RLIMIT_NOFILE)
  if fsoftlimit[0] < channel_list.size * 10:
    resource.setrlimit(resource.RLIMIT_NOFILE, (channel_list.size * 10,-1))
  #We could do this more precisely but seriously, who has 10 spikes per electrode?

  #open files for all the units
  if channel_list.size > 0:
    fout = [None]*(channel_list.max())
  else:
    fout = []
  neuw = extended_header['neural event waveform']
  #print 'fragment: warning, for debugging purposes, fixing sorted units as 4 per electrode'
  for channel in channel_list:
    if not ignore_spike_sorting:
      units_classified = neuw[channel]['number of sorted units']
    else:
      units_classified = 0
      
    #units_classified = 4
    fout[channel-1] = [None]*(units_classified+1) #0 is always the unclassified one
    for m in range(units_classified+1):
      fname = frag_dir + '/channel%02dunit%02d.bin' %(channel,m)
      fout[channel-1][m] = open(fname, 'wb')
  
  #Now just rewind and start redirecting the packets
  rewind(f, basic_header)
  
  bytes_in_data_packets = basic_header['bytes in data packets']
  eof = False
  premature_eof = False  
  nnev_counter = 0
  
  percent_packets_read = 0
  one_percent_packets = basic_header['total packets']/100
  packet_counter = one_percent_packets
  
  while not eof:
    #Read a the packet header
    header_buffer = f.read(6)
    if len(header_buffer) < 6:
      eof = True
      if len(header_buffer) > 0:
        #This means we got cut off in an odd manner
        logger.warning('fragment : premature end of file : rewinding')
        premature_eof = True
        rewind(f, basic_header)
    else:
      pi, = struct.unpack('H', header_buffer[4:])#packet id
      if not pi:
        fnonneural.write(header_buffer)
        fnonneural.write(f.read(bytes_in_data_packets - 6))
        nnev_counter += 1
      elif pi in channel_list:
        buffer = f.read(bytes_in_data_packets - 6)
        if not ignore_spike_sorting:
          sub_unit, =  struct.unpack('B', buffer[0])
        else:
          sub_unit = 0
        thisfptr = fout[pi-1][sub_unit]
        
        thisfptr.write(header_buffer)
        thisfptr.write(buffer)
        #Note that even if we ignore online spike sorting, we preserve the unit
        #identity
      else:  
        f.seek(bytes_in_data_packets - 6, 1) #skip appropriate number of bytes
  
    packet_counter -=1
    if not packet_counter:
      percent_packets_read += 1
      #print 'dump_spike_data: %d%% packets read' %(percent_packets_read)
      packet_counter = one_percent_packets
            
  logger.debug('fragment: found %d non neural packets' %(nnev_counter))
  
  return not premature_eof #return false if there was a problem

# Code to read data packets from fragmented files ------------------------------

def read_frag_nonneural_digital(frag_dir, basic_header):
  """Read the nonneural packets (packet id 0) and return the timestamps and the
  value of the digital input port.
  
  Inputs:
  frag_dir - the directory the fragmented data is in
  basic_header - from reading the nev file
  
  Ouputs:
  time_stamps - absolute times in ms (to match with lablib convention)
  codes - value of the digital input port"""
  
  #Open file
  f = open(frag_dir + '/nonneural.bin') 
  f.seek(0,2)
  file_length_in_bytes = f.tell()
  f.seek(0)#skip back to start

  Fs = float(basic_header['time stamp resolution Hz'])
  T_ms = 1000.0/Fs #ms in one clock cycle
  bytes_in_data_packets = basic_header['bytes in data packets']
  
  N = file_length_in_bytes/bytes_in_data_packets
  
  logger.debug('read_frag_nonneural_digital : expecting %d packets' %(N))
  
  time_stamp_ms = numpy.zeros(N, dtype='float32')
  codes = numpy.zeros(N, dtype='uint16')
  
  eof = False
  premature_eof = False
  counter = 0
  while not eof:
    buffer = f.read(bytes_in_data_packets)
    if len(buffer) < bytes_in_data_packets:
      eof = True
      if len(buffer) > 0:
        #This means we got cut off in an odd manner
        logger.warning('read_frag_nonneural_digital: premature end of file indicative of serious error in dump_spike_data')
        premature_eof = True
    else:
      timestamp, = struct.unpack('I', buffer[0:4])
      time_stamp_ms[counter] = timestamp * T_ms
      codes[counter], = struct.unpack('H', buffer[8:10])
      counter += 1

  return time_stamp_ms[:counter], codes[:counter]


def read_frag_unit(frag_dir, basic_header, extended_header,
                   channel = 1, 
                   unit = 0,
                   tstart_ms = 0.0,
                   tdur_ms = 10.0,
                   load_waveform = False,
                   buffer_increment_size = 1000):
  """Given channel and unit number (0 for unsorted) and the time brackets
  return us the spike data with waveform if needed.
  
  frag_dir - the directory the fragmented data is in
  basic_header - from reading the nev file
  extended_header - from the nev file
  channel - channel(pin) number
  unit - unit number 0 for unclassified
  tstart_ms - take spikes from here
  tdur_ms - for this window. If set to -1 then all the data to the end of the 
            file is read
  load_waveform - if true loads the actual spike waveform as well
  buffer_increment_size - an internal thing - we increase the size of the data
                        arrays by this amount each time we overrun. We squeeze
                        the arrays back to actual size once we are done
  """
  
  if channel < 1 or channel > 255:
    logger.warning('nev.read_frag_unit: Channel given (%d) out of range' %(channel))
    return
  
  f = open(frag_dir + '/channel%02dunit%02d.bin' %(channel,unit))
  
  Fs = float(basic_header['time stamp resolution Hz'])
  Ts = 1.0/Fs
  channel_info_dict = extended_header['neural event waveform'][channel]  
  bytes_in_data_packets = basic_header['bytes in data packets']
  
  spike_time_ms = numpy.zeros(buffer_increment_size, dtype='float32')
  if load_waveform:
    mVperLSB = channel_info_dict['nV per LSB'] * 1e-3 #gives mV
    bytes_per_waveform_sample = channel_info_dict['bytes per waveform sample']  
    if bytes_per_waveform_sample == 2:
      waveform_format = 'h' #signed short
    else:
      logger.warning('%d bytes in data packets not implemented yet' %(bytes_in_data_packets))
      return
    waveform_size = (bytes_in_data_packets - 8)/bytes_per_waveform_sample
    waveform = numpy.zeros((buffer_increment_size,waveform_size), dtype='float32')
 
  eof = False
  premature_eof = False
  counter = 0
  current_buffer_size = spike_time_ms.size
  tstop_ms = tstart_ms + tdur_ms
  while not eof:
    buffer = f.read(bytes_in_data_packets)
    if len(buffer) < bytes_in_data_packets:
      eof = True
      if len(buffer) > 0:
        #This means we got cut off in an odd manner
        logger.warning('read_frag_unit: premature end of file indicative of serious error in dump_spike_data')
        premature_eof = True
    else:
      timestamp, = struct.unpack('I', buffer[0:4])
      tm_ms = timestamp * Ts * 1000 #cerebus time stamps are in clock cycles Ts is in s
      
      if tm_ms < tstart_ms:
        continue
      
      if tdur_ms < 0 or tm_ms < tstop_ms:
        spike_time_ms[counter] = tm_ms 
        if load_waveform:
          waveform[counter,:] = \
            numpy.array(struct.unpack(waveform_size*waveform_format, buffer[8:])) * mVperLSB
        
        counter += 1
        #print 'Loaded %d spikes' %(counter)
        if counter == current_buffer_size: #Need to enlarge the buffer
          spike_time_ms = \
                numpy.concatenate((spike_time_ms, 
                                   numpy.zeros(buffer_increment_size, dtype='float32')))
          if load_waveform:
            waveform = numpy.concatenate((waveform,
                                         numpy.zeros((buffer_increment_size,waveform_size), dtype='float32')))
          current_buffer_size = spike_time_ms.size
      else:
        eof = True
        
  data = {'spike time ms': spike_time_ms[:counter]}
  if load_waveform:
    data['waveform mV'] = waveform[:counter,:]
  
  return data, not premature_eof



#Move to utility or delete -----------------------------------------------------

def total_histogram(fname = None,
                    frag_dir = None,
                    channel_list = numpy.arange(1,97),
                    bin_ms = 1,
                    ignore_spike_sorting = True):
  """Go through the fragmented files and add up all the spikes in the listed 
  channels in bins of bin_ms. This is useful when we want to identify noise
  spikes etc."""

  total_bins = 10
  bins = numpy.zeros(total_bins,dtype=int)#Just to get us started
  
  f = open(fname)
  basic_header = read_basic_header(f)
  extended_header = read_extended_header(f, basic_header)
  
  for channel in channel_list:  
    data, premature_eof = \
    read_frag_unit(frag_dir, basic_header, extended_header,
                     channel = channel, 
                     unit = 0,
                     tstart_ms = 0.0,
                     tdur_ms = -1.0,
                     load_waveform = False,
                     buffer_increment_size = 1000)
    
    this_spike_time_ms = data['spike time ms']
    spike_bin_idx = this_spike_time_ms/bin_ms
    for n_spike in range(len(this_spike_time_ms)):
      n_bin = int(spike_bin_idx[n_spike])
      if n_bin >= total_bins:
        bins.resize(n_bin+1)#fills with zeros
       # print bins.size, n_bin
        bins[n_bin] = 1
        total_bins = len(bins)
      else:
        bins[n_bin] += 1
  
  f.close()
    
  return bins
  
# Convenience functions that bundle together operations ------------------------
def frag(nev_fname, frag_dir):
  f = open(nev_fname)
  basic_header = nev.read_basic_header(f)
  extended_header = nev.read_extended_header(f, basic_header)
  nev.fragment(f, basic_header, extended_header,
               frag_dir = frag_dir)
  f.close()