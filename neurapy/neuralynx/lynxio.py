"""Functions to read the flotilla of files produced by the Neuralynx system."""

from struct import unpack as upk, pack as pk, calcsize as csize
import pylab, logging, argparse
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

def read_nse(fin):
  """Read single electrode spike record."""
  time_stamp = []
  saen = []
  cell_no = []
  features = []
  waveform = []

  hdr = read_header(fin)
  fmt = '=QII8I32h'
  sz = csize(fmt)
  while fin:
    dain = fin.read(sz)
    if len(dain) < sz:
      break
      #data.append(upk(fmt, dain))
    daup = upk(fmt, dain)
    qwTimeStamp, dwScNumber, dwCellNumber = daup[:3]
    dnParams = daup[3:12]
    snData = daup[12:]

    time_stamp.append(qwTimeStamp)
    saen.append(dwScNumber)
    cell_no.append(dwCellNumber)
    features.append(dnParams)
    waveform.append(snData)

  return {'header': hdr, 'time stamp': time_stamp, 'spike acquisition entity': saen, 'cell number': cell_no,
          'features': features, 'waveform': waveform}

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

#fin = open('/Users/kghose/Downloads/data/raw/A18.Ncs')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-10-31_12-46-18/CSC_photo.ncs')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-10-31_12-46-18/CSC_lfp1.ncs')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-12-04_11-18-55/CSC_lfp1.ncs')
#data = read_csc(fin)
#fin.close()

#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-12-04_11-18-55/Events.nev')
#data = read_nev(fin)
#fin.close()

#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-12-04_11-18-55/SE2.nse')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-12-04_11-18-55/fake1.nse')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-12-04_11-18-55/fake4.nse')
#data = read_nse(fin)
#fin.close()