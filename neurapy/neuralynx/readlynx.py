"""Functions to read the flotilla of files produced by the Neuralynx system."""

from struct import unpack as upk, calcsize as csize
import pylab, logging, argparse
logger = logging.getLogger(__name__)

def read_header(fin):
  """Standard 16 kB header."""
  return fin.read(16*1024).strip('\00')

def read_csc(fin):
  """Read a continuous record file"""
  hdr = read_header(fin)
  packet_info = []
  snSamples = []
  fmt = '=QIII512h' #= prevents padding and alignment nonsense 512 based on the Neuralynx we have
  sz = csize(fmt)
  while fin:
    dain = fin.read(sz)
    if len(dain) < sz:
      break
    #data.append(upk(fmt, dain))
    #qwTimeStamp, dwChannelNumber, dwSampleFreq, dwNumValidSamples, snSamples = upk(fmt, dain)
    daup = upk(fmt, dain)
    packet_info.append(daup[:4])
    snSamples.append(daup[4:])

  waveforms = [] #We might have several 'runs'
  tstarts = [] #Start times of the runs
  epsilon = 1e6 #We assume the system isn't dropping samples, and that we need to separate out the packets only if we
                #have have an actual gap in the record (if we stopped and started again). We set this threshold to be 1s
                #Contact NeuraLynx to see if this is reasonable
  wf = None
  end_t_us = -1
  for pak,samp in zip(packet_info, snSamples):
    start_t_us = pak[0]
    if start_t_us - end_t_us > epsilon:
      #Start a new packet
      if wf is not None:#close out the old packet
        waveforms.append(wf[:Nsamps])
      tstarts.append(start_t_us)
      Nsamps = 0
      wf = pylab.zeros(len(packet_info)*512) #based on the Neuralynx we have

    Fs = pak[2]
    N = pak[3]
    wf[Nsamps:Nsamps+N] = samp[:N]
    Nsamps += N
    end_t_us = start_t_us + 1e6*N/Fs

  if wf is not None:#close out the old packet
    waveforms.append(wf[:Nsamps])

  return {'header': hdr, 'packet': packet_info, 'waveforms': waveforms, 'tstarts': tstarts}

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

#fin = open('/Users/kghose/Downloads/data/raw/A18.Ncs')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-10-31_12-46-18/CSC_photo.ncs')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-10-31_12-46-18/CSC_lfp1.ncs')
#fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-12-04_11-18-55/CSC_lfp1.ncs')
#data = read_csc(fin)
#fin.close()

fin = open('/Users/kghose/Research/2012/Projects/Workingmemory/Data/NeuraLynx/2012-12-04_11-18-55/Events.nev')
data = read_nev(fin)
fin.close()
