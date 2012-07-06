"""
Kaushik Ghose (kaushik_ghose@med.harvard.edu)
March 2009

Contains routines that allow us to read in .NSx files as saved by the
I2S/Cyberkinetics/Blackrock/Whatelse Cerebus system
"""

import struct
#For unpacking binary data

import numpy
#I wanted to make this module independent of 'non-standard' python modules, but
#decided that numpy.array would save space and time

import os
#Because we need to create directories if needed

import resource
#to set open file limits

def read_basic_header(f, verbose = False):
  """Given a freshly opened file handle, read us the basic nsx header"""
  
  basic_header = {}
  message = ''
  f.seek(0, 2)#skip to end
  file_length_in_bytes = f.tell()
  f.seek(0)#skip back to start
  basic_header = {}
  file_type_id = basic_header['file type id'] = f.read(8)
  message += file_type_id + '\n'
  if file_type_id != 'NEURALSG':
    message += 'Not a NSx 2.1 file : %s. Handling not implemented\n' %(file_type_id)
    print message
    return
  
  label = basic_header['label'] = f.read(16)
  message += label + '\n'
  period, = struct.unpack('I', f.read(4))
  Fs = basic_header['Fs Hz'] = 30000.0/period
  message += 'Fs = %f (period = %d)\n' %(Fs,period)
  channel_count, = struct.unpack('I', f.read(4))
  basic_header['number of channels'] = channel_count
  message += '%d channels\n' %(channel_count)
  channel_id = numpy.array(struct.unpack(channel_count*'I', f.read(channel_count*4)))
  basic_header['channel ids'] = channel_id
  message += channel_id.__str__() + '\n'
  
  #This is for us, not stored in the actual file
  basic_header['bytes in header'] = 32 + channel_count*4
  basic_header['file size'] = file_length_in_bytes
  total_samples = (basic_header['file size'] - basic_header['bytes in header'])/2
  samples = total_samples / channel_count
  basic_header['samples per channel'] = samples
  basic_header['waveform bytes'] = 2
  basic_header['waveform format'] = 'h' #short signed int16

  if verbose:
    print message
      
  return basic_header

def rewind(f, basic_header):
  """Position file pointer at start of packet data"""
  bytes_in_header = basic_header['bytes in header']
  f.seek(bytes_in_header, 0) #we are now positioned at the start of the data packets

def length_of_lfp(basic_header, t_dur_ms):
  """Utility function - given the time of the lfp trace we want, return how many
  samples it will have"""
  
  Fs = float(basic_header['Fs Hz'])  
  return int(Fs * t_dur_ms/1000.0 + 0.5)

def read_channel(f, basic_header,
                 channel, 
                 tstart_ms = 0.0,
                 tdur_ms = 100.0):
  """Given channel and the time brackets return us the lfp from the .ns3 file
  directly.
  """

  Fs = float(basic_header['Fs Hz'])   
  channel_count = basic_header['number of channels']
  waveform_bytes = basic_header['waveform bytes']
  fmt = basic_header['waveform format']
  
  if tdur_ms < 0: #Read to end
    tdur_ms = 1000*basic_header['samples per channel']/Fs
    
  Nwave = int(Fs * tdur_ms/1000.0 + 0.5)
  Nstart = int(tstart_ms/1000.0 * Fs + 0.5)
  Nstop = Nstart + Nwave
  lfp = numpy.zeros(Nwave,dtype='short')
  
  offset = waveform_bytes * (channel-1) #We are assuming channels are in order
  start_byte = Nstart * waveform_bytes * channel_count + offset
  skip_bytes = waveform_bytes * (channel_count-1)
  
  rewind(f, basic_header)  
  f.seek(start_byte,1)
  for n in range(Nwave):
    lfp[n], = struct.unpack(fmt, f.read(waveform_bytes))
    f.seek(skip_bytes,1)
  
  return lfp