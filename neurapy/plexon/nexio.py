"""Module contains methods to read .nex files produced from the offline sorter.
Based on reverse engineering .nex file structure from http://www.neuroexplorer.com/code.html

"""

import struct, pylab, logging
logger = logging.getLogger(__name__)

upk = struct.unpack
sz = struct.calcsize
rd = lambda f,fmt: upk(fmt, f.read(sz(fmt))) #My first lambda, you'll thank me later

def a_neuron(f, dt, v):
  """Read a neuron from the stream."""
  fmt = str(v['N']) + 'i'
  time_stamps = pylab.array(rd(f,fmt))/dt['Header']['Freq']
  this_neuron = {
    'name': v['name'],
    'version': v['version'],
    'wire no': v['wireNo'],
    'unit no': v['unitNo'],
    'x pos': v['xPos'],
    'y pos': v['yPos'],
    'timestamps': time_stamps
  }
  dt['Neurons'].append(this_neuron)

  return dt

def an_event(f, dt, v):
  """Read an event type from the stream."""
  fmt = str(v['N']) + 'i'
  time_stamps = pylab.array(rd(f,fmt))/dt['Header']['Freq']
  this_event = {
    'name': v['name'],
    'version': v['version'],
    'timestamps': time_stamps
  }
  dt['Events'].append(this_event)

  return dt

def an_interval(f, dt, vars):
  """Read an interval from the stream."""
  logger.debug('An interval')
  return dt

def a_waveform(f, dt, v):
  """Read a waveform set from the stream."""
  fmt = str(v['N']) + 'i'
  time_stamps = pylab.array(rd(f,fmt))/dt['Header']['Freq']

  fmt = str(v['N']*v['npW']) + 'h'
  waveforms = pylab.array(rd(f,fmt)).reshape((v['N'], v['npW'])) * v['AD2mV']

  this_waveform = {
    'name': v['name'],
    'version': v['version'],
    'timestamps': time_stamps,
    'sampling freq': v['wSampF'],
    'waveforms': waveforms
  }
  dt['Waveforms'].append(this_waveform)
  return dt

def a_pop_vector(f, dt, vars):
  """Read a neuron from the stream."""
  return dt

def cont_var(f, dt, vars):
  """Read a neuron from the stream."""
  logger.debug('A continuous variable')
  return dt

def a_marker(f, dt, v):
  """Read a set of markers (which are what we dump from our experiment control software) from the stream."""
  fmt = str(v['N']) + 'i'
  time_stamps = pylab.array(rd(f,fmt))/dt['Header']['Freq']
  this_marker = {
    'name': v['name'],
    'version': v['version'],
    'timestamps': time_stamps
  }
  tmv = this_marker['channel'] = [] #Best name I could come up with
  fmt = str(v['markerLen']) + 's'
  for n in range(v['Nmarkers']):
    mk_name = f.read(64).strip('\x00')
    val = [''] * v['N']
    for m in range(v['N']):
      val[m] = f.read(v['markerLen']).strip('\x00')
    tmv.append({'name': mk_name, 'string': val})

  dt['Markers'].append(this_marker)
  #It's kind of annoying, but this is a direct copy of the matlab code, so presumably this is general enough to handle
  #all .nex files. You have to access the event codes by doing
  #dt['Markers'][0]['channel'][0]['string']

  return dt

def read_nex(fname = '../../Data/SortedNex/Space vs Object Learning-Flippe-07-22-2009-KG.nex'):
  #Python switch statement made as a dictionary. Will lead to keyValueError if we get something out of scope
  switch = {
    0: a_neuron,
    1: an_event,
    2: an_interval,
    3: a_waveform,
    4: a_pop_vector,
    5: cont_var,
    6: a_marker
  }

  f = open(fname, 'r')

  h = {}
  ftid = f.read(4)
  if ftid != 'NEX1':
    logger.error('Not a .nex file : ' + fname)
  else:
    logger.debug('Reading ' + fname)

  fmt = '=i 256s d i i i 260s' # '=' means don't align
  h['Version'], h['Comment'], h['Freq'], h['t begin'], h['t end'], nvar, dummy,  = rd(f,fmt)
  h['Comment'] = h['Comment'].strip('\x00')
  h['t begin'] = h['t begin']/h['Freq']
  h['t end'] = h['t end']/h['Freq']

  dt = {
    'Header': h,
    'Neurons': [],
    'Events': [],
    'Intervals': [],
    'Waveforms': [],
    'Population vectors': [],
    'Continuous variables': [],
    'Markers': []
  }

  for k in range(nvar):
    fmt = '=i i 64s i i i i i i d d d d i i i d 60s'
    v = {}
    type, v['version'], name, offset, v['N'], v['wireNo'], v['unitNo'], v['gain'], v['filter'], \
    v['xPos'], v['yPos'], v['wSampF'], v['AD2mV'], v['npW'], v['Nmarkers'], v['markerLen'], v['mVoffset'], dummy = rd(f,fmt)
    v['name'] = name.strip('\x00')

    this_place = f.tell()
    f.seek(offset)
    dt = switch[type](f, dt, v)
    f.seek(this_place)

  f.close()
  return dt

if __name__ == "__main__":
  dt = read_nex() #Do a test run