"""Module contains methods to read .nex files produced from the offline sorter.
Based on documents from http://www.neuroexplorer.com/code.html. The text file
HowToReadAndWriteNexFilesInCPlusPlus.txt is especially useful
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

def cont_var(f, dt, v):
  """Read a continuous variable from the stream."""
  fmt = str(v['N']) + 'i'
  time_stamps = pylab.array(rd(f,fmt))/dt['Header']['Freq']

  fmt = str(v['N']) + 'i'
  indexes = pylab.array(rd(f,fmt))

  fmt = str(v['npW']) + 'h'
  waveform = pylab.array(rd(f,fmt)) * v['AD2mV']

  waveforms = []
  st = indexes[0]
  for n in xrange(1,indexes.size):
    nd = indexes[n]
    waveforms.append(pylab.array(waveform[st:nd]))
    st = indexes[n]
  waveforms.append(pylab.array(waveform[st:]))

  this_continuous = {
    'name': v['name'],
    'version': v['version'],
    'timestamps': time_stamps,
    'indexes': indexes,
    'sampling freq': v['wSampF'],
    'waveform': waveforms
  }
  dt['Continuous'].append(this_continuous)
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
  # Each time stamp can have several fields each with an associated value
  # The neuroexplorer system only dumps one field whose value is the strobed word stored as a string

  for n in range(v['Nmarkers']):
    mk_name = f.read(64).strip('\x00')
    if mk_name in ['name', 'version', 'timestamps']:#Try to avoid a name clash (it is possible)
      mk_name = 'my' + mk_name
    val = [f.read(v['markerLen']).strip('\x00') for m in range(v['N'])]
    if v['name'] == 'Strobed':
      logger.debug('Marker name is Strobed, treating it as Plexon strobed word and converting it to a numerical array')
      val = pylab.array([int(va) for va in val])
    this_marker[mk_name] = val

  dt['Markers'].append(this_marker)
  return dt

def read_nex(fname = '../../Data/SortedNex/Space vs Object Learning-Flippe-07-22-2009-KG.nex',
             load = ['neurons', 'events', 'intervals', 'waveforms', 'popvectors', 'continuous', 'markers']):
  """Reads nex file into a standard python dictionary.
  fname - name of nex file
  load - list of strings that instruct us what to load (skipping over others). load is a list that includes
         one or more of the following:
         'neurons' - neuron time stamps,
         'events' - events strobed in,
         'intervals' - any defined intervals,
         'waveforms' - waveforms of the units,
         'popvectors' - don't know what this is
         'continuous' - any continuous channels recorded
         'markers' - markers
  """
  data_type = {
    'neurons': 0,
    'events': 1,
    'intervals': 2,
    'waveforms': 3,
    'popvectors': 4,
    'continuous': 5,
    'markers': 6
  }
  types_to_load = [data_type[dt] for dt in load]

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
    'Continuous': [],
    'Markers': []
  }

  for k in range(nvar):
    fmt = '=i i 64s i i i i i i d d d d i i i d 60s'
    v = {}
    type, v['version'], name, offset, v['N'], v['wireNo'], v['unitNo'], v['gain'], v['filter'], \
    v['xPos'], v['yPos'], v['wSampF'], v['AD2mV'], v['npW'], v['Nmarkers'], v['markerLen'], v['mVoffset'], dummy = rd(f,fmt)
    v['name'] = name.strip('\x00')

    if type in types_to_load:
      this_place = f.tell()
      f.seek(offset)
      dt = switch[type](f, dt, v)
      f.seek(this_place)

  f.close()
  return dt

if __name__ == "__main__":
  dt = read_nex() #Do a test run