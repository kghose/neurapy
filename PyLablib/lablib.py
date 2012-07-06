"""

Kaushik Ghose (kaushik_ghose@med.harvard.edu)
November 2008

Contains routines that allow us to read in a lablib data file and extract the
data in terms of a python dictionary of lists. There are also methods to
save this extracted data into a Python binary file (pickle) and load it back.

This code is known to work on lablib file versions 6.1 and 6.2

This code is written for readability and not speed. The assumption is that it
will be used to convert .dat files into python pickle files which will then be
used for data analysis

Path:
Under
/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/
Add a file called lablib.pth with just one line
/Users/kghose/Research/2008/Projects/Pylablib
or whatever your path to lablib.py is
(From http://www.python.org/doc/2.5.2/inst/search-path.html#SECTION000410000000000000000)


Example usage:

import lablib

R = lablib.LLDataFileReader(fname = "dj-2008-10-30-01.dat")
R.pickle(fname = "dj-2008-10-30-01.pkl")

#To ignore some events do
R = lablib.LLDataFileReader(fname = "dj-2008-11-10-01.dat", ignore_events = ['eyeXData','eyeYData','eyePData'])

#To force the reader to read bigendian files
R = lablib.LLDataFileReader(fname = "dj-2008-11-10-01.dat", endian = '>')

#To force littleendian
R = lablib.LLDataFileReader(fname = "dj-2008-11-10-01.dat", endian = '<')

Dictionary structure:
"Header"
    "File Name" - the name of the file
    "File Create Date" - 
    "File Create Time" - 
    "Version" - the lablib file version
    "Events" - a dictionary of event definitions, keyed by event name
    "EventsByCode" - a list of events ordered by event code in the lablib file
                     such that calling file_header["EventsByCode"][2] will 
                     return the event definition for code 2

"Experiment" - any data associated not with a single trial but with the expt
               as a whole. Found in events before the first trialStart. 
               Consists of a dictionary of event structures

"Trials" - a list of dictionaries, each of which stores data for one trial

Notes:
1. Events that come in between trials (after trialEnd and before trialStart) are
ignored.


Todo:
1. Need to handle bundling

Performance data:
File: dj-2008-11-10-01.dat
Trials: 1822
Original Lablib file:  48.1 MB
Python conversion time: 32s
Python pickling time: 26s
Pickled file size: 56.1 MB
Matlab m file: 125.9 MB
mat file uncreatable....

Without eye data:
Python conversion time: 16s
Python pickling time: 4s
Pickled file size: 4.7 MB
Matlab m file: 7.5 MB

Notes:

"""
import logging
logger = logging.getLogger('lablib')

import struct 
#Needed when we treat the buffer as a binary stream e.g. have long ints etc that
#are numbers taking up more than one byte

import sys
#Need for error messages when readLLEvent raises an exception and we are not 
#done with the file -> we have an error in the file somewhere

import copy
#Need for deepcopy of the data_dict structure - we have one for the experiment
#header and another for the trials

import cPickle
#Needed for pickling the data structure to file

import glob
#Needed for directory listing for pickle_all

# Stateless functions --------------------------------------------------------
#These functions just require a file pointer, and not any complex state 
#information
def readLLString(f):
  """Lablib convention for a string is to have a single byte with the length of
  the string followed by the string itself (not null terminated)
  
  Inputs:
  f - file object
  
  Output:
  string read
  """
  str_len = ord(f.read(1))
  return f.read(str_len)

def readLLNumeric(f, fmt):
  """We use struct.unpack to grab the binary data into a specific numeric format
  
  Inputs:
  f - file object
  fmt - format string used by struct.unpack
        See http://www.python.org/doc/2.5.2/lib/module-struct.html  
  
  Output:
  Number
  """
  return struct.unpack(fmt, f.read(struct.calcsize(fmt)))[0]

def readLLDataEventDef(f, endian = '', root = False):
  """Read in and return a dictionary corresponding to a lablib dataevent 
  definition found in the header
  
    Lablib event definitions are stored in the file in sequence of nesting. So if an event 
  is of type struct, it has 'tags' that refer to substructs. These structs are
  stored in sequence in the file and are loaded in sequence into the array
  dataDefs (in LLDataEventDef). 
  
  Say that a 'root event' X has two sub events A and B. A has two events AA and 
  AB while B has three events BA, BB and BC. LabLib loads the file such that
  the X data structure ends up with 7 LLDataDefs in dataDefs 
  A AA AB B BA BB BC
  
  The way we know which is which is because when we access A we see it has 2 
  tags and these will be the two LLDataDefs immediately following A.
  
  Python allows us to make the logical structure of the nesting more explicit.
  Here we only have one datatype - LLDataDef. Each LLDataDef can have zero or
  more child LLDataDefs and so on recursively.
   
  I have maintained John's naming conventions (camelCase) to make debugging
  easier.
  
  Notes:
  1. dataBytes = 0 no data
               > 0 fixed length data
               = -1 variable length data. data length in file immediately after
                    the event code
  """

  theEventDef = {}
  data_dict = {}
  if root:
    #Root event. These are dictionaries and have time stamps
    data_dict['Absolute Time'] = []
    data_dict['Trial Time'] = []
    #We are screwed if some one uses these exact strings in a root event
    #definition
  
  #read in this event data
  theEventDef["typeName"] = readLLString(f)
  theEventDef["dataName"]= readLLString(f) 
  theEventDef["offsetBytes"] = readLLNumeric(f, endian + 'L')
  theEventDef["elements"] = readLLNumeric(f, endian + 'l')
  theEventDef["elementBytes"] = readLLNumeric(f, endian + 'L')
  tags = readLLNumeric(f, endian + 'L')
  #Ignoring size consistency check...
  
  #If there are sub nodes    
  if theEventDef["typeName"] == "struct":
    theEventDef["Children"] = {}
    for n_tag in xrange(tags):
      childEventDef, child_data_dict = readLLDataEventDef(f, endian)
      theEventDef["Children"][childEventDef["dataName"]] = childEventDef
      data_dict[childEventDef["dataName"]] = child_data_dict
  else:
    if root:
      if theEventDef["typeName"] != 'no data':
        data_dict['Data Values'] = []
        #On root nodes where we have data we store it under this key. Again, we
        #are screwed if some one uses these as event names...
    else:
      data_dict = []
      #This is the leaf node - i.e. this is where the data is at and should be a 
      #list, not a dictionary

  return theEventDef, data_dict
  
def readLLHeader(f, endian):
  """Read in the header part of file, including date and version strings and the
  event definitions. Call this first.
  Input:
  f - file pointer
  Output:
  file_header - dictionary with keys:
    "File Name" - the name of the file
    "File Create Date" - 
    "File Create Time" - 
    "Version" - the lablib file version
    "Events" - a dictionary of event definitions, keyed by event name
    "EventsByCode" - a list of events ordered by event code in the lablib file
                     such that calling file_header["EventsByCode"][2] will 
                     return the event definition for code 2
  data_dictionary - a template dictionary to store the trial data
                     """
                     

  file_header = {}
  #Read header and find out if everything is all right
  
  data_dictionary = {}
  #The dictionary into which we put our data
  
  #See if this is the right format for us
  buffer = f.read(2)
  #Don't know what this tests for, copied from 
  #LLDataFileReader.m 945
  if (ord(buffer[0]) != 7) or (ord(buffer[1]) < 2) or (ord(buffer[1]) > 6):
    f.close()
    logger.error("Cannot parse format specifier in file")
    return file_header
  
  len_fmt_spec = ord(buffer[1]) #This is the length of the format specifier
  buffer = f.read(len_fmt_spec)
  format_version = float(buffer)
  
  if format_version < 6.1:
    f.close()
    logger.error("File is version %f, need version 6.2" %(float(buffer[2:7])))
    return file_header
  
  file_header["File Name"] = glob.os.path.basename(f.name)
  file_header["Version"] = format_version
  logger.debug("Lablib data file version %f" % format_version)
  
  #This is the python version of John's tricky code in line 974 LLDataFileReader
  N_event_types = readLLNumeric(f, endian + 'l')
  if N_event_types > 1000:
    logger.error('I read %d events which is probably an error. It is possible that the'\
    ' endian-ness (%c) is wrong' %(N_event_types, endian))
    N_event_types = 0
  else:
    logger.debug("%d Events " % N_event_types)

  file_header["Events"] = {}# - eventsDict
  file_header["EventsByCode"] = [None]*N_event_types #Lablib data file indexes events by this
  data_dictionary = {}
  for n_event_code in xrange(N_event_types):
    eventName = readLLString(f)
    dataBytes = readLLNumeric(f, endian + 'l') #How many bytes does the event have
    
    file_header["Events"][eventName], dd = \
      readLLDataEventDef(f, endian, root = True)
    data_dictionary[file_header["Events"][eventName]['dataName']] = dd
    file_header["Events"][eventName]["dataBytes"] = dataBytes
    file_header["Events"][eventName]["EventName"] = eventName 
    file_header["Events"][eventName]["EventCode"] = n_event_code
    #Trust me, we need these two (for when we are analysing event data in file)    
    file_header["EventsByCode"][n_event_code] = file_header["Events"][eventName]
  
  file_header["File Create Date"] = readLLString(f)
  file_header["File Create Time"] = readLLString(f)
  
  return file_header, data_dictionary

def append_new_data_dict_trial(data_dict):
  """Prep the data_dict (returned from readLLHeader) to receive a new trial.
  This works because dictionaries are passed by value and can will remember the
  changes made to them by the called function"""
  for key in data_dict.keys():
    if isinstance(data_dict[key], dict):
      #http://www.python.org/download/releases/2.2/descrintro/#introspection
      #we gotta recurse
      append_new_data_dict_trial(data_dict[key])
    else:
      #Append an empty list for this trial
      data_dict[key].append([])

def read_llevent_code_from_file(f, event_code_fmt):
  """This 'read ahead' is needed because what the event is often decides where
  it should go"""
  try:
    event_code = readLLNumeric(f, event_code_fmt)
  except:
    event_code = None
  return event_code

def skip_llevent_from_file(f, event_code, event_data_def_by_code, endian):
  """Skip over this event
  Inputs:
  f - file object
  event_code - the event code
  events_by_code - dictionary of event definitions indexed by event code
  """
  event_def = event_data_def_by_code[event_code]
  numBytes = event_def["dataBytes"]
  #LLDataFileReader.m:746 and 766
  if numBytes < 0:
    #Variable length data, figure it out LLDataFileReader.m:774
    numBytes = readLLNumeric(f, endian + 'L')
  f.seek(numBytes,1)
  f.seek(struct.calcsize('L'),1)
    
def read_llevent_from_file(f, event_code, event_data_def_by_code, endian):
  """Read an event from the file.
  Inputs:
  f - file object
  event_code - the event code
  events_by_code - dictionary of event definitions indexed by event code
  
  Outputs:
  absolute_time - the time stamp from data
  event_data_buffer - all the data we need from disk for this event"""
  #event_code = readLLNumeric(f, event_code_fmt)
  event_def = event_data_def_by_code[event_code]
  numBytes = event_def["dataBytes"]
  event_data_buffer = None
  #LLDataFileReader.m:746 and 766
  if numBytes < 0:
    #Variable length data, figure it out LLDataFileReader.m:774
    numBytes = readLLNumeric(f, endian + 'L')
  event_data_buffer = f.read(numBytes)
  absolute_time = readLLNumeric(f, endian + 'L')
  
  return absolute_time, event_data_buffer

def append_llevent_in_data_dict(data_dict, 
                                event_code, absolute_time, trial_time,
                                event_data_buffer,
                                event_def,
                                type_dict,
                                endian = ''):
  """Get the data and push it into the appropriate place in the data_dict.
  This is the root level where we have a few special things, like timing and
  text"""

  event_def_name = event_def['dataName']
  event = data_dict[event_def_name]
  event['Absolute Time'][-1].append(absolute_time)
  event['Trial Time'][-1].append(trial_time)

  #LLDataEventDef.m:175 'text' is treated specially 
  if event_def_name == 'text':
    event['Data Values'][-1].append(event_data_buffer)
  else:
    recurse_into_data_dict(event, event_data_buffer, event_def, type_dict, endian, root = True)

def recurse_into_data_dict(data_dict, event_data_buffer, event_def_part, type_dict, endian = '', root = False):
  """This recursive function is used to analyze nested structures. 
  It is passed an event structure (or part thereof) read from the file and the
  pertinent event definition. The job of this function is to analyze the event
  definition, extract pertinent data from the event_data_buffer, assign data
  and/or recursively call itself to analyze any children. It returns a python
  dictionary"""
  offsetBytes = event_def_part["offsetBytes"]
  if event_def_part["typeName"] != "struct":
    #No need to recurse further. Assign and return
    elementBytes = event_def_part["elementBytes"]
    if not elementBytes:
      #Carries no data, data_dict expects none
      return
    else:
      elements = event_def_part["elements"]
      if elements == -1:
        #Variable number of elements, will read to end of structure
        stopBytes = len(event_data_buffer)
        elements = (stopBytes - offsetBytes)/elementBytes
      else:
        stopBytes = offsetBytes+elements*elementBytes
      
      buf = event_data_buffer[offsetBytes:stopBytes]
      fmt = endian + str(elements) + type_dict[event_def_part["typeName"]]
      #See http://mail.python.org/pipermail/tutor/2008-March/060698.html
      #For all data, you can ask for a multitude of values to be returned
      #by repeating the format string or placing a number in front of the 
      #format code. The returned is a tuple. Ain't Python grand?
      try:
        value = list(struct.unpack(fmt, buf))
        if root:
          data_dict['Data Values'][-1].append(value)
        else:
          data_dict[-1].append(value)
      except:
        print event_def_part, fmt, offsetBytes, stopBytes
        print sys.exc_info()
  else:
    #Need to recurse according to children.
    children = event_def_part["Children"] 
    child_keys = children.keys()
    for key in child_keys:
      event_def_child = children[key]
      recurse_into_data_dict(data_dict[key], event_data_buffer[offsetBytes:], event_def_child, type_dict, endian)
      #Important: the offsets are relative, hence we have to pass in this way

  
# Stateful functions ---------------------------------------------------------
#These functions require a state, so they are written as methods of a 
#LLDataFileReader class, so that we can store the state in the instance and
#let functions use it as needed

class LLDataFileReader:
  def __init__(self, fname = None, ignore_events = None, endian = ''):
    self.initialized = False
    #LLDataEventDef.m:13
    self.type_dict = {'short': 'h',
                      'unsigned short': 'H', 
                      'long': 'l',
                      'unsigned long': 'L',
                      'double': 'd',
                      'float': 'f',
                      'string': 's', 
                      'char': 's', 
                      'boolean': 'b'}
    #This dictionary converts lablib's data type convention to struct.unpack's
    if fname is not None:
      self.read(fname, ignore_events, endian)
      
  def read(self, fname = "dj-2008-10-30-01.dat", ignore_events = None, endian = ''):
    """The returned dictionary 'file_data' is described in the doc string of 
    this module.
    ignore_events should be a list of strings"""
    
    self.ignore_events = ignore_events
    self.endian = endian
    
    self.data = {}
    data = self.data

    f = open(fname,"rb")
    
    #Read in header that stores file name, date and the dictionary of event
    #definitions
    data['Header'], data_dict_template = readLLHeader(f, self.endian)
    
    self.data['Trials'] = copy.deepcopy(data_dict_template)
    self.data['Experiment Header'] = copy.deepcopy(data_dict_template)
    self.data['Junk Events'] = copy.deepcopy(data_dict_template)
    
    #Initialize things according to the header etc.
    self.set_state()
    
    #Read in information presented before the trials start - data that 
    #represents experiment parameters etc. that don't change over the file

    data_dict = self.data['Experiment Header']
    append_new_data_dict_trial(data_dict)
    #event_code = self.read_and_append_llevent(f, data_dict)
    event_code = read_llevent_code_from_file(f, self.event_code_fmt)
    while event_code != self.trialStart_code and event_code is not None:
      self.read_and_append_llevent(f, event_code, data_dict)
      event_code = read_llevent_code_from_file(f, self.event_code_fmt)

    #Read in trial events
    data_dict = self.data['Trials']
    data_dict_jnk = self.data['Junk Events']
    append_new_data_dict_trial(data_dict_jnk)
    n_trials = 0
    trial_properly_closed = True
    while event_code is not None:
      while event_code != self.trialStart_code and event_code is not None:
        #self.read_and_append_llevent(f, event_code, data_dict_jnk)
        skip_llevent_from_file(f, event_code, self.events_by_code, self.endian)
        #We could potentially have junk in between the trials
        event_code = read_llevent_code_from_file(f, self.event_code_fmt)
      
      if event_code is not None:
        append_new_data_dict_trial(data_dict) #New trial starts
      while event_code != self.trialEnd_code and event_code is not None:
        trial_properly_closed = False
        if event_code not in self.ignore_event_codes:
          self.read_and_append_llevent(f, event_code, data_dict)
        else:
          skip_llevent_from_file(f, event_code, self.events_by_code, self.endian)
        event_code = read_llevent_code_from_file(f, self.event_code_fmt)

      if event_code == self.trialEnd_code:
        self.read_and_append_llevent(f, event_code, data_dict)
        event_code = read_llevent_code_from_file(f, self.event_code_fmt)
        trial_properly_closed = True
        n_trials += 1
    
      if not trial_properly_closed:
        logger.error("Possible premature end of file. %d trials, %d bytes read" \
        %(n_trials, f.tell()))
        break
    
    f.close()
        
  def set_state(self):
    """This function takes the file_header structure (as returned by 
    readLLHeader) and sets up a state for LLDataFileReader that enables us to
    call the other stateful methods in this class"""

    file_header = self.data["Header"]
    self.events_by_code = file_header["EventsByCode"]

    #Figure out what we have to ignore
    ev = file_header["Events"]
    self.ignore_event_codes = []
    if self.ignore_events is not None:
      for n in xrange(len(self.ignore_events)):
        if ev.has_key(self.ignore_events[n]):
          self.ignore_event_codes.append(ev[self.ignore_events[n]]['EventCode'])
        else:
          logger.warning('Told to ignore event %s which does not exist' %(self.ignore_events[n]))

    #LLDataFileReader.m:728
    #The size of the eventCode depends on the number of events
    #LLDataDoc.m:279. 
    #charCode = unsigned char,
    #shortCode = unsigned short
    #eventCode = long
    self.event_code_fmt = 'b'  
    if len(file_header["Events"]) > 0xff:
      self.event_code_fmt = 'H'
    if len(file_header["Events"]) > 0xffff:
      self.event_code_fmt = 'l'

    #Some important event codes (speeds up ops by avoiding string comparisons and
    #repeated dictionary look ups)
    self.text_code = file_header["Events"]["text"]["EventCode"]
    self.trialStart_code = file_header["Events"]["trialStart"]["EventCode"]
    self.trialEnd_code = file_header["Events"]["trialEnd"]["EventCode"]
    self.sampleZero_code = file_header["Events"]["sampleZero"]["EventCode"]
    self.spikeZero_code = file_header["Events"]["spikeZero"]["EventCode"]  
    self.sample01_code = file_header["Events"]["sample01"]["EventCode"] 
    self.sample_code = file_header["Events"]["sample"]["EventCode"]
    self.spike_code = file_header["Events"]["spike"]["EventCode"]
    self.spike0_code = file_header["Events"]["spike0"]["EventCode"]
    self.fileEnd_code = file_header["Events"]["fileEnd"]["EventCode"]
        
    #Some timing and other variables that we need as we suck in the events
    self.trialStartTime = -1
    self.sampleIntervalMS = -1
    self.kADChannels = 8 #LLIODevice.h:14
    self.lastSampleTime = [None]*self.kADChannels
    self.spikeStartTime = -1
  
    #LLDataFileReader.m:363 BundledEvents
    #We find all the events that begin with bundledEventPrefixes and don't contain
    #bundledEventStops, because these events need to be concatenated together
    #for a given trial
    bundledEventCodes = []
    bundledEventPrefixes = ("sample", "eye", "spike", "timestamp", "VBL", "vbl", "eStimData")
    bundledEventStops = ("calibration", "zero", "window", "eyeCal", "Break")
    for n in xrange(len(file_header["Events"])):
      #Note that EventName is the same as the key for ["Events"]
      if file_header["EventsByCode"][n]["EventName"].startswith(bundledEventPrefixes):
        stop_found = False
        for m in xrange(len(bundledEventStops)):
          if file_header["EventsByCode"][n]["EventName"].lower().find(bundledEventStops[m]) > -1:
            stop_found = True
            break
        if not stop_found:
          bundledEventCodes.append(n)
    self.bundledEventCodes = bundledEventCodes
    #Now ain't this Pythonic, especially compared to the original C code?
    
    #LLDataFileReader.m:165 timedEvents
    timedEventCodes = []
    for n in xrange(len(file_header["Events"])):
      if file_header["EventsByCode"][n]["EventName"] == "intervalOne":
        timedEventCodes.append(n)
      elif file_header["EventsByCode"][n]["EventName"] == "intervalTwo":
        timedEventCodes.append(n)
    self.timedEventCodes = timedEventCodes
    #These are events for which we store the times

    self.initialized = True          

  def read_and_append_llevent(self, f, event_code, data_dict):
    """This reads in an event from file and places it in the appropriate node in
    the data_dict. The event is appended to the last trial started in the
    data_dict. Remember to append_new_data_dict_trial before starting"""
    
    absolute_time, event_data_buffer = \
      read_llevent_from_file(f, event_code, self.events_by_code, self.endian)
    trial_time = self.set_timing_state(event_code, absolute_time, event_data_buffer)
    event_def = self.events_by_code[event_code]
    append_llevent_in_data_dict(data_dict, 
                                event_code, absolute_time, trial_time,
                                event_data_buffer,
                                event_def,
                                self.type_dict,
                                self.endian) 
    
    return event_code  
  
  def set_timing_state(self, event_code, absolute_time, event_data_buffer):
    """Some events cause us to reset some timers"""

    endian = self.endian
    trial_time = absolute_time
    #LLdataFileReader.m:787
    if event_code == self.trialStart_code:
      self.trialStartTime = absolute_time
    #LLdataFileReader.m:790
    if self.trialStartTime == -1:  
      trial_time = -1
    else:
      trial_time = absolute_time - self.trialStartTime
    #LLDataEventDef.m:794
    if event_code == self.sampleZero_code:
      self.sampleIntervalMS = struct.unpack(endian + 'l', event_data_buffer)[0]
      for idx in xrange(self.kADChannels):
        self.lastSampleTime[idx] = 0
    elif event_code == self.spikeZero_code:
      self.spikeStartTime = absolute_time
    elif event_code == self.sample01_code:
      trial_time = self.lastSampleTime[0]
      self.lastSampleTime[0] += self.sampleIntervalMS
      self.lastSampleTime[1] += self.sampleIntervalMS
    elif event_code == self.sample_code:
      #LLDataEventDef.m:808
      #LLStandardDataEvents.h:16 - ADData is composed of two shorts - 'channel' 
      #and 'data'. Assuming sequential storage, 'channel' comes first
      channel = struct.unpack(endian + 'hh', event_data_buffer)[0]
      trial_time = self.lastSampleTime[channel]
      self.lastSampleTime[channel] += self.sampleIntervalMS
    elif event_code == self.spike_code:
      #LLStandardDataEvents.h:11 - TimestampData is composed of a 
      #short ('channel') and long ('time')
      timeStamp = struct.unpack(endian + 'hl', event_data_buffer)[1] 
      trial_time = self.spikeStartTime - self.trialStartTime + timeStamp
    elif event_code == self.spike0_code:
      theTime = struct.unpack(endian + 'l', event_data_buffer)[0]
      trial_time = self.spikeStartTime - self.trialStartTime + theTime
    
    return trial_time  
               
  def recurse_event_data(self, thisEventPartDef, event_data_buffer):
    """This recursive function is used to analyze nested structures. 
    It is passed an event structure (or part thereof) read from the file and the
    pertinent event definition. The job of this function is to analyze the event
    definition, extract pertinent data from the event_data_buffer, assign data
    and/or recursively call itself to analyze any children. It returns a python
    dictionary"""
    offsetBytes = thisEventPartDef["offsetBytes"]
    if thisEventPartDef["typeName"] != "struct":
      #No need to recurse further. Assign and return
      elementBytes = thisEventPartDef["elementBytes"]
      if not elementBytes:
        #Carries no data
        value = None
      else:
        elements = thisEventPartDef["elements"]
        if elements == -1:
          #Variable number of elements, will read to end of structure
          stopBytes = len(event_data_buffer)
          elements = (stopBytes - offsetBytes)/elementBytes
        else:
          stopBytes = offsetBytes+elements*elementBytes
        
        buf = event_data_buffer[offsetBytes:stopBytes]
        fmt = endian + str(elements) + self.type_dict[thisEventPartDef["typeName"]]
        #See http://mail.python.org/pipermail/tutor/2008-March/060698.html
        #For all data, you can ask for a multitude of values to be returned
        #by repeating the format string or placing a number in front of the 
        #format code. The returned is a tuple. Ain't Python grand?
        try:
          value = list(struct.unpack(fmt, buf))
        except:
          logger.error(str(thisEventPartDef) + ' ' + str(fmt) + ' ' + str(offsetBytes) + ' ' + str(stopBytes))
    else:
      #Need to recurse according to children. value becomes a dictionary
      value = {}
      children = thisEventPartDef["Children"] 
      child_keys = children.keys()
      for key in child_keys:
        ch = children[key]
        value[ch["dataName"]] = \
          self.recurse_event_data(ch, event_data_buffer[offsetBytes:])
        #Important: the offsets are relative, hence we have to pass in this way
          
    return value
  
  def pickle(self, fname = 'test.pkl'):
    """Pickle the data"""
    
    f = open(fname,'wb')
    cPickle.dump(self.data, f, protocol = 2)
    f.close()
    
  def un_pickle(self, fname = 'test.pkl'):
    """Un pickle an existing file"""
    
    f = open(fname,'rb')
    self.data = cPickle.load(f)
    f.close()

#Utility functions ----------------------------------------------------------
def load_pickled(fname = 'test.pkl'):
  """Un pickle an existing file and return a LLDataFileReader structure"""
    
  R = LLDataFileReader()
  f = open(fname,'rb')
  R.data = cPickle.load(f)
  f.close()
  return R
    
def pickle(fname = '../Data/dj-2008-11-10-01.dat', 
           ignore_eye_data = True):
  """A convenient function to load a lablib data file and pickle it to speed up
  analysis.
  Inputs:
  fname - filename of lablib .dat file
  ignore_eye_data - true or false depending on if we wanna load the eye positon
                    and pupil data
  Outputs:
  R - the data structure. Pickled file is saved automatically, substituting the
      extension .dat with .pkl"""
  
  if ignore_eye_data:
    ignore_events = ['eyeXData','eyeYData','eyePData']
  else:
    ignore_events = None
  R = LLDataFileReader(fname = fname, ignore_events = ignore_events)
  if ignore_eye_data:
    foutname = fname.replace('.dat','.pkl')
  else:
    foutname = fname.replace('.dat','-eye.pkl')    
  R.pickle(fname = foutname)
  
  return R

def pickle_all(dir = '.', ignore_eye_data = True, force = False):
  """Go through a directory converting lablib data files to pickle files.
  Inputs:
  dir - where are the files located
  ignore_eye_data - if true don't bother to pickle the eye data
  force - if true reconvert files even if they have an existing pickle file"""
  
  dat_files = glob.glob(dir + '/*.dat')
  for n in xrange(len(dat_files)):
    if ignore_eye_data:
      pkl_file = dat_files[n].replace('.dat','.pkl')
    else:
      pkl_file = dat_files[n].replace('.dat','-eye.pkl')
    
    if glob.os.path.exists(pkl_file) and not force:
        continue
      
    pickle(fname = dat_files[n], ignore_eye_data = ignore_eye_data)  
    