"""Demo/test script, puts nev through its paces"""
import logging
LOG_FILENAME = 'testnev.log'
logging.basicConfig(filename=LOG_FILENAME,
                    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                    level=logging.DEBUG)
logger = logging.getLogger('TestNev')

import pylab
from optparse import OptionParser #For command line arguments

import nev

parser = OptionParser()
parser.add_option("-f", "--nevfile",
                  dest = "nevfile",
                  default = './',
                  help = 'nev file to read in [%default]')
parser.add_option("-o", "--fragdir",
                  dest = "fragdir",
                  default = './Frag',
                  help = 'directory to store the fragmented files [%default]')
(options, args) = parser.parse_args()

fragment = True
read_non_neural = False
load_waveform = False

logger.info('Reading headers')
f = open(options.nevfile)
basic_header = nev.read_basic_header(f)
extended_header = nev.read_extended_header(f, basic_header)

if fragment:
  logger.info('Fragmenting')
  nev.fragment(f, basic_header, extended_header,
               frag_dir = options.fragdir,
               channel_list = pylab.arange(1,97))

f.close()

if read_non_neural:
  logger.info('Reading non neural events')
  time_stamp_ms, codes = \
    nev.read_frag_nonneural_digital(options.fragdir, basic_header)
    
if load_waveform:
  logger.info('Loading all waveforms')  
  data, read_errors = \
    nev.read_frag_unit(options.fragdir, basic_header, extended_header,
                 channel = 1, 
                 unit = 0,
                 tstart_ms = 0.0,
                 tdur_ms = -1.0,
                 load_waveform = True)
  