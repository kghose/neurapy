"""Read in nev file from cerebus and spit out an ASCII file that can be read by
nev2lkit"""
import logging
from neurapy.cerebus import nev

LOG_FILENAME = 'nev2ascii.log'
logging.basicConfig(filename=LOG_FILENAME,
                    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                    level=logging.DEBUG)
logger = logging.getLogger('Nev2ASCII')

from optparse import OptionParser #For command line arguments

def write_header(fout, basic_header, extended_header, channel, total_channels):
  """Given a (ascii) file handle and the basic and extended headers, write the
  header for the ascii file. nev2lkit is fussy about the formatting. Specifically
  it needs the header string to be exact (it's hard coded in their src)
  http://nev2lkit.cvs.sourceforge.net/viewvc/nev2lkit/nev2lkit/lib/ascii.cpp?revision=1.11&view=markup
  """
  Fs = basic_header['neural sample resolution Hz']
  channel_info_dict = extended_header['neural event waveform'][channel]  
  bytes_in_data_packets = basic_header['bytes in data packets']  
  bytes_per_waveform_sample = channel_info_dict['bytes per waveform sample']
  waveform_size = (bytes_in_data_packets - 8)/bytes_per_waveform_sample
  fout.write('%d\n%d\n%d\n' %(int(Fs), waveform_size, total_channels))#Needs all ints
  #Set to the max number of channels possible (255)
  fout.write('Electrode.Unit\tTime stamp\tWave form\n')#Has to be exactly like this

def write_data(fout, channel, unit, data, threshold):
  spike_time_ms = data['spike time ms']
  waveform = data['waveform mV']
  for n in range(spike_time_ms.size):
    if waveform[n,:].max() < threshold:      
      fout.write('%d.%d %f ' %(channel, unit, spike_time_ms[n] * 1000))#Need spike time in s
      for m in range(waveform[n,:].size):
        fout.write('%3.2f\t' %(waveform[n,m]))
    fout.write('\n')
    
  
  
parser = OptionParser()
parser.add_option("-f", "--nevfile",
                  dest = "nevfile",
                  default = '/Users/kghose/Research/2008-20XX (Monkeys)/Data/Mstim1/DOBBY2010/Neural/20100712/afc005.nev',
                  help = 'nev file to read in [%default]')
parser.add_option("--fragdir",
                  dest = "fragdir",
                  default = '/Users/kghose/Research/2008-20XX (Monkeys)/Data/Mstim1/DOBBY2010/NeuralProcessed/NevFrag/2010-07-12-05/',
                  help = 'directory to store the fragmented files [%default]')
parser.add_option("-o", "--outdir",
                  dest = "outdir",
                  default = '/Users/kghose/Research/2008-20XX (Monkeys)/Data/Mstim1/DOBBY2010/NeuralProcessed/ASCII/2010-07-12-05/',
                  help = 'Output file [%default]')
parser.add_option("--threshold",
                  dest = "threshold",
                  default = 500,
                  type = 'float',
                  help = '[%default] Reject spikes above this as artefacts.')


(options, args) = parser.parse_args()

fragment = False
read_non_neural = False
load_waveform = True

logger.info('Reading headers')
f = open(options.nevfile)
basic_header = nev.read_basic_header(f)
extended_header = nev.read_extended_header(f, basic_header)
unit = 0#Unsorted

for channel in range(1,97):
  fout_name =  options.outdir + '/channel%02d.asc' %(channel)
  fout = open(fout_name,'w')
  write_header(fout, basic_header, extended_header, 1, 1)#each file gets a separate channel
  data, data_ok = \
    nev.read_frag_unit(options.fragdir, basic_header, extended_header,
                 channel = channel, 
                 unit = unit,
                 tstart_ms = 0.0,
                 tdur_ms = -1.0,#Read all of it
                 load_waveform = True)
  if data_ok:
    write_data(fout, 1, unit, data, options.threshold)
  fout.close()

f.close()
