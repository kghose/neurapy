"""The functions in lynxio.py allow you to extract raw binary data from the raw Neuralynx files (.nrd)
This script does a quick check of the extraction by comparing the data in the TTL events and the Events.nev file.
Tp be able to run this you should be saving the Events.nev file while recording the data.
"""

import os, argparse, logging
logger = logging.getLogger(__name__)
from neurapy.neuralynx import lynxio

def check_codes(codes, times, codes_raw, times_raw):
  """Cycle through each event in turn (ignore events with code 0 - those are recording start and do not match in the
  raw file), find the raw code that corresponds to the matching time and see if the code matches. Flag the file if
  any code does not match up and spit out the first index and time stamp where things don't match. Return True if no
  errors with the file, False otherwise"""

  cri = 0
  for n, (code, time) in enumerate(zip(codes, times)):
    time = int(time)
    if code == 0: continue
    while times_raw[cri] != time:
      #cri += int((time - int(times_raw[cri]))/30.8 + .5)
      fwd = (time - int(times_raw[cri])) >> 5
      if fwd == 0: fwd = 1
      cri += fwd
      if cri >= times_raw.size:
        logger.error('Problem with matching code times for code #{:d}: Ran out of raw data'.format(n))
        return False

    #If we get here we successfully matched up times
    if (code ^ codes_raw[cri]) & 0xff: #More compact way of testing if the lower 8 bits of the two are same
      logger.error('Times matchup, but codes do not for code #{:d}'.format(n))
      return False

  return True #If we get here, all is OK

def process_session(events_nev_fname='Events.nev', timestamps_raw_fname='timestamps.raw', ttl_raw_fname='ttl.raw'):
  """Wrapper around check_codes."""
  for fname in [events_nev_fname, timestamps_raw_fname, ttl_raw_fname]:
    if not os.path.exists(fname):
      logger.info('Missing file {:s}'.format(fname))
      return

  events = lynxio.read_nev(open(events_nev_fname))
  codes = events['packets']['nttl']
  times = events['packets']['timestamp']

  if not len(codes): return #These lead to 0 byte .raw files and screw us up

  times_raw = lynxio.read_extracted_data(timestamps_raw_fname, 'ts')
  codes_raw = lynxio.read_extracted_data(ttl_raw_fname, 'ttl')

  if not check_codes(codes, times, codes_raw, times_raw):
    logger.error('Problem with file')
    return

  logger.info('File conversion checks out')
  return

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('-e', default='Events.nev', help='Full path to Events.nev file')
  parser.add_argument('-t', default='timestamps.raw', help='Full path to timestamps raw file')
  parser.add_argument('-l', default='ttl.raw', help='Full path to ttl raw file')
  args = parser.parse_args()

  logging.basicConfig(level=logging.DEBUG)
  process_session(args.e, args.t, args.l)