"""Some methods for dealing with continuous data. We assume that the original data is in files and that they are
annoyingly large. So all the methods here work on buffered input, using memory maps.
"""
import pylab
from scipy.signal import filtfilt

def filtfiltlong(finname, foutname, fmt, b, a, buffer_len=100000, overlap_len=100, max_len=-1):
  """Use memmap and chunking to filter continuous data.
  Inputs:
    finname -
    foutname    -
    fmt         - data format eg 'int32'
    b,a         - filter coefficients
    buffer_len  - how much data to process at a time
    overlap_len - how much data do we add to the end of each chunk to smooth out filter transients
    max_len     - how many samples to process. If set to -1, processes the whole file
  Outputs:
    y           - The memmapped array pointing to the written file


  Notes on algorithm:
    1. The arrays are memmapped, so we let pylab (numpy) take care of handling large arrays
    2. The filtering is done in chunks:

    Chunking details:
                |<----- buffer------>|
    -----[------*--------------------*------]-------
         |<--------- chunk ---------------->|

    From the array of data we cut out contiguous buffers and to each buffer we add some extra overlap to make a chunk
    The first chunk starts at the first buffer and the last chunk ends at the last buffer

  """
  x = pylab.memmap(finname, dtype=fmt, mode='r')
  if max_len == -1:
    max_len = x.size
  y = pylab.memmap(foutname, dtype=fmt, mode='w+', shape=max_len)

  for buff_st_idx in xrange(0, max_len, buffer_len):
    chk_st_idx = max(0, buff_st_idx - overlap_len)
    buff_nd_idx = min(max_len, buff_st_idx + buffer_len)
    chk_nd_idx = min(x.size, buff_nd_idx + overlap_len)
    rel_st_idx = buff_st_idx - chk_st_idx
    rel_nd_idx = buff_nd_idx - chk_st_idx
    this_y_chk = filtfilt(b, a, x[chk_st_idx:chk_nd_idx])
    y[buff_st_idx:buff_nd_idx] = this_y_chk[rel_st_idx:rel_nd_idx]

  return y



