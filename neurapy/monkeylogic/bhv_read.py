#!/usr/bin/env python
"""
Reads behavioral data files created by MonkeyLogic.

Kaushik Ghose (kaushik.ghose@gmail.com)

"""
from struct import unpack as upk, calcsize as csize
import pylab, logging
logger = logging.getLogger(__name__)

def unpack(fmt,f):
  """A utiltiy function that reads the appropriate number of bytes for the
  format string passed and, if any of the elements are string, strips them of
  whitespace
  Inputs:
    fmt - format string for struct.unpack
    f - file handle
  Returns:
    tuple of read variables."""
  ans = upk(fmt, f.read(csize(fmt)))
  return tuple(an.strip() if isinstance(an, str) else an for an in ans)

def read_header(bhv, fin):
  """Reads the header of a .bhv file. Call this first.
  Inputs:
    bhv - a dictionary (can be empty)
    fin - file handle
  Returns:
    True if no errors
    False otherwise
    bhv is filled in implicitly
  """
  mn, fh, fv = unpack('=I64sd', fin)
  bhv['MagicNumber'] = mn
  bhv['FileHeader'] = fh
  bhv['FileVersion'] = fv
  logger.info('File version {:1.2f}'.format(fv))

  if fv < 1.5:
    logger.error('File ' + fin.name + 'version is ' + str(fv) + '. ' + __name__ + ' only handles versions greater than 2.72')
    return False

  st, en = unpack('=32s128s', fin)
  bhv['StartTime'] = st
  bhv['ExperimentName'] = en
  if fv > 1.5:
    bhv['Investigator'], = unpack('128s', fin)
  bhv['SubjectName'] = unpack('128s', fin)
  if fv > 2.1:
    bhv['ComputerName'], = unpack('128s', fin)
  bhv['ConditionsFile'], = unpack('128s', fin)

  nc, opc = unpack('HH', fin)
  to = [[unpack('64s', fin)[0] for n in xrange(nc)] for m in xrange(opc)]
  bhv['NumConds'] = nc
  bhv['ObjectsPerCond'] = opc
  bhv['TaskObject'] = to

  if fv >= 2.65:
    bhv['TimingFileByCond'] = [unpack('128s', fin)[0] for n in xrange(nc)]
    if fv >=2.71:
      maxblks, = unpack('B', fin)
      bbc = [[unpack('B', fin)[0] for n in xrange(nc)] for m in xrange(maxblks)]
    else:
      bbc = [unpack('B', fin)[0] for n in xrange(nc)]
    ibc = [unpack('128s', fin)[0] for n in xrange(nc)]
    bhv['BlockByCond'] = bbc
    bhv['InfoByCond'] = ibc

  ntf, = unpack('B',fin)
  tf = [unpack('128s', fin)[0] for n in xrange(ntf)]
  bhv['NumTimingFiles'] = ntf
  bhv['TimingFiles'] = tf

  el, bl, cd, bsf, csf = unpack('=64s64s64s64s64s', fin)
  bhv['ErrorLogic'] = el
  bhv['BlockLogic'] = bl
  bhv['CondLogic'] = cd
  bhv['BlockSelectFunction'] = bsf
  bhv['CondSelectFunction'] = csf
  if fv>2.0:
    bhv['VideoRefreshRate'] = unpack('d', fin)[0]
    if fv>3.0:
      bhv['ActualVideoRefreshRate'] = unpack('d', fin)[0]
    bhv['VideoBufferPages'] = unpack('H', fin)[0]

  sxr, syr, vd, ppd = unpack('=HHdd', fin)
  bhv['ScreenXresolution'] = sxr
  bhv['ScreenYresolution'] = syr
  bhv['ViewingDistance'] = vd
  bhv['PixelsPerDegree'] = ppd
  if fv > 2.01:
    ait, aif = unpack('=32sd', fin)
    bhv['AnalogInputType'] = ait
    bhv['AnalogInputFrequency'] = aif
  if fv > 2.0:
    bhv['AnalogInputDuplication'], = unpack('32s', fin)

  escm, tmatrix = unpack('32sB', fin)
  et = {}
  if tmatrix:
    et['ndims in'], et['ndims out'], et['fwdfcn'], et['invfcn'] = unpack('=HH64s64s', fin)
    tsize, = unpack('H', fin)
    tsqrt = int(tsize ** 0.5)
    et['tdata T'] = [[unpack('d', fin)[0] for n in xrange(tsqrt)] for m in xrange(tsqrt)]
    et['tdata Tinv'] = [[unpack('d', fin)[0] for n in xrange(tsqrt)] for m in xrange(tsqrt)]
  bhv['EyeSignalCalibrationMethod'] = escm
  bhv['EyeTransform'] = et


  jcm, tmatrix = unpack('32sB', fin)
  et = {}
  if tmatrix:
    et['ndims in'], et['ndims out'], et['fwdfcn'], et['invfcn'] = unpack('=HH64s64s', fin)
    tsize, = unpack('H', fin)
    tsqrt = int(tsize ** 0.5)
    et['tdata T'] = [[unpack('d', fin)[0] for n in xrange(tsqrt)] for m in xrange(tsqrt)]
    et['tdata Tinv'] = [[unpack('d', fin)[0] for n in xrange(tsqrt)] for m in xrange(tsqrt)]
  bhv['JoystickCalibrationMethod'] = jcm
  bhv['JoyTransform'] = et

  bhv['PhotoDiodePosition'], = unpack('12s', fin)

  if fv >= 1.9:
    bhv['ScreenBackgroundColor'] = unpack('3d', fin)
    bhv['EyeTraceColor'] = unpack('3d', fin)
    bhv['JoyTraceColor'] = unpack('3d', fin)

  np, = unpack('H', fin)
  pic = [{} for n in xrange(np)]
  for n in xrange(np):
    pic[n]['Name'], = unpack('128s', fin)
  for n in xrange(np):
    pic[n]['Size'] = unpack('3H', fin)
  for n in xrange(np):
    sz = pic[n]['Size']
    pic[n]['Data'] = pylab.array([[unpack(sz[0]*'B', fin) for j in xrange(sz[1])] for k in xrange(sz[2])]).T
  bhv['Stimuli'] = {
    'NumPics': np,
    'Pic': pic
  }

  if fv >= 2.5:
    nm, = unpack('H', fin)
    mov = [{} for n in xrange(nm)]
    for n in xrange(nm):
      mov[n]['Name'], = unpack('128s', fin)
    for n in xrange(nm):
      mov[n]['Size'], = unpack('3H', fin)
      mov[n]['NumFrames'], = unpack('H', fin)
    for n in xrange(nm):
      sz = mov[n]['Size']
      numframes = mov[n]['NumFrames']
      mov[n]['Data'] = []
      for fr in xrange(numframes):
        mov[n]['Data'].append(pylab.array([[unpack(sz[0]*'B', fin) for j in xrange(sz[1])] for k in xrange(sz[2])]).T)
    bhv['Stimuli']['NumMovs'] = nm
    bhv['Stimuli']['Mov'] = mov

  logger.info('Read header (' + str(fin.tell()) + ' bytes)')
  return True

def read_trials(bhv, fin):
  """Reads the meat of the .bhv file: the trials. Call this after reading the
  header.
  Inputs:
    bhv - a dictionary filled out by read_header
    fin - file handle positioned at end of header
  Returns:
    True if no errors
    False otherwise
    bhv is filled in implicitly
  """
  fv = bhv['FileVersion']

  bhv['Padding'] = unpack('1024B', fin)
  bhv['NumTrials'], = unpack('H', fin)
  nTrials = bhv['NumTrials']
  logger.info('File contains ' + str(nTrials) + ' trials')
  tno = [0]*nTrials
  atst = [None]*nTrials
  bn = [0]*nTrials
  bi = [0]*nTrials
  cn = [0]*nTrials
  te = [0]*nTrials
  cr = [0]*nTrials
  mcr = [0]*nTrials
  ncodes = [0]*nTrials
  coden = [0]*nTrials
  codet = [0]*nTrials
  xeye = [None]*nTrials
  yeye = [None]*nTrials
  xjoy = [None]*nTrials
  yjoy = [None]*nTrials
  other_analog_data = [[None]*9]*nTrials

  photodiode = [None]*nTrials
  rt = [0]*nTrials
  osr_status = [None]*nTrials
  osr_time = [0]*nTrials
  osr_data = [None]*nTrials
  rr_ront = [None]*nTrials
  rr_rofft = [None]*nTrials
  usrvars = [None]*nTrials

  for trl in xrange(nTrials):
    tno[trl], = unpack('H', fin)
    if fv > 2.2:
      numc, = unpack('B', fin)
      atst[trl] = unpack('d'*numc, fin)

    bn[trl], cn[trl], te[trl] = unpack('3H', fin)
    if fv >=2.05:
      cr[trl], = unpack('H', fin)
      if fv >= 2.72:
        mcr[trl] = unpack('H', fin)[0] if cr[trl]>0 else 0

    ncodes[trl], = unpack('H', fin)
    coden[trl] = unpack('H'*ncodes[trl], fin)
    if fv>=3.0:
      codet[trl] = unpack('I'*ncodes[trl], fin)
    else:
      codet[trl] = unpack('H'*ncodes[trl], fin)

    if fv > 1.5:
      if fv > 1.6:
        fmt='f'
      else:
        fmt='d'

      nxeyep, = unpack('I', fin)
      xeye[trl] = pylab.array(unpack(fmt*nxeyep, fin), dtype=pylab.float16)
      nyeyep, = unpack('I', fin)
      yeye[trl] = pylab.array(unpack(fmt*nyeyep, fin), dtype=pylab.float16)

      nxjoyp, = unpack('I', fin)
      xjoy[trl] = pylab.array(unpack(fmt*nxjoyp, fin), dtype=pylab.float16)
      nyjoyp, = unpack('I', fin)
      yjoy[trl] = pylab.array(unpack(fmt*nyjoyp, fin), dtype=pylab.float16)

      if fv > 2.5:
        for n in xrange(9):
          npts, = unpack('I', fin)
          other_analog_data[trl][n] = pylab.array(unpack('f'*npts, fin), dtype=pylab.float16)
      if fv >= 1.8:
        npts, = unpack('I', fin)
        photodiode[trl] = pylab.array(unpack('f'*npts, fin), dtype=pylab.float16)

      rt[trl], = unpack('h', fin)

    if fv >= 1.9:
      numstat, = unpack('I', fin)
      osr_status_trl = [None]*numstat
      osr_time_trl = [0]*numstat
      osr_data_trl = [None]*numstat
      if fv >= 2.00:
        for n in xrange(numstat):
          nb, = unpack('I', fin)
          osr_status_trl[n] = unpack('B'*nb, fin)
          osr_time_trl[n], = unpack('I', fin)
          if any(x>1 for x in osr_status_trl[n]):
            nf, = unpack('B', fin)
            osr_data_trl[n] = [None]*nf
            for fnum in xrange(nf):
              dc, = unpack('I', fin)
              osr_data_trl[n][fnum] = unpack('d'*dc, fin)
      else:
        for n in xrange(numstat):
          nb, = unpack('I', fin)
          osr_status_trl[n] = unpack('B'*nb, fin)#This might need to be masked in some way
          osr_time_trl[n], = unpack('I', fin)

      osr_status[trl] = osr_status_trl
      osr_time[trl] = osr_time_trl
      osr_data[trl] = osr_data_trl

    if fv >= 1.95:
      nrew, = unpack('I', fin)
      rr_ront[trl] = unpack('I'*nrew, fin)
      rr_rofft[trl] = unpack('I'*nrew, fin)

    if fv >= 2.7:
      nuv, = unpack('B', fin)
      usrvars[trl] = ['',[None]*nuv]
      for n in xrange(nuv):
        varv = None
        varn, type = unpack('32sc', fin)
        if type == 'd':
          lenv, = unpack('B', fin)
          varv = unpack('d'*lenv, fin)
        elif type == 'c':
          varv = unpack('128s', fin)
        usrvars[trl][n][0] = varn
        usrvars[trl][n][1] = varv

  bhv['TrialNumber'] = tno
  bhv['AbsoluteTrialStartTime'] = atst
  bhv['BlockNumber'] = bn
  bhv['BlockIndex'] = bi
  bhv['ConditionNumber'] = cn
  bhv['TrialError'] = te
  bhv['CycleRate'] = cr
  bhv['MinCycleRate'] = mcr
  bhv['NumCodes'] = ncodes
  bhv['CodeNumbers'] = coden
  bhv['CodeTimes'] = codet
  bhv['XEye'] = xeye
  bhv['YEye'] = yeye
  bhv['XJoy'] = xjoy
  bhv['YJoy'] = yjoy
  bhv['OtherAnalogData'] = other_analog_data
  bhv['PhotoDiode'] = photodiode
  bhv['ReactionTime'] = rt
  bhv['ObjectStatusRecord'] = {
    'Status': osr_status,
    'Time': osr_time,
    'Data': osr_data
  }
  bhv['RewardRecord'] = {
    'RewardOnTime': rr_ront,
    'RewardOffTime': rr_rofft
  }
  bhv['UserVars'] = usrvars

  logger.info('Finished reading trials (' + str(fin.tell()) + ' bytes)')
  return True

def read_footer(bhv, fin):
  """Reads the tail of the .bhv file: the summary. Call this after reading the
  trials.
  Inputs:
    bhv - a dictionary filled out by read_header
    fin - file handle positioned at end of all the trials
  Returns:
    True if no errors
    False otherwise
    bhv is filled in implicitly
  """
  try:
    nbcu, = unpack('H', fin)
    bhv['CodeNumbersUsed'] = unpack('H'*nbcu, fin)
    bhv['CodeNamesUsed'] = unpack('64s'*nbcu, fin)

    if bhv['FileVersion'] > 2.05:
      numf, = unpack('H', fin)
      vc = {}
      for n in xrange(numf):
        fn, cnt = unpack('64sH', fin)
        vc[fn] = {
          'Trial': unpack('H'*cnt, fin),
          'Value': unpack('d'*cnt, fin)
        }
      bhv['VariableChanges'] = vc

    bhv['FinishTime'] = unpack('32s', fin)
    logger.info('Finished reading file (' + str(fin.tell()) + ' bytes)')
    return True
  except:
    logger.error('Error reading file footer')
    return False



def read_bhv(fname = '../SampleData/WMHU-MJT-06-04-2012.bhv'):
  """Reads a behavioral file. TODO: partial recovery of files.
  Input:
    fname - name of the file
  Output:
    bhv - dictionary with the behavioral data
  """
  logger.info('Opening ' + fname)
  with open(fname, "rb") as fin:
    bhv = {}
    if read_header(bhv, fin):
      read_trials(bhv, fin)
      read_footer(bhv, fin)

  return bhv