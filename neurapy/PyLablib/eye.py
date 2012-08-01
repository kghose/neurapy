"""This module is part of the pylablib package. It contains methods to analyze
eye position data stored in lablib files"""

import logging
import pylab

logger = logging.getLogger('Eye')

# Some diagnostic functions
def eye_sample_insert_interval(R):
  tt = R.data['Trials']['eyeXData']['Trial Time']  
  n_trials = len(tt)
  d_esii = pylab.array([],dtype=float)
  for tr in range(n_trials):
    d_esii = pylab.concatenate((d_esii,pylab.diff(tt[tr])))

  return d_esii

def eye_sample_count_per_packet(R):
  dv = R.data['Trials']['eyeXData']['Data Values']  
  n_trials = len(dv)
  p_cnt = []
  for tr in range(n_trials):
    n_packets = len(dv[tr])
    for pak in range(n_packets):
      p_cnt.append(len(dv[tr][pak]))

  return pylab.array(p_cnt,dtype=float)

def eye_calibrations(R):
  """Return us the eye calibration transform for each trial.
  Input:
   R - lablib data file from lablib.py
  
  Output:
   M - tr x 2 x 2 array
      m11, m12
      m21, m22
   C - tr x 2 array
      tx, ty
   """
  cal = R.data['Trials']['eyeCalibrationData']['cal']
  tX = cal['tX']
  tY = cal['tY']
  m21 = cal['m21']
  m22 = cal['m22']
  m11 = cal['m11']
  m12 = cal['m12']
  n_trials = len(tX)
  M = pylab.zeros((n_trials,2,2), dtype=float)
  C = pylab.zeros((n_trials,2), dtype=float)
  '''John uses col,row when refering to the calibration matrix. So his
  notation for the matrix is  
     m11, m21.
     m12  m22
  Hence th transposition in the code below.'''
  for tr in range(n_trials):
    M[tr,0,0] = m11[tr][0][0]
    M[tr,0,1] = m21[tr][0][0]
    M[tr,1,0] = m12[tr][0][0]
    M[tr,1,1] = m22[tr][0][0]
    C[tr,0] = tX[tr][0][0]
    C[tr,1] = tY[tr][0][0]
    
  return M,C

def fixation_window(R):
  """Return the coordinates of the fixation point and its size.
  Output is a n x 4 array . Each row is a trial columns are
  0 - fp x
  1 - fp y
  2 - width
  3 - height
  Replaces fix_window"""
  fx_x = pylab.array(R.data['Trials']['fixWindowData']['windowDeg']['origin']['x'])[:,0,0]
  fx_y = pylab.array(R.data['Trials']['fixWindowData']['windowDeg']['origin']['y'])[:,0,0]
  fx_w = pylab.array(R.data['Trials']['fixWindowData']['windowDeg']['size']['width'])[:,0,0]
  fx_h = pylab.array(R.data['Trials']['fixWindowData']['windowDeg']['size']['height'])[:,0,0]
  fix_w = pylab.zeros((fx_x.size,4),dtype=float)
  fix_w[:,0] = fx_x + fx_w/2.0
  fix_w[:,1] = fx_y + fx_h/2.0
  fix_w[:,2] = fx_w
  fix_w[:,3] = fx_h
      
  return fix_w
  
  
    
def fix_window(R):
  """Return the fixation window for each trial."""
  fx_x = R.data['Trials']['fixWindowData']['windowDeg']['origin']['x']
  fx_y = R.data['Trials']['fixWindowData']['windowDeg']['origin']['y']
  fx_w = R.data['Trials']['fixWindowData']['windowDeg']['size']['width']
  fx_h = R.data['Trials']['fixWindowData']['windowDeg']['size']['height']
  n_trials = len(fx_x)
  fix_w = pylab.zeros((n_trials,5,2),dtype=float)
  for tr in range(n_trials):
    x = fx_x[tr][0][0]
    y = fx_y[tr][0][0]
    w = fx_w[tr][0][0]
    h = fx_h[tr][0][0]
    fix_w[tr,:,0] = pylab.array([x, x+w, x+w, x, x])
    fix_w[tr,:,1] = pylab.array([y, y, y+h, y+h, y])
    
  return fix_w

def raw_eye_xy(R, old_format = None):
  """Return the raw x,y data points without calibration for diagnostic purposes."""
  if dat_format is None:
    dvx = R.data['Trials']['eyeXData']['Data Values']  
    dvy = R.data['Trials']['eyeYData']['Data Values']
  elif dat_format == 'fixate':
    #Annoying old format
    dvxy = R.data['Trials']['eyeData']['Data Values']
    n_trials = len(dvxy)
    dvx = []
    dvy = []
    for tr in range(n_trials):
      dvx_tr = []
      dvy_tr = []
      for pak in range(len(dvxy[tr])):
        dvx_tr.append(dvxy[tr][pak][::2])
        dvy_tr.append(dvxy[tr][pak][1::2])
      dvx.append(dvx_tr)
      dvy.append(dvy_tr)

  n_trials = len(dvx)
  
  all_x = []
  all_y = []
  for tr in range(n_trials):
    raw_x = []
    raw_y = []
    for pak in range(len(dvx[tr])):
      raw_x += dvx[tr][pak]
    for pak in range(len(dvy[tr])):
      raw_y += dvy[tr][pak]
      
    if len(dvx[tr]) != len(dvy[tr]):
      print 'eye.eye_xy: trial %d : unequal number of packets x=%d y=%d' %(tr, len(dvx[tr]), len(dvy[tr]))
    
    uncal_x = pylab.array(raw_x, dtype=float)
    uncal_y = pylab.array(raw_y, dtype=float)
    
    if uncal_x.size != uncal_y.size:
      print 'eye.eye_xy: trial %d : unequal number of samples x=%d y=%d' %(tr, uncal_x.size, uncal_y.size)
      minnsamp = min(uncal_x.size, uncal_y.size)
      uncal_x = uncal_x[:minnsamp]
      uncal_y = uncal_y[:minnsamp]
    
    all_x.append(pylab.array(uncal_x, dtype=float))
    all_y.append(pylab.array(uncal_y, dtype=float))
                 
  return all_x, all_y
  
  
def eye_xy(R, M, C, old_format = False):
  """Give us eye position data for each trial.
  Inputs:
    R - lablib data structure from reader 
    M - m matrix
    C - c matrix
    old_format - if False then we load the usual MstimPair dat files where the
                 eye xy is stored in eyeXData and eyeYData
                 if True then we look for it in eyeData
  Outputs:
    all_x - list of arrays of the eyeposition for each trial
    all_y - list of arrays of the eyeposition for each trial
  """
  if not old_format:
    dvx = R.data['Trials']['eyeXData']['Data Values']  
    dvy = R.data['Trials']['eyeYData']['Data Values']
  else:
    #Annoying old format
    dvxy = R.data['Trials']['eyeData']['Data Values']
    n_trials = len(dvxy)
    dvx = []
    dvy = []
    for tr in range(n_trials):
      dvx_tr = []
      dvy_tr = []
      for pak in range(len(dvxy[tr])):
        dvx_tr.append(dvxy[tr][pak][::2])
        dvy_tr.append(dvxy[tr][pak][1::2])
      dvx.append(dvx_tr)
      dvy.append(dvy_tr)

  n_trials = len(dvx)
  
  all_x = []
  all_y = []
  for tr in range(n_trials):
    raw_x = []
    raw_y = []
    for pak in range(len(dvx[tr])):
      raw_x += dvx[tr][pak]
    for pak in range(len(dvy[tr])):
      raw_y += dvy[tr][pak]
      
    if len(dvx[tr]) != len(dvy[tr]):
      print 'eye.eye_xy: trial %d : unequal number of packets x=%d y=%d' %(tr, len(dvx[tr]), len(dvy[tr]))
    
    uncal_x = pylab.array(raw_x, dtype=float)
    uncal_y = pylab.array(raw_y, dtype=float)
    
    if uncal_x.size != uncal_y.size:
      print 'eye.eye_xy: trial %d : unequal number of samples x=%d y=%d' %(tr, uncal_x.size, uncal_y.size)
      minnsamp = min(uncal_x.size, uncal_y.size)
      uncal_x = uncal_x[:minnsamp]
      uncal_y = uncal_y[:minnsamp]

    #The formula that John uses for the calibration is (counter to matrix notation)
    # x = m_11 x + m_21 y + tx
    # y = m_12 x + m_22 y + ty 
    # Note the back-diagonal terms are flipped!
    #x = M[tr,0,0]*uncal_x + M[tr,0,1]*uncal_y + C[tr,0]#WRONG CODE!
    #y = M[tr,1,0]*uncal_x + M[tr,1,1]*uncal_y + C[tr,1]
    x = M[tr,0,0]*uncal_x + M[tr,1,0]*uncal_y + C[tr,0]
    y = M[tr,0,1]*uncal_x + M[tr,1,1]*uncal_y + C[tr,1]
      
    all_x.append(pylab.array(x, dtype=float))
    all_y.append(pylab.array(y, dtype=float))
                 
  return all_x, all_y

def eye_xy_selected(all_x, all_y, trial_no, start_ms, stop_ms, f_samp = 200.0):
  """For the given trial give us the eye samples between the start_ms and stop_ms
  which are equal sized arrays.
  start_ms, stop_ms - trial times (from lablib)
  
  Returns-
   x, y - list of arrays for the eye pos
   mx, my - array of mean eye pos"""

  x = [None] * start_ms.size
  y = [None] * start_ms.size
  mx = pylab.zeros(start_ms.size,dtype=float)
  my = pylab.zeros(start_ms.size,dtype=float)
  
  for n in range(start_ms.size):
    start_idx = int(f_samp * start_ms[n]/1000.0)
    stop_idx = int(f_samp * stop_ms[n]/1000.0)
    x[n] = all_x[trial_no][start_idx:stop_idx]
    mx[n] = x[n].mean()
    y[n] = all_y[trial_no][start_idx:stop_idx]
    my[n] = y[n].mean()
  
  return x, y, mx, my

def fixation_box_dwell_times(R):
  """Return us the fixation box dwell start and end times. We assume that 
  'fixate' indicates the entry of the eye into the fixation box, and that 
  'saccade' indicates the exit of the eye from the box. If 'fixate'
  is empty, no fixation occured (-1 is filled) if 'saccade' is empty, the 
  fixation was not broken until the trial stopped (-1 is filled)."""
  
#  if R.data['Trials'].has_key('fixate'):
  f_on = R.data['Trials']['fixate']['Trial Time']
  if R.data['Trials'].has_key('saccade'):  
    f_off = R.data['Trials']['saccade']['Trial Time']
  else:
    f_off = [[] for n in range(len(f_on))]
  
#  else:
#    f_on = R.data['Trials']['fixOn']['Trial Time']
#    f_off = R.data['Trials']['fixOff']['Trial Time']
  
  n_trials = len(f_on)
  dwell_times = pylab.zeros((n_trials,2)) # Start and stop
  for tr in range(n_trials):
    st = f_on[tr]
    if len(st) == 1:
      st = st[0]
    else:
      st = -1
    nd =f_off[tr]
    if len(nd) == 1:
      nd = nd[0]
    else:
      nd = -1
    dwell_times[tr,:] = [st, nd] 
  return dwell_times
  
def fixation_box_samples(all_x, all_y, fix_w, dwell_times, f_samp = 200.0):
  """Collect all x and ys for all trials for when the eye is within the fixation
  box."""
  n_trials = len(all_x)
  in_fix_box_x = pylab.array([],dtype=float)
  in_fix_box_y = pylab.array([],dtype=float)
  for tr in range(n_trials):
    if dwell_times[tr,0] >= 0:
      # We got a fixation
      start_idx = int(f_samp * dwell_times[tr,0]/1000.0)
      end_idx = -1
      if dwell_times[tr,1] >= 0:
        end_idx = int(f_samp * dwell_times[tr,1]/1000.0) - 5
      in_fix_box_x = pylab.concatenate((in_fix_box_x, all_x[tr][start_idx:end_idx]))
      in_fix_box_y = pylab.concatenate((in_fix_box_y, all_y[tr][start_idx:end_idx]))
  return in_fix_box_x, in_fix_box_y    
  
  
def plot_eye_pos(trial_no, all_x, all_y, fix_w, dwell_times, f_samp = 200.0):
  """."""
  n = trial_no
  pylab.plot(all_x[n], all_y[n], c='gray')
  pylab.plot(fix_w[n,:,0], fix_w[n,:,1],'k:')
  if dwell_times[n,0] >= 0:
    # We got a fixation
    start_idx = int(f_samp * dwell_times[n,0]/1000.0)
    end_idx = -1
    if dwell_times[n,1] >= 0:
      end_idx = int(f_samp * dwell_times[n,1]/1000.0) - 5
      pylab.plot(all_x[n][start_idx:end_idx], all_y[n][start_idx:end_idx], 'k')
  
  pylab.axis('scaled')