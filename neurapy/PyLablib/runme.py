import matplotlib

matplotlib.use("Agg")

import pylab

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
  """Return us the eye calibration transform for each trial."""
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
  for tr in range(n_trials):
    M[tr,0,0] = m11[tr][0][0]
    M[tr,0,1] = m12[tr][0][0]
    M[tr,1,0] = m21[tr][0][0]
    M[tr,1,1] = m22[tr][0][0]
    C[tr,0] = tX[tr][0][0]
    C[tr,1] = tY[tr][0][0]
    
  return M,C
    
def fix_window(R):
  """Return the fixation window."""
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
    
def eye_xy(R, M, C):
  """Give us eye position data for each trial."""
  dvx = R.data['Trials']['eyeXData']['Data Values']  
  dvy = R.data['Trials']['eyeYData']['Data Values']
  n_trials = len(dvx)
  
  all_x = []
  all_y = []
  for tr in range(n_trials):
    n_packets = len(dvx[tr])
    raw_x = []
    raw_y = []
    for pak in range(n_packets):
      raw_x += dvx[tr][pak]
      raw_y += dvy[tr][pak]
    
    uncal_x = pylab.array(raw_x, dtype=float)
    uncal_y = pylab.array(raw_y, dtype=float)
    
    x = M[tr,0,0]*uncal_x + M[tr,0,1]*uncal_y + C[tr,0]
    y = M[tr,1,0]*uncal_x + M[tr,1,1]*uncal_y + C[tr,1]
     
    all_x.append(pylab.array(x, dtype=float))
    all_y.append(pylab.array(y, dtype=float))
                 
  return all_x, all_y

def eye_xy_selected(all_x, all_y, trial_no, start_ms, stop_ms, f_samp = 200.0):
  """For the given trial give us the eye samples between the start_ms and stop_ms
  which are equal sized arrays.
  start_ms, stop_ms - trial times (from lablib)"""

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
  f_on = R.data['Trials']['fixate']['Trial Time']
  f_off = R.data['Trials']['saccade']['Trial Time']
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
  pylab.plot(fix_w[20,:,0], fix_w[20,:,1],'k:')
  if dwell_times[n,0] >= 0:
    # We got a fixation
    start_idx = int(f_samp * dwell_times[n,0]/1000.0)
    end_idx = -1
    if dwell_times[n,1] >= 0:
      end_idx = int(f_samp * dwell_times[n,1]/1000.0) - 5
      pylab.plot(all_x[n][start_idx:end_idx], all_y[n][start_idx:end_idx], 'k')
  
  pylab.axis('scaled')
    
#R = lablib.LLDataFileReader(fname = "/Users/kghose/Research/2008-20XX (Monkeys)/Data/Mstim1/DJ2009/Dat/dj-map-2009-06-12-03.dat")

#d_esii = eye_sample_insert_interval(R)
#p_cnt = eye_sample_count_per_packet(R)
#
#N_esii, be = pylab.histogram(d_esii, bins = pylab.linspace(10,30,101), new=True)
#bc_esii = (be[:-1] + be[1:])/2.0
#
#N_escpp, be = pylab.histogram(d_esii, bins = pylab.linspace(10,30,101), new=True)
#bc_escpp = (be[:-1] + be[1:])/2.0

#M, C = eye_calibrations(R)
#
#all_x, all_y = eye_xy(R, M, C)
#fix_w = fix_window(R)
#dwell_times = fixation_box_dwell_times(R)
#
#in_fix_box_x, in_fix_box_y = fixation_box_samples(all_x, all_y, fix_w, dwell_times, f_samp = 200.0)
#
#H, xedges, yedges = pylab.histogram2d(in_fix_box_x,in_fix_box_y, 
#                                      bins=[pylab.linspace(-1.5,1.5,30), pylab.linspace(-1.5,1.5,30)], normed=True)
#
#X,Y = pylab.meshgrid(xedges, yedges)
#pylab.figure(figsize=(3,3))
#pylab.pcolor(X, Y, H, cmap=pylab.cm.gray_r)
#pylab.axis('scaled')
#pylab.axis([-1.5,1.5,-1.5,1.5])
#pylab.savefig('eyepos.png')


trial_no = 1
start_ms = pylab.array([100, 310])
stop_ms = pylab.array([200, 400])
x, y, mx, my = eye_xy_selected(all_x, all_y, trial_no, start_ms, stop_ms, f_samp = 200.0)
pylab.plot(all_x[trial_no], all_y[trial_no], c='gray')
for n in range(start_ms.size):
  pylab.plot(x[n], y[n],'k.-')
pylab.plot(mx, my,'yo')
  
pylab.axis('scaled')

#plot_eye_pos(40, all_x, all_y, fix_w, dwell_times, f_samp = 200.0)
#pylab.plot(bc, N/float(N.sum()), 'k.-')

