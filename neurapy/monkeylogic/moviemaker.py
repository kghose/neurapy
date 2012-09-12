"""This module (that can be run as a script) reads in a .bhv file and lets you
make movies of the monkey's eye position behavior. Movies of single trials or a
range of trials can be made. Sequentially numbered frames are dumped into a
temporary directory and ffmpeg is used to stitch them together into a movie.

Right now, the movie does not handle TTL objects and movies
"""
import matplotlib
matplotlib.use("Agg") #Don't need to see the frames
import pylab, bhv_read as brd, logging, tempfile, re, argparse, os, subprocess
logger = logging.getLogger(__name__)

def parse_task_object_data(bhv):
  """Convert all the objects into image data and parse their initial positions."""
  obj_data = bhv['Stimuli']['Pic'] #Only handling pics now
  obj_r = re.compile("(\w+)\(") #Regexp to find task object description
  args_r = re.compile("([-.\w]+)[,\)]")#Regexp to extract arguments
  to = bhv['TaskObject']

  objects = []
  initial_pos = []
  for n in xrange(len(to)):
    oname = obj_r.findall(to[n][0])[0]
    if oname == 'fix':
      odata = pylab.ones((5,5,3),dtype=float)#Arbitrary square for FP
      args = args_r.findall(to[n][0])
      p = [float(p) for p in args]
    elif oname =='pic':
      args = args_r.findall(to[n][0])
      picname = args[0] #First one is object name
      p = [float(p) for p in args[1:]]
      for oidx in xrange(len(obj_data)):
        if obj_data[oidx]['Name'] == picname:
          odata = obj_data[oidx]['Data']/255.0 #matplotlib needs [0,1]
          break
    else:
      odata = pylab.zeros((4,4,3))
      logger.error('Could not find object')

    objects.append(odata)
    initial_pos.append(p)

  return objects, pylab.array(initial_pos)

def prepare_trial(bhv, trl, options):
  """From the bhv file extract all the information necessary for plotting a
  trial.

  1. Extract the trial object image data
  2. Interpret the plotting options
  3. Create a timeline setting out what objects must be shown for each frame
     and where they are to be shown including eye position
  4. Package all this data into a dictionary that can be passed on to subsequent
     functions
  """

  scr_col = bhv['ScreenBackgroundColor']
  ppd = bhv['PixelsPerDegree']
  scr_size = (bhv['ScreenXresolution']/ppd, bhv['ScreenYresolution']/ppd)

  osr_status = bhv['ObjectStatusRecord']['Status'][trl]
  osr_data = bhv['ObjectStatusRecord']['Data'][trl]
  osr_time = bhv['ObjectStatusRecord']['Time'][trl]
  rew = [bhv['RewardRecord']['RewardOnTime'][trl], bhv['RewardRecord']['RewardOffTime'][trl]]

  eyex = bhv['XEye'][trl]
  eyey = bhv['YEye'][trl]
  fs = bhv['AnalogInputFrequency']

  objects, pos = parse_task_object_data(bhv)
  ocount = len(objects)
  object_size = pylab.zeros((ocount,2))
  for n in xrange(ocount):
    object_size[n,0] = objects[n].shape[0]/ppd
    object_size[n,1] = objects[n].shape[1]/ppd

  tstep = options['tstep'] #ms between frames
  tstart = bhv['CodeTimes'][trl][2]
  tend = bhv['CodeTimes'][trl][-3]
  fcount = int(tend/tstep)
  logger.debug(str(fcount) + ' frames')

  tframe = pylab.arange(fcount)*tstep #in ms
  eyex_i = pylab.interp(tframe, 1000*pylab.arange(eyex.size)/fs, eyex)#1000x to get ms
  eyey_i = pylab.interp(tframe, 1000*pylab.arange(eyey.size)/fs, eyey)

  frame_data = pylab.zeros((fcount,1+ocount,3))
  #frame_data tells us, for each frame, which objects are visible
  #and where they are located
  #index 0 - frame number
  #index 1 - object number 0=eye, 1 ... ocount = objects
  #index 2 -> [0] - visible (1) or not (0)
  #           [1,2] - xy position in degrees

  cur_ev = 0 #This keeps track of which is the next event in the object status record we are hoing to hit

  #Fill in the eye data
  frame_data[:,0,0] = 1 #Eye is always present
  frame_data[:,0,1] = eyex_i
  frame_data[:,0,2] = eyey_i

  for fr in xrange(fcount):
    while tframe[fr] >= osr_time[cur_ev]:
      for n,ojst in enumerate(osr_status[cur_ev]):
        if ojst == 0:
          frame_data[fr:,n+1,0] = 0 #Switch it off
        elif ojst == 1:
          frame_data[fr:,n+1,0] = 1 #Switch it on
          frame_data[fr:,n+1,1:] = pos[n,:]
        elif ojst == 2: #Move it
          pos[n,:] = osr_data[cur_ev][0]
          frame_data[fr:,n+1,1:] = pos[n,:]
      cur_ev += 1
      if cur_ev >= len(osr_time):
        break

    if cur_ev >= len(osr_time):
      break

  movie_data = {
    'speed': options['speed'],
    'screen color': scr_col,
    'screen size': scr_size,
    'pixels per degree': ppd,
    'frame data': frame_data,
    'tstart': tstart,
    'tframe': tframe,
    'objects': objects,
    'object size': object_size
  }

  return movie_data


def single_frame(movie_data, frame_no):
  """Do the dirty work of converting the movie data into a plot of a single frame."""
  speed = movie_data['speed']
  t_ms = movie_data['tframe'][frame_no] - movie_data['tstart']
  scr_col = movie_data['screen color']
  sx, sy = movie_data['screen size']
  objects = movie_data['objects']
  osize = movie_data['object size']

  pylab.plot([-sx/2, sx/2, sx/2, -sx/2, -sx/2],[sy/2, sy/2, -sy/2, -sy/2, sy/2],'k')
  ax = pylab.gca()
  fd = movie_data['frame data'][frame_no, :, :]
  for n in xrange(fd.shape[0]-1,0,-1): #Need to go backwards to ensure proper z-stack for plotting. Objects earlier in the conditions file obscure later objects
    if fd[n,0]:#Show this image
      cx = fd[n,1]
      cy = fd[n,2]
      ext = [cx - osize[n-1,0]/2, cx + osize[n-1,0]/2,
             -cy + osize[n-1,1]/2, -cy - osize[n-1,1]/2]#l,r,b,t
      pylab.imshow(objects[n-1], extent=ext, interpolation='none')

  pylab.plot(fd[0,1], -fd[0,2], 'w.')  #Plot eye position

  pylab.plot([8, 10, 10], [10, 10, 8], 'y-') #Scale bar
  pylab.text(9, 9, r'$2^{o}$', size=6, horizontalalignment='center', verticalalignment='center',) #Scale bar text

  pylab.text(-0.99*sx/2,0.99*sy/2, '{:04.0f} ms  (x{:1.2f})'.format(t_ms, speed), name='mono', size=7)#Frame info

  pylab.axis('scaled')
  pylab.setp(ax, 'xticks', [], 'yticks', [], 'axis_bgcolor', scr_col, 'ylim',[sy/2, -sy/2], 'xlim', [-sx/2, sx/2])#Ensures reversed y-axis uniformly (otherwise images will flip y-axis w/o warning)


def play(movie_data):
  """Cycle through all the frames, save them as needed."""
  tdir = tempfile.mkdtemp(suffix='_ml_momkr') #Temporary directory for frames
  logger.debug('Frames are being placed in ' + tdir)

  for fr in xrange(movie_data['tframe'].size):
    pylab.figure()
    single_frame(movie_data, frame_no=fr)
    fname = tdir + '/frame{:05d}.png'.format(fr)
    pylab.savefig(fname)
    pylab.close()

  ffmpeg_command = ['ffmpeg','-i', tdir + '/frame%05d.png', '-vcodec', 'libx264', '-x264opts', 'keyint=123:min-keyint=20', '-r', str(options['fps']), '-an', '-y', '-f', 'avi', options['movie name']]
  logger.debug(ffmpeg_command)
  logger.debug(subprocess.call(ffmpeg_command))

if __name__ == "__main__":
  logger.setLevel(logging.DEBUG)
  logger.addHandler(logging.StreamHandler())

  parser = argparse.ArgumentParser()
  parser.add_argument('-f','--file', help=".bhv file full path")
  parser.add_argument('-t', '--trial', help="Trial number", default=0, type=int)
  parser.add_argument('-x', '--speed', help="Speed multiplier of movie relative to actual speed", default=1.0, type=float)
  parser.add_argument('--fps', help="fps of video", default=25.0, type=float)
  parser.add_argument('-v','--verbose', action="store_true", default=False, help="Print logger messages")
  args = parser.parse_args()

  if args.verbose:
    level = logging.DEBUG
  else:
    level = logging.ERROR
  logging.basicConfig(level=level)

  options = {
    'speed': args.speed,
    'tstep': (1000* args.speed)/args.fps,
    'movie name': os.path.basename(args.file)[:-4] + '_t{:04d}_{:1.1f}.avi'.format(args.trial, args.speed),
    #'fps': (1000.0/args.tstep) * args.speed
    'fps': args.fps
  }

  bhv = brd.read_bhv(fname = args.file)
  movie_data = prepare_trial(bhv, args.trial - 1, options)
  logger.debug('{:d} frames {:0.2f} ms'.format(movie_data['tframe'].size, movie_data['tframe'][-1]))
  play(movie_data)