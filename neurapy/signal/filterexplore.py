"""This code explores the creation of bandpass filters using the scipy.signal module. It shows us how
changing the constraints affects the order of the filter and how well its response can be fitted to the requested one.
The program also prints the a and b coefficients to the command line, so we can use the designed filter in our programs.
The coefficients are also available as vf.a and vf.b in the workspace after exiting the program.

References:

http://docs.scipy.org/doc/scipy/reference/signal.html

http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirdesign.html#scipy.signal.iirdesign

http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html

http://matplotlib.org/users/event_handling.html

Notes:

The whole drag thing is implemented as follows: when we click on the figure, if we are within a small tolerance of the line, we figure our which point we are on, and set up self.being_dragged with the index of the point. We have some constraints we impose. FOr example, the end points can't be dragged and therefore don't respond to clicks. When we move the mouse pointer, if we are in the middle of a drag, the dragged point is moved and out design_line is redrawn. The drag needs to obey certain constraints to keep the design meaningful. When we release the mouse pointer we call the update_design function to compute and draw the new designed filter
"""
import matplotlib
matplotlib.use('macosx')
import pylab, scipy.signal as ss
freqz = ss.freqz

class VisualizeFilter:
  def __init__(self):
    self.being_dragged = None
    self.wp = [.2, .4]
    self.ws = [.1, .6]
    self.gpass = .01
    self.gstop = 30
    self.design_line = None
    self.ftype = 'butter'
    self.setup_main_screen()
    self.connect()
    self.update_design()


  def setup_main_screen(self):
    self.fig = pylab.figure(figsize=(8,8))
    pylab.subplots_adjust(top=.93, bottom=0.1, left=.15, right=.97, hspace=.01, wspace=.15)
    pylab.suptitle('Filter explore')
    self.ax = pylab.subplot2grid((2, 5), (0, 0), colspan=4) #Amplitude response
    self.ax2 = pylab.subplot2grid((2, 5), (1, 0), colspan=4) #Phase response

    meax = self.menuax = pylab.subplot2grid((2, 5), (0, 4))
    menu_names = ['butterworth', 'chebyshev I', 'chebyshev II', 'elliptic', 'bessel']
    self.filter_types = ['butter', 'cheby1', 'cheby2', 'ellip', 'bessel']
    self.menu_h = []
    for n, item in enumerate(menu_names):
      self.menu_h.append(meax.text(0.1, n, item, picker=5))
    self.menu_h[0].set_weight('bold')
    pylab.setp(meax, 'xticks', [], 'yticks', [], 'xlim',[0,1], 'ylim', [-.5, len(menu_names)-.5])

  def connect(self):
    self.cidpick = self.fig.canvas.mpl_connect(
      'pick_event', self.on_pick)
    self.cidrelease = self.fig.canvas.mpl_connect(
      'button_release_event', self.on_release)
    self.cidmotion = self.fig.canvas.mpl_connect(
      'motion_notify_event', self.on_motion)

  def update_design(self):
    ax = self.ax
    ax.cla()
    ax2 = self.ax2
    ax2.cla()

    wp = self.wp
    ws = self.ws
    gpass = self.gpass
    gstop = self.gstop
    b, a = ss.iirdesign(wp, ws, gpass, gstop, ftype=self.ftype, output='ba')
    self.a = a
    self.b = b
    #b = [1,2]; a = [1,2]
    #Print this on command line so we can use it in our programs
    print 'b = ', pylab.array_repr(b)
    print 'a = ', pylab.array_repr(a)

    my_w = pylab.logspace(pylab.log10(.1*self.ws[0]), 0.0, num=512)
    #import pdb;pdb.set_trace()
    w, h = freqz(b, a, worN=my_w*pylab.pi)
    gp = 10**(-gpass/20.)#Go from db to regular
    gs = 10**(-gstop/20.)
    self.design_line, = ax.plot([.1*self.ws[0], self.ws[0], wp[0], wp[1], ws[1], 1.0], [gs, gs, gp, gp, gs, gs], 'ko:', lw=2, picker=5)
    ax.semilogx(w/pylab.pi, pylab.absolute(h),lw=2)
    ax.text(.5,1.0, '{:d}/{:d}'.format(len(b), len(a)))
    pylab.setp(ax, 'xlim', [.1*self.ws[0], 1.2], 'ylim', [-.1, max(1.1,1.1*pylab.absolute(h).max())], 'xticklabels', [])

    ax2.semilogx(w/pylab.pi, pylab.unwrap(pylab.angle(h)),lw=2)
    pylab.setp(ax2, 'xlim', [.1*self.ws[0], 1.2])
    ax2.set_xlabel('Normalized frequency')

    pylab.draw()

  def process_menu(self, event):
    """If we get a pick event and we think it's a menu item, we run this."""
    ea = event.artist
    for n, item in enumerate(self.menu_h):
      if item == ea:
        item.set_weight('bold')
        self.ftype = self.filter_types[n]
      else:
        item.set_weight('normal')
    pylab.draw()
    self.update_design()
    pylab.draw()

  def on_pick(self, event):
    """When we click on the figure and hit either the line or the menu items this gets called."""
    if event.artist != self.design_line: #Might be a menu item
      self.process_menu(event)
      return
    n = len(event.ind)
    if not n: return
    N = event.ind[0]
    if N < 1 or N > 4: return #We don't move end-points
    self.being_dragged = event.ind[0]


  def on_motion(self, event):
    """Move the design line around as we drag the control points."""
    if event.inaxes != self.ax: return
    if self.being_dragged is None: return
    N = self.being_dragged
    line = self.design_line
    x = event.xdata
    y = event.ydata
    xy = line.get_data()
    #Now apply constraints
    if N==1 or N==4: #We are moving stop-band
      if (N==1 and x < xy[0][2]) or \
          (N==4 and x > xy[0][3] and x < xy[0][5]):
        xy[0][N] = x #Move the frequency
      if N==1:
        xy[0][0] = .1*x
      if y > 0 and y < xy[1][2]: #Need to be less than pass band
        xy[1][0] = y
        xy[1][1] = y
        xy[1][4] = y #Move the gain jointly
        xy[1][5] = y

    if N==2 or N==3: #We are moving pass-band
      if (N==2 and x > xy[0][1] and x < xy[0][3]) or \
          (N==3 and x < xy[0][4] and x > xy[0][2]):
        xy[0][N] = x #Move the frequency
      if y > xy[1][1] and y < 1: #Need to be more than stop band
        xy[1][2] = y
        xy[1][3] = y #Move the gain jointly

    line.set_data(xy)
    pylab.draw()

  def on_release(self, event):
    """When we release the mouse, if we were dragging the line, recompute the filter and draw it."""
    if self.being_dragged is None: return

    self.being_dragged = None
    xy = self.design_line.get_data()
    self.wp = [xy[0][2], xy[0][3]]
    self.ws = [xy[0][1], xy[0][4]]
    self.gpass = -20*pylab.log10(xy[1][2])
    self.gstop = -20*pylab.log10(xy[1][1])
    self.update_design()

vf = VisualizeFilter()
pylab.show()