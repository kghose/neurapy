"""Some statistical utilities."""
import pylab

def boot_p(pc, nsamp, bootstraps=2000):
  """Given a probability value p and sample size n, return us bootstraps
  number of probability values obtained by random resampling based on p."""
  r = pylab.rand(nsamp, bootstraps)
  z = pylab.zeros((nsamp, bootstraps))
  idx = pylab.find(r < pc)
  z.flat[idx] = 1
  booted_p = z.mean(axis=0)
  return booted_p


def bin_confint(pc, nsamp, ci = .05, bootstraps=2000):
  """Shortcut to computing confidence intervals on bernoulli trials (like
  percent correct).

  Inputs:
    pc - list/array of percent corrects
    nsamp - number of trials used to obtain each pc
    ci - confidence level (e.g. 0.01, 0.05)
    bootstraps - number of bootstraps to use

  Output:
    3xN array - first row is median (should be approximately same as pc)
                last two rows are lower and upper ci as expected by pylab.errorbar
  """
  def one_ci(pc, nsamp, ci, bootstraps):
    booted_p = boot_p(pc, nsamp, bootstraps)
#    r = pylab.rand(nsamp, bootstraps)
#    z = pylab.zeros((nsamp, bootstraps))
#    idx = pylab.find(r < pc)
#    z.flat[idx] = 1
#    booted_p = z.mean(axis=0)
    booted_p.sort()

    p = pylab.median(booted_p)
    idx_lo = int(bootstraps * ci/2.0)
    idx_hi = int(bootstraps * (1.0-ci/2))

    return p, p-booted_p[idx_lo], booted_p[idx_hi]-p

  return pylab.array([one_ci(p, ns, ci, bootstraps) for p,ns in zip(pc, nsamp)]).T

def boot_confint(val, ci = .05, bootstraps=2000):
  """Full blown bootstrapping for arbitrary variable types.

  val - list of lists/arrays
  ci - confidence interval
  bootstraps - bootstraps to use
  """
  def one_ci(v, ci, bootstraps):
    v = pylab.array(v)
    v = pylab.ma.masked_array(v,pylab.isnan(v)).compressed()
    if v.size == 0:
      return pylab.nan, 0, 0 #Nothing to compute

    r = pylab.randint(v.size, size=(v.size, bootstraps))
    booted_samp = pylab.array([pylab.median(v[r[:,n]]) for n in xrange(bootstraps)])
    booted_samp.sort()

    med = pylab.median(booted_samp)
    idx_lo = int(bootstraps * ci/2.0)
    idx_hi = int(bootstraps * (1.0-ci/2))

    return med, med-booted_samp[idx_lo], booted_samp[idx_hi]-med

  return pylab.array([one_ci(v, ci, bootstraps) for v in val]).T

def boot_curvefit(x,y,fit, p0, ci = .05, bootstraps=2000):
  """use of bootstrapping to perform curve fitting.
  Inputs:
    x - x values
    y - corresponding y values
    fit - a packaged fitting function
    p0 - intial parameter list that fit will use

    fit should be a function of the form
      p1 = fit(x, y, p0)
    with p1 being the optimized parameter vector

  Outputs:
    ci - 3xn array (n = number of parameters: median, low_ci, high_ci)
    booted_p - an bxn array of parameter values (b = number of bootstraps)

  An example fit function is:

  def fit(x, y, p0):
    func = lambda p, t: p[0]*pylab.exp(-t/abs(p[1])) + p[2]
    errfunc = lambda p, t, y: func(p, t) - y
    p1, success = optimize.leastsq(errfunc, p0, args=(t, y))
    return p1

  """

  p0 = pylab.array(p0) #Make it an array in case it isn't one
  idx = pylab.randint(x.size, size=(x.size, bootstraps))
  booted_p = pylab.zeros((p0.size, bootstraps))
  for n in xrange(bootstraps):
    booted_p[:,n] = fit(x[idx[:,n]], y[idx[:,n]], p0)

  p_ci = pylab.zeros((3, p0.size))
  for p in xrange(p0.size):
    booted_samp = pylab.sort(booted_p[p])
    med = pylab.median(booted_samp)
    idx_lo = int(bootstraps * ci/2.0)
    idx_hi = int(bootstraps * (1.0-ci/2))
    p_ci[:,p] = [med, med-booted_samp[idx_lo], booted_samp[idx_hi]-med]

  return p_ci, booted_p