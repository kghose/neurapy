"""Some statistical utilities."""
import pylab, cPickle, os
this_dir, this_filename = os.path.split(__file__) #Needed for the lookup table
#see http://stackoverflow.com/questions/779495/python-access-data-in-package-subdirectory
ci_table_fname = os.path.join(this_dir, "ci_table.pkl")
ci_table = cPickle.load(open(ci_table_fname,'rb'))

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
    pc - array (get back several cis) or single value (get back one ci) of percent corrects
    nsamp - number of trials used to obtain each pc
    ci - confidence level (e.g. 0.01, 0.05)
    bootstraps - number of bootstraps to use

  Output:
    3xN array - first row is median (should be approximately same as pc)
                last two rows are lower and upper ci as expected by pylab.errorbar
  """
  def one_ci(pc, nsamp, ci, bootstraps):
    booted_p = boot_p(pc, nsamp, bootstraps)
    booted_p.sort()

    p = pylab.median(booted_p)
    idx_lo = int(bootstraps * ci/2.0)
    idx_hi = int(bootstraps * (1.0-ci/2))

    return p, p-booted_p[idx_lo], booted_p[idx_hi]-p

  #Need to make it user friendly here to handle array/single numbers intelligently
  if pylab.isscalar(pc):
    pc = pylab.array([pc])
    nsamp = pylab.array([nsamp])
  return pylab.array([one_ci(p, ns, ci, bootstraps) for p,ns in zip(pc, nsamp)]).T

def bin_confint_lookup(pc, nsamp, ci = .05):
  """Return the confidence interval from the lookup table.
  Inputs:
    pc - array (get back several cis) or single value (get back one ci) of percent corrects
    nsamp - number of trials used to obtain each pc
    ci - confidence level (e.g. 0.01, 0.05)
    bootstraps - number of bootstraps to use
    use_table - if true then use a precomputed table instead of doing the bootstraps

  Output:
    3xN array - first row is pc
                last two rows are lower and upper ci as expected by pylab.errorbar
  """
  points = ci_table['points']
  values_lo = ci_table['values_lo']
  values_high = ci_table['values_high']

  from scipy.interpolate import griddata
  if pylab.isscalar(pc):
    pc = pylab.array([pc])
    nsamp = pylab.array([nsamp])
  ci_a = pylab.ones(pc.size)*ci
  xi = pylab.array((pc,nsamp,ci_a)).T

  low_ci = griddata(points, values_lo, xi, method='linear')
  high_ci = griddata(points, values_high, xi, method='linear')

  return pylab.array((pc,low_ci,high_ci))


def compute_section(inputs):
  pc = inputs[0]
  nsamp = inputs[1]
  ci = inputs[2]
  points = pylab.zeros((pc.size,3))
  points[:,0] = pc
  points[:,1] = nsamp
  points[:,2] = ci
  this_ci = bin_confint(pc, pylab.ones(pc.size)*nsamp, ci = ci, bootstraps=2000)
  values_lo = this_ci[1,:]
  values_high = this_ci[2,:]
  return points, values_lo, values_high

def generate_ci_table():
  """Generate the ci table that bin_confint uses."""

  pc = pylab.arange(1,step=.1)
  nsamp = 2**pylab.arange(1,10)
  ci = pylab.array([0.01, 0.05, 0.1])
  points = pylab.zeros((pc.size*nsamp.size*ci.size,3)) #pc, nsamp, ci
  values_lo = pylab.zeros(pc.size*nsamp.size*ci.size)
  values_high = pylab.zeros(pc.size*nsamp.size*ci.size)

  from multiprocessing import Pool
  inputs = []
  for i in xrange(ci.size):
    for j in xrange(nsamp.size):
      inputs.append([pc, nsamp[j], ci[i]])

  pool = Pool(processes=10)
  outputs = pool.map(compute_section, inputs)

  for i in xrange(ci.size):
    for j in xrange(nsamp.size):
      idx = (i*nsamp.size + j)*pc.size
      points[idx:idx+pc.size], values_lo[idx:idx+pc.size], values_high[idx:idx+pc.size] = outputs.pop()

  data = {
    'pc': pc,
    'N': nsamp,
    'ci': ci,
    'points': points,
    'values_lo': values_lo,
    'values_high': values_high
  }

  cPickle.dump(data, open(ci_table_fname,'wb'), protocol=-1)
  return data

def boot_confint(val, ci = .05, bootstraps=2000):
  """Full blown bootstrapping for arbitrary variable types.

  val - list of arrays (get back several cis) or array (get back one ci)
  ci - confidence interval level (1-confidence level)
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

  #Need to make it user friendly here to handle list of lists intelligently
  if val.__class__ == list:#We are sending in a list of lists/arrays
    return pylab.array([one_ci(v, ci, bootstraps) for v in val]).T
  else:
    return one_ci(val, ci, bootstraps)

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
  if bootstraps > 1:
    idx = pylab.randint(x.size, size=(x.size, bootstraps))
  else:
    idx = pylab.zeros((x.size,1),dtype=int)
    idx[:,0] = pylab.arange(x.size)
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

def binary_phi(b1,b2):
  """Given two binary variables b1 and b2 return the ."""
  n00 = pylab.find(b1==b2==0).size
  n01 = pylab.find(b1==0 & b2==1).size
  n10 = pylab.find(b1==1 & b2==0).size
  n11 = pylab.find(b1==b2==1).size

  n1s = n10+n11
  n0s = n00+n01
  ns0 = n00+n10
  ns1 = n01+n11

  return (n11*n00 - n10*n01)/(n1s*n0s*ns0*ns1)**.5

if __name__ == '__main__':
  data = generate_ci_table()