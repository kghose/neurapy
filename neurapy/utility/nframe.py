"""
nframe

In behaving neurophysiological studies we generally present the subject with trials where several events occur, and we are interested in neural responses that occur within some time window relative to those events.

nframe is a framework for processing such data. nframe is based on the pandas DataFrame object. It simply encourages you to organize your data in the form of a two dimensional table (described below) and offers different methods for processing the data.

nframe expects you to organize your behavioral data into a 2d table (a DataFrame) where each row is a trial and each column represents the time of a particular behavioral epoch you are interested in. If you then point nframe to a file containing spike time stamps it will read in the time stamps and chop them up into time bins and put them in a new DataFrame containing the neural data (spike counts)
"""
import pandas as pd, pylab, multiprocessing as mp
import logging
logger = logging.getLogger(__name__)

def epoch_spike_count(df, nnames=[], epochs=[], bin_edges=[]):
  """
  This function chops up the spike train into bins around the chosen epochs and saves the spike counts.

  Inputs:
    df          - the trials DataFrame e.g. trials = session_store['trials']
    nnames      - names of the neurons (list).
    epochs      - names of the columns (these columns will be treated as containing time stamps for the event)
    bin_edges   - list of arrays. Each array is a set of bin edges that demarcates the

  The function expects that df[nnames[0]] will return an array containing spike timestamps and df[epochs[0]] will return a single number that tells us when the event in epochs[0] occured.

  Outputs:
    DataFrame with hierarchical column names.
      The column names are (neuronname, epoch, 0) (neuronname, epoch, 1) etc.

  e.g

  run loadsess.py -f='2013-03-11_15-59-25'
  sc_df = nframe.epoch_spike_count(trials, nnames=['n130311s3c1u1'], epochs=['fixstart','s1on'], bin_edges=[pylab.arange(-.2,.9,.01)]*2)
  pylab.plot(sc_df.n130311s3c1u1.fixstart[trials.trial_type==1].mean())
  pylab.plot(sc_df.n130311s3c1u1.fixstart[trials.trial_type==2].mean())
  """
  sc_df = []
  for nrn in nnames:
    logger.debug('Processing {:s}'.format(nrn))
    if nrn in df.columns:
      ts = df[nrn]#We expect this column to contain spike timestamps
      for epoch,be in zip(epochs, bin_edges):
        this_ts = ts - df[epoch]
        sc = pylab.array([[pylab.nan]*(be.size-1)]*len(df.index))
        #If a row is null then it means that there is no spike data for that trial. This does not mean the neuron did not spike during that period, it means that we did not hold the neuron during that period
        notnulls = pd.notnull(ts) & pd.notnull(df[epoch])
        for n in xrange(len(df.index)):
          if notnulls[n]:
            for m in xrange(be.size-1):
              sc[n,m] = pylab.find((be[m] <= this_ts[n]) & (this_ts[n] < be[m+1])).size
        col_tuples = [(nrn, epoch, n) for n in xrange(be.size-1)]
        col_index = pd.MultiIndex.from_tuples(col_tuples, names=['neuron', 'epoch', 'bin'])
        sc_df.append(pd.DataFrame(sc, columns=col_index, index=df.index))
    else:
      logger.error('No such neuron in file {:s}'.format(nrn))
  return pd.concat(sc_df,axis=1)
