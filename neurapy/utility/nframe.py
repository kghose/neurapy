"""
nframe

In behaving neurophysiological studies we generally present the subject with trials where several events occur, and we
are often interested in neural responses that occur within some time window relative to those events. We wish to group
trials based on some trial parameter and compare aggregate neural responses within and across those grouped trials.

nframe is a framework for processing such data, built around the pandas DataFrame object. Pandas was chosen for
* Convenient column and row indexing, including hierarchical indexes and tab completion in ipython (makes it easy to do
interactive data exploration)
* Convenient saving to disk
* Easy interoperability with numpy

nframe simply encourages you to organize your data in two dimensional tables and provides a bunch of filters that
allow you to process this organized data.

Basic behavioral and trial data is stored in a session dataframe.

Session dataframe
-----------------

Each row is a trial and each column corresponds to some trial parameter or event time.

In addition, when you have neural data, the spike times/lfp trace corresponding to each trial are stored as a numpy array
under columns corresponding to the neuron identity, eg. n1, n2 etc.

e.g.
index  | trial_no | correct | stimulus | on_time | off_time | n1              |n2           |
---------------------------------------------------------------------------------------------
 ...   |    1     |   0     |   s1     | .1      |    .3    | [23,45,67 ...]  | [.....]

Here
- index is an index we invent, possibly adding a sesssion code to distinguish this file from others we have.
- trial_no, correct and stimulus are parameters from our trial that we are interested in.
- on_time and off_time are timestamps corresponding to events we are interested in
- n1,n2 are neural data.
  - These can be spike times or continuous data like LFPs
  - If data do not exist for a trial (the session can run longer than there is neural data) then we should fill with a Null
  - We section the data according to some rule, e.g. the data start at on_time and end at off_time

In general it is advantageous to name the columns starting with an alphabet since that permits tab completion in Ipython

"""
import pandas as pd, pylab
import logging
logger = logging.getLogger(__name__)

def epoch_spike_count(df, nnames=[], epochs=[], epoch_names=[], bin_edges=[]):
  """
  This function chops up the spike train into bins around the chosen epochs and saves the spike counts.

  Inputs:
    df          - a session data frame
    nnames      - names of the neurons (list).
    epochs      - names of the columns (these columns will be treated as containing time stamps for the event)
    epoch_names - names of the epochs for the binned data
    bin_edges   - list of arrays. Each array is a set of bin edges that demarcates the bins

  The function expects that df[nnames[0]] will return an array containing spike timestamps and df[epochs[0]] will return a single number that tells us when the event in epochs[0] occured.

  Outputs:
    DataFrame with hierarchical column names.
      The column names are (neuronname, epoch, 0) (neuronname, epoch, 1) etc.

    neuron  | n1                          |n2                             |
    epoch   | e1           |e2            |e1              |e2            |
    bin     | 0 | 1 ....   | 0 | 1 ...    | 0 | 1 ....                    |
    ----------------------------------------------------------------------

    This structure simply makes it more convenient to handle the data for epoched based analysis, since we now have the spikes
    organized as 2D arrays, each element of which is an integer [0,1,2...] indicating how many spikes are present in that bin
    for that trial.

    It also makes accessing epoch based data very convenient, e.g

    df.n1.e2  will return us a n x m array (n trials x m time slices) corresponding to the response of neuron n1 during
    epoch e1. The peri-event spike count plot is simply
    pylab.plot(df.n1.e2.sum())
  """
  def valid_neurons(df, names):
    """The session data frame is differently structured, with no multiindex columns, so we have to treat this is a
    little different."""
    for nrn in nnames:
      if nrn not in df.columns.unique():
        logger.error('No such neuron ({:s}) in file'.format(nrn))
        continue
      else:
        yield nrn

  sc_df = []
  for nrn in valid_neurons(df, nnames):
    logger.debug('Processing {:s}'.format(nrn))
    ts = df[nrn]#We expect this column to contain spike timestamps
    for epoch,epoch_name,be in zip(epochs, epoch_names, bin_edges):
      this_ts = ts - df[epoch]
      sc = pylab.array([[pylab.nan]*(be.size-1)]*len(df.index))
      #If a row is null then it means that there is no spike data for that trial. This does not mean the neuron did not spike during that period, it means that we did not hold the neuron during that period
      notnulls = pd.notnull(ts) & pd.notnull(df[epoch])
      for n in xrange(len(df.index)):
        if notnulls[n]:
          for m in xrange(be.size-1):
            sc[n,m] = pylab.find((be[m] <= this_ts[n]) & (this_ts[n] < be[m+1])).size
      col_tuples = [(nrn, epoch_name, n) for n in xrange(be.size-1)]
      col_index = pd.MultiIndex.from_tuples(col_tuples, names=['neuron', 'epoch', 'bin'])
      sc_df.append(pd.DataFrame(sc, columns=col_index, index=df.index))
  return pd.concat(sc_df,axis=1)

def remove_baseline(df, baseline_epoch_name):
  """
  Pick an epoch as the baseline condition. Get out a nframe with the mean of the baseline subtracted from all the other
  epochs for that neuron

  Inputs:
    df                   - data frame obtained from epoch_spike_count
    nnames               - names of the neurons (list).
    baseline_epoch_name  - the epoch we want to use as a baseline

  Outputs:
    new_df               - identical format, but with the baseline for each neuron subtracted out
  """
  new_df = df.copy()
  epochs = df.columns.get_level_values('epoch').unique()
  if baseline_epoch_name not in epochs:
    logger.error('requested baseline epoch ({:s}) does not exist'.format(baseline_epoch_name))
    return new_df
  for nrn in df.columns.get_level_values('neuron').unique():
    new_df[nrn] -= df[nrn][baseline_epoch_name].mean().mean()
  return  new_df

def epoch_window_average(df, nname, epochs=[], epoch_names=[], bin_edges=[]):
  """
  Pass in a spike count data frame
  """