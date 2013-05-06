"""
Dumps trial data from bhv file into a csv file. The csv file structure is as follows:

trial_no, block, condition, result, rt

To convert to standalone macos script use
python ~/bin/pyinstaller-2.0/pyinstaller.py -F -c '/Users/kghose/Research/2008-20XX (Monkeys)/Software/NeuraPy/neurapy/monkeylogic/convert_bhv.py'
"""

import matplotlib
matplotlib.use('macosx')
import argparse, glob, os, bhv_read as brd, logging
logger = logging.getLogger(__name__)


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('-f','--file', help='Convert a single file')
  parser.add_argument('-d', '--dir', help='Convert whole directory')
  parser.add_argument('-v','--verbose', action="store_true", default=False, help="Print logger messages")

  args = parser.parse_args()
  if args.verbose:
    level = logging.DEBUG
  else:
    level = logging.INFO
  logging.basicConfig(level=level)

  if args.file is not None:
    file_list = [args.file]
  elif args.dir is not None:
    file_list = glob.glob(os.path.join(args.dir, '*.bhv'))
  else:
    parser.print_help()
    file_list = []

  for file in file_list:
    print 'Converting {:s}'.format(file)
    bhv = brd.read_bhv(fname = file)

    fout = file.replace('.bhv','.csv')
    import csv
    with open(fout, 'wb') as csvfile:
      print 'Writing {:s}'.format(fout)
      writer = csv.writer(csvfile)
      writer.writerow(['Trial no', 'Block no', 'Condition no', 'Trial error', 'reaction time'])
      for tn,blk,cond,res,rt in zip(bhv['TrialNumber'], bhv['BlockNumber'], bhv['ConditionNumber'], bhv['TrialError'], bhv['ReactionTime']):
        writer.writerow([tn, blk, cond, res, rt])
