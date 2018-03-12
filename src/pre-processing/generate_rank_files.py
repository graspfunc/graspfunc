#!/usr/bin/python

import os
import sys
import re

# For generating *.rank files for all proteins

# USAGE:
#  $ python __FILE__.py <rank-files.dat>

# This assumes that there is a rank_files.dat in the
# current directory containing a list of the original
# POOL-generated rank files.

def main():
  for l in open(sys.argv[1], 'r'):
    if l == "": continue
    generate_rank_file(l.split('\n')[0], l[3])


def generate_rank_file(fl, n):
  f = open(fl, 'r')
  lines = f.readlines()
  f.close()

  newdir="."
  if not os.path.isdir(newdir):
    os.mkdir(newdir)
  f = open(newdir+'/'+fl.split('/')[-1].split('.')[0]+'.rank', 'w')
  
  for l in lines:
    #if count > 200:
    #    break
    strs = re.split("\s+", l)
    f.write(strs[1] + " " + strs[2] + "\n")
  f.close()

main()
