#!/usr/bin/python

import os
import sys

# For generating *.sc files for all proteins.

# NOTE: This requires qhull files as input.

# USAGE:
#  $ python __FILE__.py <pdb-files.dat>

def main():
  pdbFile = open(sys.argv[1], 'r')
  pdbs = pdbFile.readlines();
  pdbFile.close()
  cleanUp = lambda x: x.split('\n')[0]
  qhullify = lambda x: x.replace('.pdb', '.qhull')
  pdbs = map(cleanUp, pdbs)
  qhulls = map(qhullify, pdbs)
  for q in qhulls:
    if os.path.isfile(q) and q[-6:]=='.qhull':
      print 'Processing: {0}'.format(q)
      generate_sc(os.path.dirname(q), os.path.basename(q))

def generate_sc(qhullPath, qhull):
  qhullPath = "."
  qhullFile = qhullPath + '/' + qhull
  command = "qhull TI \""\
      + qhullFile + "\" d QJ Ft | grep -v \"\\.\" | cut -d\' \' -f2,3,4,5 >\"" + qhullPath + "/" + qhull[:-6] + ".sc\""
  print "Command: ", command
  os.system(command)

main()
