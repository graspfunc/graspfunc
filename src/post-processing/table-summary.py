#!/usr/bin/python2.7

import sys

if len(sys.argv) <= 2:
  sys.exit("Usage:  {0}".format(sys.argv[0]) + \
           " <results_intra_group> <results_inter_group> [verbose=<0/1/False/True>]")

matchType = 'tetrahedra'
powerLaw = 2.0
verbose = False
epsilon = 1e-2
for arg in sys.argv:
  if arg.startswith('matchType='):
    matchType = arg[len('matchType='):].lower()
  if arg.startswith('powerLaw='):
    powerLaw = float(arg[len('powerLaw='):])
  if arg.startswith('verbose='):
    verboseString = arg[len('verbose='):]
    if verboseString.isdigit():
        verbose = bool(int(verboseString))
    if verboseString.capitalize() == 'True':
        verbose = True

# ====
# Read the proteins and matching scores from the file

# A match is a 5-tuple: [prot1, prot2, matchSize, grp1, grp2]

matches = []

def readMatches(filen):
  global verbose
  global matches
  parity = 0
  linenum = 1
  for line in open(filen):
    if verbose:
      print filen + ': ' + str(linenum) + ':' + line
      linenum += 1
    if line.startswith('#'):  # ignore as comment character
      continue
    if "../input" not in line and not line.startswith("["):
      continue
    if parity == 0: # starting new match
      match = range(5)
    if parity in [0,2]:
      line = line[line.index('../input'):]
      assert(line.startswith('../input'))
      match[parity/2] = line.split('/')[-1].split(':')[0]
      match[parity/2+3] = line.split('/')[-2][1:]
      if (match[parity/2+3] not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"]):
        print line
        print filen
        print linenum
        print match[parity/2+3]
        assert(False)
      # Optionally change name to:  subgroup:protein
      match[parity/2] = match[parity/2+3] + ':' + match[parity/2]
    if parity == 4:
      assert(line.startswith('['))
      if matchType == 'tetrahedra':
          match[2] = int(line.split(':')[1].split(',')[0])
      elif matchType == 'blosum62':
          match[2] = float(line.split(':')[1].split(',')[1])
      else:
          sys.exit("<matchType> must be 'tetrahedra' (default) or 'blusomu62'")
      matches += [match]
    parity = (parity+1)%5

readMatches(sys.argv[1])
readMatches(sys.argv[2])
print "len(matches): ", len(matches)

proteins = []
for match in matches:
  proteins.append(match[0])
  proteins.append(match[1])

# This would be nice, but it's not order-preserving
#   proteins = list(set(proteins))
# So, I'll use a O(n^2) algorithm for now:
for i in range(len(proteins)):
  if proteins[i] in proteins[:i]:
    proteins[i] = ""
proteins = [p for p in proteins if p]
assert(len(proteins) == len(set(proteins)))

print "len(proteins): ", len(proteins)
# If this fails, it's because the files of matches have some duplicates:
if len(proteins)**2 - len(proteins) > 2*len(matches):
  print "**** WARNING:  Some protein matches missing"
assert(len(proteins)**2 - len(proteins) >= 2*len(matches))

# ====
# Now print the tables

print "\nINTRA-GROUP (matching score for all protein pairs within a group):"
for m in matches[:]:  # Will modify original list
     matches.append([m[1], m[0], m[2], m[4], m[3]])
proteins = {}
for p in set([m[0] for m in matches]):
    proteins[p] = ""
for m in matches:
    if not proteins[m[0]]:
        proteins[m[0]] = m[3]
for p in proteins:
    matches.append([p, p, 0, proteins[p], proteins[p]])
groups = set([m[3] for m in matches])
for g in sorted(groups):
    g_matches = [m for m in matches if m[3] == m[4] == g]
    # assert([m for m in matches if m[0] == g_matches[1][1] if m[1] == g_matches[1][0]])
    # assert([m for m in g_matches if m[0] == g_matches[1][1] if m[1] == g_matches[1][0]])
    table = {}
    for m in g_matches:
        if not table.has_key(m[0]):
            table[m[0]] = {}
        table[m[0]][m[1]] = m[2]
    print "*******", table.keys()
    for row in table.keys():
        print "  " + row + " ",
        for col in table.keys():
            if col == row:
                print '  *',
            elif table[row].has_key(col):
                print "{0:>3}".format(table[row][col]),
            else:
                print '  -',
        print ""

print "\nINTER-GROUP (maximum matching score for all protein pairs between two groups):"
maxMatch = {}
for g1 in sorted(groups):
    maxMatch[g1] = {}
    for g2 in sorted(groups):
        maxMatch[g1][g2] = 0
for m in matches:
    if m[2] > maxMatch[m[3]][m[4]]:
        maxMatch[m[3]][m[4]] = m[2]
print "***", sorted(maxMatch.keys())
for g1 in sorted(maxMatch.keys()):
    print "  " + g1 + ":",
    for g2 in sorted(maxMatch[g1].keys()):
        if g1 == g2 and maxMatch[g1][g2] == 0:
                print "{0:>3}".format('*'),
        else:
                print "{0:>3}".format(maxMatch[g1][g2]),
    print ""
