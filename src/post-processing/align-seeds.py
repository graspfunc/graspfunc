#!/usr/bin/python

# USAGE:
# maxRankedResidue should normally agree with first parameter given to:
#    mains_XXX.py:matchHighlyRankedResidues()
# Second parameter of above function is distCutoff (in Angstroms)
# STRATEGY:  Choose small distCutoff (e.g.:  1 Angstrom), and
#     run mains_seeds.py,
#     and progressively raise distCutoff and re-run this alignment file.
#     Then take the union of those alignments.  (As distCutoff grows,
#     mains_XXX.py will tend to produce more matches of maximum length,
#     and will then take an intersection.  Typically, the intersection will
#     be either a superset or subset of the previous match, although
#     a disagreement (neither super- nor subset) could be possible.)
#     If we get a superset, use it.  If it's a subset, discard the match.

import sys
import re

if len(sys.argv) <= 2:
  sys.exit("Usage:  {0}".format(sys.argv[0]) + \
           " <results_intra_group> <results_inter_group> [maxRankedResidue=<24>] [verbose=<0/1/False/True>]")

maxRankedResidue = 24
minMatchSizeUsed = 5
maxColsPerMergedGroup = 120
displayRank = False

matchType = 'tetrahedra'
verbose = False
epsilon = 1e-2
for arg in sys.argv:
  if arg.startswith('matchType='):
    matchType = arg[len('matchType='):].lower()
  if arg.startswith('maxRankedResidue='):
    maxRankedResidue = float(arg[len('maxRankedResidue='):])
  if arg.startswith('verbose='):
    verboseString = arg[len('verbose='):]
    if verboseString.isdigit():
        verbose = bool(int(verboseString))
    if verboseString.capitalize() == 'True':
        verbose = True

# ====
# Read the proteins and matching scores from the file

# A match is a 5-tuple: [prot1, prot2, matchSize, grp1, grp2]



# Set line for testing.
# line = "[5]: 7, [1/D/64, 3/Q/35, 4/D/13, 10/G/37, 11/I/89, 17/G/196, 19/A/11], [4/D/59, 3/Q/30, 1/D/8, 20/G/32, 13/I/106, 8/G/187, 5/A/6]"
def align(line):
  expr = re.sub(r'([0-9]+/[A-Z]/[0-9]+)', r'"\1"', line)
  expr2 = eval(expr.split(' ', 2)[2])
  for lst in expr2:
    # 3/E/125 -> E125* ; 73/E/125 -> E125
    for i in range(len(lst)):
      if int(lst[i].split('/')[0]) <= maxRankedResidue:
        highRank = '*'
      else:
        highRank = ''
      if displayRank:
        lst[i] = ''.join(lst[i].rsplit('/',1)) + highRank
      else:
        lst[i] = ''.join(lst[i].split('/')[1:]) + highRank
  return zip(expr2[0], expr2[1])

matches = []

def readMatches(file):
  global verbose
  global matches
  parity = 0
  linenum = 1
  for line in open(file):
    if verbose:
      print str(linenum)+':', line
      linenum += 1
    if line.startswith('#'):  # ignore as comment character
      continue
    if "../input" not in line and not line.startswith("["):
      continue
    if parity == 0: # starting new match
      match = range(6)
    if parity in [0,2]:
      line = line[line.index('../input'):]
      assert line.startswith('../input')
      match[parity/2] = line.split('/')[-1].split(':')[0]
      match[parity/2+3] = line.split('/')[-2][1:]
      if (match[parity/2+3] not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"]):
        print line
        print file
        assert(False)
      # Optionally change name to:  subgroup:protein
      match[parity/2] = match[parity/2+3] + ':' + match[parity/2]
    if parity == 4:
      assert line.startswith('[')
      if matchType == 'tetrahedra':
          match[2] = int(line.split(':')[1].split(',')[0])
      elif matchType == 'blosum62':
          match[2] = float(line.split(':')[1].split(',')[1])
      else:
          sys.exit("<matchType> must be 'tetrahedra' (default) or 'blusomu62'")
      match[5] = align(line)
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
assert len(proteins) == len(set(proteins))

print "len(proteins): ", len(proteins)
# If this fails, it's because the files of matches have some duplicates:
if len(proteins)**2 - len(proteins) > 2*len(matches):
  print "**** WARNING:  Some protein matches missing"
assert len(proteins)**2 - len(proteins) >= 2*len(matches)

# ====
# Now sort the largest alignments first.
# Use heuristic similar to TOPOMAX??
matches.sort(None, lambda match: match[2], True)

# Data structures (... means zero or more like previous example):
#  alignedProteins = [prot1 ...]
#  alignedProteinGroups = [alignedProteinGroup1  ...]
#    alignedProteinGroup = [prot1 ...]
#  alignedResidueGroups = [alignedResidueGroup1 ...]
#    alignedResidueGroup = [row1 ...]
#      row = [prot1 ...]
#        where res1 is of form "prot:annotatedResidue" and "prot" is a
#        protein in similar position in alignedProteinGroups
alignedProteins = []
alignedProteinGroups = []
alignedResidueGroups = []
conflictResOccurred = False
for match in [m for m in matches if m[2] >= minMatchSizeUsed]:
  # Guarantee that alignedProteinGroups has no duplicates:
  print "DEBUG:  alignedProteins:", alignedProteins
  print "DEBUG:  alignedProteinGroups:", alignedProteinGroups
  tmp = reduce(lambda x,y: x+y, alignedProteinGroups, [])
  assert len(tmp) == len(set(tmp))
  print "DEBUG: starting new match: ", match
  alignedProteinGroup0 = None
  alignedProteinGroup1 = None
  alignedResidueGroup0 = None
  alignedResidueGroup1 = None
  if match[0] in alignedProteins and match[1] in alignedProteins:
    for a in alignedProteinGroups:
      if match[0] in a:
        alignedProteinGroup0 = a
        alignedResidueGroup0=alignedResidueGroups[alignedProteinGroups.index(a)]
      if match[1] in a:
        alignedProteinGroup1 = a
        alignedResidueGroup1=alignedResidueGroups[alignedProteinGroups.index(a)]
  if match[0] in alignedProteins and match[1] in alignedProteins and \
     alignedProteinGroup0 != alignedProteinGroup1:
    assert(match[0] in alignedProteinGroup0 and
           match[1] in alignedProteinGroup1)
    # match[0] and match[1] are in different groups
    resPair = match[5][0]
    res0 = match[0]+':'+resPair[0]
    res1 = match[1]+':'+resPair[1]
    print "**** Merging two groups:", \
      res0.rsplit(':',1)[0], res1.rsplit(':',1)[0]
    if len(alignedResidueGroup0) + len(alignedResidueGroup1) \
       > maxColsPerMergedGroup:
      print "****  Will not merge; Total size would be", \
        str(len(alignedResidueGroup0)), '+', str(len(alignedResidueGroup1)), \
        '>', str(maxColsPerMergedGroup)
      continue
    alignedProteinGroup0.extend(alignedProteinGroup1)
    del alignedProteinGroups[alignedProteinGroups.index(alignedProteinGroup1)]
    alignedResidueGroup0.extend(alignedResidueGroup1)
    del alignedResidueGroups[alignedResidueGroups.index(alignedResidueGroup1)]
    alignedProteinGroup1 = alignedProteinGroup0
    alignedResidueGroup1 = alignedResidueGroup0
    assert(match[0] in alignedProteinGroup0 and
           match[1] in alignedProteinGroup0)
    # Done combining groups:  Will fall through to case of extending alignment
  if match[0] in alignedProteins and match[1] in alignedProteins:
    # match[0] and match[1] are in same protein group
    assert alignedProteinGroup0 and alignedProteinGroup0==alignedProteinGroup1
    print "Already aligned: ", match[0], ',', match[1], "(extending alignment)"
    # Will detect and print if there is a conflicting alignment
    conflictRes = []
    for resPair in match[5]:
      matchRes0 = match[0]+':'+resPair[0]
      matchRes1 = match[1]+':'+resPair[1]
      col0 = col1 = []
      for col in alignedResidueGroup0:
        if matchRes0 in col:
          col0 = col
        if matchRes1 in col:
          col1 = col
      if col0 and col0 == col1:
        continue  # Do nothing; col0 already has match.

      # Record any conflicts (residues from same protein in same column).
      colProteins = [elt.rsplit(':',1)[0] for elt in col0+col1]
      if not col0 and match[0] in colProteins: # match[0] protein already found
        conflictRes.append((matchRes0,'-->',col1))
        continue  # Conflict:  skip this column
      if not col1 and match[1] in colProteins: # match[1] protein already found
        conflictRes.append((matchRes1,'-->',col0))
        continue  # Conflict:  skip this column
      if len(colProteins) != len(set(colProteins)):  # if duplicate proteins
        conflictRes.append(((matchRes0,matchRes1),': ',col0,col1))
        continue  # Conflict:  skip this column

      # If col0 and col0 == col1, do nothing.  col0 already contains match.
      if col0 and not col1:
        col0.append(matchRes1)
      elif col1 and not col0:
        col1.append(matchRes0)
      elif not col0 and not col1:  # Else create a new column
        print "  DEBUG:  Creating a new column:", [matchRes0, matchRes1]
        alignedResidueGroup0.append([matchRes0,matchRes1])
      elif col0 != col1:
        assert col0 and col1
        # Test if col0 and col1 have a protein in common
        # No proteins in common, from above.  It's safe to merge the columns.
        col0.extend(col1)
        del alignedResidueGroup0[alignedResidueGroup0.index(col1)]
      else:
        assert False, "Internal error"
    if conflictRes:
      conflictResOccurred = True
      print "  \n  **********************************************"
      print "  **** Conflicting alignment: ", conflictRes
      print "  **********************************************\n"
    continue # Finished match for which match[0] and match[1]
             #   in same protein group.  Continue to next match
  assert match[0] not in alignedProteins or match[1] not in alignedProteins
  if match[0] not in alignedProteins and match[1] not in alignedProteins:
    print "Protein", match[0], "is beginning of new group"
    alignedProteins.append(match[0])
    alignedProteinGroups.append([match[0]])
    alignedResidueGroups.append(
                          [[match[0]+':'+resPair[0]] for resPair in match[5]] )
    # Fall through to next case, below.
  # Now, one protein is in data structure, and one is not.
  assert (match[0] in alignedProteins) ^ (match[1] in alignedProteins)  # xor
  if match[0] in alignedProteins:
    consensusProtein = 0
  if match[1] in alignedProteins:
    consensusProtein = 1
  print "DEBUG:  ", match[consensusProtein], match[1-consensusProtein]
  alignedResidueGroup = None
  for a in alignedProteinGroups:
    if match[consensusProtein] in a:
      alignedProteins.append(match[1-consensusProtein])
      a.append(match[1-consensusProtein])
      alignedResidueGroup = alignedResidueGroups[alignedProteinGroups.index(a)]
  for resPair in match[5]:
    consensusRes = match[consensusProtein]+':'+resPair[consensusProtein]
    for resMatch in alignedResidueGroup:
      if consensusRes in resMatch:
        resMatch.append(
                    match[1-consensusProtein]+':'+resPair[1-consensusProtein] )

print "DEBUG: alignedProteins:", alignedProteins
print "DEBUG: alignedProteinGroups:", alignedProteinGroups

# Sort so that aligned residues will be in order of increasing sequence number
for group in alignedResidueGroups:
  if displayRank:
    group.sort(None, lambda x: int(x[0].split('/')[-1][1:].rstrip('*')))
  else:
    group.sort(None, lambda x: int(x[0].split(':')[-1][1:].rstrip('*')))

# reorderCols:  Input:  cols of a protein group; Output:  reorder cols
# Algorithm:  Reorder cols so that sequence numbers increasing in first row,
#    and if col has no first row, then find a non-gap entry of col, and a prev
#    col that has both the non-gap entry and an earlier row with a gap.  Then
#    the difference in sequence numbers for the non-gap entry between the
#    current col and previous col allows one to calibrate and compute a
#    presumed sequence number for the earlier gap row of the current col.
# reorderTheCol:  When col and prev. col found as above, reorder
def seqNumOfResidue(residue):
  if '/' in residue:
    residue = residue.split('/')[1]
  # residue now of form: resType (one letter) + seq num
  return int(residue[1:].rstrip('*'))
def moveColEarlier(alignedRows, origCol, col):
  assert col < origCol
  for row in alignedRows:
    colValue = row[origCol]
    del row[origCol]
    row.insert(col, colValue)
def reorderColsByFirstResidue(rowsOfAlignedProteinGroup):
  row = rowsOfAlignedProteinGroup[0]
  lastSeqNum = 0
  for i in range(1, len(row)):  # row[0] would be protein name
    if row[i] != '--':
      if lastSeqNum > seqNumOfResidue(row[i]):
        for j in range(1, i-1):
          if row[j] != '--' and \
             seqNumOfResidue(row[j]) > seqNumOfResidue(row[i]):
            moveColEarlier(rowsOfAlignedProteinGroup, i, j)
            break
      else:
        lastSeqNum = seqNumOfResidue(row[i])
def reorderTheCol(alignedRows, col, prevCol, gapRowMatch):
  assert gapRowMatch
  gapRowMatch = []  # Start over, in case caller wasn't complete.
  for row in alignedRows:
    if row[col] == '--' and row[prevCol] != '--':
      gapRowMatch.append(row)
    elif not gapRowMatch and row[col] != '--' and row[prevCol] == '--':
      pass  # Maybe caller couldn't find non-gap prevCol for this row.
    elif row[col] != '--' and row[prevCol] != '--':
      if seqNumOfResidue(row[col]) < seqNumOfResidue(row[prevCol]):
        # Caller will try to reorder again from the top after this.
        moveColEarlier(alignedRows, col, prevCol)
        return True
      elif gapRowMatch:
        # Have prevCol w/ residues in gapRow[prevCol] and row[prevCol], row[col]
        # Find i with prevCol < i < col, w/ residue in gapRow[i], and test it.
        knownIncr = seqNumOfResidue(row[col]) - seqNumOfResidue(row[prevCol])
        for gapRow in gapRowMatch:
          for i in range(prevCol+1, col):
            if gapRow[i] != '--' and seqNumOfResidue(gapRow[i]) > \
                                 seqNumOfResidue(gapRow[prevCol]) + knownIncr \
                 and (row[i] == '--'
                      or seqNumOfResidue(row[i]) >= seqNumOfResidue(row[col])):
              moveColEarlier(alignedRows, col, i)
              return True
      return False
    else:
      assert row[col] == '--' or gapRowMatch # Give up if this fails.
  return False
def reorderCols(alignedRows, alignedProteinGroup):
  proteinOrder = {}
  for i in range(len(alignedProteinGroup)): # dict to remember protein ordering
    proteinOrder[alignedProteinGroup[i]] = i
  for row in alignedRows:
    assert len(row) == len(alignedRows[0])
  needsReordering = True
  while needsReordering:
    needsReordering = False
    for col in range(1, len(alignedRows[0])):
      skipCol = False
      if alignedRows[0][col] == '--':
        for prevCol in range(col-1, 1, -1):
          gapRowMatch = []
          for row in alignedRows:
            if row[col] == '--' and row[prevCol] != '--':
              gapRowMatch.append(row)
            elif row[col] != '--' and row[prevCol] == '--':
              break # Wait until we find a row w/ non-gaps for col and prevCol
            elif row[col] != '--' and row[prevCol] != '--':
              if seqNumOfResidue(row[col]) < seqNumOfResidue(row[prevCol]):
                moveColEarlier(alignedRows, col, prevCol)
                needsReordering = True
                # print "DEBUG:  REORDERED", row[prevCol:col+1]
                break
              if gapRowMatch:
                assert seqNumOfResidue(row[col]) > seqNumOfResidue(row[prevCol])
                if reorderTheCol(alignedRows, col, prevCol, gapRowMatch):
                  needsReordering = True
                  # print "DEBUG:  REORDERED", row[prevCol:col+1]
                  break
                else:
                  break
              else:
                assert not gapRowMatch
                skipCol = True
                break
          if needsReordering or skipCol:
            break # break out of "for prevCol in ..."
      if needsReordering:
        break  # break out of "for col in ..."

csvFile = open('multipleAlignmentOct07.csv', 'w')
csvFile.write(',"Position of aligned residues",,,,,,,,,,,,,,,,,,,,,,,' + '\n')
csvFile.write(
  '"PDB ID",' + ','.join(map(str, range(1,maxColsPerMergedGroup+1))) + '\n')
for i in range(len(alignedProteinGroups)):
  alignedProteinGroup = alignedProteinGroups[i]
  alignedResidueGroup = alignedResidueGroups[i]
  rows = [[protein] for protein in alignedProteinGroup]
  # A residue is:  grp:protein:resSeqNum[*] ; e.g., 3:1XO6:E12*
  
  for residues in alignedResidueGroup:
    resEntries = [res.rsplit(':', 1) for res in residues]
    rowLen = len(rows[0])
    # TODO:
    # MUST SORT resEntries ACCORDING TO SEQ NUM FOR FIRST PROTEIN PRESENT
    # Find first protein in rows also in resEntries[i], and sequence number
    #  for that entry in resEntries[i] is the key on which to sort.
    for res in resEntries:
      for row in rows:
        if res[0] == row[0]:
          row.append(res[1])
    for row in rows:
      if len(row) == rowLen:
        row.append('--')
    for row in rows:
      assert len(rows[0]) == len(rows[1]), \
             ("TWO RESIDUES FROM SAME PROTEIN?", map(len,rows), residues)

  rowsOfAlignedProteinGroup = []
  for protein in alignedProteinGroups[i]: # Retain proteins in order matched
    for row in rows:
      if row[0] == protein:
        rowsOfAlignedProteinGroup.append(row)
  reorderColsByFirstResidue(rowsOfAlignedProteinGroup)
  reorderCols(rowsOfAlignedProteinGroup, alignedProteinGroup)
  lastProteinGroup = None
  for protein in alignedProteinGroups[i]: # Retain proteins in order matched
    for row in rows:
      if row[0] == protein:
        if lastProteinGroup and row[0].split(':')[0] != lastProteinGroup:
          # alignedProteinGroups[i] could have aligned proteins from
          #  more than one group, as determined by hand.  Add small separator.
          csvFile.write('"---"\n')
        csvFile.write('"' + '","'.join(row) + '"\n')
        lastProteinGroup = row[0].split(':')[0]
  csvFile.write('"--- --- ---"\n')
csvFile.close()

if conflictResOccurred:
  print "  \n  ****************************************************************"
  print "  **** NOTE:  Some conflicts occurred in generating the alignment."
  print "  ****************************************************************\n"

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
    # assert [m for m in matches if m[0] == g_matches[1][1] if m[1] == g_matches[1][0]]
    # assert [m for m in g_matches if m[0] == g_matches[1][1] if m[1] == g_matches[1][0]]
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
