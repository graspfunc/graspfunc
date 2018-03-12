#!/usr/bin/python

# For generating *.qhull files from *.pdb files for all proteins.

#   For a protein, a qhull file is a list of residues with their
# coordinates. It is required for calculating the lengths of sides,
# areas of tetrahedra, etc., and is also used as an input to qhull
# for generating the protein's *.sc file.

# USAGE (for standalone):
#  $ python __FILE__.py <protein.pdb> <protein.rank>

def generate_qhull_file(pdbFile, rankFile):
  pdbData = []
  with open(pdbFile,"r") as f:
    lines = [elem for elem in f.read().split('\n') if elem]
  for line in lines:
    pdbData.append(line.split())
  
  rank = []
  with open(rankFile, "r") as g:
    lines = [elem for elem in g.read().split('\n') if elem]
  for line in lines:
    rank.append(line.split())
  
  h_name = pdbFile[:-4]+'.qhull'
  h = open(h_name, "w")
  strp = "3\n"+str(len(lines))+"\n"
  h.write(strp)
  tempJ = []
  for i in range(len(lines)):
    try:
      a = rank[i][1]
    except Exception:
      import pdb
      pdb.set_trace()
  #   print a
    leng = 0
    found = False
    for j in pdbData:
      leng += 1
      if len(j)>6:
        #if j[0] == 'ATOM' and j[2] == 'CA' and a in j[4]: # for g1/1ug9, g8/3AJX
        #if j[0] == 'ATOM' and j[2] == 'CA' and a == j[5]: # for other proteins
        if (j[0] == 'ATOM' and j[2] == 'CA' and j[4] == 'A' + str(a)):
          tempJ = j
          found = True
          h.write('{0:>6}  {1:>6}  {2:>6}\n'.format(j[5], j[6], j[7])) # for g8/3AJX
        elif (j[0] == 'ATOM' and j[2] == 'CA' and j[4] == 'A' and j[5] == a): # for g9/3eu8
          tempJ = j
          found = True
          h.write('{0:>6}  {1:>6}  {2:>6}\n'.format(j[6], j[7], j[8])) # for other proteins
          #pdb_rank.append(j)
          #if len(tempJ) > 0:
            #assert(not(tempJ[5] == j[5])) # verify no duplicates
            # if (tempJ[5] == j[5]): # verify no duplicates; this is true for g6/1PII: residue #370 (ALYS and BLYS), and residue #324 (APHE and BPHE)
            #     import pdb
            #     pdb.set_trace()
          #h.write('{0:>6}  {1:>6}  {2:>6}\n'.format(j[5], j[6], j[7])) # for g8/3AJX
    #assert(found)  # for each residue we should write a line to the qhull file
    if not found:
      print "Warning: Residue #" + str(a) + " (at pos. " + str(i) + ") Not found!"
      #import pdb
      #pdb.set_trace()
  h.close()

# if this is being run standalone
if __name__ == '__main__':
    import sys
    generate_qhull_file(sys.argv[1], sys.argv[2])
