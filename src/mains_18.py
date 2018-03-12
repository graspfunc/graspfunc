#!/usr/bin/python

try:
    import networkx as nx
    # import matplotlib.pyplot as plt
    # import pygraphviz
    import sys
    import math
    import signal
    import os
    import re
    import pdb
    import pickle
    import itertools
except:
    raise

#####################################################################
#  Heuristics used:  (Please keep ths up to date.)
# 1. Triangulation:
#    A. Triangulation of the entire protein is done
#    B. Only the num-protein-residues most highly ranked residues are kept.
#         50 is recommended as a good choice.
#    C. Any vertices adjacent to those highly ranked residues are added back in.
#    D. Only those tetrahedra are kept whose vertices are a subset of the above.
# 2. The tetrahedra are used as seed pairs if each has a vertex contained
#      among the top-n-seed-pairs most highly ranked residues.
#      3 or 4 is recommended as a good choice.
# 3. There is an optional heuristic to insist that one tetrahedra of the pair
#      contains the most highly ranked vertex.  This is commented out here.
# 4. A seed pair of triangulations is rejected if the parity does not match.
# 5. We should replace the residue comparison test.  Reject if:
#      |r1-r2| > 20 for r1 or r2 <= 24; and |r1-r2| > 10 for r1 or r2 < 10;
#      and |r1-r2| > 2 for r1 or r2 < 3
#      Worse test (commented out) rejects if:  max(tmp1,tmp2) > 2*min(tmp1,tmp2)
#          where tmp1 and tmp2 are the larger of vertex rank and 10.
#      This test now added to get_blosum_score (bonus/penalty), but no reject
# 6. If vertex 1,3,3,4,5 is present (5 most highly ranked residues), then the
#      chemical similarity types must match, based on aminoAcidSimilarityList3.
#      (See NOTE A, below.)
# 7. If the volume of a tetrahedron doesn't match, reject.
#      Using hardwired value of 14.4.  Should be parameter.
# 8. If the sum of lengths of edges of a tetrahedron doesn't match, reject.
#      Using hardwired value of 9.6.  Should be parameter.
# 9. Search for best BLOSUM62 score (with normalization) + 5 for each residue
#      of rank 1 matched, 3 for rank 2, 1 for rank 3.  (NOTE: matching a rank 1
#      residue from each protein is worth 10.)  This should be merged into
#      get_blosum_match_score() routine.
# 10. If depth > 10, and Blosum score 10 less than best, give up (heuristic-8)
# 11. Setting top-n-seed-pairs = 4 on cmd line and using num-protein-residues=25
#      seems better for this version. (See NOTE B, below.)
# 12. If the distance between the first vertex and the latest vertex for
#      protein1 differs by more than 5 from that for the matching vertices
#      of protein2, then reject those matching vertices.  Also apply it
#      between the middle vertex in visitedVerticesStack1 and latest vertex.
#      [ Two residues on backbone have typical distance 3.7.  The choice of 5,
#        is somewhat arbitrary.  Ideally, we would also consider distances
#        between other pairs of residues.  But this would be inefficient. ]
# FUTURE:
# A. We should reject a tetrahedron if the matching faces don't agree
#      on whether it is an exterior face of the protein.
# B. Lots of opportunities to speed up the code!
# C. For other protein families, should scale num-protein-residues with
#      size of protein.
# D. When two highly ranked residues of a protein are matched, see if original
#      geometric distances between residues in protein1 and protein2 agree.
#      If they disagree, reject matching both residues.
# NOTE:
# A. POOL matchings have found chemical property almost always conserved in
#      multiple alignment.  This is a key to clustering proteins successfully.
# B. Reducing num-protein-residues, or allowing larger wait-time, is useful
#      to avoid timeouts.  A timeout will stop searching a given branch,
#      and often miss a good match, because it was spending too much time
#      searching a smaller sub-branch of that same branch.
#####################################################################

# g: Vertex -> {Tetrahedron Index: Tetrahedon includes Vertex and Index in list of allTetrahedra}
g_1_dict = []
g_2_dict = []

stopped_once = 0

# Maximum wait time for the search function to return
waitMax = 0

# TODO:  Include geometric information
#   Tetrahedra match only if volumes, sumoflengths, thickThin (nbrs on backbone), and exterior face all match
#   volume of tetrahedra is determinant of 3 edges viewed as vectors;  Arbitrary epsilon parameter for rejecting mismatched volue, sumoflength
# We need to check that highly ranked residues correspond (POOL method); This is a big difference from Topofit

# Indexed array of residues;  index is reside # for protein1; value is matching residue # for protein2
vertex1Match = [] # initialized later; residues are numbered from 1, with vertex1Match[0] unused
# Stack of residue #'s from protein1 that have already been matched.
matchedIndices = []
visitedVertices2 = [] # initialized later
visitedVerticesStack1 = []

# RENAME THIS TO g_numTetrahedraVisited TO MAKE CLEAR IT'S GLOBAL.
numTetrahedraVisited = 0

bestNumTetrahedraMatched = 0
bestVisitedVerticesStack1 = []
bestVertex1Match = [] # This should be called bestVertex1Match

visitedTetraStack_1 = [] # list of matching tetrahedra
visitedTetraStack_2 = []

bestMatchingTetra_1 = []
bestMatchingTetra_2 = []

allTetrahedra_1 = []
allTetrahedra_2 = []
allTetrahedraIndex_1 = {}
allTetrahedraIndex_2 = {}
allTetrahedraAsLists_1 = []
allTetrahedraAsLists_2 = []
g_1_face_adj_vertices = {}
g_2_face_adj_vertices = {}

residueNum_1 = [] # list of residue numbers from the rank file
residueNum_2 = []

residueNames_1 = [] # list of residue names from the rank file
residueNames_2 = []

volList_1 = []  # list of volumes of tetrahedra; volList[i] corresponds to the volume(allTetrahedra[i])
volList_2 = []

sideLengthsList_1 = [] # list of sum of side lengths of tetrahedra
sideLengthsList_2 = []

# residueCoords[i] --> 3-D coordinates of the i-th residue from the rankings file
residueCoords_1 = []
residueCoords_2 = []

# The BLOSUM62 matrix
blosumMatrix = []

totalBlosumScore = 0 # blosum62 score for a particular path
bestTotalBlosumScore = 0 # keeps track of the lowest blosum62 score

# See http://en.wikipedia.org/wiki/Amino_acid for explanation (incl. B/Z/X/*)
aminoAcidList = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
"GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
"THR", "TRP", "TYR", "VAL", "B", "Z", "X", "*"]
# Traditionally, sometimes "J" added after Z, or B/Z/J/X omitted
# These are the letters used in BLOSUM63:
#  B = N/D (ASN/ASP); Z = E/Q (GLU/GLN); J = I/L (ILE/LEU); X = unkonwn; *=gap
aminoAcid1LetterList = list("ARNDCQEGHILKMFPSTWYVBZX*")
# Similarity groups for values of 1 to 3, according to CSM.txt:
aminoAcidSimilarityList1 = ["A", "RHK", "NDQE", "C", "G", "ILV", "M", "FWY", "P", "ST"]
# NOTE:  N/E have a match score of only 1, but all others in NDQE have
#        a match score of at least 2.
aminoAcidSimilarityList2 = ["A", "RK", "NDQE", "C", "G", "H", "ILV", "M", "FWY", "P", "ST"]
aminoAcidSimilarityList3 = ["A", "R", "NQ", "DE", "C", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

def annotateResidue(vertex, residueNames, residueNum):
    name = aminoAcid1LetterList[aminoAcidList.index(residueNames[vertex])]
    num = residueNum[vertex]
    return str(vertex) + '/' + name + '/' + str(num)
def annotateResidues(vertices, residueNames, residueNum):
    return "[" + ", ".join([annotateResidue(v, residueNames, residueNum)
                            for v in vertices]) + "]"
# NOTE:  Order of args matters.  Protein1 must come first.
def residuesAreSimilar(vertex_1, vertex_2):
    global residueNames_1, residueNames_2
    global aminoAcidSimilarityList3
    similarityList = aminoAcidSimilarityList3
    a1 = aminoAcid1LetterList[ aminoAcidList.index(residueNames_1[vertex_1]) ]
    a2 = aminoAcid1LetterList[ aminoAcidList.index(residueNames_2[vertex_2]) ]
    sim1 = [s for s in similarityList if a1 in s]
    sim2 = [s for s in similarityList if a2 in s]
    return sim1 == sim2

def read_blosum(filename):
    blosum62 = []

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        # ignore all the comments
        if line[0] == "#":
            continue
        data = re.split("\s+", line)
        # ignore if it's the first line of the matrix
        if data[1] == "A":
            continue
        blosum62.append([int(data[v]) for v in range(1, len(data) - 1)])
        #print len(data), data[24]
    return blosum62

def normalize_blosum(blosumMatrix, normalizingFactor):
    rows = len(blosumMatrix)
    cols = len(blosumMatrix)

    assert(rows == cols)

    for r in range(rows):
        for c in range(cols):
            blosumMatrix[r][c] = float(blosumMatrix[r][c] / normalizingFactor)

    return blosumMatrix

def sum_of_sides(v1, v2, v3, v4):
    vec1 = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]] # v1 - v2
    vec2 = [v2[0] - v3[0], v2[1] - v3[1], v2[2] - v3[2]] # v2 - v3
    vec3 = [v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]] # v3 - v1
    vec4 = [v4[0] - v1[0], v4[1] - v1[1], v4[2] - v1[2]] # v4 - v1
    vec5 = [v4[0] - v2[0], v4[1] - v2[1], v4[2] - v2[2]] # v4 - v2
    vec6 = [v4[0] - v3[0], v4[1] - v3[1], v4[2] - v3[2]] # v4 - v3

    len1 = math.sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]) # |vec1|
    len2 = math.sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]) # |vec2|
    len3 = math.sqrt(vec3[0]*vec3[0] + vec3[1]*vec3[1] + vec3[2]*vec3[2]) # |vec3|
    len4 = math.sqrt(vec4[0]*vec4[0] + vec4[1]*vec4[1] + vec4[2]*vec4[2]) # |vec4|
    len5 = math.sqrt(vec5[0]*vec5[0] + vec5[1]*vec5[1] + vec5[2]*vec5[2]) # |vec5|
    len6 = math.sqrt(vec6[0]*vec6[0] + vec6[1]*vec6[1] + vec6[2]*vec6[2]) # |vec6|

    sumOfLengths = len1 + len2 + len3 + len4 + len5 + len6

    return sumOfLengths

def calc_tetra_signed_vol(v1, v2, v3, v4):
    vec1 = [v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]] # v1 - v2
    vec2 = [v2[0] - v3[0], v2[1] - v3[1], v2[2] - v3[2]] # v2 - v3
    vec3 = [v3[0] - v4[0], v3[1] - v4[1], v3[2] - v4[2]] # v3 - v4

    # (v1 - v2).((v2 - v3)x(v3 - v4))
    det = vec1[0]*vec2[1]*vec3[2] + vec1[1]*vec2[2]*vec3[0] + vec1[2]*vec2[0]*vec3[1] - vec1[2]*vec2[1]*vec3[0] - vec1[1]*vec2[0]*vec3[2] - vec1[0]*vec2[2]*vec3[1]

    # Divide by 6 because this is a tetrahedron, not a 3-D parallelepiped
    return det/6.0

def distance(v1, v2):
    return math.sqrt(sum(map(lambda x,y: (x-y)**2, v1, v2)))

def distance_1(v1, v2):
    return distance(residueCoords_1[v1], residueCoords_1[v2])

def distance_2(v1, v2):
    return distance(residueCoords_2[v1], residueCoords_2[v2])

def get_blosum_match_score(blosumMatrix, v1, v2):
    poolScore = 0
    if residuesAreSimilar(v1, v2):
        if max(v1,v2) < 10 and residueNames_1[v1] == residueNames_2[v2]:
            poolScore += 15 / normalizingFactor
        elif max(v1,v2) < 10:
            poolScore += 5 / normalizingFactor
        elif max(v1,v2) < 20:
            poolScore += 5 / normalizingFactor
        elif max(v1,v2) < 30:
            poolScore += 2 / normalizingFactor
    else:  # else residues are not similar
        # RG: Change this so that we keep the matched pair,
        #       even if both the highly ranked residues differ
        if max(v1,v2) < 25:
            poolScore += 2 / normalizingFactor
        # Discourage matching of a pair of dissimilar residues
        #   only one of which is highly ranked
        elif min(v1, v2) < 25:
            # NOTE: This also matches when v1, v2 < 25 but
            #         that is caught by the if clause above
            poolScore -= 1 / normalizingFactor
        else:
            poolScore -= 5 / normalizingFactor # This is really bad.
    # poolScore *= 3
    #ifdef heuristic-1
    if v1 <= 24 or v2 <= 24 and abs(v1 - v2) > 20:
            poolScore -= 8 / normalizingFactor
    elif v1 < 10 or v2 < 10 and abs(v1 - v2) > 10:
            poolScore -= 8 / normalizingFactor
    elif v1 < 3 or v2 < 3 and abs(v1 - v2) > 1:
            poolScore -= 8 / normalizingFactor
    else:
        poolScore += 5 / normalizingFactor
    #endif
    # SHOULD MOVE THIS LOGIC INTO ABOVE
    if v1 <= 3:
        poolScore += 12 * ([999999999, 5, 3, 1])[v1] / normalizingFactor
    if v2 <= 3:
        poolScore += 12 * ([999999999, 5, 3, 1])[v2] / normalizingFactor

    blosumScore = 0
    aminoAcid_1 = residueNames_1[v1]
    aminoAcid_2 = residueNames_2[v2]
    row = aminoAcidList.index(aminoAcid_1)
    col = aminoAcidList.index(aminoAcid_2)
    blosumScore = blosumMatrix[row][col] - 5.0 / normalizingFactor  # Add bias of 0.5
    return poolScore + blosumScore

def read_pickle_tetra(filen, g_dict, allTetrahedra, allTetrahedraIndex, allTetrahedraAsLists, g_face_adj_vertices, residueNum, residueNames, residueCoords, volList, sideLengthsList):

    rankFileName = filen + ".rank"
    dtFileName = filen + ".sc"
    qhullFileName = filen + ".qhull"
    #ifdef heuristic-9
    #dtFileName = filen + ".scnew"
    #qhullFileName = filen + ".qhullnew"
    #endif

    f = open(qhullFileName, 'r')
    lines = f.readlines()
    f.close()

    # get the coordinates of each residue from the qhull file
    residueCoords.append([])
    for i in range(0, len(lines) - 2):
        residueCoords.append([])
        numbers = re.split('\s+', lines[i+2].strip()) # line #2 onwards contain the vertex data
        residueCoords[i + 1].append(float(numbers[0]))
        residueCoords[i + 1].append(float(numbers[1]))
        residueCoords[i + 1].append(float(numbers[2]))

    #print dtFileName, rankFileName
    f = open(dtFileName, 'r')
    g = nx.Graph()
    lines = f.readlines()
    f.close()

    # get the tetrahedra from the sc file
    tetras = []
    for i in range(0, len(lines) - 2):
        tetras.append([])
        numbers = lines[i+2].split(' ') # line #2 onwards contain the tetrahedra data
        tetras[i].append(int(numbers[0]) + 1) # Add 1 as the vertices are [0, 199]
        tetras[i].append(int(numbers[1]) + 1)
        tetras[i].append(int(numbers[2]) + 1)
        tetras[i].append(int(numbers[3]) + 1)


    f = open(rankFileName, 'r')
    lines = f.readlines()
    f.close()

    # get the residue names and numbers from the rank file
    idx = 1
    residueNum.append(0)
    residueNames.append("")
    for l in lines:
        data = l.split(' ')
        residueNum.append(int(data[0]))
        residueNames.append(data[1][0:3])
        idx += 1

    g_dict.append(list())
    for i in range(1, idx):
        g_dict.append(list())  #(g.number_of_nodes() + 1)

    for i in range(len(tetras)):
        r1 = tetras[i][0]
        r2 = tetras[i][1]
        r3 = tetras[i][2]
        r4 = tetras[i][3]

        # No residue should be 0
        assert(r1 != 0 and r2 != 0 and r3 != 0 and r4 != 0)

        #strp = str(r1) + " " + str(r2) + " " + str(r3) + " " + str(r4)
        #print strp

        # add edges in the graph for each vertex pair in the tetrahedron
        g.add_edge(r1, r2)
        g.add_edge(r2, r3)
        g.add_edge(r3, r1)
        g.add_edge(r4, r1)
        g.add_edge(r4, r2)
        g.add_edge(r4, r3)

        # add the tetrahedron the list of all tetrahedra
        s = set([r1, r2, r3, r4])
        allTetrahedra.append(s)
        allTetrahedraIndex[frozenset(s)] = i
        for v in s:
            face = frozenset(s-set([v]))
            if g_face_adj_vertices.has_key(face):
                g_face_adj_vertices[face].append(v)
            else:
                g_face_adj_vertices[face] = [v]

        # add the corresponding volume and the sum of lengths to the global map
        signedVol = calc_tetra_signed_vol(residueCoords[r1], residueCoords[r2],
                                          residueCoords[r3], residueCoords[r4])
        sumOfLengths = sum_of_sides(residueCoords[r1], residueCoords[r2],
                                    residueCoords[r3], residueCoords[r4])

        # Record vertices of tetrahedron in ordering for positive orientation
        if signedVol > 0:
            allTetrahedraAsLists.append([r1, r2, r3, r4])
        else:
            allTetrahedraAsLists.append([r2, r1, r3, r4])
        assert( calc_tetra_signed_vol(*[residueCoords[r] for r in allTetrahedraAsLists[-1]]) > 0 )
        volList.append(abs(signedVol))
        sideLengthsList.append(sumOfLengths)
        # add the tetrahedron index to the Vertex x Tetrahedra map
        assert(i not in g_dict[r1]) # sanity check
        g_dict[r1].append(i)
        assert(i not in g_dict[r2])
        g_dict[r2].append(i)
        assert(i not in g_dict[r3])
        g_dict[r3].append(i)
        assert(i not in g_dict[r4])
        g_dict[r4].append(i)

    print filen + ": " + str(g.number_of_nodes()) + ", " + str(g.number_of_edges()) + ", " + str(nx.degree_histogram(g))
    return g

def compare_top(g, g_dict, allTetrahedra, top, filen):
    try:
        #print str(g.number_of_nodes()) + ", " + str(g.number_of_edges()) + ", " + str(nx.degree_histogram(g))
        # stores all the vertices that must not be deleted
        topVerticesSet = set()
        for i in g.nodes():  # defined in networkx
            if i <= top:
                for k in xrange(len(g_dict[i]) - 1, -1, -1): # for all the adjacent tetrahedra to i
                    tetrahedron = allTetrahedra[g_dict[i][k]]
                    for v in tetrahedron:
                        topVerticesSet.add(v)
        # delete the rest of the vertices and the corresponding tetrahedra
        for i in g.nodes():  # defined in networkx
            if i > top:
                if i in topVerticesSet:
                    continue
                g.remove_node(i)
                # iterate over all tetrahedra adjacent to a vertex
                for k in xrange(len(g_dict[i]) - 1, -1, -1): # for all the adjacent tetrahedra to i
                    # g_dict[i][k] --> Index of the k-th (k \in [0, N)) tetrahedron adjacent to i-th vertex
                    tetrahedron = allTetrahedra[g_dict[i][k]]
                    for v in tetrahedron:
                        if v == i:
                            continue
                        else:
                            # remove the tetrahedron index from g_dict[v]
                            g_dict[v].remove(g_dict[i][k])
                    # remove the tetrahedron index from g_dict[i]
                    g_dict[i].remove(g_dict[i][k])
    except:
        print "error: " + str(i)
        raise
    print filen + ": " + str(g.number_of_nodes()) + ", " + str(g.number_of_edges()) + ", " + str(nx.degree_histogram(g))
    return g

# Function is redundant because of the check in do_dfs_tetra()
def verify_matching():
    global visitedTetraStack_1
    global visitedTetraStack_2
    temp_s = list()
    strp = ""
    flag = 0
    for i in range(len(visitedTetraStack_1) - 0):
        if not(visitedTetraStack_1[i] in temp_s):
            temp_s.append(visitedTetraStack_1[i])
        else:
            flag = 1
            strp = strp + ", " + str(i)
    if flag == 1:
        print strp
        print len(temp_s)
        return False
    else:
        return True

def make_graph(matching_l, fname):
    g = nx.Graph()
    for tetra_s in matching_l:
        tetra_l = list(tetra_s)
        g.add_edge(tetra_l[0], tetra_l[1])
        g.add_edge(tetra_l[1], tetra_l[2])
        g.add_edge(tetra_l[2], tetra_l[0])
        g.add_edge(tetra_l[3], tetra_l[0])
        g.add_edge(tetra_l[3], tetra_l[1])
        g.add_edge(tetra_l[3], tetra_l[2])
    nx.write_dot(g, fname)
    print nx.degree_histogram(g)

def parityTetrahedron(tetra):
    assert(type(tetra) is list or type(tetra) is tuple)
    parity = True
    # Note that v_i != v_j below.  So '>' suffices for parity test.
    for j in range(len(tetra)):
        v_j = tetra[j]
        for v_i in tetra[:j]:
            if v_i > v_j:
                parity = not parity
    return parity

def handler(signum, frame):
    global waitMax
    print "Took more than " + str(waitMax) + " seconds."
    raise Exception("end of time")

def newTetrahedron(g, allTetrahedraIndex, g_face_adj_vertices, oldtetra, oldface):
    newvertex = None
    assert(1 <= len(g_face_adj_vertices[frozenset(oldface)]) <= 2)
    for v in g_face_adj_vertices[frozenset(oldface)]:
        # If v is higher than num-protein-residues and not a nbr, then skip it.
        if v not in oldtetra and g.has_node(v):
            newvertex = v
            newtetra = oldface | set([newvertex])
    if not newvertex:
        return (None, -1, 0)

    tetraIndex = allTetrahedraIndex[frozenset(newtetra)]
    return (newtetra, tetraIndex, newvertex)

def do_dfs_tetra(g_1, g_2, oldtetra_1, oldtetra_2, oldface_1, oldface_2, seedsList_1, seedsList_2):

    global vertex1Match, matchedIndices, visitedVertices2, numTetrahedraVisited, bestNumTetrahedraMatched, bestVertex1Match
    global visitedTetraStack_1, visitedTetraStack_2, bestMatchingTetra_1, bestMatchingTetra_2, visitedVerticesStack1
    global residueNum_1, residueNum_2, residueNames_1, residueNames_2, blosumMatrix, totalBlosumScore, bestTotalBlosumScore

    (newtetra_1, tetraIndex_1, newvertex_1) = \
      newTetrahedron(g_1, allTetrahedraIndex_1, g_1_face_adj_vertices, oldtetra_1, oldface_1)
    if newtetra_1 == None:
        return

    #ifdef heuristic-8
    # If Blosum score too far below maximum, give up.
    if len(visitedVerticesStack1) > 10 and \
       totalBlosumScore + 10 < bestTotalBlosumScore:
        return
    #endif

    (newtetra_2, tetraIndex_2, newvertex_2) = \
      newTetrahedron(g_2, allTetrahedraIndex_2, g_2_face_adj_vertices, oldtetra_2, oldface_2)
    if newtetra_2 == None:
        return

    #ifdef heuristic-7
    # If matching faces disagree on whether they are exterior faces, then return
    # NOT CURREBTLY USED -- (NOT YET SHOWN USEFUL)
    # for v in oldface_1:
    #     newface_1 = frozenset(oldface_1 - set([v]) | set([newvertex_1]))
    #     newface_2 = frozenset([vertex1Match[v2] for v2 in oldface_1 - set([v])]
    #                           + [newvertex_2])
    #     if len(g_1_face_adj_vertices[newface_1]) != \
    #        len(g_2_face_adj_vertices[newface_2]):
    #         return
    #endif

    #ifdef heuristic-6
    assert(newvertex_1 and newvertex_2)
    # RG: if both the residues are highly ranked they should match
    if newvertex_1 < 11 or newvertex_2 < 11:
        if not residuesAreSimilar(newvertex_1, newvertex_2):
            return # require highly ranked match residues be similar
            # pass # Could use get_blosum_match_score instead
    #endif

    if newtetra_1 in visitedTetraStack_1 or newtetra_2 in visitedTetraStack_2:
        # This implies we must have come back a full circle
        return

    # Check that the new vertices match
    if vertex1Match[newvertex_1] and vertex1Match[newvertex_1] != newvertex_2:
        return
    if not vertex1Match[newvertex_1] and visitedVertices2[newvertex_2]:
        return


    #ifdef heuristic-1
    # tmp1 = max(newvertex_1, 10)
    # tmp2 = max(newvertex_2, 10)
    # if max(tmp1, tmp2) > 2 * min(tmp1, tmp2):
    #     return
    # Never return here.  Use get_blosum_match_score instead.
    if newvertex_1 <= 24 or newvertex_2 <= 24:
        if False and abs(newvertex_1 - newvertex_2) > 20:
            return
    if newvertex_1 < 10 or newvertex_2 < 10:
        if False and abs(newvertex_1 - newvertex_2) > 10:
            return
    if newvertex_1 < 3 or newvertex_2 < 3:
        if False and abs(newvertex_1 - newvertex_2) > 1:
            return
    #endif

    #ifdef heuristic-2
    #for v in oldface_1:
    #    isBackbone1 = abs(residueNum_1[v] - residueNum_1[newvertex_1]) == 1
    #    isBackbone2 = abs(residueNum_2[vertex1Match[v]] - residueNum_2[newvertex_2]) == 1
    #    if isBackbone1 != isBackbone2:
    #        return
    #endif

    #ifdef heuristic-3
    # Matching volumes, sum of side lengths
    if (volList_1[tetraIndex_1] - volList_2[tetraIndex_2]) > 14.4:
        #print "Rejecting based on Volume: " + str(volList_1[tetraIndex_1]) + ", " + str(volList_2[tetraIndex_2])
        return
    #endif
    #ifdef heuristic-4
    if abs(sideLengthsList_1[tetraIndex_1] - sideLengthsList_2[tetraIndex_2]) > 9.6:
        #print "Rejecting based on Sum of Side Lengths" + str(sideLengthsList_1[tetraIndex_1]) + ", " + str(sideLengthsList_2[tetraIndex_2])
        return
    #endif
    #ifdef heuristic-10
    # RG: Distances from all seed-pairs should match
    for (seed1, seed2) in zip(seedsList_1, seedsList_2):
      if abs(distance_1(newvertex_1, seed1)) - abs(distance_2(newvertex_2, seed2)) > 2:
        return
    #endif

    # =========================================================================
    # We are now committing the new tetrahedron and new vertex; adjust the score
    # newvertex_1 and newvertex_2 are residues to be matched when extending
    if not vertex1Match[newvertex_1]:
        addedNewVertex = newvertex_1
        vertex1Match[newvertex_1] = newvertex_2
        visitedVerticesStack1.append(newvertex_1)
        visitedVertices2[newvertex_2] = True
    else:
        addedNewVertex = None
    visitedTetraStack_1.append(newtetra_1)
    visitedTetraStack_2.append(newtetra_2)
    numTetrahedraVisited  += 1
    #print newtetra_1, newtetra_2, new_f1, new_f2
    assert(len(visitedTetraStack_1) == len(visitedTetraStack_2))
    assert(len(visitedTetraStack_2) == numTetrahedraVisited)

    #   blosum62 scoring
    matchScore = 0
    # HEURISTIC:  ADDED "True or" BECAUSE IT WORKS.  FIGURE OUT WHY, AND RE-WRITE.
    if addedNewVertex:
        totalBlosumScore += \
            get_blosum_match_score(blosumMatrix, newvertex_1, newvertex_2)

    # bestTotalScore = numTetrahedraVisited + totalBlosumScore

    # Record best match only after ensuring that this really is a new tetrahedron
    # bestScoreCriterion = (numTetrahedraVisited > bestNumTetrahedraMatched)
    bestScoreCriterion = (totalBlosumScore > bestTotalBlosumScore)
    if bestScoreCriterion:
        # bestVertex1Match = vertex1Match.deepCopy()
        bestNumTetrahedraMatched = numTetrahedraVisited
        bestTotalBlosumScore = totalBlosumScore
        bestVisitedVerticesStack1[:] = visitedVerticesStack1
        bestVertex1Match[:] = vertex1Match
        bestMatchingTetra_1[:] = visitedTetraStack_1
        bestMatchingTetra_2[:] = visitedTetraStack_2
        assert(set([v for l in visitedTetraStack_1 for v in l]) == set(visitedVerticesStack1))

    # =========================================================================
    # We will now search for a new tetrahedron and a new face to explore
    # Compute the faces of the new tetrahedra
    newTetraVertexList_1 = list(newtetra_1)
    newTetraVertexList_2 = list(newtetra_2)

    origfaces_1 = list([set([newTetraVertexList_1[0], newTetraVertexList_1[1],
                             newTetraVertexList_1[2]]), set([newTetraVertexList_1[0],
                                                           newTetraVertexList_1[1], newTetraVertexList_1[3]]),
                        set([newTetraVertexList_1[2], newTetraVertexList_1[1],
                             newTetraVertexList_1[3]]), set([newTetraVertexList_1[0],
                                                           newTetraVertexList_1[2], newTetraVertexList_1[3]])])
    # We only want to explore new faces
    if oldface_1 in origfaces_1:
        origfaces_1.remove(oldface_1)


    faces_2 = list([set([newTetraVertexList_2[0], newTetraVertexList_2[1],
                         newTetraVertexList_2[2]]), set([newTetraVertexList_2[0],
                                                       newTetraVertexList_2[1], newTetraVertexList_2[3]]),
                    set([newTetraVertexList_2[2], newTetraVertexList_2[1],
                         newTetraVertexList_2[3]]), set([newTetraVertexList_2[0],
                                                       newTetraVertexList_2[2], newTetraVertexList_2[3]])])
    # ISN'T THIS A DUPLICATE OF ABOVE CODE?  IF SO, DELETE.
    # We only want to explore new faces
    if oldface_2 in faces_2:
        faces_2.remove(oldface_2)

    # Next, try all possible new faces
    # try all sequences for origfaces_1
    # Convert generator to list for ease of debugging
    face_seqs_1 = list(itertools.permutations(origfaces_1))

    for faces_1 in face_seqs_1:
        # At end of trying all faces in faces_1,
        #   we might be able to skip other seqs
        verticesFoundInFaceSeqGraph1 = []
        verticesFoundInFaceSeqGraph2 = []

        origNumTetraVisited = numTetrahedraVisited
        origBlosumScore = totalBlosumScore
        numVisitedVertices1 = len(visitedVerticesStack1)
        #print visitedVerticesStack1
        global stopped_once
        for new_f1 in faces_1:
            if new_f1 == oldface_1:
                continue
            new_f2 = [vertex1Match[v] for v in new_f1]
            # ABOVE LINE SHOULD ALLOW US TO DIRECTLY CALL do_dfs_tetra
            # VERIFY THAT REMOVING "for new_f2" DOES NOT CHANGE ANSWER.
            for new_f2 in faces_2:
                if new_f2 == oldface_2:
                    continue

                commonEdge1 = new_f1 & oldface_1
                commonEdge2 = new_f2 & oldface_2
                if set([vertex1Match[v] for v in commonEdge1]) != set(commonEdge2):
                    continue

                # The recursive call will update bestTotalBlosumScore
                #   if this search branch (new_f1) improves on best so far.
                do_dfs_tetra(g_1, g_2, newtetra_1, newtetra_2, new_f1, new_f2, seedsList_1, seedsList_2)

        # Backtrack and try a different face sequence
        numTetrahedraVisited  = origNumTetraVisited
        totalBlosumScore = origBlosumScore

        # Optimization possible:  if all branches starting at one face
        #   always match a second face of the tetrahedron, then the order
        #   of trying the two faces doesn't matter.  (Similarly for 3 faces.)

	# BUG:  This finds only the vertices in a maximal match, instead of
	#   all vertices.  It will sometimes return early and not explore
	#   non-trivial face sequences.
        verticesFoundInFaceSeqGraph1 += \
            visitedVerticesStack1[numVisitedVertices1:]
        for v in visitedVerticesStack1[numVisitedVertices1:]:
            verticesFoundInFaceSeqGraph2.append(vertex1Match[v])

        for v in visitedVerticesStack1[numVisitedVertices1:]:
            visitedVertices2[vertex1Match[v]] = False
            vertex1Match[v] = 0  # 0 means it's unset
        del visitedVerticesStack1[numVisitedVertices1:]
        del visitedTetraStack_1[numTetrahedraVisited:]
        del visitedTetraStack_2[numTetrahedraVisited:]

        if len(verticesFoundInFaceSeqGraph1) == \
            len(set(verticesFoundInFaceSeqGraph1)) and \
            len(verticesFoundInFaceSeqGraph2) == \
            len(set(verticesFoundInFaceSeqGraph2)):
              # Vertices found from each face are distinct.
              # So, no need to try a different sequence
              break
    # print "END: visitedVerticesStack1:", visitedVerticesStack1
    # print "END: len(visitedVerticesStack1):", len(visitedVerticesStack1)

    # pops corresponding to the appends before
    if addedNewVertex:
        assert(addedNewVertex == visitedVerticesStack1[-1])
        visitedVerticesStack1.pop()
        visitedVertices2[vertex1Match[addedNewVertex]] = False
        vertex1Match[addedNewVertex] = 0 # 0 means it's unset
    visitedTetraStack_1.pop()
    visitedTetraStack_2.pop()
    numTetrahedraVisited -= 1
    return


# find a vertex which is not in the tetra_1 and find a vertex which is not in tetra_2
# if vertex is not visited, add it to the list; recurse on the new tetra_1' and tetra_2'
# if vertex is visited, add the tetra to matching;



#g_1 = nx.Graph()
#g_1.add_nodes_from([1, 2, 3, 4])
#g_1.add_edges_from([(1,2), (2,3), (3,1), (4,1), (4,2), (4,3), (5,1), (5,2), (5,4), (5,6), (5,7), (5,8), (6,7), (7,8), (8,6)])
#g_1.add_edges_from([(1, 2), (2, 3), (3, 1), (2, 4)])
#g_1.add_edges_from([(1,2), (2,3), (3,4), (4,1), (5,6), (6,7), (7,8), (8,5), (1, 5), (2, 6), (3, 7), (8, 4)])

#g_2 = nx.Graph()
#g_2.add_edges_from([(1,2), (2,3), (3,1), (4,1), (4,2), (4,3), (5,1), (5,2), (5,4), (5,6), (5,7), (5,8), (6,7), (7,8), (8,6)])
#g_2.add_edges_from([(1,2), (2,3), (3,4), (4,1), (5,6), (6,7), (7,8), (8,5), (1, 5), (2, 6), (3, 7), (8, 4)])
#g_2.add_nodes_from(['A', 'B', 'C', 'D'])
#g_2.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'A')])


# If this is a module, it needs special handling
if __name__ != '__main__':
    if '-c' in sys.argv:
      del sys.argv[sys.argv.index('-c')]
    print __name__ + ' being loaded as a module.'

# The parsing of the command line, expects sys.argv[0] to be this file.
if not sys.argv[0].endswith('.py'):
    sys.argv = ['-m '+__name__+'.so'] + sys.argv

if sys.argv[0] == 'profile.py' or '-waittorun' in sys.argv:
  waitToRun = True
  if '-waittorun' in sys.argv:
    del sys.argv[sys.argv.index('-waittorun')]
else:
  waitToRun = False


if len(sys.argv) < 6:
    sys.exit("Usage: python {0} </path/to/protein_1.sc> </path/to/protein_2.sc> <top-n-seed-pairs> <num-protein-residues> <wait-time>  [</path/to/blosum> <normalizing-factor>]".format(sys.argv[0]))
if not os.path.exists(sys.argv[1]+".sc"):
    sys.exit("File {0} not found!".format(sys.argv[1]+".sc"))
if not os.path.exists(sys.argv[2]+".sc"):
    sys.exit("File {0} not found!".format(sys.argv[2]+".sc"))
if not os.path.exists(sys.argv[6]):
    sys.exit("File {0} not found!".format(sys.argv[6]))

# read the file, create graph, and build dictionary
g_1 = read_pickle_tetra(sys.argv[1], g_1_dict, allTetrahedra_1, allTetrahedraIndex_1, allTetrahedraAsLists_1, g_1_face_adj_vertices, residueNum_1, residueNames_1, residueCoords_1, volList_1, sideLengthsList_1)
g_1 = compare_top(g_1, g_1_dict, allTetrahedra_1, int(sys.argv[4]), sys.argv[1])

g_2 = read_pickle_tetra(sys.argv[2], g_2_dict, allTetrahedra_2, allTetrahedraIndex_2, allTetrahedraAsLists_2, g_2_face_adj_vertices, residueNum_2, residueNames_2, residueCoords_2, volList_2, sideLengthsList_2)
g_2 = compare_top(g_2, g_2_dict, allTetrahedra_2, int(sys.argv[4]), sys.argv[2])

# Indexed array of residues;  index is reside # for protein1; value is matching residue # for protein2
vertex1Match = [0 for vi in range(len(residueNum_1))]  # residues are numbered from 1, with vertex1Match[0] unused
# Stack of residue #'s from protein1 that have already been matched.
matchedIndices = []
visitedVertices2 = [False for v in range(len(residueNum_2))]
visitedVerticesStack1 = []

top = int(sys.argv[3])

waitMax = int(sys.argv[5])

normalizingFactor = float(sys.argv[7])

# read blosum matrix from file; expects blosum62 format
blosumMatrix = read_blosum(sys.argv[6])

# normalize the blosum matrix
blosumMatrix = normalize_blosum(blosumMatrix, normalizingFactor)

maxReportedBlosumScore = -999999

# reset the data for a new seed pair
def reset_data():
  global bestTotalBlosumScore, maxReportedBlosumScore
  global bestNumTetrahedraMatched, maxReportedTetraMatched
  (maxReportedTetraMatched, bestNumTetrahedraMatched) = (0, 0)
  (bestTotalBlosumScore, maxReportedBlosumScore) = (0, -999999)

def print_data():
  global bestTotalBlosumScore, maxReportedBlosumScore
  global bestNumTetrahedraMatched, maxReportedTetraMatched
  print maxReportedTetraMatched, bestNumTetrahedraMatched
  print bestTotalBlosumScore, maxReportedBlosumScore

# Checks if a residuePair (\in residuePairsList) can be added to the resultsList
#   and returns the index in resultsList where it should be inserted
def check_residue_pair(residuePair, residuePairsList, resultsList, lengthsList, consMatch):
  (x_m, y_m) = residuePair
  i = 0

  # Residue-pair shouldn't contradict with the one in consMatch
  for (x, y) in consMatch:
    if (x_m == x and not(y_m == y)) or (not(x_m == x) and (y_m == y)):
      print "########## MISMATCH: ", residuePair, x, y
      return (i, False)

  # Iterate over the existing resultsList
  for (r1, r2) in resultsList:
    # No need to add if already added
    if r1 == x_m and r2 == y_m:
      return (i, False)
    # If one of the residues is common...
    elif (r1 == x_m) and not(r2 == y_m):
      if len(residuePairsList) <= lengthsList[i]:
        return (i, False)
      else:
        # ...add only when this list is bigger than the last one
        return (i, True)
    elif not(r1 == x_m) and (r2 == y_m):
      if len(residuePairsList) <= lengthsList[i]:
        return (i, False)
      else:
        return (i, True)
    else:
      i += 1
  return (i, True)

def calc_matches_union(allTetraMatches, consMatch):
  print "*************** Calculating the union ***************"
  unionList = []
  lengthsList = []
  for matchedPairsLists in allTetraMatches:
    # NOTE: There is only one list
    for matchedPairsList in matchedPairsLists:
      # Iterate over all the residue pairs in the list
      for matchedPair in matchedPairsList:
        # Check if the pair can be added, and at what position
        (idx, canAdd) = check_residue_pair(matchedPair, matchedPairsList,\
                                           unionList, lengthsList, consMatch)
        if canAdd:
          n = len(unionList)
          if not(n) or idx == (n):
            unionList.append(matchedPair)
            lengthsList.append(len(matchedPairsList))
          else:
            print "DEBUG:", idx, matchedPair, unionList, n
            unionList[idx] = matchedPair
            lengthsList[idx] = len(matchedPairsList)
  return unionList

def add_matching_2(g_1, g_2, rank1List, rank2List, consMatch):

    global vertex1Match, matchedIndices, visitedVertices2, numTetrahedraVisited, bestNumTetrahedraMatched, bestVertex1Match, bestVisitedVerticesStack1
    global visitedTetraStack_1, visitedTetraStack_2, bestMatchingTetra_1, bestMatchingTetra_2
    global totalBlosumScore, bestTotalBlosumScore, maxReportedBlosumScore

    allTetraMatches = []
    lastTetraMatches = []

    #print g_1.degree()
    #print nx.degree_histogram(g_2) #g_1.degree()
    nodes_1 = rank1List
    nodes_2 = rank2List
    matching_temp2 = []
    matching_temp = []
    maxReportedTetraMatched = 0
    maxReportedBlosumScore = -999999 # less than min. possible Blosum score
    continueToNextSeedPair = False
    for (i, j) in zip(nodes_1, nodes_2): # for each node j in g_2
      #print str(i) + ": " + str(len(g_1_dict[i])) + "-- " + str(len(g_1[i]))
      #print str(j) + ": " + str(len(g_2_dict[j])) + "-- " + str(len(g_2[j]))
      # Heuristic:  One of the 2 tetrahedra should have rank 1 vertex.
      #if i != 1 and j != 1:
      #    continue

      reset_data()
      allTetraMatches.append(lastTetraMatches)
      # for all the adjacent tetrahedra to i
      for k in range(len(g_1_dict[i]) - 1, -1, -1):
        # for all the adjacent tetrahedra to j
        for l in range(len(g_2_dict[j]) - 1, -1, -1):
          continueToNextSeedPair = False

          # get the k-th tetrahedron adjacent to i
          tetra_1 = allTetrahedra_1[g_1_dict[i][k]]
          # get the l-th tetrahedron adjacent to j
          tetra_2 = allTetrahedra_2[g_2_dict[j][l]]

          # EXPENSIVE OPERATION
          tetraIndex_1 = allTetrahedraIndex_1[frozenset(tetra_1)]
          tetraIndex_2 = allTetrahedraIndex_2[frozenset(tetra_2)]
          #ifdef heuristic-3
          # Matching volumes, sum of side lengths
          if (volList_1[tetraIndex_1] - volList_2[tetraIndex_2]) > 14.4:
            #print "Rejecting based on Volume: " + str(volList_1[tetraIndex_1]) + ", " + str(volList_2[tetraIndex_2])
            continue
          #endif
          #ifdef heuristic-4
          if abs(sideLengthsList_1[tetraIndex_1] - sideLengthsList_2[tetraIndex_2]) > 9.6:
            #print "Rejecting based on Sum of Side Lengths: " + str(sideLengthsList_1[tetraIndex_1]) + ", " + str(sideLengthsList_2[tetraIndex_2])
            continue
          #endif

          # The price we pay for prematurely converting to sets. :-)
          tetraAsLists_1 = next(t for t in allTetrahedraAsLists_1 \
              if set(t) == tetra_1)
          tetraAsLists_2 = next(t for t in allTetrahedraAsLists_2 \
              if set(t) == tetra_2)

          # sanity check
          assert(len(tetra_1) == 4 and len(tetra_2) == 4)

          # Compute the faces of the new tetrahedra
          tetraVertexList_1 = list(tetra_1)
          tetraVertexList_2 = list(tetra_2)

          faces_1 = list([set([tetraVertexList_1[0], tetraVertexList_1[1],
            tetraVertexList_1[2]]), set([tetraVertexList_1[0],
              tetraVertexList_1[1], tetraVertexList_1[3]]),
            set([tetraVertexList_1[2], tetraVertexList_1[1],
              tetraVertexList_1[3]]), set([tetraVertexList_1[0],
                tetraVertexList_1[2], tetraVertexList_1[3]])])

          vertex_seqs2 = list(itertools.permutations(tetraVertexList_2))
          #ifdef heuristic-5
          # Parity computation based on this line of code from below:
          #   vertex1Match[tetraVertexList_1[idx]] = vertex_seq2[idx]
          # Vertex ordering of orig tetras chosen for positive orient.
          # An even re-ordering of vertices keeps positive orient.
          parity1_orig = parityTetrahedron(tetraAsLists_1)
          parity1_now = parityTetrahedron(tetraVertexList_1[:4])
          parity2_orig = parityTetrahedron(tetraAsLists_2)
          #endif
          for vertex_seq2 in vertex_seqs2:
            #ifdef heuristic-5
            parity2_now = parityTetrahedron(vertex_seq2[:4])
            if (parity1_orig == parity1_now) != \
                (parity2_orig == parity2_now):
                  continue
            #endif

            #ifdef heuristic-6
            skipThis = False
            assert(len(tetraVertexList_1) == len(vertex_seq2) == 4)
            for (v1, v2) in zip(tetraVertexList_1, vertex_seq2):
              # RG: if both the residues are highly ranked they should match
                if v1 < 11 or v2 < 11:
                  if not residuesAreSimilar(v1, v2):
                    skipThis = True
                    # pass # Could use get_blosum_match_score instead
            if skipThis:
              continue
            #endif

            #   blosum62 scoring
            # Reset Blosum to 0, and begin new initial seed pair.
            totalBlosumScore = 0
            for (v1, v2) in zip(tetraVertexList_1, vertex_seq2):
              totalBlosumScore += \
                  get_blosum_match_score(blosumMatrix, v1, v2)

            bestNumTetrahedraMatched = 0
            bestVisitedVerticesStack1[:] = visitedVerticesStack1
            bestVertex1Match[:] = vertex1Match

            visitedTetraStack_1.append(tetra_1)
            visitedTetraStack_2.append(tetra_2)
            numTetrahedraVisited  += 1
            assert(visitedVerticesStack1 == [])
            # Actual vertex list ordering determined by code below:
            for idx in range(4):
              visitedVerticesStack1.append(tetraVertexList_1[idx])
              vertex1Match[tetraVertexList_1[idx]] = vertex_seq2[idx]
              visitedVertices2[vertex_seq2[idx]] = True
            # As vertex_seq2 varies, faces_1 will include all
            #  orderings of the four faces of tetra_1.
            # Next, try all possible faces
            for f1 in faces_1:
              if continueToNextSeedPair:
                continue

              # Choose f2 such that vertices of f1 and f2 match
              f2 = set([vertex1Match[v] for v in f1])

              #ifdef heuristic-1
              skipThis = False
              for v in tetra_1:
                # tmp1 = max(v, 10)
                # tmp2 = max(vertex1Match[v], 10)
                # if max(tmp1,tmp2) > 2 * min(tmp1,tmp2):
                #     skipThis = True
                if v <= 24 or vertex1Match[v] <= 24:
                  if abs(v-vertex1Match[v]) > 20:
                    skipThis = True
                if v < 10 or vertex1Match[v] < 10:
                  if abs(v-vertex1Match[v]) > 10:
                    skipThis = True
                if v < 3 or vertex1Match[v] < 3:
                  if abs(v-vertex1Match[v]) > 1:
                    skipThis = True
              skipThis = False # Use get_blosum_match_score instead
              if skipThis:
                continue
              #endif

              try:
                global waitMax
                # register a signal handler of SIGALRM
                signal.signal(signal.SIGALRM, handler)
                # We've now matched the seed tetrahedron pair:
                bestNumTetrahedraMatched = \
                    max(bestNumTetrahedraMatched, 1)
                # set signal to be fired off in waitMax seconds
                signal.alarm(waitMax)
                do_dfs_tetra(g_1, g_2, tetra_1, tetra_2, f1, f2, nodes_1, nodes_2)
                #print visitedTetraStack_1
                # if the function returned, cancel the signal
                signal.alarm(0)
                #bestScoreCriterion = \
                    #  (bestNumTetrahedraMatched > maxReportedTetraMatched)
                bestScoreCriterion = \
                    bestTotalBlosumScore > maxReportedBlosumScore
                if bestScoreCriterion:
                  # print the usual
                  maxReportedTetraMatched = bestNumTetrahedraMatched
                  maxReportedBlosumScore = \
                      bestTotalBlosumScore
                  lastTetraMatches = []
                  lastTetraMatches.append(zip(\
                  [annotateResidue(v, residueNames_1, residueNum_1)\
                                      for v in bestVisitedVerticesStack1],\
                  [annotateResidue(bestVertex1Match[v], residueNames_2,\
                                                          residueNum_2)\
                                      for v in  bestVisitedVerticesStack1]))
                  print "[" + str(i) + ", " + str(j) + "]: " + str(bestNumTetrahedraMatched) + ", " + str(bestTotalBlosumScore) + ", " + annotateResidues(bestVisitedVerticesStack1, residueNames_1, residueNum_1) + ", " + annotateResidues([bestVertex1Match[idx] for idx in bestVisitedVerticesStack1], residueNames_2, residueNum_2)
                  verify_matching()
              except Exception as exc:
                x,  = exc.args
                if x == ("end of time"):
                  # print the usual if it took too long
                  # On interrupt, must take max to find larger
                  bestScoreCriterion = \
                      bestTotalBlosumScore > maxReportedBlosumScore
                  if bestScoreCriterion:
                    # print the usual
                      maxReportedTetraMatched = max(bestNumTetrahedraMatched,
                          maxReportedTetraMatched)
                      maxReportedBlosumScore = max(bestTotalBlosumScore,
                          maxReportedBlosumScore)
                      lastTetraMatches = []
                      lastTetraMatches.append(zip(\
                      [annotateResidue(v, residueNames_1, residueNum_1)\
                                          for v in bestVisitedVerticesStack1],\
                      [annotateResidue(bestVertex1Match[v], residueNames_2,\
                                                              residueNum_2)\
                                          for v in  bestVisitedVerticesStack1]))
                      print "[" + str(i) + ", " + str(j) + "]: " + str(maxReportedTetraMatched) + ", " + str(maxReportedBlosumScore) + ", " + annotateResidues(bestVisitedVerticesStack1, residueNames_1, residueNum_1) + ", " + annotateResidues([bestVertex1Match[idx] for idx in bestVisitedVerticesStack1], residueNames_2, residueNum_2)
                  for idx in visitedVerticesStack1[4:]:
                    vertex1Match[idx] = 0 # unset the values
                  del visitedVerticesStack1[4:]
                  del visitedTetraStack_1[1:]
                  del visitedTetraStack_2[1:]
                  numTetrahedraVisited  = 1

                  # continueToNextSeedPair = True
                  continue
                  # return # cancel the run
                else:
                  raise

            for idx in range(4):
              vertex1Match[tetraVertexList_1[idx]] = 0 # unset the values
              visitedVerticesStack1.pop()
            visitedTetraStack_1.pop()
            visitedTetraStack_2.pop()
            numTetrahedraVisited  -= 1
            assert(len(visitedVerticesStack1) == 0)
            assert(len(visitedTetraStack_1) == 0)
    #make_graph(bestMatchingTetra_1, 'test_1.dot')
    #make_graph(bestMatchingTetra_2, 'test_2.dot')
    if maxReportedBlosumScore <= -999999:  # same value as when initialized
      # Print the answer now.  We never did it earlier.
      print "[" + str(0) + ", " + str(0) + "]: " + str(0) + ", " + str(bestTotalBlosumScore) + ", " + str(bestVisitedVerticesStack1) + ", " + str([bestVertex1Match[idx] for idx in bestVisitedVerticesStack1])
    allTetraMatches.append(lastTetraMatches)
    return calc_matches_union(allTetraMatches, consMatch)

def matchHighlyRankedResiduesNextLevel(numResidues, distCutoff,
                                       level, oldres1, rank1, rank2):
  global deepestLevel, count, allMatches
  for newres1 in range(oldres1+1,numResidues+1):
    rank1.append(newres1)
    rank2.append(-1)
    for newres2 in range(1,numResidues+1):
      if newres2 in rank2:
        continue
      rank2[-1] = newres2
      if not residuesAreSimilar(rank1[-1], rank2[-1]):
        continue
      distancesMatch = True
      for v1,v2 in zip(rank1[:-1],rank2[:-1]):
        if abs(distance_1(v1,newres1)-distance_2(v2,newres2)) > distCutoff:
          distancesMatch = False
          break
      if not distancesMatch:
        continue
      if deepestLevel < level:
        deepestLevel = level 
        count = 0
        allMatches = []
      if deepestLevel <= level:
        count += 1
        allMatches.append(zip(
          [annotateResidue(v, residueNames_1, residueNum_1) for v in rank1],
          [annotateResidue(v, residueNames_2, residueNum_2) for v in rank2]))
        print '['+str(count)+']: '+str(level)+', '+ \
              annotateResidues(rank1, residueNames_1, residueNum_1) \
               + ', ' + \
              annotateResidues(rank2, residueNames_2, residueNum_2)
      matchHighlyRankedResiduesNextLevel(numResidues, distCutoff,
                                         level+1, newres1, rank1, rank2)
    rank2.pop()
    rank1.pop()

def matchHighlyRankedResidues(numResidues, distCutoff):
  global deepestLevel, count, allMatches
  deepestLevel = count = 1
  allMatches = []
  print '[0]: 0, [], []'  # For automated tables, we need default result.
  rank1 = []
  rank2 = []
  i1 = 1
  while i1 < numResidues+1:
    rank1 = [i1]
    for j1 in range(1,numResidues+1):
      rank2 = [j1]
      # TODO: Allow dissimilar highly ranked residues??
      if not residuesAreSimilar(rank1[-1], rank2[-1]):
        continue

      i2 = i1+1
      while i2 < numResidues+1:
        rank1.append(i2)
        rank2.append(-1)
        for j2 in range(1,numResidues+1):
          if j2 in rank2:
            continue
          rank2[-1] = j2
          if not residuesAreSimilar(rank1[-1], rank2[-1]):
            continue
          if abs(distance_1(i1,i2)-distance_2(j1,j2)) > distCutoff:
            continue

          i3 = i2+1
          while i3 < numResidues+1:
            rank1.append(i3)
            rank2.append(-1)
            for j3 in range(1,numResidues+1):
              if j3 in rank2:
                continue
              rank2[-1] = j3
              if not residuesAreSimilar(rank1[-1], rank2[-1]):
                continue
              if abs(distance_1(i1,i3)-distance_2(j1,j3)) > distCutoff:
                continue
              if abs(distance_1(i2,i3)-distance_2(j2,j3)) > distCutoff:
                continue

              for i4 in range(i3+1,numResidues+1):
                rank1.append(i4)
                rank2.append(-1)
                for j4 in range(1,numResidues+1):
                  if j4 in rank2:
                    continue
                  rank2[-1] = j4
                  if not residuesAreSimilar(rank1[-1], rank2[-1]):
                    continue
                  distancesMatch = True
                  for v1,v2 in zip(rank1[:-1],rank2[:-1]):
                    if abs(distance_1(v1,i4)-distance_2(v2,j4)) > distCutoff:
                      distancesMatch = False
                      break
                  if not distancesMatch:
                    continue
                  # We now have 4 residues.  Test if parities agree.
                  if parityTetrahedron((j1,j2,j3,j4)) != \
                     parityTetrahedron((j1,j2,j3,j4)):
                    continue
                  level = 4  # corresponding to i4
                  if deepestLevel < level:
                    deepestLevel = level 
                    count = 0
                    allMatches = []
                  if deepestLevel <= level:
                    count += 1
                    allMatches.append(zip(
        [annotateResidue(v, residueNames_1, residueNum_1) for v in rank1],
        [annotateResidue(v, residueNames_2, residueNum_2) for v in rank2]))
                    print '['+str(count)+']: '+str(level)+', '+ \
                          annotateResidues(rank1, residueNames_1, residueNum_1)\
                           + ', ' + \
                          annotateResidues(rank2, residueNames_2, residueNum_2)

                  matchHighlyRankedResiduesNextLevel(numResidues, distCutoff,
                                                     5, i4, rank1, rank2)

                rank2.pop()
                rank1.pop()
            rank2.pop()
            rank1.pop()
            i3 += 1
        rank2.pop()
        rank1.pop()
        i2 += 1
    i1 += 1

if waitToRun:
  # Define run() command to be called interactively;  Useful if this is a module
  def run():
    matchHighlyRankedResidues(24, 2)
    #(p1_consensus_match, p2_consensus_match) = get_consensus_match(rank1List, rank2List)
    #add_matching_2(g_1, g_2, p1_consensus_match, p2_consensus_match)
else:
  matchHighlyRankedResidues(24, 2)
  if len(allMatches) > 1:
    print "# *** Intersection of all matches, for maximum length match"
  if not len(allMatches):
    print '[0]: 0, [], []'  # For automated tables, we need default result.
    sys.exit(0)
  consMatch = allMatches[0]
  for match in allMatches[1:]:
    # like set intersection, but order-preserving
    consMatch = [m for m in match if m in consMatch]
  consMatch1 =  "[" + ", ".join([x for (x,y) in consMatch]) + "]"
  consMatch2 =  "[" + ", ".join([y for (x,y) in consMatch]) + "]"
  if len(allMatches) > 1:
    print '[1]: '+str(len(consMatch))+', '+consMatch1+', '+consMatch2
  p1_consensus_match = [int(x.split('/')[0]) for (x,y) in consMatch]
  p2_consensus_match = [int(y.split('/')[0]) for (x,y) in consMatch]
  allTetraMatches = add_matching_2(g_1, g_2, p1_consensus_match, p2_consensus_match, consMatch)
  print "****************** Consensus match ********************"
  consMatch1 =  "[" + ", ".join([x for (x,y) in allTetraMatches]) + "]"
  consMatch2 =  "[" + ", ".join([y for (x,y) in allTetraMatches]) + "]"
  print "[1]: " + str(len(allTetraMatches)) + ", " + consMatch1+", "+consMatch2
