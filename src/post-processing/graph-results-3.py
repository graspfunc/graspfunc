#!/usr/bin/python3

import graph_tool.all as gt
from graph_tool.all import *
import sys

if len(sys.argv) <= 2:
  sys.exit("Usage:  {0}".format(sys.argv[0]) + \
           " <results_inter_group> <results_intra_group>" + \
	   " [matchType=<tetrahedra|blosum62; default: tetrahedra>]" + \
	   " [saveFile=<graph.pdf>] [powerLaw=<number; default: 2]" + \
	   " [highlightGroups=<comma-separated list of groups>]" + \
           " [verbose=<bool>; default: false]")

# Set this to filename to save display graph.
matchType = 'tetrahedra'
saveGraphicsFile = None
powerLaw = 2.0
verbose = False
epsilon = 1e-2
highlightGroups = []
for arg in sys.argv:
  if arg.startswith('matchType='):
    matchType = arg[len('matchType='):].lower()
  if arg.startswith('saveFile='):
    saveGraphicsFile = arg[len('saveFile='):]
  if arg.startswith('highlightGroups='):
    highlightGroups = arg[len('highlightGroups='):].split(',')
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

matches = list()

def readMatches(filen):
  global verbose
  global matches
  parity = 0
  for line in open(filen):
    if verbose:
      print(line)
    if line.startswith('#'):  # ignore as comment character
      continue
    if "../input" not in line and not line.startswith("["):
      continue
    if parity == 0: # starting new match
      match = list(range(5))
    if parity in [0,2]:
      line = line[line.index('../input'):]
      assert(line.startswith('../input'))
      prot = int(parity/2)
      group = int(parity/2)+3
      match[prot] = (line.split('/')[-1]).split(':')[0]
      match[group] = line.split('/')[-2][1:]
      #assert(match[parity/2+3] in "12345678")
      if (match[group] not in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"]):
        print(line)
        print(filen)
        assert(False)
      # Optionally change name to:  subgroup:protein
      match[prot] = match[group] + ':' + match[prot]
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
print("len(matches): ", len(matches))

proteins = []
for match in matches:
  proteins.append(match[0])
  proteins.append(match[1])

print(len(list(set(proteins))))
# This would be nice, but it's not order-preserving
# proteins = list(set(proteins))
# So, I'll use a O(n^2) algorithm for now:
for i in range(len(proteins)):
  if proteins[i] in proteins[:i]:
    proteins[i] = ""
proteins = [p for p in proteins if p]

print("len(proteins): ", len(proteins))

# ==
# This is the new (or correct) mapping
# Supergroup, Group, PDB IDs, Biochemical Function
# 1 1 1pii(N), 1i4n, 2c3z IGPS
# 1 2 1pii(C), 1lbm PRAI
# 2 3 1geq, 1qop, 1rd5, 1xc4 TrpA
# 2 4 1rpx, 2fli, 1h1y, 1tqj, 3ovp RPE
# 3 5 1qo2, 1vzw, 2y85 HisA
# 3 6 1thf, 1h5y, 1ox6 HisF
# 3 7 1dbt, 1dv7, 1dqw, 1l2u, 2za1, 3qw3, 3l0k OMPDC
# 3 8 1xbv, 3exr KGPDC
# 3 9 3ajx, P42405 HPS
def mapOldProteinToNew(name):

  pdbId = name.split(":")[1].lower()
  origGroup = name.split(":")[0]
  newGroup = origGroup

  if pdbId in ["1i4n", "2c3z"]:
    newGroup = "1"
  elif pdbId in ["1lbm"]:
    newGroup = "3"
  elif pdbId in ["1pii"]:
    if origGroup == "5":
      newGroup = "1"
    else:
      newGroup = "3"
  elif pdbId in ["1geq", "1qop", "1rd5", "1xc4"]:
    newGroup = "2"
  elif pdbId in ["1rpx", "2fli", "1h1y", "1tqj", "3ovp"]:
    newGroup = "6"
  elif pdbId in ["1qo2", "1vzw", "2y85"]:
    newGroup = "4"
  elif pdbId in ["1thf", "1h5y", "1ox6"]:
    newGroup = "5"
  elif pdbId in ["1dbt", "1dv7", "1dqw", "1l2u", "2za1", "3qw3", "3l0k"]:
    newGroup = "7"
  elif pdbId in ["1xbv", "3exr"]:
    newGroup = "8"
  elif pdbId in ["3ajx", "p42405"]:
    newGroup = "9"
    if pdbId in "p42405":
      pdbId = "HPS1"
  elif origGroup in [str(x) for x in range(9, 14)]:
    newGroup = str(int(origGroup) + 1)
  print(name + " --> " + newGroup + ":" + pdbId.upper())
  return newGroup + ":" + pdbId.upper()

# ====
# Now build the graph

g = Graph()
g.add_vertex(len(proteins))
  
shape = g.new_vertex_property("string")
size = g.new_vertex_property("int")
name = g.new_vertex_property("string")
color = g.new_vertex_property("string")

g.vertex_properties["shape"] = shape
g.vertex_properties["size"] = size
g.vertex_properties["name"] = name
g.vertex_properties["color"] = color

proteinVertex = {}
proteinsIter = (protein for protein in proteins)
for vertex in g.vertices():
  shape[vertex] = "circle"
  p = proteinsIter.__next__()
  # XXX: Enable this for RPBB
  name[vertex] = mapOldProteinToNew(p.upper())
  # XXX: Enable this for all other superfamilies
  #name[vertex] = p.upper()
  if highlightGroups and any(s in name[vertex] for s in highlightGroups):
    color[vertex] = 'green'
  else:
    color[vertex] = 'cyan'
  # name[vertex] = p[0] ; assert(name[vertex] in "12345678")
  proteinVertex[p] = vertex

# Instead of using raw capacity (matchSize) for edges, we will use a
#   power of the capacity:  pow(capacity, n).
# The idea is that at a very large n, the min-cut problem reduces to adding
#   edges in order of capacity, and stopping just before forming a connected component.
# If we can guarantee that intra-group edges are always larger than inter-group
#   edges, then this always works.  Since there will always be some inter-group
#   edges that are too large, we reduce the power n.  Then a single edge between
#   groups does not destroy our clustering.
# If the power n is too small, the chosen cut is around the single source or target.
# We input the expected ratio between inter- and intra-edge matches.
# We also input the expectted size of a group.  The result is the desired power.
def capacityPower (graphSize, interIntraRatio, groupSize, margin=2):
  # Normalize to assume all intra-group edges have size 1, and inter: interIntraRatio
  # (groupSize - 1) = (graphSize - groupSize) * interIntraRatio^n
  # Multiply right hand size by margin, so intra-group capacity larger than inter-group
  return (log(groupSize - 1) - log(margin) - log(graphSize - groupSize)) \
	  / log(interIntraRatio)

capacity = g.new_edge_property("float")
g.edge_properties["capacity"] = capacity
pen_width = g.new_edge_property("float")
g.edge_properties["pen_width"] = pen_width

# The graph will have two directed edges between a pair of vertices
#   for the sake of Boykov-Komogorov below.
g.set_directed(True)
maxCapacity = max([m[2] for m in matches])
print( "maximum match:", maxCapacity)
intraMatches = [m[2] for m in matches if m[3] == m[4]]
medianMatches = [ min(sorted(intraMatches[int(len(intraMatches)/2):])),
                  max(sorted(intraMatches[:int((len(intraMatches)+1)/2)])) ]
if medianMatches[0] == 0 and medianMatches[1] == 0:
  medianMatches[0] = 5
  medianMatches[1] = 5  # For GH superfamily
medianIntraCapacity = sum(medianMatches)/2
print("medianIntraGroupMatch", medianIntraCapacity)
for match in matches:
  assert(match[0] != match[1])
  edge1 = g.add_edge(proteinVertex[match[0]], proteinVertex[match[1]])
  power = powerLaw
  # Absolute value of capacity determines attractive force in graphing for sfdp_layout.
  capacity[edge1] = 3.0*pow(match[2],power)/pow(maxCapacity,power)
  pen_width[edge1] = 3.0 * min(capacity[edge1], 3*pow(medianIntraCapacity,power))
  # max-flow, min-cut routines below don't like edge capaicty of zero or near zero
  if capacity[edge1] < epsilon:
    capacity[edge1] = 0
  # capacity[edge1] = 2.0*float(match[2])/maxCapacity
  edge2 = g.add_edge(proteinVertex[match[1]], proteinVertex[match[0]])
  capacity[edge2] = capacity[edge1]
  pen_width[edge2] = pen_width[edge1]

for e in g.edges():
  if e.source() == e.target():
    print("self-edge:", e.source)

# ====
# Now partition the graph using max flow.

if True:
  # Really should choose src and target with min flow between them (far apart).
  src = target = None
  for i in range(len(proteins)):
    # if '1L2U' in proteins[i]:
    if '1:' in proteins[i]:
      src = proteinVertex[proteins[i]]
    # if '1RD5' in proteins[i]:
    if '4:' in proteins[i]:
      target = proteinVertex[proteins[i]]
  assert(src and target)
else:

  # Find (src,target) pair having the smallest "max flow" between them.
  # There are much more efficient algorithms, but it's not a bottleneck right now.
  min_pair = [sum(capacity.a),None, None]
  for src in g.vertices():
    for target in g.vertices():
      if src != target:
        res = gt.boykov_kolmogorov_max_flow(g, src, target, capacity)
        minCut, part = gt.min_st_cut(g, src, res)
        # res.a = res.a - capacity.a
        # max_flow = sum(res[e] for e in target.in_edges())
        if minCut > 0 and minCut < min_pair[0] and sum(capacity[e] for e in src.out_edges()) > minCut + epsilon and sum(capacity[e] for e in target.in_edges()) > minCut + epsilon:
          print("Found better pair:", minCut, name[src], name[target])
          print("minCut: ", minCut, "; sum1/src: ", sum(capacity[e] for e in src.out_edges()), "sum2/tgt: ", sum(capacity[e] for e in target.in_edges()), "; part: ")
          print([name[v] for v in g.vertices() if part[v] == 0], [name[v] for v in g.vertices() if part[v] == 1])
          min_pair = [minCut, src, target]
  src, target = min_pair[1], min_pair[2]

#shape[src] = shape[target] = "square"
# Note, must be at least large enough for font_size
for v in g.vertices():
  size[v] = 10
size[src] = size[target] = 80

# Can also get cut of undirected graph, but usually cuts into singleton group
# view = gt.GraphView(g, directed=False)
# mc,part = min_cut(view, capacity)
# print "Min cut of undirected graph:", mc
# print "Min cut parttion:", [name[v] for v in view.vertices() if part[v] == 0], \
#       "      ", [name[v] for v in view.vertices() if part[v] == 1]

g.set_directed(True)  # required for max_flow routine below.
# Compute residual capacity after max flow.  This partitions the graph.
# Consider also using:  push_relabel_max_flow  O(V^3)
# Boykov-Kolmogorov is: O(EV^2 C), but often performs much better.
res = gt.boykov_kolmogorov_max_flow(g, src, target, capacity)
#mc, part = gt.min_st_cut(g, src, capacity, res)
part = gt.min_st_cut(g, src, capacity, res)
mc = sum([capacity[e] - res[e] for e in g.edges() if part[e.source()] != part[e.target()]])
print(mc)

part1 = [name[v] for v in g.vertices() if part[v] == 0]
part2 = [name[v] for v in g.vertices() if part[v] == 1]

print("src:", name[src], "; target:", name[target])
print("partition1:", part1)
print("partition2:", part2)

# ====
# Now display the graph

pos = gt.sfdp_layout(g, eweight=capacity)

# Commented out.  (Displaying double edges makes the graph too cluttered.)
# gt.graph_draw(g, pos=pos, edge_pen_width=capacity, vertex_fill_color=part,
#              vertex_text=name, vertex_font_size=9,
#              output=saveGraphicsFile)

forwardEdge = g.new_edge_property("bool")
g.edge_properties["forwardEdge"] = forwardEdge
for i in range(g.num_vertices()):
  for j in range(i):
    if g.edge(g.vertex(i), g.vertex(j)):
      forwardEdge[g.edge(g.vertex(i), g.vertex(j))] = True

view = gt.GraphView(g, efilt=forwardEdge, directed=False)
view.set_directed(False)

pos = gt.sfdp_layout(view, eweight=capacity)

# gt.prop_to_size(capacity,power=3.0,log=True)

gt.graph_draw(view, pos=pos, edge_pen_width=pen_width, vertex_fill_color=color,
    	      vertex_shape=shape,
    	      vertex_size=size,
              # vertex_text=g.vertex_index, vertex_font_size=18,
              output_size=(1916, 1037),
              vertex_text=name, vertex_font_size=14,
              output=saveGraphicsFile,
              fit_view=True)

# gt.graphviz_draw(view,
#                  pos=pos,
#                  vsize=0.3,
#                  pin=True,
#                  vcolor="cyan",
#                  penwidth=pen_width,
#                  overlap=False,
#                  vprops={"label":name},
#                  output=saveGraphicsFile)
