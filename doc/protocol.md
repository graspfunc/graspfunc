# OVERVIEW:

For each pair of proteins from the superfamily, there are many ways in which
the residues from the two proteins could be matched.

  1.  We search through all possible matches, rejecting some matches
       early due to a cutoff criterion.

  2.  Among the remaining matches, we assign a score for each match.

  3.  We take the match with the highest score for the given pair of
       proteins, and we then count the total number of residues that
       were matched for that highest-score match.  (If no match was
       found for the given pair of proteins, the number of residues
       matched is zero.)

  4.  We then set a threshold in terms of the above number of residues matched
       for each pair of proteins.  Any pair of proteins with a
       higher number of residues matched is considered to come from
       the same subgroup.

  5.  We apply the conclusions of step 4 transitively:  If the match (A,B)
       is above the threshold, and if the match (B,C) is above the threshold,
       then we also consider proteins A and C to be in the same subgroup.
       We do this, even if the match (A,C) produces a low number of residues
       matched.

------------------------------------------------------------------------------

## Details of protocol

The rest of this discussion is concerned only with matching a pair of proteins,
P~1~ and P~2~.  For this pair, we are concerned with three issues:

  1.  The conditions for a "cutoff" (rejection of any further matches in one
       branch of the search tree.

  2.  The criteria for the score of a protein match.

  3.  The manner of choosing the threshold for the number of residues matched.
       (Above that threshold, two proteins are considered in the same subgroup.)

*NOTE THAT WHILE THIS IS USED TO PRODUCE A GOOD SPREADSHEET FOR PROTEINS
WITHIN A SINGLE SUBGROUP, WE ARE STILL WORKING ON HOW TO EXTEND THAT TO A
MATCH ACROSS PROTEIN SUBGROUPS.  WE WILL ADD THAT TO THE PROTOCOL LATER.*

DELAUNAY:

  This is the pre-processing step. Triangulation of the protein structure
  is done in this step.

  In this preprocessing stage, each protein in isolation is converted
  from its PDB file to a partitioning of 3-dimensional space into a set
  of tetrahedra.  Each vertex, edge, or face of a tetrahedron can be
  (and usually is) shared among more than one tetrahedron.  Each vertex
  corresponds to the alpha carbon of a residue of the protein.

  The particular partitioning used is a "Delaunay triangulation".
  A Delaunay triangulation has the special property of forming a unique set
  of tetrahedra such that no additional vertices (alpha carbons) lie within
  the circumscribed sphere around a given tetrahedron.  This restriction
  also implies that vertices (residues) connected by an edge are usually
  close to each other.

PAIR_WISE_MATCHING:

  After the pre-processing step, pairwise matching of two proteins in a
  superfamily is done.

  There are two parts to the algorithm:

   - SEED_PAIRS: This discovers pairs of matching residues -

   - DELAUNAY_EXTENSION: This uses the matching residue pairs, and the set of
                          tetrahedra from triangulation to extend the match
                          size.

#### DELAUNAY (Preprocessing):

  The first step is to triangulate the entire protein structure. The
  alpha carbons of the residues on the protein's backbone are used for
  the triangulation.  This gives us a list of tetrahedra (and triangles)
  for the protein. Only the tetrahedra which are in the vicinity of
  the active site are retained for matching in the later stages of the
  algorithm.

  Steps:

   A. Triangulation of the entire protein is done using the "qhull" software
       (http://www.qhull.org/).

   B. Only the 50 most highly ranked residues are kept.

   C. Any vertices adjacent to those highly ranked residues are added back in.

   D. Only those tetrahedra are kept whose vertices are a subset of the above.

   E. Return the set of tetrahedra corresponding to the triangulation


#### PAIR_WISE_MATCHING: 

  The next step is to obtain a list of matching residue pairs from two proteins
  in the superfamily. This pair-wise matching is done in two steps.

  S~1~ = (r~1~, r~2~, ...) = Ordered sequence of POOL-ranked residues of P~1~

  S~2~ = (r~1~, r~2~, ...) = Ordered sequence of POOL-ranked residues of P~2~

  1. SEED_PAIRS: 

    The aim is discover a list of matching highly ranked residues from
    P~1~ and P~2~. This is a list of pairs (r~i~, r~j~) such that r~i~ is
    in S~1~, and r~j~ is in S~2~. The relative distance between any two
    r~i~'s from this list is within 2 Angstrom's of the relative distance
    between the two corresponding r~j~'s. At the end of the algorithm, the
    largest possible list of seed pairs is returned.

    Steps: 

     A. For each seed pair, (r~i~, r~j~):

      1. If the residues are not similar based on the chemical similarity
          matrix reject the seed pair.

      2. If the distance(r~i~, r'~i~) - distance(r~j~, r'~j~) is greater
          than the cutoff, reject the pair. Here, (r~i~, r~j~) is the seed
          pair that we are trying to match, and (r'~i~, r'~j~) is an already
          matched seed pair. 

      3. Store the seed pair (r~i~, r~j~).

     B. Return list of matching residue pairs

  2. DELAUNAY_EXTENSION:

    Residue pairs from step SEED_PAIRS are used as seeds in this part of the
    algorithm to extend the match size. The aim is to discover a list of pairs
    of matching tetrahedra (one from each protein).

    Steps:

     A. For each seed pair (r~i~, r~j~):

      1. Initialize the matching score to 0.

      2. For each pair of tetrahedra, (t~i~, t~j~), such that t~i~ has r~i~
          as one of its vertices, and t~j~ has r~j~ as one of its vertices.

        2.1. If the parities (orientation) of t~i~ and t~j~ do not match,
              reject the pair.

        2.2. If volumes of t~i~ and t~j~ differ by more than 14.4, reject
              the pair.

        2.3. If the sums of lengths of edges t~i~ and t~j~ differ by more
              than 9.6, reject the pair.

        2.4. For each matched vertex pair, (v~i~, v~j~), of t~i~ and t~j~:

          + If v~1~ or v~2~ is among the top-11 POOL-ranked resides in P~1~
             and P~2~, respectively, and v~1~ is not chemically similar to
             v~2~, reject the pair of tetrahedra and move to the next one.

          + If v~1~ or v~2~ is among the top-24 POOL ranked residues in their
             respective proteins, and |POOL-rank(v~1~) - POOL-rank(v~2~)|
             > 24, reject this pair of tetrahedra.

          + If v~1~ or v~2~ is among the top-10 POOL ranked residues in their
             respective proteins, and |POOL-rank(v~1~) - POOL-rank(v~2~)|
             > 10, reject this pair of tetrahedra.

          + If v~1~ or v~2~ is among the top-3 POOL ranked residues in their
             respective proteins, and |POOL-rank(v~1~) - POOL-rank(v~2~)|
             > 1, reject this pair of tetrahedra.

          + Calculate a matching score based on BLOSUM-62 matrix, and
             store it.

      3. At this point we have pair of matching tetrahedra. Using this
          pair, the algorithm extends the list by doing a depth-first
          search for matching tetrahedra. The algorithm sequentially
          searches in four directions given by the four pairs of
          matching faces of the tetrahedra pair.

        + For each new discovered pair of tetrahedra, apply the same
           rejection criteria (described above), calculate the total
           matching score, and add it to the existing score along this
           search path.

        + Use the new matched pair of tetrahedra and repeat the step 3.

      4. Record this path (with matching residue pairs) and the matching score.

     B. Return the list of matching residue pairs with the best matching score.


------------------------------------------------------------------------------

<!--

###  Heuristics used:

  1. Triangulation:
     A. Triangulation of the entire protein is done
     B. Only the num-protein-residues most highly ranked residues are kept.
          50 is recommended as a good choice.
     C. Any vertices adjacent to those highly ranked residues are added back in.
     D. Only those tetrahedra are kept whose vertices are a subset of the above.
  2. The tetrahedra are used as seed pairs if each has a vertex contained
       among the top-n-seed-pairs most highly ranked residues.
       3 or 4 is recommended as a good choice.
  3. There is an optional heuristic to insist that one tetrahedra of the pair
       contains the most highly ranked vertex.  This is commented out here.
  4. A seed pair of triangulations is rejected if the parity does not match.
  5. We should replace the residue comparison test.  Reject if:
       |r1-r2| > 20 for r1 or r2 <= 24; and |r1-r2| > 10 for r1 or r2 < 10;
       and |r1-r2| > 2 for r1 or r2 < 3
       Worse test (commented out) rejects if:  max(tmp1,tmp2) > 2*min(tmp1,tmp2)
           where tmp1 and tmp2 are the larger of vertex rank and 10.
       This test now added to get_blosum_score (bonus/penalty), but no reject
  6. If vertex 1,2,3,4,5 is present (5 most highly ranked residues), then the
       chemical similarity types must match, based on aminoAcidSimilarityList3.
       (See NOTE A, below.)
  7. If the volume of a tetrahedron doesn't match, reject.
       Using hardwired value of 14.4.  Should be parameter.
  8. If the sum of lengths of edges of a tetrahedron doesn't match, reject.
       Using hardwired value of 9.6.  Should be parameter.
  9. Search for best BLOSUM62 score (with normalization) + 5 for each residue
       of rank 1 matched, 3 for rank 2, 1 for rank 3.  (NOTE: matching a rank 1
       residue from each protein is worth 10.)  This should be merged into
       get_blosum_match_score() routine.
  10. If depth > 10, and Blosum score 10 less than best, give up (heuristic-8)
  11. Setting top-n-seed-pairs = 4 on cmd line and using num-protein-residues=25
       seems better for this version. (See NOTE B, below.)
  12. If the distance between the first vertex and the latest vertex for
       protein1 differs by more than 5 from that for the matching vertices
       of protein2, then reject those matching vertices.  Also apply it
       between the middle vertex in visitedVerticesStack1 and latest vertex.
       [ Two residues on backbone have typical distance 3.7.  The choice of 5,
         is somewhat arbitrary.  Ideally, we would also consider distances
         between other pairs of residues.  But this would be inefficient. ]
### FUTURE:
  A. We should reject a tetrahedron if the matching faces don't agree
       on whether it is an exterior face of the protein.
  B. Lots of opportunities to speed up the code!
  C. For other protein families, should scale num-protein-residues with
       size of protein.
  D. When two highly ranked residues of a protein are matched, see if original
       geometric distances between residues in protein1 and protein2 agree.
       If they disagree, reject matching both residues.
  NOTE:
  A. POOL matchings have found chemical property almost always conserved in
       multiple alignment.  This is a key to clustering proteins successfully.
  B. Reducing num-protein-residues, or allowing larger wait-time, is useful
       to avoid timeouts.  A timeout will stop searching a given branch,
       and often miss a good match, because it was spending too much time
       searching a smaller sub-branch of that same branch.

-->
