# GRASP-Func Protocol

## Accessing the repository

 1. On repository page, click Clone or download > Use SSH > copy

        $ git clone git@github.com:graspfunc/graspfunc.git

    Now all the files that are currently on GitHub will be added to your home directory

## Preprocessing

The following set of steps are used to prepare the input files that
the main program can process. Also, see the [document](doc/pre-processing-steps.txt).

 1. Make input folder (name input-"family")

 2. In input-"family", add raw files (.pdb and .ranks -- TC or TIC)

 3. Make text file with list of PDBs (pfiles.txt) and rank (rfiles.txt) in input-"family"

 4. Make sure list matches exactly with original files. One possible ways is to do the following:

        $ ls -l *.pdb > pdb-files.txt
        $ ls -l *rank* > rank-files.txt

    Edit the two files in any text editor so they only have the file names. Also, ensure that the
    filenames in the two text files are in the same order.
 
 5. Generate input files for the Qhull program.

        $ src/pre-processing/generate_qhull_files.py pfiles.txt rfiles.txt

    This one works on many proteins; calls generate_qhull_file.py for each protein in the given list

 6. Extract the ranks of residues from the POOL-generated rank files.

        $ src/pre-processing/generate_rank_files.py rfiles.txt
 
 7. Run the Qhull program to generate Delauney Tesselation of the protein structure.

        $ src/pre-processing/generate_sc_files.py pfiles.txt

 8.  Copy all files from the above steps and add them into the newly created folder.

        $ cp *.qhull *.rank *.sc input-"family"/

## Processing

 1. Run `command-line.sh`

        $ src/command-line.sh input-"family"/

## Post Processing (Analysis)

 1. To see a summary of results, execute the following command.

        $ src/post-processing/table-summary.py results_intra.txt results_inter.txt

 2. To generate an alignment csv file, execute the following command.

        $ src/post-processing/align-seeds.py results_intra.txt results_inter.txt

 3. To see a graph of matches, execute the following command.

        $ src/post-processing/graph_results.py results_intra.txt results_inter.txt
    
    Note when running this command on a remote server, this requires X forwarding.
    (`ssh -X`).

## Pre-requisites

 1. Python 2.7

 2. [Qhull](http://www.qhull.org/)

 3. [NetworkX](https://networkx.github.io/)

 4. [Graph-Tool](https://graph-tool.skewed.de/)
