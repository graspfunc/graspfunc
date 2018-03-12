#!/bin/bash

# this assumes that the input files are arranged in subdirectories: g1, g2, ...
DIRS=$1/g* # enable this for all groups
# DIRS=$1  # enable this for just one group

for d in $DIRS
do
    echo "***********************************************************************************"
    echo $d
    echo "***********************************************************************************"
    for i in 4
    do
        FILES1=($d/*.sc)
        echo "#############################################################################"
        echo $i
        echo "#############################################################################"
        for ((j=0; j < ${#FILES1[@]}; j++))
        do
            for ((k=j+1; k < ${#FILES1[@]}; k++))
            do
                file1=`dirname ${FILES1[$j]}`/`basename ${FILES1[$j]} .sc`
                file2=`dirname ${FILES1[$k]}`/`basename ${FILES1[$k]} .sc`
                echo ============== $j $k
                echo python2.7 mains_18.py $file1 $file2 $i 25 10 $d/../matrix/BLOSUM62.txt 12
                echo "#############################################################################"
                python2.7 mains_18.py $file1 $file2 $i 25 10 $d/../matrix/BLOSUM62.txt 12 | tee intra.temp
                head -4 intra.temp >> $2
                tail -1 intra.temp >> $2
            done
        done
    done
done

rm intra.temp
