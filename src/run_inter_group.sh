#!/bin/bash

DIRS=($1/g*)

for ((p=0; p < ${#DIRS[@]}; p++))
do
    for ((q=p+1; q < ${#DIRS[@]}; q++))
    do
        d=${DIRS[$p]}
        c=${DIRS[$q]}
        echo "***********************************************************************************"
        echo $d $c
        echo "***********************************************************************************"
        for i in 4
        do
            FILES1=($d/*.sc)
            FILES2=($c/*.sc)
            echo "#############################################################################"
            echo $i
            echo "#############################################################################"
            for ((j=0; j < ${#FILES1[@]}; j++))
            do
                for ((k=0; k < ${#FILES2[@]}; k++))
                do
                    file1=`dirname ${FILES1[$j]}`/`basename ${FILES1[$j]} .sc`
                    file2=`dirname ${FILES2[$k]}`/`basename ${FILES2[$k]} .sc`
                    echo ============ $j $k
                    echo python2.7 mains_18.py $file1 $file2 $i 25 10 $d/../matrix/BLOSUM62.txt 12
                    echo "#############################################################################"
                    python2.7 mains_18.py $file1 $file2 $i 25 10 $d/../matrix/BLOSUM62.txt 12 | tee inter.temp
                    head -4 inter.temp >> $2
                    tail -1 inter.temp >> $2
                done
            done
        done
    done
done

rm inter.temp
