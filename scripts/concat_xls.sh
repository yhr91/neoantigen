#!/bin/bash

ALLELE=$1 

cat ./$ALLELE"_0.xls" > out0
mv ./$ALLELE"_0.xls" ./$ALLELE"_0.xls_temp"

# Run netMHCpan on each split file
find . -name "*_?.xls" -print0| sort -zV | xargs --null -n 1 tail -n +3 > out1
find . -name "*_??.xls" -print0| sort -zV | xargs --null -n 1 tail -n +3 > out2

# mv ./$ALLELE"_0.xls_temp" ./$ALLELE"_0.xls"

cat out0 out1 out2 > $ALLELE"_all.xls"
