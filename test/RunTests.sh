#!/bin/bash
# First test set, this is just a start. There will be more tests and hopefully a better testing system...
TESTS="cyclic4 cyclic5 katsura7"
P=${0%/*}; # Binary should be in the same path at this file
INPUT="input/"
COMPARE="gb/"
for TEST in ${TESTS}; 
do
for PROC in `seq 1 4`; do
R="`${P}/test-f4.bin ${INPUT}/${TEST}.txt ${PROC} 0 1 | diff -q - ${COMPARE}/${TEST}.txt`"
if [ $? -eq 1 ]; then
	echo -e "$TEST ($PROC processors)\t\t\t\t\e[01;31mfailed\e[00;30m:\t${R}"
else
	echo -e "$TEST ($PROC processors)\t\t\t\t\e[01;32mpassed\e[00;30m"
fi;
done
done
