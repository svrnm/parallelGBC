#!/bin/bash
FAILED="\033[1;31mfailed\033[0m"
PASSED="\033[0;32mpassed\033[0m"

ACOUNT=0;
FCOUNT=0;
declare -i FCOUNT;
declare -i ACOUNT;

function failed() {
	echo -e ${FAILED}
	FCOUNT=$FCOUNT+1;
}

function passed() {
	echo -e ${PASSED}
}
for c in 1 2 4;
	do
	echo -e "\nRunning tests with \033[1;34m${c} core(s)\033[0m:"
	for f in gb/*;
	do
		ACOUNT=$ACOUNT+1;
		i=input/${f##"gb/"};
		echo -en "${f##"gb/"} ... ";
		./test/test-f4.bin $i $c 0 1 | diff -q - $f >> /dev/null && passed || failed
	done;
done;

if [ $FCOUNT -gt 0 ]
then
echo -e "\n\033[1;31m${FCOUNT} of ${ACOUNT} tests failed, so something went wrong.\033[0m"
else
echo -e "\n\033[1;32mAll tests passed!\033[0m";
fi
