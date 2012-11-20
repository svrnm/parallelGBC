#!/bin/bash

##
# This file provides a test suite for the F4 implementation used in
# parallelGBC. For all files in gb/ the script looks up the matching
# file in input/ and computes the groebner basis and compares the
# result with the pre-computed expected result.
#
######
#
# parallelGBC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# parallelGBC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with parallelGBC.  If not, see <http://www.gnu.org/licenses/>.


# Message text if computation fails
FAILED="\033[1;31mfailed\033[0m"
# Message text if computation passes
PASSED="\033[0;32mpassed\033[0m"

# Counter for the number of executed tests
ACOUNT=0;
declare -i ACOUNT;

# Counter for the number of failed tests.
FCOUNT=0;
declare -i FCOUNT;

# Fail function: Print error message and count the error
function failed() {
	echo -e ${FAILED}
	FCOUNT=$FCOUNT+1;
}

# Success function: Print success message
function passed() {
	echo -e ${PASSED}
}

# For 1 to 4 processors do ...
for c in 1 2 4;
	do
	echo -e "\nRunning tests with \033[1;34m${c} core(s)\033[0m:"
	# For all files in gb/ do ...
	for f in gb/*;
	do
		# Count this run
		ACOUNT=$ACOUNT+1;
		# Set the path to the input file
		i=input/${f##"gb/"};
		# Output the input file name
		echo -en "${f##"gb/"} ... ";
		# Run the test
		./test/test-f4.bin $i $c 0 1 | diff -q - $f >> /dev/null && passed || failed
	done;
done;

# If not all tests passed print a statistic how many tests failed.
if [ $FCOUNT -gt 0 ]
then
echo -e "\n\033[1;31m${FCOUNT} of ${ACOUNT} tests failed, so something went wrong.\033[0m"
else
echo -e "\n\033[1;32mAll tests passed!\033[0m";
fi
