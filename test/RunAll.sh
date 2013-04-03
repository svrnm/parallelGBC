#!/bin/bash
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


# Counter for the number of executed tests
ACOUNT=0;
declare -i ACOUNT;

# For 1 to 4 processors do ...
c=$1;
verbosity=$2
simplify=$3
sugar=$4
echo -e "\nRunning tests with \033[1;34m${c} core(s)\033[0m:"
# For all files in gb/ do ...
for i in input/*;
do
	# Count this run
	ACOUNT=$ACOUNT+1;
	# Output the input file name
	echo -en "${i##"input/"} ... ";
	# Run the test
	./test/test-f4.bin $i $c $verbosity 0 1024 $simplify $sugar
done;

# If not all tests passed print a statistic how many tests failed.
if [ $FCOUNT -gt 0 ]
then
	echo -e "\n\033[1;31m${FCOUNT} of ${ACOUNT} tests failed, so something went wrong.\033[0m"
else
	echo -e "\n\033[1;32mAll tests passed!\033[0m";
fi
