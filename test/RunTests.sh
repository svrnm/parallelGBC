#!/bin/bash
set -o pipefail

##
# This file provides a test suite for the F4 implementation used in
# parallelGBC. For all files in gb/ the script looks up the matching
# file in input/ and computes the groebner basis and compares the
# result with the pre-computed expected result on MIN_C (first CORE_LIST value) only.
# Higher thread counts stress parallel F4; with VERIFY_GB=1, Buchberger is checked on MAX_C.
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
	echo -e "${FAILED}"
	FCOUNT=$FCOUNT+1;
}

# Success function: Print success message
function passed() {
	echo -e "${PASSED}"
}

# Optional Buchberger criterion verification (set VERIFY_GB=1 to enable).
# When enabled, verifies computed GB on the highest processor count in CORE_LIST (parallel F4 + parallel S-pair check).
# Verification can be expensive on larger benchmarks; use VERIFY_TIMEOUT (seconds) to cap runtime per case.
# VERIFY_MAX_GB: skip Buchberger check when the expected |G| (comma-separated in gb/*.txt) exceeds this.
# 0 = no limit. CI sets a cap (e.g. 100); use 0 to verify all sizes (can be very slow).
VERIFY=${VERIFY_GB:-0};
VERIFY_PROGRESS=${VERIFY_PROGRESS:-1};
VERIFY_TIMEOUT=${VERIFY_TIMEOUT:-120};
VERIFY_MAX_GB=${VERIFY_MAX_GB:-0};
TIMEOUT_BIN=$(command -v gtimeout || command -v timeout || true)

# Processor counts for F4 (Buchberger verify runs on MAX_C only). Override e.g. CORE_LIST="1 2 4".
CORE_LIST=${CORE_LIST:-"1 8"}
MIN_C=$(echo "$CORE_LIST" | awk '{print $1}')
MAX_C=$(echo "$CORE_LIST" | awk '{print $NF}')

# test-f4 executable (default: CMake build tree). Override if needed.
TEST_F4_BIN=${TEST_F4_BIN:-build/test-f4}

for c in $CORE_LIST;
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
		# Run the test (sort polynomials for comparison - GB element order can vary with parallel execution)
		# When VERIFY_GB=1 and c is MAX_C, verify Buchberger criterion (unless |G| too large)
		VERIFY_THIS=0
		if [ "$VERIFY" != "0" ] && [ "$c" = "$MAX_C" ]; then
			VERIFY_THIS=1
			if [ "$VERIFY_MAX_GB" != "0" ]; then
				EXPECTED_GB_N=$(head -1 "$f" | awk -F', ' '{print NF}')
				if [ "$EXPECTED_GB_N" -gt "$VERIFY_MAX_GB" ]; then
					VERIFY_THIS=0
					echo -e "\033[1;33mverify skipped (expected |G|=${EXPECTED_GB_N} > VERIFY_MAX_GB=${VERIFY_MAX_GB})\033[0m"
				fi
			fi
		fi
		if [ "$VERIFY_THIS" = "1" ]; then
			if [ -n "$TIMEOUT_BIN" ]; then
				"$TIMEOUT_BIN" "${VERIFY_TIMEOUT}s" "$TEST_F4_BIN" $i $c 0 1 1024 0 1 1 "$VERIFY_PROGRESS" | awk -F', ' '{for(i=1;i<=NF;i++) print $i}' | sort > /tmp/pgbc_actual.$$
			else
				perl -e 'my $t=shift; my $pid=fork(); exit 125 unless defined $pid; if($pid==0){ exec @ARGV or exit 127; } my $deadline=time+$t; while(1){ my $r=waitpid($pid, 1); if($r==$pid){ exit($? >> 8); } if(time >= $deadline){ kill 9, $pid; waitpid($pid, 0); exit(124); } select undef,undef,undef,0.1; }' "$VERIFY_TIMEOUT" "$TEST_F4_BIN" $i $c 0 1 1024 0 1 1 "$VERIFY_PROGRESS" | awk -F', ' '{for(i=1;i<=NF;i++) print $i}' | sort > /tmp/pgbc_actual.$$
			fi
			CMD_STATUS=${PIPESTATUS[0]}
			if [ $CMD_STATUS -eq 124 ] || [ $CMD_STATUS -eq 137 ] || [ $CMD_STATUS -eq 142 ]; then
				echo -e "\033[1;33mverify timeout (${VERIFY_TIMEOUT}s)\033[0m";
				failed
				rm -f /tmp/pgbc_actual.$$ /tmp/pgbc_expected.$$
				continue
			fi
			if [ $CMD_STATUS -ne 0 ]; then
				failed
				rm -f /tmp/pgbc_actual.$$ /tmp/pgbc_expected.$$
				continue
			fi
		else
			"$TEST_F4_BIN" $i $c 0 1 | awk -F', ' '{for(i=1;i<=NF;i++) print $i}' | sort > /tmp/pgbc_actual.$$
			CMD_STATUS=${PIPESTATUS[0]}
			if [ $CMD_STATUS -ne 0 ]; then
				failed
				rm -f /tmp/pgbc_actual.$$
				continue
			fi
		fi
		# Reference gb/*.txt matches ApCoCoA-style output; parallel F4 may differ as a set while still being a GB.
		if [ "$c" = "$MIN_C" ]; then
			awk -F', ' '{for(i=1;i<=NF;i++) print $i}' $f | sort > /tmp/pgbc_expected.$$
			diff -q /tmp/pgbc_actual.$$ /tmp/pgbc_expected.$$ >> /dev/null && passed || failed
			rm -f /tmp/pgbc_actual.$$ /tmp/pgbc_expected.$$
		else
			passed
			rm -f /tmp/pgbc_actual.$$
		fi
	done;
done;

# If not all tests passed print a statistic how many tests failed.
if [ $FCOUNT -gt 0 ]
then
echo -e "\n\033[1;31m${FCOUNT} of ${ACOUNT} tests failed, so something went wrong.\033[0m"
else
echo -e "\n\033[1;32mAll tests passed!\033[0m";
fi
