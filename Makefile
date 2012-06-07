# Makefile for parallelGBC 0.2 and later
# 
# Just execute "make" to compile the tool. If you want to configure
# compiler options, please change Makefile.rules. If you need more
# informations, please read on.
#
# This file is part of parallelGBC, a parallel groebner basis computation tool.
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


# Run through all Makefile options
all: test

# Delete everything and rebuild all files
fresh:
	$(MAKE) clean;
	$(MAKE) all;

# The Makefile.rules includes all compiler flags and generic MAKE rules.
# Last is the reason why "all" and "fresh" are above.
include Makefile.rules

# Run all tests in the test directory
run: test
	$(MAKE) -C test run

# Compile all files in test
.PHONY: test
test: library
	$(MAKE) -C test all

# Build the library
library: src
	ar rcs lib/libf4.a src/*.o

# Compile all src files
.PHONY: src
src: 
	$(MAKE) -C src all

# Clean the whole project
clean:
	$(MAKE) -C src objclean
	$(MAKE) -C test objclean
	rm -f lib/*.a

# Do primitive checks
check: test
	./test/RunTests.sh
