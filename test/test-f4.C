/*
 *  Example and test file for parallelGBC. To use it execute
 *
 *  	# ./test-f4 input.txt <NUM_OF_PROCS>
 *
 *  where input.txt is a file providing polynomials in one line, such as
 * 
 *  	x[1]+x[2]+x[3], x[1]*x[2]+x[1]*x[3]+x[2]*x[3], x[1]*x[2]*x[3]-1
 *
 *  and <NUM_OF_PROCS> is the number of processors you want to use in parallel.
 *
 *  If you want to use the library for your own code, please look below for an 
 *  usage example!
 *
 *  This file is part of parallelGBC, a parallel groebner basis computation tool.
 *
 *  parallelGBC is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  parallelGBC is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with parallelGBC.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "../include/F4.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/regex.hpp> 

using namespace boost;
using namespace std;

int main(int argc, char* argv[]) {
	// Input stuff, example below ...
	
	// Is there a file provided?
	if(argc < 2) {
		cerr << "Please provide a file and a optional number of threads.\n";
		exit(-1);
	}
	
	// If the second parameter provides the number of
	// threads, use it, if not use the default value.
	int threads = __F4_ALGORITHM_THREADS;
	if(argc > 2) {
		istringstream( argv[2] ) >> threads;
	}
	int verbosity = 0;
	if(argc > 3) {
		istringstream( argv[3] ) >> verbosity;
	}
	// Read the provided input file. Example still below.
	fstream filestr (argv[1], fstream::in);
	std::string s,t;
	t = "";
	if(filestr.is_open()) {
		while(filestr >> s) {
			t += s;
		}
	} else {
		cerr << "Could not open file\n";
		exit(-1);
	}
	// Count the indeterminants automaticly
	boost::regex expression("x\\[(\\d*)\\]");
	std::vector<std::string> numbers;
	std::copy(boost::sregex_token_iterator(s.begin(), s.end(), expression, 1), boost::sregex_token_iterator(),   std::back_inserter(numbers));
	degreeType max = 1, val;
	for(size_t i = 0; i < numbers.size(); i++) {
		istringstream ( numbers[i] ) >> val;
		if( val > max ) max = val;
	}
	filestr.close();

	// [[[ EXAMPLE ]]] //
	// Use the following example as guidance if you want to use it for your own code:
	
	// 1. Create a term ordering for 'max' indeterminants, e.g. if you use x_1...x_n, the parameter is n.
	TOrdering* o = new DegRevLexOrdering(max);
	// 2. Create a power product monoid for the terms. Pay attention that ordering and monoid match.
	TMonoid m(max);
	// 3. Create a coefficient field.
	CoeffField* cf = new CoeffField(32003);
	// 4. Read in the polynomials from string 't'. The second parameter is the power product monoid.
	vector<Polynomial> list = Polynomial::createList(t, m);
	// 5. Before you can compute the groebner basis, you have to order your polynomials by term
	// ordering and have to bring in the coefficients to your coefficient field. Finally you have
	// to normalize your polynomials. Remark: This step will be merged into f4(...) in a later release,
	// doing everything twice shouldn't harm.
	for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::order), _1, o));
	for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::bringIn), _1, cf, false));
	// Create the f4 computer.
	F4 f4;
	// Compute the groebner basis for the polynomials in 'list' with 'threads' threads/processors 
	vector<Polynomial> result = f4(list, o, cf, threads, verbosity);
	// Return the size of the groebner basis
	//cout << "Size of GB:\t" << result.size() << "\n";
	for(size_t i = 0; i < result.size(); i++) {
		if(i > 0) {
			cout << ", ";
		}
		cout << result[i];
	}
	cout << "\n";
	// Clean up your memory
	delete o;
	delete cf;
}


