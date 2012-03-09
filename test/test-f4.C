#include "../include/F4.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/regex.hpp> 

using namespace boost;
using namespace std;

int main(int argc, char* argv[]) {
	if(argc < 2) {
		cerr << "Please provide a file.\n";
		exit(-1);
	}
	int threads = __F4_ALGORITHM_THREADS;
	if(argc > 2) {
		istringstream( argv[2] ) >> threads;
	}
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
	boost::regex expression("x\\[(\\d*)\\]");
	std::vector<std::string> numbers;
	std::copy(boost::sregex_token_iterator(s.begin(), s.end(), expression, 1), boost::sregex_token_iterator(),   std::back_inserter(numbers));
	degreeType max = 1, val;
	for(size_t i = 0; i < numbers.size(); i++) {
		istringstream ( numbers[i] ) >> val;
		if( val > max ) max = val;
	}
	cout << "Found " << max << "indeterminants\n";
	TOrdering* o = new DegRevLexOrdering(max);
	TMonoid m(max);
	filestr.close();
	std::cout << "Testing ...\n";
	CoeffField* cf = new CoeffField(32003);
	cout << t << "\n";
	vector<Polynomial> list = Polynomial::createList(t, m);
	for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::order), _1, o));
	for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::bringIn), _1, cf, false));
	for_each(list.begin(), list.end(), bind(mem_fn(&Polynomial::normalize), _1, cf));
	F4 f4;
	
	vector<Polynomial> result = f4(list, o, cf, threads);
	cout << result.size() << "\n";
	delete o;
	delete cf;
}


