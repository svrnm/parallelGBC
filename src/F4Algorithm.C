/*
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
#include "../include/F4Algorithm.H"
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <sstream>

using namespace std;
using namespace tbb;

namespace parallelGBC {

void F4::updatePairs(vector<Polynomial>& polys) 
{
	// Setup timer to measuere how long the 'update' takes.
	double timer = seconds();
	// The index of the current new element in the groebner basis.
	size_t t = groebnerBasis.size();
	// Since the groebner basis grows, we store the current size
	// of the groebner basis in a variable. 
	size_t is = groebnerBasis.size();
	// Iterate over all new polynomials, which have been found in the last
	// run of reduce(...)
	for(size_t i = 0; i < polys.size(); i++)
	{
		Polynomial& h = polys[i];
		
		// True if the current polynomial h will be inserted into the groebner basis
		bool insertIntoG = true;
		
		// Check if h should be inserted: h will not be inserted into g if there
		// is already a (new!) element in the groebner basis, which has a leading
		// term that divides the leading term of h
		// Note: An old element of the groebner basis means that this element was
		// already known before the last call of reduce, so h can't have a leading
		// term which is divisible by the leading term of the old element, because
		// if this would be the case, there would have been a reduction polynomial
		// reducing this leading term.
		for(size_t i = is; insertIntoG && i < groebnerBasis.size(); i++)
		{
			if(inGroebnerBasis[i] && h.LT().isDivisibleBy(groebnerBasis[i].LT())) 
			{
				insertIntoG = false;
			}
		}

		//Check the criteria only if h will be inserted into the groebner basis
		if(insertIntoG)
		{   
			// Do the first criterium: 
			// Cancel in P all pairs (i,j) which satisfy T(i,j) = T(i,j,t), T(i,t) != T(i,j) != T(j,t) [ B_t(i,j) ]
			F4PairSet P1(pairs.key_comp());
			for(set<F4Pair>::iterator it = pairs.begin(); it != pairs.end(); it++) 
			{
				if( !it->LCM.isDivisibleBy(h.LT())  || h.lcmLT(groebnerBasis[it->i]) == it->LCM  ||  h.lcmLT(groebnerBasis[it->j]) == it->LCM  ) { 
					P1.insert( *it );
				}
			}   
			swap(pairs, P1);

			// Do the second criterium:
			// Cancel in D1 each (i,t) for which a (j,t) exists s.t. T(i,t) is a proper multiple of T(j,t) [ M(i,t) ]
			vector<bool> D1(inGroebnerBasis.begin(), inGroebnerBasis.end());
			for(size_t i = 0; i < D1.size(); i++) 
			{
				for(size_t j = i+1; D1[i] && j < D1.size(); j++)
				{
					if(D1[j])				
					{
						Term a = h.lcmLT(groebnerBasis[i]);
						Term b = h.lcmLT(groebnerBasis[j]);
						if(a != b) {
							if(a.isDivisibleBy(b))
							{
								D1[i] = false;
							}
							if(b.isDivisibleBy(a)) 
							{
								D1[j] = false;
							}
						}
					}
				}
			}

			// Do the thrd criterium:
			// In each nonvoid subset { (j,t) | T(j,t) = tau } ...
			// Attention P2 is not a multiset, so each element is unique.
			set<F4Pair, F4Pair::comparator> P2(pairs.key_comp());
			for(size_t i = 0; i < D1.size(); i++)
			{
				if(D1[i])
				{
					Term LCM = groebnerBasis[i].lcmLT(h);
					// Create a new pair: LCM, i, t, marked, sugar degree
					F4Pair newpair( LCM, i, t, LCM == groebnerBasis[i].LT().mul(h.LT()), max(groebnerBasis[i].sugar() - groebnerBasis[i].LT().deg(), h.sugar() - h.LT().deg()) + LCM.deg() );
					pair<F4PairSet::iterator,bool> ret;
					ret = P2.insert( newpair );
					// If there is a marked pair for the given LCM, store this,
					// since all pairs with this LCM will be deleted.
					if(newpair.marked && ret.second)
					{
						P2.erase(ret.first);
						P2.insert(newpair);
					}
				}
			}

			// Finally delete all (i,t) with T(i)T(j) = T(i,t).
			for(set<F4Pair, F4Pair::comparator>::iterator it = P2.begin(); it != P2.end(); it++)
			{   
				if(!it->marked)
				{ 
					pairs.insert(*it);
				}  
			}

			// Check all old elements of the groebner basis if the current element
			// divides the leading term, so the old element is reducible and can
			// removed from the result set.
			for(size_t j = 0; j < groebnerBasis.size(); j++)
			{   
				if(inGroebnerBasis[j] && groebnerBasis[j].LT().isDivisibleBy(h.LT()))
				{   
					inGroebnerBasis[j] = false;
				}
			}
			// Insert h into the groebner basis
			groebnerBasis.push_back( h );
			inGroebnerBasis.push_back( insertIntoG );
			t++;
		}
	}
	updateTime += seconds() - timer;
}

void F4::gauss()
{
	// Iterate over all elements in the matrix matrix which
	// aren't having a leading term. These are the rows
	// in the upper part which have an odd index.
	for(size_t i = 0; i < upper/2; i++)
	{
		// BEGIN: Find the fimatrixt non zero entry in the current row, which will be the pivot element
		// Store the index in 'p' and the value in 'factor'
		size_t p = 0;
		bool found = false;
		coeffType factor = 0;
		for(p = 0; !found && p < matrix[i].size(); p++) {
			found = matrix[i][p] != 0;
			factor = matrix[i][p];
		}
		p--;
		// END: Find the first non zero entry

		// If no pivot was found, the current row is empty, store this information,
		// so this row is not relevant for the result.
		empty[i] = !found;
		if(found) {
			// We need a helper row, which stores the logarithms of all entries in
			// the current row matrix[i] since the mulSub operation expects that the 
			// operator row is given by the logarithms of the entries.
			coeffRow logRow(matrix[i].size(), 0);
			// Normalize the entries of matrix[i] if factor is not 1
			if(factor != 1) {
				factor = field->inv(factor);
				for(size_t j = p; j < matrix[i].size(); j++) {
					matrix[i][j] = field->mul(matrix[i][j], factor);
					logRow[j] = field->getFactor(matrix[i][j]);
				}
			} else {
				for(size_t j = p; j < matrix[i].size(); j++) {
					logRow[j] = field->getFactor(matrix[i][j]);
				}
			}
			
			// Iterate over all rows behind and before the current index
			// TODO: Do this with tbb:parallel_for()?
			#pragma omp parallel for num_threads( threads )
			for(size_t j = 1; j < upper/2; j++)
			{
				// Compute the relative index
				size_t k = (i+j)%(upper/2);
				// If the row matrix[k] has an entry at p, reduce this entry.
				if(matrix[k][p] != 0) {
					// The prefix are all entries before p.
					size_t prefix = (p/field->pad)*field->pad;
					// Reduce matrix[k] by matrix[i] multiplied by matrix[k][p].
					field->mulSub(matrix[k], logRow, matrix[k][p], prefix, matrix[k].size());
				}
			}
		}
	}

}


void F4::pReduceRange(coeffMatrix& rs, vector<size_t>& prefixes, vector<size_t>& suffixes, size_t i, tbb::blocked_range<size_t>& range) {
	// Iterate ofer the given range of operations.
	for(size_t j = range.begin(); j < range.end(); j++)
	{
		size_t target = ops[i].target( j );
		// Get the target row
		// Subtract from the target row the operator row multiplied with the factor. The prefixes and the suffixes
		// for the operator row are precomputed
		field->mulSub(rs[target], rs[ops[i].oper( j )], ops[i].factor( j ), prefixes[ ops[i].oper( j ) ],suffixes[ ops[i].oper( j ) ]);
		// Reduce the dependencies of the target by one
		deps[ target ]--;
		// If the current target is fully reduced (= has no dependencies), is not empty
		// and is not in a row, which will never be a target before gauss(), compute the
		// prefixes and suffixes and convert the row into log(row).
		if(deps[ target ] == 0 && rs[ target ].size() > 0 && (target > upper || target % 2 == 0 )) {
			size_t prefix, suffix;
			// Find the 0-padding of the target row at the beginning
			for(prefix = 0; prefix < rs[ target ].size() && rs[ target ][prefix] == 0; prefix++);
			prefixes[ target ] = ( prefix/field->pad )*field->pad;
			// Find the 0-padding of the target row at the end
			for(suffix = rs[target].size()-1; suffix > 0 && rs[target][suffix] == 0; suffix--) {
//				cout << suffix << " SUFFIX\n";
			}
			suffixes[target] = ( (suffix+field->pad)/field->pad )*field->pad;
			// Convert each entry of the target row into its logarithm value.
			for(size_t j = prefixes[target]; j < suffixes[target]; j++) {
				rs[target][j] = field->getFactor(rs[target][j]);
			}
		}
	}
}

void F4::pReduce()
{
	// The 'pad' is the padding which is required by the SSE vectors,
	// so even if a prefix may be some elements longer it has to be
	// cut down to a value which is divisible by pad
	size_t pad = field->pad;

	size_t end;
	size_t s = (( terms.size()+reduceBlockSize-1 )/ reduceBlockSize ) * reduceBlockSize;

	matrix.assign(upper/2, coeffRow());
	// Preset the 0-paddings of all rows to zero


	for(size_t start = 0; start < s; start+=reduceBlockSize) {
		if(start+reduceBlockSize < terms.size()) {
			end = start+reduceBlockSize;
		} else {
			end = terms.size();
		}

		coeffMatrix rs(rowCount, coeffRow(reduceBlockSize, 0) );
		vector<size_t> savedDeps(deps.begin(), deps.end());

		std::vector<size_t> prefixes(rs.size(), 0);
		std::vector<size_t> suffixes(rs.size(), 0);
		//cout << start << " - " << end << "\n";


		tbb::parallel_for(blocked_range<size_t>(start, end), F4SetupDenseRow(*this, rs, start));

		// Iterate over all matrix rows
		for(size_t i = 0; i < rs.size(); i++){
			// If the current row is already full reduced (=has no dependencies), is not
			// empty and is not in a row, which will never be atarget before gauss(),
			// compute the prefixes and suffixes and convert the row into log(row)
			// TODO: Merge this and the same procedure above into a function
			if(deps[i] == 0 && rs[i].size() > 0 && (i > upper || i % 2 == 0 )) {
				size_t prefix, suffix;
				// Find the 0-paddings of the row at the beginning
				for(prefix = 0; prefix < rs[i].size() && rs[i][prefix] == 0; prefix++);
				prefixes[i] = ( (prefix)/pad ) * pad; // prefix is the first non zero entry
				// Find the 0-padding of thw row at the end
				for(suffix = rs[i].size()-1; suffix >= prefix && rs[i][suffix] == 0; suffix--);
				suffixes[i] = ( (suffix+pad-1+1)/pad) *pad;
				// Convert each entry of the row into its logarithm value
				for(size_t j = prefixes[i]; j < suffixes[i]; j++) {
					rs[i][j] = field->getFactor(rs[i][j]);
				}
			}
		}

		// Iterate over all operation blocks
		for(size_t i = 0; i < ops.size(); i++) {
			// Split the given set of operations for parallel reduction
			tbb::parallel_for(blocked_range<size_t>(0, ops[i].size()), F4PReduceRange(*this, rs, prefixes, suffixes, i));
		}

		//for(size_t i = 0; i < upper; i++) {
		for(size_t i = 1, j = 0; i < upper; i+=2, j++) {
			matrix[j].insert(matrix[j].end(), rs[i].begin(), rs[i].end()); // copy rows to matrix;
		}
		std::copy(savedDeps.begin(), savedDeps.end(), deps.begin());
		/*std::fill(rs.begin(),rs.end(), coeffRow(reduceBlockSize, 0) );
		std::fill(prefixes.begin(),prefixes.end(),0);
		std::fill(suffixes.begin(),suffixes.end(),0);*/
	}
	
	
	
	rightSide.clear();
}

void F4::setupRow(Polynomial& current, Term& ir, size_t i, tbb::blocked_range<size_t>& range) 
{
	for(size_t j = range.begin(); j != range.end(); j++) {
		coeffType coeff = current.coeff(j);
		Term t = ir.mul(current.term(j));
		bool wontFound = termsUnordered.count(t) > 0;	
		bool found = false;
		// If there is not yet a pivot for t, try to create one
		if(!wontFound) {
			found = pivots.count(t) > 0;
			if(!found)
			{
				size_t element = 0;
				for(size_t k = 0; !found && k < groebnerBasis.size(); k++) 
				{
					if(inGroebnerBasis[k] && t.isDivisibleBy(groebnerBasis[k].LT())) {
						found = true;
						element = k;
					}
				}
				if(found) {
					tbb::concurrent_vector<std::pair<size_t, Term> >::iterator ret = rows.push_back(make_pair(element, t));
					pivots.insert(make_pair(t, std::distance(rows.begin(), ret)));
				}
			}
		}

		// Eliminate if possible
		if(found) {
			pivotOps[t].push_back( make_pair(i, coeff) );
		} else {
			if(!wontFound) {
				uint32_t temp = termCounter.fetch_and_increment();
				termsUnordered.insert( make_pair(t,temp) );
			}
			// Column -> (Entry, Row)
			rightSide[ termsUnordered[t] ].push_back( make_pair(coeff, i) );
		}
	}
}


void F4::setupDenseRow(coeffMatrix& rs, size_t offset, tbb::blocked_range<size_t>& range)
{
	for(size_t i = range.begin(); i != range.end(); i++) {
		for(size_t j = 0; j < rightSide[i].size(); j++) {
			rs[ rightSide[i][j].second ][i-offset] = rightSide[i][j].first;
		}
		rightSide[i].clear();
	}
}

void F4::prepare()
{
	double timer = seconds();
	// SELECTION
	vector<F4Pair> tmp(pairs.begin(), pairs.end());
	sort(tmp.begin(), tmp.end(), F4Pair::sugarComparator());
	currentDegree = tmp.begin()->sugar;
	if(verbosity & 16) {
		*out << "Sugar degree:\t" <<currentDegree << "\n";
	}
	size_t index;

	for(index = 0; index < tmp.size() && tmp[index].sugar == currentDegree; index++)
	{
		rows.push_back(make_pair(tmp[index].i, tmp[index].LCM));
		rows.push_back(make_pair(tmp[index].j, tmp[index].LCM));
		pivots.insert(make_pair(tmp[index].LCM, 2*index));
	}
	pairs.clear();
	pairs.insert(tmp.begin() + index, tmp.end());
	tmp.clear();
	// SELECTION END

	upper = 2*index;
	//rightSide.assign(rows.size(), tbb::concurrent_vector<pair<coeffType, uint32_t> >() );
	//double testtimer = 0;
	
	for(size_t i = 0; i < rows.size(); i++) 
	{
		Polynomial& current = groebnerBasis[ rows[i].first ];
		Term ir = rows[i].second.div(current.LT());
		// For the pivot rows (even rows and lower part) we start at 1 
		rightSide.grow_to_at_least( termsUnordered.size() + current.size() );
		tbb::parallel_for(blocked_range<size_t>((i > upper || i % 2 == 0 ? 1 : 0), current.size()), F4SetupRow(*this, current, ir, i));
	}

	// Clear unneeded data structures.
	rowCount = rows.size();
	rows.clear();
	pivotsOrdered.insert(pivots.begin(), pivots.end());
	pivots.clear();

	terms.insert(termsUnordered.begin(), termsUnordered.end());
	termsUnordered.clear();
	termCounter = 0;

	if(verbosity & 64) {
			*out << "Matrix (r x c):\t" << rowCount << " x " << terms.size() << "+" << pivotsOrdered.size() << "\n";
			size_t counter = 0;
			for(size_t i = 0; i < rightSide.size(); i++) {
				counter += rightSide[i].size();
			}
			*out << "RS density:\t" << ((double)counter / (double)(rowCount * terms.size())  ) << "\n";
	}

	ops.push_back( F4Operations() );
	deps.assign(rowCount, 0);
	size_t oCounter = 0;
#if PGBC_SORTING == 1 
	vector<size_t> l(rowCount,0);
	for(map<Term, uint32_t, Term::comparator>::reverse_iterator it = pivotsOrdered.rbegin(); it != pivotsOrdered.rend(); it++) 
	{
		uint32_t o = it->second;
		vector<pair<uint32_t, coeffType> >& entries = pivotOps[it->first];
		for(size_t i = 0; i < entries.size(); i++)
		{
			uint32_t t = entries[i].first;
			if(l[ o ] > l[t]) {
				l[t] = l[o];
			}
			oCounter++;
			ops[ l[t] ].push_back( t,o,entries[i].second );
			//oCount++;
			deps[t]++; // target t has deps[t] operations do be done before it can be used as operator
			l[t]++; // one operation per level, attention this also affects the following if-statements

			if(l[t] >= ops.size()) {
				ops.push_back( F4Operations() );
			}
		}
		pivotOps[it->first].clear();
	}
#else
	size_t l = 0;
	for(map<Term, uint32_t, Term::comparator>::reverse_iterator it = pivotsOrdered.rbegin(); it != pivotsOrdered.rend(); it++) 
	{
		uint32_t o = it->second;
		vector<pair<uint32_t, coeffType> >& entries = pivotOps[it->first];
		for(size_t i = 0; i < entries.size(); i++)
		{
			size_t t = entries[i].first;
			oCounter++;
			ops[ l ].push_back( t, o, entries[i].second );
			deps[t]++;
		}
		l++;
		ops.push_back( F4Operations() );
	}
#endif
		if(verbosity & 64) {
			*out << "Operations:\t" << oCounter << "\n";
			*out << "in levels:\t" << ops.size() << "\n";
			*out << "Op. Density:\t" << ( (double)oCounter /  (double)(rowCount * pivotsOrdered.size()) ) << "\n";
		}
	pivotOps.clear();
	pivotsOrdered.clear();
	ops.pop_back();
	
	prepareTime += seconds() - timer;
}

void F4::reduce(vector<Polynomial>& polys)
{
	prepare();

	// ELIMINATE
	double timer = seconds();
	pReduce();
	
	
	ops.clear();
	//timer = seconds();

	empty.assign(upper/2, false); // too large, FIX?
	//double gt = seconds();	
	
	coeffRow temp;
	// Finally sort the resulting matrix using terms. The previous intermediate matrix was not in in term ordering
	size_t nCounter = 0;
	for(size_t i = 0; i < upper/2; i++) {
		temp.assign(matrix[i].size(), 0);
		size_t j = 0;
		for(map<Term, uint32_t, Term::comparator>::iterator it = terms.begin(); it != terms.end(); it++) {
			if(matrix[i][it->second] != 0) {
				nCounter++;
			}
			temp[j] = matrix[i][it->second];
			j++;
		}
		matrix[i].swap(temp);
	}

	if(verbosity & 64) {
		cout << "Final Matrix:\t" << (upper/2) << "x" << matrix[0].size() << "\n";
		cout << "Entries:\t" << nCounter << "\n";
		cout << "Density:\t" << ((double) nCounter / (double)( (upper/2) * matrix[0].size() ) ) << "\n";
	}
	
	gauss();
	
	reductionTime += seconds()-timer;
	if(verbosity & 32) {
		*out << "Red. step (s):\t" << seconds()-timer << "\n";
	}
	for(size_t i = 0; i < upper/2; i++)
	{
		if(!empty[i])
		{
			Polynomial p(currentDegree);
			size_t j = 0;
			for(map<Term, uint32_t, Term::comparator>::iterator it = terms.begin(); it != terms.end(); it++) 
			{
				if(matrix[i][j] != 0)
				{
					p.push_back(matrix[i][j], it->first);
				}
				j++;
			}
			if(verbosity & 128) {
				*out << "NEW POLY:\t" << p << "\n";
			}
			polys.push_back( p );
		}
	}
	// Reset matrix.
	empty.clear();
	terms.clear();
	matrix.clear();
}

vector<Polynomial> F4::operator()(vector<Polynomial>& generators, const TOrdering* o, CoeffField* field, int threads, int verbosity, std::ostream& output)
{
	double start = seconds();
	this->field = field;
	this->threads = threads;
	this->O = o;
	this->verbosity = verbosity;
	this->out = &output;
	this->reduceBlockSize = 2048;

	termCounter = 0;
	pairs = F4PairSet( F4Pair::comparator(O) );
	terms = map<Term, uint32_t, Term::comparator>( Term::comparator(O, true) );
	pivotsOrdered = map<Term, uint32_t, Term::comparator>( Term::comparator(O, true) );

	updateTime = 0;
	prepareTime = 0;
	reductionTime = 0;

	tbb::task_scheduler_init init(threads);

	sort(generators.begin(), generators.end(), Polynomial::comparator(o, true));		

	//normalize
	for_each(generators.begin(), generators.end(), bind(mem_fn(&Polynomial::normalize), _1, field));

	updatePairs(generators);

	while( !pairs.empty() ) {
		vector<Polynomial> polys;
		reduce(polys);
		if(!polys.empty()) {
			updatePairs(polys);
		}
	}

	if(verbosity & 2) {
		*out << "Reduction (s): \t" << reductionTime << "\n";
	}
	if(verbosity & 4) {
		*out << "Prepare (s): \t" << prepareTime << "\n";
	}
	if(verbosity & 8) {
		*out << "Update (s): \t" << updateTime << "\n";
	}
	vector<Polynomial> result;
	for(size_t i = 0; i < groebnerBasis.size(); i++)
	{
		if(inGroebnerBasis[i])
			result.push_back(groebnerBasis[i]);
	}
	if(verbosity & 1) {
		*out << "Runtime (s):\t" << seconds() - start << "\n";
	}

	/* Enable this to collect timings
	*out << threads << ";";
	*out << seconds() - start << ";";
	*out << reductionTime << ";";
	*out << updateTime << ";";
	*out << prepareTime << "\n";
	*/

	return result;
}

}
