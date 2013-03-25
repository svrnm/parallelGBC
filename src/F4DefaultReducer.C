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
#include "../include/F4DefaultReducer.H"
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

	void F4DefaultReducer::gauss()
	{
		// Iterate over all elements in the matrix matrix which
		// aren't having a leading term. These are the rows
		// in the upper part which have an odd index.
		for(size_t i = 0; i < upper/2; i++)
		{
			// BEGIN: Find the first matrix with non zero entry in the current row, which will be the pivot element
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
					factor = f4->field->inv(factor);
					for(size_t j = p; j < matrix[i].size(); j++) {
						matrix[i][j] = f4->field->mul(matrix[i][j], factor);
						logRow[j] = f4->field->getFactor(matrix[i][j]);
					}
				} else {
					for(size_t j = p; j < matrix[i].size(); j++) {
						logRow[j] = f4->field->getFactor(matrix[i][j]);
					}
				}

				// Iterate over all rows behind and before the current index
				// TODO: Do this with tbb:parallel_for()?
#pragma omp parallel for num_threads( f4->threads )
				for(size_t j = 1; j < upper/2; j++)
				{
					// Compute the relative index
					size_t k = (i+j)%(upper/2);
					// If the row matrix[k] has an entry at p, reduce this entry.
					if(matrix[k][p] != 0) {
						// The prefix are all entries before p.
						size_t prefix = (p/f4->field->pad)*f4->field->pad;
						// Reduce matrix[k] by matrix[i] multiplied by matrix[k][p].
						f4->field->mulSub(matrix[k], logRow, matrix[k][p], prefix, matrix[k].size());
					}
				}
			}
		}

	}


	void F4DefaultReducer::prepareOperator(coeffRow& row, size_t index, size_t& prefix, size_t& suffix, size_t offset) {
		for(prefix = 0; prefix < row.size() && row[prefix] == 0; prefix++);
		prefix = ( prefix/f4->field->pad )*f4->field->pad;
		for(suffix = row.size()-1; suffix >= prefix && row[suffix] == 0; suffix--);
		suffix = ( (suffix+f4->field->pad)/f4->field->pad )*f4->field->pad;
		for(size_t j = prefix; j < suffix; j++) {
			if(row[j] != 0) {
				//rightSide[ j ].push_back( make_pair(row[j], index) );
				row[j] = f4->field->getFactor(row[j]);
			}
		}
	}

	void F4DefaultReducer::pReduceRange(coeffMatrix& rs, vector<size_t>& prefixes, vector<size_t>& suffixes, size_t i, size_t offset, tbb::blocked_range<size_t>& range) {
		// Iterate over the given range of operations.
		for(size_t j = range.begin(); j < range.end(); j++)
		{
			size_t target = ops[i].target( j );
			size_t oper = ops[i].oper( j );
			// Get the target row
			// Subtract from the target row the operator row multiplied with the factor. The prefixes and the suffixes
			// for the operator row are precomputed
			f4->field->mulSub(rs[target], rs[oper], ops[i].factor( j ), prefixes[oper],suffixes[oper]);
			// Reduce the dependencies of the target by one
			deps[ target ]--;
			// If the current target is fully reduced (= has no dependencies), is not empty
			// and is not in a row, which will never be a target before gauss(), compute the
			// prefixes and suffixes and convert the row into log(row).
			if(deps[ target ] == 0 && rs[ target ].size() > 0 && (target > upper || target % 2 == 0 )) {
				prepareOperator(rs[target], target, prefixes[target], suffixes[target], offset);
			}
		}
	}

	void F4DefaultReducer::pReduce()
	{
		size_t end;
		size_t s = (( terms.size()+reduceBlockSize-1 )/ reduceBlockSize ) * reduceBlockSize;

		matrix.assign(upper/2, coeffRow());

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

			tbb::parallel_for(blocked_range<size_t>(start, end), F4SetupDenseRow(*this, rs, start));

			// Iterate over all matrix rows
			for(size_t i = 0; i < rs.size(); i++){
				// If the current row is already full reduced (=has no dependencies), is not
				// empty and is not in a row, which will never be atarget before gauss(),
				// compute the prefixes and suffixes and convert the row into log(row)
				// TODO: Merge this and the same procedure above into a function
				if(deps[i] == 0 && rs[i].size() > 0 && (i > upper || i % 2 == 0 )) {
					prepareOperator(rs[i], i, prefixes[i], suffixes[i], start);
				}
			}

			// Iterate over all operation blocks
			for(size_t i = 0; i < ops.size(); i++) {
				// Split the given set of operations for parallel reduction
				tbb::parallel_for(blocked_range<size_t>(0, ops[i].size()), F4PReduceRange(*this, rs, prefixes, suffixes, i, start));
			}

			//for(size_t i = 0; i < upper; i++) {
			for(size_t i = 1, j = 0; i < upper; i+=2, j++) {
				matrix[j].insert(matrix[j].end(), rs[i].begin(), rs[i].end()); // copy rows to matrix;
			}
			std::copy(savedDeps.begin(), savedDeps.end(), deps.begin());
		}



		}

		void F4DefaultReducer::setupRow(Polynomial& current, Term& ir, size_t i, tbb::blocked_range<size_t>& range) 
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
						for(size_t k = 0; !found && k < f4->groebnerBasis.size(); k++) 
						{
							if(f4->inGroebnerBasis[k] && t.isDivisibleBy(f4->groebnerBasis[k].LT())) {
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


		void F4DefaultReducer::setupDenseRow(coeffMatrix& rs, size_t offset, tbb::blocked_range<size_t>& range)
		{
			for(size_t i = range.begin(); i != range.end(); i++) {
				for(size_t j = 0; j < rightSide[i].size(); j++) {
					rs[ rightSide[i][j].second ][i-offset] = rightSide[i][j].first;
				}
				rightSide[i].clear();
			}
		}

		void F4DefaultReducer::prepare()
		{
			double timer = F4Logger::seconds();
			// SELECTION

			// SELECTION END

			upper *= 2;
			//rightSide.assign(rows.size(), tbb::concurrent_vector<pair<coeffType, uint32_t> >() );
			//double testtimer = 0;

			for(size_t i = 0; i < rows.size(); i++) 
			{
				Polynomial current = f4->groebnerBasis[ rows[i].first ];
				Term ir = rows[i].second.div(current.LT());
				rowOrigin.push_back( make_pair( rows[i].first, ir ) );
				rightSide.grow_to_at_least( termsUnordered.size() + current.size() );
				tbb::parallel_for(blocked_range<size_t>((i > upper || i % 2 == 0 ? 1 : 0), current.size()), F4SetupRow(*this, current, ir, i));
			}
			rowCount = rows.size();
			rows.clear();
			pivotsOrdered.insert(pivots.begin(), pivots.end());
			pivots.clear();

			terms.insert(termsUnordered.begin(), termsUnordered.end());
			termsUnordered.clear();
			termCounter = 0;


			if(f4->log->verbosity & 64) {
				*(f4->log->out) << "Matrix (r x c):\t" << rowCount << " x " << terms.size() << "+" << pivotsOrdered.size() << "\n";
				size_t counter = 0;
				for(size_t i = 0; i < rightSide.size(); i++) {
					counter += rightSide[i].size();
				}
				*(f4->log->out) << "RS density:\t" << ((double)counter / (double)(rowCount * terms.size())  ) << "\n";
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
			if(f4->log->verbosity & 64) {
				*(f4->log->out) << "Operations:\t" << oCounter << "\n";
				*(f4->log->out) << "in levels:\t" << ops.size() << "\n";
				*(f4->log->out) << "Op. Density:\t" << ( (double)oCounter /  (double)(rowCount * pivotsOrdered.size()) ) << "\n";
			}
			pivotOps.clear();
			pivotsOrdered.clear();
			ops.pop_back();

			f4->log->prepareTime += F4Logger::seconds() - timer;
		}

		void F4DefaultReducer::reduce(vector<Polynomial>& polys, degreeType currentDegree)
		{

			prepare();

			// ELIMINATE
			double timer = F4Logger::seconds();
			pReduce();
			ops.clear();
			empty.assign(upper/2, false);

			coeffRow temp;

			// Finally sort the resulting matrix using terms. The previous intermediate matrix was not in in term ordering
			termMapping.assign(terms.size(), 0);
			size_t c = 0;
			for(map<Term, uint32_t, Term::comparator>::iterator it = terms.begin(); it != terms.end(); it++, c++) {
				termMapping[it->second] = c; 
			}

			size_t nCounter = 0;
			for(size_t i = 0; i < upper/2; i++) {
				temp.assign(matrix[i].size(), 0);
				for(size_t j = 0; j < termMapping.size(); j++) {
					if(matrix[i][j] != 0) {
						nCounter++;
					}
					temp[termMapping[j]] = matrix[i][j];
				}
				matrix[i].swap(temp);
			}

			if(f4->log->verbosity & 64) {
				(*f4->log->out) << "Final Matrix:\t" << (upper/2) << "x" << matrix[0].size() << "\n";
				(*f4->log->out) << "Entries:\t" << nCounter << "\n";
				(*f4->log->out) << "Density:\t" << ((double) nCounter / (double)( (upper/2) * matrix[0].size() ) ) << "\n";
			}

			gauss();

			f4->log->reductionTime += F4Logger::seconds()-timer;
			if(f4->log->verbosity & 32) {
				*(f4->log->out) << "Red. step (s):\t" << F4Logger::seconds()-timer << "\n";
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
					if(f4->log->verbosity & 128) {
						*(f4->log->out) << p << "\n";
					}
					polys.push_back( p );
				}
			}
			if(f4->log->verbosity & 64) {
				*(f4->log->out) << "Polys:\t" << polys.size() << "\n";
			}

			timer = F4Logger::seconds();
			// Reset matrix.
			rowOrigin.clear();
			empty.clear();
			terms.clear();
			matrix.clear();
			rightSide.clear();
			upper = 0;
		}

		void F4DefaultReducer::addSPolynomial(size_t i, size_t j, Term& lcm) {
			rows.push_back(make_pair(i, lcm));
			rows.push_back(make_pair(j, lcm));
			pivots.insert(make_pair(lcm, 2*upper));
			upper++;
		}
	}
