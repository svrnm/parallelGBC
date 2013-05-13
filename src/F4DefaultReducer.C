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
using namespace boost;

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
				newPivots.push_back(make_pair(p,i));
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
				if(doSimplify > 0) {
					savedRows[index].push_back( make_pair(j + offset, row[j]) );
				}
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

#if PGBC_WITH_MPI == 1
		size_t columnCount = rightSide.size();
#else
		size_t columnCount = terms.size();
#endif

		size_t s = (( columnCount+reduceBlockSize-1 )/ reduceBlockSize ) * reduceBlockSize;

		matrix.assign(upper/2, coeffRow());

		for(size_t start = 0; start < s; start+=reduceBlockSize) {
			size_t last = reduceBlockSize;
			if(start+reduceBlockSize < columnCount) {
				end = start+reduceBlockSize;
			} else {
				end = columnCount;
				last = end - start;
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

			for(size_t i = 1, j = 0; i < upper; i+=2, j++) {
				matrix[j].insert(matrix[j].end(), rs[i].begin(), rs[i].begin() + last); // copy rows to matrix;
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
#if PGBC_WITH_MPI == 1
					int dst = termsUnordered[t] % f4->world.size();
					toSend[dst].push_back(  make_pair( termsUnordered[t]/f4->world.size(), make_pair(coeff, (uint32_t)i)) );
#else
					rightSide[ termsUnordered[t] ].push_back( make_pair(coeff, i) );
#endif		
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
			upper *= 2;

#if PGBC_WITH_MPI == 1
			int status = 0;
			std::vector<std::pair<uint32_t, std::pair<coeffType, uint32_t> > > ownPart;
#endif
			for(size_t i = 0; i < rows.size(); i++) 
			{
#if PGBC_WITH_MPI == 1
				toSend.assign(f4->world.size(), tbb::concurrent_vector<std::pair<uint32_t, std::pair<coeffType, uint32_t> > >());
#endif
				Polynomial current = f4->groebnerBasis[ rows[i].first ];
				Term ir = rows[i].second.div(current.LT());
				if(doSimplify == 2) {
					rowOriginDB.push_back( make_pair( rows[i].first, ir ) );
					std::pair<Term, Polynomial> s = simplifyDB->search(rows[i].first, ir);
					if(s.first != ir) {
						ir = s.first;
						current = s.second;
					} 
				} else if(doSimplify == 1) {
						simplify->search(ir, current);
						rowOrigin.push_back( make_pair( ir,current ) );
				} 
#if PGBC_WITH_MPI == 0
				rightSide.grow_to_at_least( termsUnordered.size() + current.size() );
#endif
				tbb::parallel_for(blocked_range<size_t>((i > upper || i % 2 == 0 ? 1 : 0), current.size()), F4SetupRow(*this, current, ir, i));
#if PGBC_WITH_MPI == 1
				double mpiTimer = F4Logger::seconds();
				for(size_t k = 0; k < toSend.size(); k++) {
					toSendCopy.push_back( std::vector<std::pair<uint32_t, std::pair<coeffType, uint32_t> > >(toSend[k].begin(), toSend[k].end()) );
				}  	
				mpi::scatter(f4->world, toSendCopy, ownPart, 0);
				for(size_t k = 0; k < ownPart.size(); k++) {
					rightSide.grow_to_at_least( ownPart[k].first+1 );
					rightSide[ ownPart[k].first ].push_back( ownPart[k].second );
				}
				toSend.clear();
				toSendCopy.clear();
				ownPart.clear();
				if(i == rows.size() - 1) {
					status = 1;
				}
				mpi::broadcast(f4->world, status, 0);
				f4->log->mpiTime += F4Logger::seconds() - mpiTimer;
#endif
			}
			rowCount = rows.size();
			rows.clear();
			pivotsOrdered.insert(pivots.begin(), pivots.end());
			pivots.clear();

			terms.insert(termsUnordered.begin(), termsUnordered.end());
			termsUnordered.clear();

#if PGBC_WITH_MPI == 1
  mpi::broadcast(f4->world, upper, 0);
  mpi::broadcast(f4->world, rowCount, 0);
#endif

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

			if(doSimplify > 0) {
				savedRows.assign(rowCount, vector<pair<uint32_t, coeffType> >());
      }

			f4->log->prepareTime += F4Logger::seconds() - timer;
		}

		void F4DefaultReducer::reduce(vector<Polynomial>& polys, degreeType currentDegree)
		{
#if PGBC_WITH_MPI == 1
			int status = 1;
			mpi::broadcast(f4->world, status, 0);
#endif
			
			prepare();
#if PGBC_WITH_MPI == 1
			size_t opCount = ops.size();
			mpi::broadcast(f4->world, opCount, 0);
			for(size_t i = 0; i < ops.size(); i++) {
				mpi::broadcast(f4->world, ops[i], 0);
			}
			mpi::broadcast(f4->world, deps, 0);
#endif

			// ELIMINATE
			double timer = F4Logger::seconds();
			pReduce();

#if PGBC_WITH_MPI == 1
			double mpiTimer = F4Logger::seconds();
			if(f4->world.size() > 1) {
				vector<coeffMatrix> gatheredMatrix;
				mpi::gather(f4->world, matrix, gatheredMatrix, 0);
				// Reconstructing matrix
				coeffMatrix m(upper/2, coeffRow(terms.size(), 0));
				for(size_t i = 0; i < gatheredMatrix.size(); i++) {
					for(size_t j = 0; j < gatheredMatrix[i].size(); j++) {
						size_t n = gatheredMatrix[i][j].size();
						// Local matrix may be to large since the setup of rightSide may be to large
						if(n > (termCounter+f4->world.size()-1) / f4->world.size()) {
							n = (termCounter+f4->world.size()-1) / f4->world.size();
						}
						for(size_t k = 0 ; k < n; k++) {
							m[j][k * f4->world.size() + i] = gatheredMatrix[i][j][k];
						}
					}
				}
				matrix.swap(m);
				gatheredMatrix.clear();

				if(doSimplify > 0) {
					vector<vector<pair<uint32_t, coeffType > > > s(rowCount, vector<pair<uint32_t, coeffType > >());
					vector<vector<vector<pair<uint32_t, coeffType> > > > gatheredSavedRows; 
					mpi::gather(f4->world, savedRows, gatheredSavedRows, 0);
					for(size_t i = 0; i < gatheredSavedRows.size(); i++) {
						for(size_t j = 0; j < gatheredSavedRows[i].size(); j++) {
							for(size_t k = 0; k < gatheredSavedRows[i][j].size(); k++) {
								s[j].push_back( make_pair(gatheredSavedRows[i][j][k].first * f4->world.size() + i, gatheredSavedRows[i][j][k].second) );
							}
						}
					}
					savedRows.swap(s);
					gatheredSavedRows.clear();
				}
			}
			f4->log->mpiTime += F4Logger::seconds()-mpiTimer;
#endif

			termCounter = 0;
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
				size_t aligned = (( matrix[i].size()+reduceBlockSize-1 )/ reduceBlockSize ) * reduceBlockSize;
				temp.assign(aligned, 0);
				for(size_t j = 0; j < termMapping.size(); j++) {
					if(matrix[i][j] != 0) {
						nCounter++;
					}
					temp[termMapping[j]] = matrix[i][j];
				}
				matrix[i].swap(temp);
			}

			if(f4->log->verbosity & 64) {
				(*f4->log->out) << "Final Matrix:\t" << (upper/2) << "x" << terms.size() << "\n";
				(*f4->log->out) << "Entries:\t" << nCounter << "\n";
				(*f4->log->out) << "Density:\t" << ((double) nCounter / (double)( (upper/2) * terms.size() ) ) << "\n";
			}

			gauss();

			f4->log->reductionTime += F4Logger::seconds()-timer;
			if(f4->log->verbosity & 32) {
				*(f4->log->out) << "Red. step (s):\t" << F4Logger::seconds()-timer << "\n";
			}

			size_t t = f4->groebnerBasis.size();

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
					Term one = p.LT().getOne();
					if(doSimplify == 2) {
						simplifyDB->insert(t, one, p); // this breaks the independency from Reducer and Algorithm...
					}
					t++;
				}
			}
			if(f4->log->verbosity & 64) {
				*(f4->log->out) << "Polys:\t" << polys.size() << "\n";
			}


			timer = F4Logger::seconds();
			if(doSimplify == 2) {
				#pragma omp parallel for num_threads ( f4->threads )
        for(size_t i = 0; i < savedRows.size(); i++) {
          if(!savedRows[i].empty() && simplifyDB->check(rowOriginDB[i].first, rowOriginDB[i].second) > savedRows[i].size()) {
            Polynomial p(currentDegree);
						std::vector<coeffType> tmp(terms.size(), 0);
						for(size_t j = 0; j < savedRows[i].size(); j++) {
							tmp[termMapping[ savedRows[i][j].first  ]] = savedRows[i][j].second;
						}

						// POST Reduce
#if PGBC_POST_REDUCE == 1
						for(size_t k = 0; k < newPivots.size(); k++) {
							uint32_t pivot = newPivots[k].first;
							uint32_t index = newPivots[k].second;
							if(tmp[pivot] != 0) {
								coeffRow logRow(matrix[index].size(), 0);
								for(size_t j = pivot; j < matrix[index].size(); j++) {
									logRow[j] = f4->field->getFactor(matrix[index][j]);
								}
								size_t prefix = (pivot/f4->field->pad)*f4->field->pad;
								f4->field->mulSub(tmp, logRow, tmp[pivot], prefix, tmp.size());
							}
						}
#endif



            p.push_back(1, f4->groebnerBasis[ rowOriginDB[i].first ].LT().mul(rowOriginDB[i].second)); 
						size_t j = 0;
						for(map<Term, uint32_t, Term::comparator>::iterator it = terms.begin(); it != terms.end(); it++, j++) {
							if(tmp[j] != 0) { 
								p.push_back(tmp[j], it->first); 
							}
						}

            simplifyDB->insert(rowOriginDB[i].first, rowOriginDB[i].second, p); 
          }   
        }
			} else if(doSimplify == 1) {
				#pragma omp parallel for num_threads ( f4->threads )
				for(size_t i = 0; i < savedRows.size(); i++) {
					if(!savedRows[i].empty()) {
						Polynomial p(currentDegree);
						std::vector<coeffType> tmp(terms.size(), 0);
						for(size_t j = 0; j < savedRows[i].size(); j++) {
							tmp[termMapping[ savedRows[i][j].first  ]] = savedRows[i][j].second;
						}

#if PGBC_POST_REDUCE == 1
						for(size_t k = 0; k < newPivots.size(); k++) {
							uint32_t pivot = newPivots[k].first;
							uint32_t index = newPivots[k].second;
							if(tmp[pivot] != 0) {
								coeffRow logRow(matrix[index].size(), 0);
								for(size_t j = pivot; j < matrix[index].size(); j++) {
									logRow[j] = f4->field->getFactor(matrix[index][j]);
								}
								size_t prefix = (pivot/f4->field->pad)*f4->field->pad;
								f4->field->mulSub(tmp, logRow, tmp[pivot], prefix, tmp.size());
							}
						}
#endif

						p.push_back(1, rowOrigin[i].second.LT().mul(rowOrigin[i].first));
						size_t j = 0;
						for(map<Term, uint32_t, Term::comparator>::iterator it = terms.begin(); it != terms.end(); it++, j++) {
							if(tmp[j] != 0) { 
								p.push_back(tmp[j], it->first); 
							}
						}
						simplify->insert(rowOrigin[i].first, rowOrigin[i].second, p);
					}
				}
			}
			f4->log->simplifyTime += F4Logger::seconds()-timer;

			// Reset matrix.
			rowOrigin.clear();
			rowOriginDB.clear();
			empty.clear();
			terms.clear();
			matrix.clear();
			rightSide.clear();
			upper = 0;
			newPivots.clear();
		}

		void F4DefaultReducer::addSPolynomial(size_t i, size_t j, Term& lcm) {
			rows.push_back(make_pair(i, lcm));
			rows.push_back(make_pair(j, lcm));
			pivots.insert(make_pair(lcm, 2*upper));
			upper++;
		}

		void F4DefaultReducer::finish() {
#if PGBC_WITH_MPI == 1
			int status = 2;
			mpi::broadcast(f4->world, status, 0);
#endif
		}

		void F4DefaultReducer::client() {
#if PGBC_WITH_MPI == 1
			std::vector<std::pair<uint32_t, std::pair<coeffType, uint32_t> > > ownPart;
			int status;
			mpi::broadcast(f4->world, status, 0); 	
			while(status != 2) {
				do {
					mpi::scatter(f4->world, toSendCopy, ownPart, 0);
					for(size_t i = 0; i < ownPart.size(); i++) {
						rightSide.grow_to_at_least( ownPart[i].first+1 );
						rightSide[ ownPart[i].first ].push_back( ownPart[i].second );
					}   
					ownPart.clear();
					mpi::broadcast(f4->world, status, 0); 
				} while( status == 0 );
				mpi::broadcast(f4->world, upper, 0); 
				mpi::broadcast(f4->world, rowCount, 0); 
	
				
				if(doSimplify > 0) {
					savedRows.assign(rowCount, vector<pair<uint32_t, coeffType> >());
				}

				size_t opCount = 0;

				mpi::broadcast(f4->world, opCount, 0);

				ops.assign(opCount, F4Operations() );
				for(size_t i = 0; i < ops.size(); i++) {
					mpi::broadcast(f4->world, ops[i], 0); 
				}
				mpi::broadcast(f4->world, deps, 0); 
				if(!rightSide.empty()) {
					pReduce();
				}   
				ops.clear();
				deps.clear();
				mpi::gather(f4->world, matrix, 0); 
				if(doSimplify > 0) {
					mpi::gather(f4->world, savedRows, 0);
				}
				matrix.clear();
				mpi::broadcast(f4->world, status, 0); 	
			}
#endif
		}
	}
