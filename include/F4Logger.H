/**
 **********************
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
#ifndef F4_LOGGER_H
#define F4_LOGGER_H
#include <sys/time.h>
#include <iostream>
#include <sys/types.h>
#include <sys/sysinfo.h>

namespace parallelGBC {

	class F4Logger {
		public:
			/**
			 * Verbosity, which can be changed during runtime, nothing which should
			 * influence performance:
			 *
			 * 1 - Overall Runtime of the computation
			 * 2 - Sum of reduction time
			 * 4 - Sum of matrix setup/prepare time
			 * 8 - Sum of update time
			 * 16 - Print sugar degree during reduction step
			 * 32 - Print time of reduction step
			 * 64 - Print matrix size during reduction step 
			 *
			 * Add all debug informations to the verbosity flag of the ()-operator below, to get these
			 * details, e.g. to get all informations use '127' and to get only statistics at the end
			 * use '15'
			 */
			int verbosity;
			/**
			 * The output stream for the verbose messages. 
			 */
			std::ostream* out;

			/**
			 * Stores the time of all reduction steps
			 */
			double reductionTime;
			/**
			 * Stores the time of all matrix preparation steps
			 */
			double prepareTime;
			/**
			 * Stores the time of all critcial pairs update steps
			 */
			double updateTime;
			/**
			 * Stores the time of all simplify update steps
			 */
			double simplifyTime;
			/**
			 * Stores the time which is overhead by using MPI.
			 */
			double mpiTime;

			F4Logger(int verbosity, std::ostream* out) : verbosity(verbosity), out(out) {
				reductionTime = 0;
				prepareTime = 0;
				updateTime = 0;
				simplifyTime = 0;
			}

			static double seconds()
			{    
				struct timeval tv;
				gettimeofday(&tv, NULL);
				return (double)tv.tv_sec + ((double)tv.tv_usec)/1000000.0;
			}



	};
}
#endif
