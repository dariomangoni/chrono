//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2010 Alessandro Tasora
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//


/// Class for interfacing with Mumps

#ifndef CHMUMPSENGINE_H
#define CHMUMPSENGINE_H

///////////////////////////////////////////////////
//
//   ChMumpsEngine.h
//
//   Use this header if you want to exploit
//	 MUMPS Library from Chrono::Engine programs.
//
//   HEADER file for CHRONO,
//   Multibody dynamics engine
//
// ------------------------------------------------
//             www.deltaknowledge.com
// ------------------------------------------------
///////////////////////////////////////////////////

#include "chrono_mumps/ChApiMumps.h"
#include "chrono_mumps/ChCOOMatrix.h"

#include <dmumps_c.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define USE_COMM_WORLD -987654

/* macro s.t. indices match documentation */
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]
#define INFO(I) info[(I)-1]
#define INFOG(I) infog[(I)-1]
#define RINFO(I) rinfo[(I)-1]
#define RINFOG(I) rinfog[(I)-1]




namespace chrono
{

	class ChApiMumps ChMumpsEngine
	{
	public:
        enum mumps_SYM
        {
            UNSYMMETRIC = 0,
            SYMMETRIC_POSDEF = 1,
            SYMMETRIC_GENERAL = 2
        };

        enum mumps_JOB
        {
            INIT = -1,
            END = -2,
            ANALYZE = 1,
            FACTORIZE = 2,
            SOLVE = 3,
            ANALYZE_FACTORIZE = 4,
            FACTORIZE_SOLVE = 5,
            COMPLETE = 6
        };

	    explicit ChMumpsEngine(mumps_SYM symmetry = UNSYMMETRIC, int mumps_mpi_comm = -987654, int activate_this_node = 1);
		virtual ~ChMumpsEngine();

		void SetProblem(const ChCOOMatrix& Z, const ChMatrix<>& rhs);
		void SetMatrix(const ChCOOMatrix& Z);

        /// Set the right-hand side vector.
        /// Note that it is the caller's responsibility to ensure that the size is appropriate.
        void SetRhsVector(const ChMatrix<>& b);
        void SetRhsVector(double* b);

        int MumpsCall(mumps_JOB job_call);
 
		void PrintINFOG();

        void SetNullPivotDetection(bool val, double threshold = 0)
        {
            mumps_id.ICNTL(24) = val; // activates null pivot detection
            mumps_id.ICNTL(25) = 0; // tries to compute one of the many solutions of AX = B
            mumps_id.CNTL(5) = 1e20; // fixation value
            mumps_id.CNTL(3) = threshold; // pivot threshold
        }

        double GetCNTL(int parnum) { return mumps_id.CNTL(parnum); }
        void SetCNTL(int parnum, int parvalue) { mumps_id.CNTL(parnum) = parvalue; }

		int GetICNTL(int parnum){ return mumps_id.ICNTL(parnum); }
		void SetICNTL(int parnum, int parvalue){ mumps_id.ICNTL(parnum) = parvalue; }

		int GetINFO(int parnum){ return mumps_id.INFO(parnum); }
		int GetINFOG(int parnum){ return mumps_id.INFOG(parnum); }

		double GetRINFO(int parnum){ return mumps_id.RINFO(parnum); }
		double GetRINFOG(int parnum){ return mumps_id.RINFOG(parnum); }

        DMUMPS_STRUC_C& GetMumpsStruc() { return mumps_id; }

        /// Classes that wraps some useful functions in OpenMP
        /// (in case no OpenMP is used, it defaults to dummy functions
        /// that do nothing)
#ifdef _OPENMP

        /// Sets the number of threads in subsequent parallel
        /// regions, unless overridden by a 'num_threads' clause
        static void SetNumThreads(int mth) { omp_set_num_threads(mth); }

        /// Returns the number of threads in the parallel region.
        static int GetNumThreads() { return omp_get_num_threads(); }

        /// Returns the thread number of the thread executing
        /// within its thread team.
        static int GetThreadNum() { return omp_get_thread_num(); }

        /// Returns the number of available processors on this machine
        static int GetNumProcs() { return omp_get_num_procs(); }

        /// Returns the max.number of threads that would be used
        /// by default if num_threads not specified. This is the same
        /// number as GetNumProcs() on most OMP implementations.
        static int GetMaxThreads() { return omp_get_max_threads(); }

#else
        static void SetNumThreads(int mth) {}
        static int GetNumThreads() { return 1; }
        static int GetThreadNum() { return 0; }
        static int GetNumProcs() { return 1; }
        static int GetMaxThreads() { return 1; }
#endif

	private:
		DMUMPS_STRUC_C mumps_id;
		int myid = 0;
		int ierr = 0;
	};
} // end of namespace chrono


#endif