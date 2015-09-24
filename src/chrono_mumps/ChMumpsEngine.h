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
//   Use this header if you want to exploit Intel®
//	 MUMPS Library from Chrono::Engine programs.
//
//   HEADER file for CHRONO,
//  Multibody dynamics engine
//
// ------------------------------------------------
//             www.deltaknowledge.com
// ------------------------------------------------
///////////////////////////////////////////////////

#include "chrono_mumps/ChApiMumps.h"
#include "chrono_mumps/ChCOOMatrix.h"

#include <dmumps_c.h>
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

/* macro s.t. indices match documentation */
#define ICNTL(I) icntl[(I)-1]
#define INFO(I) info[(I)-1]
#define INFOG(I) infog[(I)-1]
#define RINFO(I) rinfo[(I)-1]
#define RINFOG(I) rinfog[(I)-1]




namespace chrono
{

	class ChApiMumps ChMumpsEngine
	{
	public:
		ChMumpsEngine(){};
		virtual ~ChMumpsEngine();

		void SetProblem(const ChCOOMatrix& Z, const ChMatrix<>& rhs);
		void SetMatrix(const ChCOOMatrix& Z);

        /// Set the right-hand side vector.
        /// Note that it is the caller's responsibility to ensure that the size is appropriate.
        void SetRhsVector(const ChMatrix<>& b);
        void SetRhsVector(double* b);


		void Initialize();
		int MumpsCall(int job_call = 6);
 
		void PrintINFOG();

		int GetICNTL(int parnum){ return mumps_id.ICNTL(parnum); }
		void SetICNTL(int parnum, int parvalue){ mumps_id.ICNTL(parnum) = parvalue; }

		int GetINFO(int parnum){ return mumps_id.INFO(parnum); }
		int GetINFOG(int parnum){ return mumps_id.INFOG(parnum); }
		double GetRINFO(int parnum){ return mumps_id.RINFO(parnum); }
		double GetRINFOG(int parnum){ return mumps_id.RINFOG(parnum); }

	private:
		DMUMPS_STRUC_C mumps_id;
		int n = 0;
		int nz = 0;
		int myid = 0;
		int ierr = 0;
	};
} // end of namespace chrono


#endif