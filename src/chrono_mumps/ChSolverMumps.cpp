#include "ChSolverMumps.h"


namespace chrono
{
    bool ChSolverMumps::Setup(ChSystemDescriptor& sysd)
    {
        m_timer_setup_assembly.start();

        sysd.ConvertToMatrixForm(&m_mat, nullptr);
        m_dim = m_mat.GetNumRows();

        // Allow the matrix to be compressed.
        bool change = m_mat.Compress();

        // Set current matrix in the MKL engine.
        m_engine.SetMatrix(m_mat);

        m_timer_setup_assembly.stop();

        // Perform the factorization with the Pardiso sparse direct solver.
        m_timer_setup_solvercall.start();
        auto mumps_message = m_engine.MumpsCall(ChMumpsEngine::mumps_JOB::ANALYZE_FACTORIZE);
        m_timer_setup_solvercall.stop();

        m_setup_call++;

        if (verbose) {
            GetLog() << " Mumps setup n = " << m_dim << "  nnz = " << m_mat.GetNNZ() << "\n";
            GetLog() << "  assembly: " << m_timer_setup_assembly.GetTimeSecondsIntermediate() << "s" <<
                        "  solver_call: " << m_timer_setup_solvercall.GetTimeSecondsIntermediate() << "\n";
        }

        if (mumps_message != 0) {
            m_engine.PrintINFOG();
            return false;
        }

        return true;
    }

    double ChSolverMumps::Solve(ChSystemDescriptor& sysd) ///< system description with constraints and variables
	{
        m_timer_solve_assembly.start();
		sysd.ConvertToMatrixForm(nullptr, &m_rhs);
        m_timer_solve_assembly.stop();

        m_timer_solve_solvercall.start();
		m_engine.SetRhsVector(m_rhs);
		m_engine.MumpsCall(ChMumpsEngine::mumps_JOB::SOLVE);
        m_timer_solve_solvercall.stop();

        m_solve_call++;
		printf("\nCall: %d\n", m_solve_call);

		sysd.FromVectorToUnknowns(m_rhs);

		return 0.0;
	}

} // namespace chrono