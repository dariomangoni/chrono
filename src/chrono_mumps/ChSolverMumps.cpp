#include "ChSolverMumps.h"


namespace chrono
{
    //bool ChSolverMumps::Setup(ChSystemDescriptor& sysd)
    //{
    //    m_timer_setup_assembly.start();

    //    // Calculate problem size at first call.
    //    if (m_setup_call == 0) {
    //        m_dim = sysd.CountActiveVariables() + sysd.CountActiveConstraints();
    //    }


    //    // If an NNZ value for the underlying matrix was specified, perform an initial resizing, *before*
    //    // a call to ChSystemDescriptor::ConvertToMatrixForm(), to allow for possible size optimizations.
    //    // Otherwise, do this only at the first call, using the default sparsity fill-in.
    //    if (m_nnz != 0) {
    //        m_mat.Reset(m_dim, m_dim, m_nnz);
    //    }
    //    else
    //        if (m_setup_call == 0) {
    //            m_mat.Reset(m_dim, m_dim, static_cast<int>(m_dim * (m_dim * SPM_DEF_FULLNESS)));
    //        }

    //    sysd.ConvertToMatrixForm(&m_mat, nullptr);



    //    // Allow the matrix to be compressed.
    //    bool change = m_mat.Compress();


    //    // Set current matrix in the MKL engine.
    //    m_engine.SetMatrix(m_mat);


    //    m_timer_setup_assembly.stop();

    //    // Perform the factorization with the Pardiso sparse direct solver.
    //    m_timer_setup_solvercall.start();
    //    int pardiso_message_phase12 = m_engine.MumpsCall(1);
    //    m_timer_setup_solvercall.stop();

    //    m_setup_call++;

    //    if (verbose) {
    //        GetLog() << " MKL setup n = " << m_dim << "  nnz = " << m_mat.GetNNZ() << "\n";
    //        GetLog() << "  assembly: " << m_timer_setup_assembly.GetTimeSecondsIntermediate() << "s" <<
    //            "  solver_call: " << m_timer_setup_solvercall.GetTimeSecondsIntermediate() << "\n";
    //    }

    //    if (pardiso_message_phase12 != 0) {
    //        GetLog() << "Pardiso analyze+reorder+factorize error code = " << pardiso_message_phase12 << "\n";
    //        return false;
    //    }

    //    return true;
    //}

    double ChSolverMumps::Solve(ChSystemDescriptor& sysd) ///< system description with constraints and variables
	{
        m_timer_solve_assembly.start();
		sysd.ConvertToMatrixForm(&m_mat, &m_rhs);
        m_timer_solve_assembly.stop();

		m_engine.Initialize();
        m_timer_solve_solvercall.start();
        m_mat.Compress();
		m_engine.SetProblem(m_mat, m_rhs);
		m_engine.MumpsCall();
        m_timer_solve_solvercall.stop();

        m_solve_call++;
		printf("\nCall: %d\n", m_solve_call);

		sysd.FromVectorToUnknowns(m_rhs);

		return 0.0;
	}

} // namespace chrono