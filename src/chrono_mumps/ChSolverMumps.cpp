#include "ChSolverMumps.h"
#include "chrono/core/ChCSR3Matrix.h"


namespace chrono
{
    bool ChSolverMumps::Setup(ChSystemDescriptor& sysd)
    {
        m_timer_setup_assembly.start();

        // Calculate problem size at first call.
        if (m_setup_call == 0) {
            m_dim = sysd.CountActiveVariables() + sysd.CountActiveConstraints();
        }

        // Let the matrix acquire the information about ChSystem
        if (m_force_sparsity_pattern_update)
        {
            m_force_sparsity_pattern_update = false;

            ChSparsityPatternLearner sparsity_learner(m_dim, m_dim, true);
            sysd.ConvertToMatrixForm(&sparsity_learner, nullptr);
            m_mat.LoadSparsityPattern(sparsity_learner);
        }
        else
        {
            // If an NNZ value for the underlying matrix was specified, perform an initial resizing, *before*
            // a call to ChSystemDescriptor::ConvertToMatrixForm(), to allow for possible size optimizations.
            // Otherwise, do this only at the first call, using the default sparsity fill-in.
            if (m_nnz != 0) {
                m_mat.Reset(m_dim, m_dim, m_nnz);
            }
            else
                if (m_setup_call == 0) {
                    m_mat.Reset(m_dim, m_dim, static_cast<int>(m_dim * (m_dim * SPM_DEF_FULLNESS)));
                }
        }

        sysd.ConvertToMatrixForm(&m_mat, nullptr);

        m_dim = m_mat.GetNumRows();

        // Allow the matrix to be compressed.
        m_mat.Compress();

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