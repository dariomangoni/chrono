#ifndef CHSOLVERMR_H
#define CHSOLVERMR_H
#include "solver/ChSolver.h"
#include "ChModelReduction.h"
#include "core/ChTimer.h"

namespace chrono
{
    class ChSolverMR: public ChSolver
    {
    public:
        ChSolverMR(ChSystem& system) : m_system(system) {}
        ~ChSolverMR() override {}

        bool SolveRequiresMatrix() const override { return false; }

        bool Setup(ChSystemDescriptor& sysd) override {

            if (m_setup_call == 0)
            {
                m_timer_setup_assembly.start();

                // Calculate problem size at first call.
                if (m_setup_call == 0) {
                    m_dim = sysd.CountActiveVariables() + sysd.CountActiveConstraints();
                    eig_vector_col.Reset(m_dim, 1);
                }

                // Let the matrix acquire the information about ChSystem
                if (m_force_sparsity_pattern_update)
                {
                    m_force_sparsity_pattern_update = false;

                    ChSparsityPatternLearner sparsity_learnerA(m_dim, m_dim, true);
                    m_system.GetStiffnessMatrix(&sparsity_learnerA);
                    matK.LoadSparsityPattern(sparsity_learnerA);

                    ChSparsityPatternLearner sparsity_learnerB(m_dim, m_dim, true);
                    m_system.GetMassMatrix(&sparsity_learnerB);
                    matM.LoadSparsityPattern(sparsity_learnerB);
                }
                else
                {
                    // If an NNZ value for the underlying matrix was specified, perform an initial resizing, *before*
                    // a call to ChSystemDescriptor::ConvertToMatrixForm(), to allow for possible size optimizations.
                    // Otherwise, do this only at the first call, using the default sparsity fill-in.
                    if (m_nnz != 0) {
                        matK.Reset(m_dim, m_dim, m_nnz);
                        matM.Reset(m_dim, m_dim, m_nnz);
                    }
                    else
                        if (m_setup_call == 0) {
                            matK.Reset(m_dim, m_dim, static_cast<int>(m_dim * (m_dim * SPM_DEF_FULLNESS)));
                            matM.Reset(m_dim, m_dim, static_cast<int>(m_dim * (m_dim * SPM_DEF_FULLNESS)));
                        }
                }

                //sysd.ConvertToMatrixForm(&matK, &matM, true); //TODO: delete this functio in ChSystemDescriptor

                m_system.GetStiffnessMatrix(&matK);
                m_system.GetMassMatrix(&matM);

                // Allow the matrix to be compressed. TODO: it should be done inside the eigen solver...
                matK.Compress();
                matM.Compress();

                m_timer_setup_assembly.stop();

                //matK.VerifyMatrix();
                //matM.VerifyMatrix();

                //matK.ExportToDatFile("C:/K", 6);
                //matM.ExportToDatFile("C:/M", 6);

                // Perform the factorization with the Pardiso sparse direct solver.
                m_timer_setup_solvercall.start();
                chrono::ChSymGEigsSolver eig_solver(matK, matM, eig_val, eig_vect);
                eig_solver.compute(num_computed_nodes == 0 ? m_dim - 1 : std::min(num_computed_nodes, m_dim - 1));
                m_timer_setup_solvercall.stop();
            }
            

            m_setup_call++;

            return true;
        }

        double Solve(ChSystemDescriptor& sysd) override { // Assemble the problem right-hand side vector.

            m_solve_call++;


            //// Scatter solution vector to the system descriptor. TODO: if the selected_mode doesn't change we can avoid to reload the vector every time
            //m_timer_solve_assembly.start();
            //eig_vector_col.PasteClippedMatrix(&eig_vect, 0, selected_mode-1, eig_vect.GetRows(), 1, 0, 0);
            //sysd.FromVectorToUnknowns(eig_vector_col);
            //m_timer_solve_assembly.stop();

            // Scatter solution vector to the system descriptor. TODO: if the selected_mode doesn't change we can avoid to reload the vector every time
            m_timer_solve_assembly.start();
            eig_vector_col.PasteClippedMatrix(&eig_vect, 0, selected_mode-1, eig_vect.GetRows(), 1, 0, 0);
            //sysd.FromVectorToUnknowns(eig_vector_col);
            m_timer_solve_assembly.stop();

            return 0.0f;
        }


        void SetVisualizedMode(int mode_to_visualize)
        {
            selected_mode = mode_to_visualize;
        }

        void SetComputedModes(int modes_to_be_extracted)
        {
            num_computed_nodes = modes_to_be_extracted;
        }

        void SetSystem(ChSystem& system) { m_system = system; }

        /// Get a handle to the underlying matrix.
        ChCSR3Matrix& GetMatrixA() { return matK; }
        ChCSR3Matrix& GetMatrixB() { return matM; }

        /// Enable/disable locking the sparsity pattern (default: false).
        /// If \a val is set to true, then the sparsity pattern of the problem matrix is assumed
        /// to be unchanged from call to call.
        void SetSparsityPatternLock(bool val)
        {
            m_lock = val;
            matK.SetSparsityPatternLock(m_lock);
            matM.SetSparsityPatternLock(m_lock);
        }

        /// Call an update of the sparsiy pattern on the underlying matrix.
        /// It is used to inform the solver (and the underlying matrices) that the sparsity pattern is changed.
        /// It is suggested to call this function just after the construction of the solver.
        void ForceSparsityPatternUpdate(bool val = true)
        {
            m_force_sparsity_pattern_update = val;
        }

        /// Reset timers for internal phases in Solve and Setup.
        void ResetTimers() {
            m_timer_setup_assembly.reset();
            m_timer_setup_solvercall.reset();
            m_timer_solve_assembly.reset();
            m_timer_solve_solvercall.reset();
        }

        /// Get cumulative time for assembly operations in Solve phase.
        double GetTimeSolve_Assembly() const { return m_timer_solve_assembly(); }
        /// Get cumulative time for Pardiso calls in Solve phase.
        double GetTimeSolve_SolverCall() const { return m_timer_solve_solvercall(); }
        /// Get cumulative time for assembly operations in Setup phase.
        double GetTimeSetup_Assembly() const { return m_timer_setup_assembly(); }
        /// Get cumulative time for Pardiso calls in Setup phase.
        double GetTimeSetup_SolverCall() const { return m_timer_setup_solvercall(); }


        /// Method to allow serialization of transient data to archives.
        virtual void ArchiveOUT(ChArchiveOut& marchive) override {
            // version number
            marchive.VersionWrite(1);
            // serialize parent class
            ChSolver::ArchiveOUT(marchive);
            // serialize all member data:
            marchive << CHNVP(m_lock);
        }

        /// Method to allow de serialization of transient data from archives.
        virtual void ArchiveIN(ChArchiveIn& marchive) override {
            // version number
            int version = marchive.VersionRead();
            // deserialize parent class
            ChSolver::ArchiveIN(marchive);
            // stream in all member data:
            marchive >> CHNVP(m_lock);
        }

        

    private:
        ChCSR3Matrix matK, matM;
        ChMatrixDynamic<double> eig_val, eig_vect, eig_vector_col;
        ChSymGEigsSolver eig_solver = { matK, matM, eig_val, eig_vect };
        ChSystem& m_system;
        

        int m_dim = 0;                      ///< problem size
        int m_nnz = 0;                      ///< user-supplied estimate of NNZ
        int m_solve_call = 0;               ///< counter for calls to Solve
        int m_setup_call = 0;               ///< counter for calls to Setup
        int num_computed_nodes = 0;         ///< number of modes to extract (0: all)
        int selected_mode = 1;              ///< mode visualized
        double modes_magnification = 1;     ///< amplify mode amplitude

        bool m_lock = false;              ///< is the matrix sparsity pattern locked?
        bool m_force_sparsity_pattern_update = false; ///< is the sparsity pattern changed compared to last call?

        ChTimer<> m_timer_setup_assembly;  ///< timer for matrix assembly
        ChTimer<> m_timer_setup_solvercall;   ///< timer for factorization
        ChTimer<> m_timer_solve_assembly;  ///< timer for RHS assembly
        ChTimer<> m_timer_solve_solvercall;   ///< timer for solution

    };
}



#endif

