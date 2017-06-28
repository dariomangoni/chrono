// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Dario Mangoni
// =============================================================================
//
// Interior-Point Header File
//
// =============================================================================

#ifndef CHIPSOLVER_H
#define CHIPSOLVER_H

#include "ChApiInteriorPoint.h"
#include "chrono/solver/ChSystemDescriptor.h"
#include "chrono/solver/ChSolver.h"


#ifdef CHRONO_MUMPS
#include "chrono_mumps/ChCOOMatrix.h"
#include "chrono_mumps/ChMumpsEngine.h"
#endif

// Interior point methdon based on Numerical Optimization by Nocedal, Wright
// minimize 0.5*vT*H*v + vT*v while Ax>=b (16.54 pag.480)
// WARNING: FOR THE MOMENT THE CONSTRAINTS MUST BE INEQUALITIES
// Further references: (all pages number refers to [1] if not otherwise specified)
// [1] Nocedal&Wright - Numerical Optimization, 2nd edition
// [2] D'Apuzzo et al. - Starting-point strategies for an infeasible potential reduction method
// [3] Mangoni D., Tasora A. - Solving Unilateral Contact Problems in Multibody Dynamics using a Primal-Dual Interior Point Method
// [3] Mangoni D., Tasora A., R. Garziera - A Primal-Dual Predictor-Corrector Interior Point Method for Non-Smooth Contact Dynamics
// [4] Meszaros - Steplengths in interior-point algorithms of quadratic programming
// [5] Vanderberghe - The CVXOPT linear and quadratic cone solvers

// Symbol conversion table from [1] to [2]
// [2]    | [1]    |
//  z     |  y     |
//  y     | lambda |
//  Q     |  H     |
// lambda |  -     |
//  b     |  b     |
//  s     |  -     |
//  u     |  -     |
//  v     |  -     |
//  d     |  -     |
//  t     |  -     |
//  H     |  -     |

// Symbol conversion table from [1] to Chrono
//                          | Chr | [1]
// Mass/Stiffness           |  H  |  H
// Acceleration/DSpeed      |  q  |  v
// Constraints              | Cq  |  EqCon + IneqCon
// Forces(internal)         |  l  |  lambda
// Forces(external)         |  f  |  -c
// Constr. compliance       |  ?  |  E (+ o - ?)
// Slack (contact distance) |  c  |  y
// Constraint rhs           |  b  |  -b

// KKT conditions [4]
// H*v - EqCon^T*gamma - IneqCon^T*lambda + c = 0; (dual)
// EqCon*v - b_eq = 0; (primal)
// IneqCon*v - y - b_ineq = 0; (primal)
// y.*lambda = 0
// y>=0
// lambda>=0

// In order to run, the ChInteriorPoint solver needs:
// - ChCOOMatrix BigMat (so it needs H, IneqCon, EqCon)
// - IPrhs_t rhs
// everything else is computed from these two elements.


namespace chrono {

/** \class ChInteriorPoint
\brief ChInteriorPoint is a class that implements an Interior-Point Primal-Dual solver for QP convex programming.

The solver is experimental;
- no friction yet
*/

class ChApiInteriorPoint ChInteriorPoint : public ChSolver {

  public:
    enum class IP_KKT_SOLUTION_METHOD { STANDARD, AUGMENTED, NORMAL };
    enum class IP_STARTING_POINT_METHOD { STP1, STP2, NOCEDAL, NOCEDAL_WS, VANDERBERGHE };

  private:
	int n = 0;  ///< size of #v, #H, #IneqCon columns, #EqCon columns
    int m_eq = 0;  ///< size of #gamma, #EqCon rows
    int m_ineq = 0;  ///< size of #lambda, #y, #IneqCon rows
    int solver_call = 0;
    int iteration_count = 0;
    int iteration_count_max = 50;

    const bool EQUAL_STEP_LENGTH = true;
    const bool ADAPTIVE_ETA = false;
    const bool ONLY_PREDICT = false;
    bool warm_start_broken = false;
    bool warm_start = true;
	bool leverage_symmetry = false;

    ChTimer<> ip_timer_solver_solvercall;
    ChTimer<> ip_timer_solve_assembly;
    int ip_solver_call = 0;
    int iteration_count_tot = 0;

    IP_KKT_SOLUTION_METHOD KKT_solve_method = IP_KKT_SOLUTION_METHOD::AUGMENTED;
    IP_STARTING_POINT_METHOD starting_point_method = IP_STARTING_POINT_METHOD::NOCEDAL;

    // Problem matrices and vectors
    ChCOOMatrix BigMat;	///< Global sparse matrix
    ChCOOMatrix E;  ///< Compliance matrix

    // Known terms
    struct IPrhs_t {
		IPrhs_t(){}
		IPrhs_t(int n, int m_eq, int m_ineq) { Resize(n, m_eq, m_ineq); }
		void Resize(int n, int m_eq, int m_ineq)
		{
			c.Resize(n, 1);
			b.Resize(m_eq + m_ineq, 1);
			b_eq.Resize(m_eq, 1);
			b_ineq.Resize(m_ineq, 1);
		}
		ChMatrixDynamic<double> c;  ///< forces (is '-f' in chrono)
        ChMatrixDynamic<double> b;  ///< rhs of constraints (is '-b' in chrono)
		ChMatrixDynamic<double> b_eq;  ///< rhs of constraints (is '-b' in chrono)
		ChMatrixDynamic<double> b_ineq;  ///< rhs of constraints (is '-b' in chrono)
    } rhs;



    // Variables
    struct IPvariables_t {
		IPvariables_t(){}
		IPvariables_t(int n, int m_eq, int m_ineq) { Resize(n, m_eq, m_ineq); }
		void Resize(int n, int m_eq, int m_ineq)
		{
			v.Resize(n, 1);
			y.Resize(m_ineq, 1);
			gamma.Resize(m_eq, 1);
			lambda.Resize(m_ineq, 1);
		}
        ChMatrixDynamic<double> v;    ///< DeltaSpeed/Acceleration ('q' in chrono)
        ChMatrixDynamic<double> y;    ///< Slack variable/Contact points distance ('c' in chrono)
        ChMatrixDynamic<double> gamma;    ///< Lagrangian multipliers/contact|constraint forces for bilateral constraints ('l' in chrono)
        ChMatrixDynamic<double> lambda;  ///< Lagrangian multipliers/contact|constraint forces for unilateral constraints ('l' in chrono)
    } var;

    // Residuals
    struct IPresidual_nnorm_t {
        double rp_gamma_nnorm = 1e-10;
        double rp_lambda_nnorm = 1e-10;
        double rd_nnorm = 1e-10;
        double mu = 1e-9;
    } res_nnorm_tol;


    struct IPresidual_t {
		ChMatrixDynamic<double> rd;    ///< Residual about dual variables (i.e. violation of constraints equations); rd = H*v - AT*lambda + c.
		ChMatrixDynamic<double> rp_gamma;    ///< Residual about primal variables regarding unilateral constraints (i.e. violation if dynamic equation of motion); rp_lambda = IneqCon*v - y - b.
        ChMatrixDynamic<double> rp_lambda;   ///< Residual about primal variables regarding bilateral constraints (i.e. violation if dynamic equation of motion); rp_lambda = IneqCon*v - y - b.
        double mu = 0;  ///< complementarity measure

		IPresidual_t(){}
		IPresidual_t(int n, int m_eq, int m_ineq) { Resize(n, m_eq, m_ineq); }
		void Resize(int n, int m_eq, int m_ineq)
		{
			rd.Resize(n, 1);
			rp_gamma.Resize(m_eq, 1);
			rp_lambda.Resize(m_ineq, 1);
		}

        bool operator<(const IPresidual_t& other) const { return rp_lambda.NormTwo() < other.rp_lambda.NormTwo() && rp_gamma.NormTwo() < other.rp_gamma.NormTwo() && rd.NormTwo() < other.rd.NormTwo() && mu < other.mu; }
        bool operator<=(const IPresidual_t& other) const { return rp_lambda.NormTwo() <= other.rp_lambda.NormTwo() && rp_gamma.NormTwo() <= other.rp_gamma.NormTwo() && rd.NormTwo() <= other.rd.NormTwo() && mu <= other.mu; }
        bool operator>(const IPresidual_t& other) const { return !(*this <= other); }
        bool operator>=(const IPresidual_t& other) const { return !(*this < other); }

        bool operator<(const IPresidual_nnorm_t& other) const { return rp_lambda.NormTwo() < other.rp_lambda_nnorm * rp_lambda.GetRows() && rp_gamma.NormTwo() < other.rp_gamma_nnorm * rp_gamma.GetRows() && rd.NormTwo() < other.rd_nnorm * rd.GetRows() && mu < other.mu; }
        bool operator<=(const IPresidual_nnorm_t& other) const { return rp_lambda.NormTwo() <= other.rp_lambda_nnorm * rp_lambda.GetRows() && rp_gamma.NormTwo() <= other.rp_gamma_nnorm * rp_gamma.GetRows() && rd.NormTwo() <= other.rd_nnorm * rd.GetRows() && mu <= other.mu; }
        bool operator>(const IPresidual_nnorm_t& other) const { return !(*this <= other); }
        bool operator>=(const IPresidual_nnorm_t& other) const { return !(*this < other); }

    } res;

    // Temporaries used in different functions
    mutable ChMatrixDynamic<double> vectn;  // temporary variable that has always size (#n,1)
    mutable ChMatrixDynamic<double> vectm_ineq;  // temporary variable that has always size (#m_eq,1)
    mutable ChMatrixDynamic<double> vectm_eq;  // temporary variable that has always size (#m_ineq,1)
    mutable ChMatrixDynamic<double> sol_chrono;  // intermediate file to inject the IP solution into Chrono used in adapt_to_Chrono()


    // MUMPS engine
    ChMumpsEngine mumps_engine;

    // IP specific functions
    void iterate();  ///< Perform an IP iteration; returns \e true if exit conditions are met.
    void setup_system_matrix(const IPvariables_t& vars);
    void makeNewtonStep(IPvariables_t& Dvar_unknown, ChMatrix<>& rhs, const IPresidual_t& residuals);
    void set_starting_point(IP_STARTING_POINT_METHOD start_point_method, int n_old = 0, int m_eq_old = 0, int m_ineq_old = 0);
    static double find_Newton_step_length(const ChMatrix<double>& vect, const ChMatrix<double>& Dvect, double tau = 1);
    void find_Newton_step_length(const IPvariables_t& vars, const IPvariables_t& Dvars, double tau, double& alfa_prim, double& alfa_dual) const;
    double evaluate_objective_function() const;  ///< Evaluate the objective function i.e. 0.5*vT*H*v + vT*v.

    // Auxiliary
    void reset_internal_dimensions(int n_old, int m_eq_new, int m_ineq_new, int m_ineq_full);
    ChMatrix<>& adapt_to_Chrono(ChSystemDescriptor& sysd, ChMatrix<>& solution_vect) const;
    void residual_fullupdate(IPresidual_t& residuals, const IPvariables_t& variables) const;
    void make_positive_definite();  ///< Change EQ|IneqCon^T to -Eq|IneqCon^T in the current system matrix.
    void multiplyIneqCon(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
	void multiplyEqCon(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
	void multiplyNegIneqConT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
	void multiplyNegEqConT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
	void multiplyH(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;

    // Debug
    std::ofstream logfile_stream;
    std::string logfile_name{ "interior_point_log" };
    bool print_history = false;
    void LoadProblem();




  public:
    ChInteriorPoint();
    ~ChInteriorPoint();
    double Solve(ChSystemDescriptor& sysd) override;

    bool SolveRequiresMatrix() const override { return true; }

    // Auxiliary
    /// Set the Karush–Kuhn–Tucker problem form that will be used to solve the IP problem. Change it before starting the solver.
    void SetKKTSolutionMethod(IP_KKT_SOLUTION_METHOD qp_solve_type_selection) { KKT_solve_method = qp_solve_type_selection; }

    /// Set the Karush–Kuhn–Tucker problem form that will be used to solve the IP problem. Change it before starting the solver.
    void SetStartingPointMethod(IP_STARTING_POINT_METHOD starting_point_method_in) { starting_point_method = starting_point_method_in; }

    /// Set the maximum number of iterations after which the iteration loop will be stopped.
    void SetMaxIterations(int max_iter) { iteration_count_max = max_iter; }

    /// Set the tolerance over the residual of the \e primal unilateral variables (i.e. violation of constraints equations).
    void SetPrimalUnilateralResidualTolerance(double rp_tol) { res_nnorm_tol.rp_lambda_nnorm = rp_tol; }
    /// Set the tolerance over the residual of the \e primal bilateral variables (i.e. violation of constraints equations).
    void SetPrimalBilateralResidualTolerance(double rp_tol) { res_nnorm_tol.rp_lambda_nnorm = rp_tol; }
    /// Set the tolerance over the residual of the \e dual variables (i.e. stationarity of the solution).
    void SetDualResidualTolerance(double rd_tol) { res_nnorm_tol.rd_nnorm = rd_tol; }
    /// Set the tolerance over the residual of complementarity measure (i.e. violation of orthogonality of forces and contact points distance)
    void SetComplementarityMeasureTolerance(double complementarity_tol) { res_nnorm_tol.mu = complementarity_tol; }

    /// Set the null pivot detection for the kernel solver.
    void SetNullPivotDetection(bool on_off, double threshold) { mumps_engine.SetNullPivotDetection(on_off, threshold); }

    /// Leverage matrix symmetry
	void SetUseSymmetry(bool val);

	/// Get cumulative time for assembly operations in Solve phase.
	double GetTimeSolve_Assembly() const { return ip_timer_solve_assembly(); }
	/// Get cumulative time for Pardiso calls in Solve phase.
	double GetTimeSolve_SolverCall() const { return ip_timer_solver_solvercall(); }

    // Test
    void DumpProblem(std::string suffix = "");
    void DumpIPStatus(std::string suffix = "") const;
    void RecordHistory(bool on_off, std::string file_name = "interior_point_log");
	void PrintIPStatus() const;
    int GetSolverCalls() const { return solver_call; }
    int GetIPSolverCalls() const { return ip_solver_call; } ///< it may differ from GetSolverCalls() in case no inequality constraint is found
    int GetIPIterations() const { return iteration_count_tot; } ///< total number of iterations of the IP
    int GetMassMatrixDimension() const { return n; }
    int GetUnilateralConstraintsMatrixRows() const { return m_ineq; }
    int GetBilateralConstraintsMatrixRows() const { return m_eq; }
    void Solve(const ChCOOMatrix& normal_mat, const ChMatrix<double>& rhs_b, const ChMatrix<double>& rhs_c, ChMatrixDynamic<double>& var_x, ChMatrixDynamic<double>& var_y, ChMatrixDynamic<double>& var_gamma, ChMatrixDynamic<double>& var_lambda);
};

}  // end of namespace chrono

#endif