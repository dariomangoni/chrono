// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
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
// minimize 0.5*xT*G*x + xT*x while Ax>=b (16.54 pag.480)
// WARNING: FOR THE MOMENT THE CONSTRAINTS MUST BE INEQUALITIES
// Further references: (all pages number refers to [1] if not otherwise specified)
// [1] Nocedal&Wright - Numerical Optimization, 2nd edition
// [2] D'Apuzzo et al. - Starting-point strategies for an infeasible potential reduction method
// [3] Mangoni D., Tasora A. - Solving Unilateral Contact Problems in Multibody Dynamics using a Primal-Dual Interior Point Method
// [4] Meszaros - Steplengths in interior-point algorithms of quadratic programming

// Symbol conversion table from [1] to [2]
// [2] | [1]
//  z  |  y
//  y  | lam
//  Q  |  G
// lam |  -
//  b  |  b
//  s  |  -
//  u  |  -
//  v  |  -
//  d  |  -
//  t  |  -
//  G  |  -

// Symbol conversion table from [1] to Chrono
//                          | Chr | [1]
// Mass/Stiffness           |  H  |  G
// Acceleration/DSpeed      |  q  |  x
// Constraints              | Cq  |  A
// Forces(internal)         |  l  |  lam
// Forces(external)         |  f  |  -c
// Constr. compliance       |  ?  |  E (+ o - ?)
// Slack (contact distance) |  c  |  y
// Constraint rhs           |  b  |  -b

// KKT conditions (16.55 pag.481)
// G*x-AT*lam+c = 0; (dual)
// A*x-y-b = 0; (primal)
// y.*lam = 0 (mixed)
// y>=0
// lam>=0

// In order to run the ChInteriorPoint solver needs:
// - ChCOOMatrix BigMat
// - IPrhs_t rhs
// everything else is computed from this two elements.


namespace chrono {

/** \class ChInteriorPoint
\brief ChInteriorPoint is a class that implements an Interior-Point Primal-Dual solver for QP convex programming.

The solver is experimental;
- no friction yet
*/

class ChApiInteriorPoint ChInteriorPoint : public ChSolver {

  public:
    enum class IP_KKT_SOLUTION_METHOD { STANDARD, AUGMENTED, NORMAL };
    enum class IP_STARTING_POINT_METHOD { STP1, STP2, NOCEDAL, NOCEDAL_WS };

  private:
    int m = 0;  // size of #lam, #y, A rows
    int n = 0;  // size of #x, G, A columns
    int solver_call = 0;
    int iteration_count = 0;
    int iteration_count_max = 50;

    const bool EQUAL_STEP_LENGTH = true;
    const bool ADAPTIVE_ETA = true;
    const bool ONLY_PREDICT = false;
    bool warm_start_broken = false;
    bool warm_start = true;

    ChTimer<> ip_timer;
    int ip_solver_call = 0;
    int iteration_count_tot = 0;

    IP_KKT_SOLUTION_METHOD KKT_solve_method = IP_KKT_SOLUTION_METHOD::AUGMENTED;
    IP_STARTING_POINT_METHOD starting_point_method = IP_STARTING_POINT_METHOD::NOCEDAL;

    // Problem matrices and vectors
    ChCOOMatrix BigMat;
    ChCOOMatrix SmallMat;
    ChCOOMatrix E;  // compliance matrix

    // Known terms
    struct IPrhs_t {
        ChMatrixDynamic<double> b;  ///< rhs of constraints (is '-b' in chrono)
        ChMatrixDynamic<double> c;  ///< forces (is '-f' in chrono)
    } rhs;



    // Variables
    struct IPvariables_t {
        ChMatrixDynamic<double> x;    ///< DeltaSpeed/Acceleration ('q' in chrono)
        ChMatrixDynamic<double> y;    ///< Slack variable/Contact points distance ('c' in chrono)
        ChMatrixDynamic<double> lam;  ///< Lagrangian multipliers/contact|constraint forces ('l' in chrono)
    } var;

    // Residuals
    struct IPresidual_nnorm_t {
        double rp_nnorm = 1e-10;
        double rd_nnorm = 1e-10;
        double mu = 1e-9;
    } res_nnorm_tol;


    struct IPresidual_t {
        ChMatrixDynamic<double> rp;    ///< Residual about primal variables (i.e. violation if dynamic equation of motion); rp = A*x - y - b.
        ChMatrixDynamic<double> rd;    ///< Residual about dual variables (i.e. violation of constraints equations); rd = G*x - AT*lam + c.
        ChMatrixDynamic<double> rpd;    ///< Residual about primal-dual variables (only for #IP_KKT_SOLUTION_METHOD#NORMAL mode)
        double mu = 0;  ///< complementarity measure

        bool operator<(const IPresidual_t& other) const { return rp.NormTwo() < other.rp.NormTwo() && rd.NormTwo() < other.rd.NormTwo() && mu < other.mu; }
        bool operator<=(const IPresidual_t& other) const { return rp.NormTwo() <= other.rp.NormTwo() && rd.NormTwo() <= other.rd.NormTwo() && mu <= other.mu; }
        bool operator>(const IPresidual_t& other) const { return !(*this <= other); }
        bool operator>=(const IPresidual_t& other) const { return !(*this < other); }

        bool operator<(const IPresidual_nnorm_t& other) const { return rp.NormTwo() < other.rp_nnorm * rp.GetRows() && rd.NormTwo() < other.rd_nnorm * rd.GetRows() && mu < other.mu; }
        bool operator<=(const IPresidual_nnorm_t& other) const { return rp.NormTwo() <= other.rp_nnorm * rp.GetRows() && rd.NormTwo() <= other.rd_nnorm * rd.GetRows() && mu <= other.mu; }
        bool operator>(const IPresidual_nnorm_t& other) const { return !(*this <= other); }
        bool operator>=(const IPresidual_nnorm_t& other) const { return !(*this < other); }

    } res;

    // Temporaries used in different functions
    mutable ChMatrixDynamic<double> vectn;  // temporary variable that has always size (#n,1)
    mutable ChMatrixDynamic<double> vectm;  // temporary variable that has always size (#m,1)
    mutable ChMatrixDynamic<double> sol_chrono;  // intermediate file to inject the IP solution into Chrono used in adapt_to_Chrono()


    // MUMPS engine
    ChMumpsEngine mumps_engine;

    // IP specific functions
    void iterate();  ///< Perform an IP iteration; returns \e true if exit conditions are met.
    void setup_system_matrix(const IPvariables_t& vars);
    void makeNewtonStep(IPvariables_t& Dvar_unknown, ChMatrix<>& rhs, const IPresidual_t& residuals);
    void set_starting_point(IP_STARTING_POINT_METHOD start_point_method, int n_old = 0, int m_old = 0);
    static double find_Newton_step_length(const ChMatrix<double>& vect, const ChMatrix<double>& Dvect, double tau = 1);
    void find_Newton_step_length(const IPvariables_t& vars, const IPvariables_t& Dvars, double tau, double& alfa_prim, double& alfa_dual) const;
    double evaluate_objective_function() const;  ///< Evaluate the objective function i.e. 0.5*xT*G*x + xT*x.

    // Auxiliary
    void reset_internal_dimensions(int n_old, int m_old);
    ChMatrix<>& adapt_to_Chrono(ChMatrix<>& solution_vect) const;
    void residual_fullupdate(IPresidual_t& residuals, const IPvariables_t& variables) const;
    void make_positive_definite();  ///< Change A^T to -A^T in the current system matrix.
    void multiplyA(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
    void multiplyNegAT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
    void multiplyG(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;

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

    /// Set the tolerance over the residual of the \a primal variables (i.e. violation of constraints equations).
    void SetPrimalResidualTolerance(double rp_tol) { res_nnorm_tol.rp_nnorm = rp_tol; }

    /// Set the tolerance over the residual of the \a dual variables (i.e. stationarity of the solution).
    void SetDualResidualTolerance(double rd_tol) { res_nnorm_tol.rd_nnorm = rd_tol; }

    /// Set the tolerance over the residual of complementarity measure (i.e. violation of orthogonality of forces and contact points distance)
    void SetComplementarityMeasureTolerance(double complementarity_tol) { res_nnorm_tol.mu = complementarity_tol; }

    /// Set the null pivot detection for the kernel solver.
    void SetNullPivotDetection(bool on_off, double threshold) { mumps_engine.SetNullPivotDetection(on_off, threshold); }

    // Test
    void DumpProblem(std::string suffix = "");
    void DumpIPStatus(std::string suffix = "") const;
    void RecordHistory(bool on_off, std::string file_name = "interior_point_log");
    int GetSolverCalls() const { return solver_call; }
    int GetIPSolverCalls() const { return ip_solver_call; }
    int GetIPIterations() const { return iteration_count_tot; }
    int GetMassMatrixDimension() const { return m; }
    int GetJacobianMatrixRows() const { return n; }
    double GetIPTimer() const { return ip_timer(); }
    //void Solve(const ChSparseMatrix& Q, const ChSparseMatrix& A, const ChMatrix<double>& rhs_b, const ChMatrix<double>& rhs_c, ChMatrix<double>& var_x, ChMatrix<double>& var_y, ChMatrix<double>& var_lam );
    void Solve(const ChCOOMatrix& normal_mat, const ChMatrix<double>& rhs_b, const ChMatrix<double>& rhs_c, ChMatrixDynamic<double>& var_x, ChMatrixDynamic<double>& var_y, ChMatrixDynamic<double>& var_lam);
};

}  // end of namespace chrono

#endif