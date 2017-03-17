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

#ifndef CHIPENGINE_H
#define CHIPENGINE_H

#include "ChApiInteriorPoint.h"
#include "chrono/solver/ChSystemDescriptor.h"
#include "chrono/solver/ChSolver.h"

//#define VTK_PLOT
#ifdef VTK_PLOT
#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkPlotPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkAxis.h>
#endif

#ifdef CHRONO_MUMPS
#include "chrono_mumps/ChCOOMatrix.h"
#include "chrono_mumps/ChMumpsEngine.h"
#endif

// Interior point methdon based on Numerical Optimization by Nocedal, Wright
// minimize 0.5*xT*G*x + xT*x while Ax>=b (16.54 pag.480)
// WARNING: FOR THE MOMENT THE CONSTRAINTS MUST BE INEQUALITIES
// Further references: (all pages number refers to [1] if not otherwise specified)
// [1] Nocedal&Wright, Numerical Optimization 2nd edition
// [2] D'Apuzzo et al., Starting-point strategies for an infeasible potential reduction method
// [3] Mangoni D., Tasora A., Solving Unilateral Contact Problems in Multibody Dynamics using a Primal-Dual Interior
// Point Method

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


namespace chrono {

/// Class for Interior-Point method
/// for QP convex programming

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

    const bool EQUAL_STEP_LENGTH = false;
    const bool ADAPTIVE_ETA = true;
    const bool ONLY_PREDICT = false;
    bool warm_start_broken = false;
    bool warm_start = true;

    double rp_nnorm_tol = 1e-8;
    double rd_nnorm_tol = 1e-8;
    double mu_tol = 1e-8;

#ifdef VTK_PLOT

    // vtk handlers
    vtkSmartPointer<vtkTable> table;
    vtkSmartPointer<vtkContextView> view;
    vtkSmartPointer<vtkIntArray> arr_call;
    vtkSmartPointer<vtkDoubleArray> arr_rpnnorm;
    vtkSmartPointer<vtkDoubleArray> arr_rdnnorm;
    vtkSmartPointer<vtkDoubleArray> arr_mu;
    vtkSmartPointer<vtkChartXY> chart;
#endif

    IP_KKT_SOLUTION_METHOD KKT_solve_method = IP_KKT_SOLUTION_METHOD::AUGMENTED;

    struct IPvariables_t {
        ChMatrixDynamic<double> x;    ///< DeltaSpeed/Acceleration ('q' in chrono)
        ChMatrixDynamic<double> y;    ///< Slack variable/Contact points distance ('c' in chrono)
        ChMatrixDynamic<double> lam;  ///< Lagrangian multipliers/contact|constraint forces ('l' in chrono)
    } var, Dvar;

    struct IPresidual_t {
        ChMatrixDynamic<double> rp;    ///< Residual about primal variables (i.e. violation if dynamic equation of motion); rp = A*x - y - b.
        ChMatrixDynamic<double> rd;    ///< Residual about dual variables (i.e. violation of constraints equations); rd = G*x - AT*lam + c.
        ChMatrixDynamic<double> rpd;    ///< Residual about primal-dual variables (only for #IP_KKT_SOLUTION_METHOD#NORMAL mode)
        double mu = 0;  ///< complementarity measure

        bool operator<(const IPresidual_t& other) const { return rp.NormTwo() < other.rp.NormTwo() && rd.NormTwo() < other.rd.NormTwo() && mu < other.mu; }
        bool operator<=(const IPresidual_t& other) const { return rp.NormTwo() <= other.rp.NormTwo() && rd.NormTwo() <= other.rd.NormTwo() && mu <= other.mu; }
        bool operator>(const IPresidual_t& other) const { return !(*this <= other); }
        bool operator>=(const IPresidual_t& other) const { return !(*this < other); }

    } res;

    struct IPresidual_nnorm_t {
        double rp_nnorm = 0;
        double rd_nnorm = 0;
        double mu = 0;
        IPresidual_nnorm_t(double rp_nnorm_in, double rd_nnorm_in, double mu_in) : rp_nnorm(rp_nnorm_in), rd_nnorm(rd_nnorm_in), mu(mu_in) {}
        IPresidual_nnorm_t() {}
        bool operator<(const IPresidual_nnorm_t& other) const { return rp_nnorm < other.rp_nnorm && rd_nnorm < other.rd_nnorm && mu < other.mu; }
        bool operator<=(const IPresidual_nnorm_t& other) const { return rp_nnorm <= other.rp_nnorm && rd_nnorm <= other.rd_nnorm && mu <= other.mu; }
        bool operator>(const IPresidual_nnorm_t& other) const { return !(*this <= other); }
        bool operator>=(const IPresidual_nnorm_t& other) const { return !(*this < other); }
        void update_residual_status(IPresidual_t& res)
        {
            rp_nnorm = res.rp.NormTwo() / res.rp.GetRows();
            rd_nnorm = res.rd.NormTwo() / res.rd.GetRows();
            mu = res.mu;
        }
    } res_nnorm_tol{1e-8, 1e-8, 1e-8};

    struct IPrhs_t {
        ChMatrixDynamic<double> b;  ///< rhs of constraints (is '-b' in chrono)
        ChMatrixDynamic<double> c;  ///< forces (is '-f' in chrono)
    } rhs;

    // Variables
    // ChMatrixDynamic<double> x; // DeltaSpeed/Acceleration ('q' in chrono)
    // ChMatrixDynamic<double> y; // Slack variable/Contact points distance ('c' in chrono)
    // ChMatrixDynamic<double> lam; // Lagrangian multipliers/contact|constraint forces ('l' in chrono)
    // ChMatrixDynamic<double> b; // rhs of constraints (is '-b' in chrono)
    // ChMatrixDynamic<double> c; // forces (is '-f' in chrono)

    // Temporaries used in iterate() function
    // ChMatrixDynamic<double> x_pred, y_pred, lam_pred;
    // ChMatrixDynamic<double> x_corr, y_corr, lam_corr;
    // ChMatrixDynamic<double> Dx, Dy, Dlam;
    //      ChMatrixDynamic<double> Dx_pre, Dy_pre, Dlam_pre;
    ChMatrixDynamic<double> vectn;  // temporary variable that has always size (#n,1)
    ChMatrixDynamic<double> vectm;  // temporary variable that has always size (#m,1)

    // Residuals
    // ChMatrixDynamic<double> rp; ///< Residual about primal variables (i.e. violation if dynamic equation of motion);
    // rp = A*x - y - b.  ChMatrixDynamic<double> rd; ///< Residual about dual variables (i.e. violation of constraints
    // equations); rd = G*x - AT*lam + c.  ChMatrixDynamic<double> rpd; ///< Residual about primal-dual variables (only
    // for #IP_KKT_SOLUTION_METHOD#NORMAL mode)  double mu = 0; // complementarity measure
    //      double rp_nnorm = 0;
    //      double rd_nnorm = 0;

    // problem matrices and vectors
    ChMatrixDynamic<double> rhs_sol;
    ChMatrixDynamic<double> sol_chrono;  // intermediate file to inject the IP solution into Chrono used in adapt_to_Chrono()
    ChCOOMatrix BigMat;
    ChCOOMatrix SmallMat;
    ChCOOMatrix E;  // compliance matrix

    // MUMPS engine
    ChMumpsEngine mumps_engine;

    // IP specific functions
    IPresidual_nnorm_t& iterate();  ///< Perform an IP iteration; returns \e true if exit conditions are met.
    void KKTsolve(double sigma, bool apply_correction);
    void set_starting_point(IP_STARTING_POINT_METHOD start_point_method, int n_old = 0, int m_old = 0);
    static double find_Newton_step_length(const ChMatrix<double>& vect, const ChMatrix<double>& Dvect, double tau = 1);
    double evaluate_objective_function();  ///< Evaluate the objective function i.e. 0.5*xT*G*x + xT*x.

    // Auxiliary
    void reset_dimensions(int n_old, int m_old);
    ChMatrix<>& adapt_to_Chrono(ChMatrix<>& solution_vect) const;
    void residual_fullupdate();  ///< Update #rp, #rd, and #mu from current #x, #y, #lam and the current system matrix.
    void make_positive_definite();  ///< Change A^T to -A^T in the current system matrix.
    void multiplyA(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
    void multiplyNegAT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;
    void multiplyG(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const;

    // Debug
    std::ofstream history_file;
    bool print_history;
    void LoadProblem();

  public:
    ChInteriorPoint();
    ~ChInteriorPoint();
    double Solve(ChSystemDescriptor& sysd) override;

    bool SolveRequiresMatrix() const override { return true; }

    // Auxiliary
    /// Set the Karush–Kuhn–Tucker problem form that will be used to solve the IP problem.
    void SetKKTSolutionMethod(IP_KKT_SOLUTION_METHOD qp_solve_type_selection) {
        KKT_solve_method = qp_solve_type_selection;
    }

    /// Set the maximum number of iterations after which the iteration loop will be stopped.
    void SetMaxIterations(int max_iter) { iteration_count_max = max_iter; }

    /// Set the tolerance over the residual of the primal variables (i.e. violation if dynamic equation of motion).
    void SetPrimalResidualTolerance(double rp_tol) { rp_nnorm_tol = rp_tol; }

    /// Set the tolerance over the residual of the dual variables (i.e. violation of constraints equations).
    void SetDualResidualTolerance(double rd_tol) { rd_nnorm_tol = rd_tol; }

    /// Set the tolerance over the residual of complementarity measure (i.e. violation of orthogonality of forces and
    /// contact points distance)
    void SetComplementarityMeasureTolerance(double complementarity_tol) { mu_tol = complementarity_tol; }

    // Test
    void DumpProblem(std::string suffix = "");
    void DumpIPStatus(std::string suffix = "") const;
    void RecordHistory(bool on_off, std::string filepath = "history_file.txt");
};

}  // end of namespace chrono

#endif