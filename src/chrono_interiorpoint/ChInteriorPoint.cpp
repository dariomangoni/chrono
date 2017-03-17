#include "ChInteriorPoint.h"
#include <algorithm>
#include "chrono_postprocess/ChVTK.h"

//#define DEBUG_MODE
#define SKIP_CONTACTS_UV true
#define ADD_COMPLIANCE false
#define REUSE_OLD_SOLUTIONS false

namespace chrono {

void PrintMatrix(ChMatrix<>& matrice) {
    for (auto i = 0; i < matrice.GetRows(); i++) {
        for (auto j = 0; j < matrice.GetColumns(); j++) {
            printf("%.1f ", matrice.GetElement(i, j));
        }
        printf("\n");
    }
}

void PrintCSR3Matrix(ChSparseMatrix& matrice) {
    for (auto i = 0; i < matrice.GetNumRows(); i++) {
        for (auto j = 0; j < matrice.GetNumColumns(); j++) {
            printf("%.1f ", matrice.GetElement(i, j));
        }
        printf("\n");
    }
}

template <class matrix>
void ExportArrayToFile(matrix mat, std::string filepath, int precision = 12) {
    std::ofstream ofile;
    ofile.open(filepath);
    ofile << std::scientific << std::setprecision(precision);

    for (auto row_sel = 0; row_sel < mat.GetRows(); row_sel++) {
        for (auto col_sel = 0; col_sel < mat.GetColumns(); col_sel++) {
            ofile << mat.GetElement(row_sel, col_sel);
        }

        ofile << std::endl;
    }

    ofile.close();
}

void ImportArrayFromFile(ChMatrix<>& output_mat, std::string filename) {
    std::ifstream my_file;
    my_file.open(filename);

    double temp;
    int row_sel = -1;
    for (row_sel = 0; row_sel < output_mat.GetRows(); row_sel++) {
        my_file >> temp;
        output_mat.SetElement(row_sel, 0, temp);
    }
    my_file.close();
}

double ChInteriorPoint::Solve(ChSystemDescriptor& sysd) {
    solver_call++;
    verbose = true;

    // The problem to be solved is loaded into the main matrix that will be used to solve the various step of the IP
    // method The initial guess is modified in order to be feasible The residuals are computed

    /********** Load system **********/
    // TODO: dimensions generally change at each call; 'm' for sure, but maybe 'n' does not?
    auto n_old = n;
    auto m_old = m;
    reset_dimensions(sysd.CountActiveVariables(), sysd.CountActiveConstraints(false, SKIP_CONTACTS_UV));

    // Load system matrix in 'BigMat', 'rhs_sol', 'b' and 'c'
    // Convert to different formats:
    // format = 0; (used throughout Chrono, but not here)
    //             | M  Cq'|*| q|-| f|=|0|
    //             | Cq  E | |-l| |-b| |c|
    //
    // format = 1; | M  Cq'|
    //             | Cq  0 |
    //
    // format = 2; | M   0  Cq'|
    //             | Cq  0   0 |
    //             | 0   0   0 |

    switch (KKT_solve_method) {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
            sysd.ConvertToMatrixForm(&BigMat, nullptr, false, SKIP_CONTACTS_UV, 2);
            make_positive_definite();
            break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
            sysd.ConvertToMatrixForm(&BigMat, nullptr, false, SKIP_CONTACTS_UV, 1);
            make_positive_definite();
            break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
            sysd.ConvertToMatrixForm(&SmallMat, &BigMat, nullptr, nullptr, nullptr, nullptr, false, true);
            break;
    }

    sysd.ConvertToMatrixForm(nullptr, nullptr, nullptr, &rhs.c, &rhs.b, nullptr, false, SKIP_CONTACTS_UV);  // load f->c and b->b
    rhs.c.MatrScale(-1);                         // adapt to InteriorPoint convention
    rhs.b.MatrScale(-1);                         // adapt to InteriorPoint convention

    /********** Check if system has constraints **********/
    if (m == 0)  // if no constraints
    {
        // Fill 'rhs_sol' with just Chrono's 'f' i.e. IP's '-c'
        for (auto row_sel = 0; row_sel < n; row_sel++)
            rhs_sol.SetElement(row_sel, 0, -rhs.c.GetElement(row_sel, 0));

        // Solve the KKT system
        BigMat.Compress();
        mumps_engine.SetProblem(BigMat, rhs_sol);
        if (mumps_engine.MumpsCall(ChMumpsEngine::COMPLETE))
            mumps_engine.PrintINFOG();

        if (verbose && mumps_engine.GetRINFOG(6) > 1e-6)
            std::cout << "MUMPS scaled residual: " << mumps_engine.GetRINFOG(6) << std::endl;

        residual_fullupdate();

        sysd.FromVectorToUnknowns(rhs_sol);

        // Export variable so that can be used in the next iteration as starting point
        //var.x.Resize(n, 1);
        //var.x = rhs_sol;

        if (verbose)
            std::cout << "IP call: " << solver_call << "; No constraints." << std::endl;

        return res.rd.NormTwo() / n;
    }

    /********* The system DOES have constraints! Start Interior Point ********/

    if (ADD_COMPLIANCE && m > 0) {
        sysd.ConvertToMatrixForm(nullptr, nullptr, &E, nullptr, nullptr, nullptr, false, SKIP_CONTACTS_UV);
        E *= -1;
    }

    //DumpProblem();

    set_starting_point(IP_STARTING_POINT_METHOD::NOCEDAL, n_old, m_old);

    for (iteration_count = 1; iteration_count < iteration_count_max && iterate() > res_nnorm_tol; iteration_count++) {
    }

    if (verbose)
        std::cout << "IP call: " << solver_call << "; iter: " << iteration_count << "/" << iteration_count_max
                  << std::endl;

    // Scatter the solution into the Chrono environment
    sysd.FromVectorToUnknowns(adapt_to_Chrono(sol_chrono));

    return 0.0;
}

// Iterating function
// output: (x, y, lam) are computed
// (res.rp, res.rd, res.mu) are updated based on most recent (x, y, lam)
// (res.rp, res.rd, res.mu, x, y, lam) are taken as they are
ChInteriorPoint::IPresidual_nnorm_t& ChInteriorPoint::iterate() {
/*********************************************************************************/
/***************************** Prediction Phase **********************************/
/*********************************************************************************/

    IPresidual_nnorm_t res_nnorm;
    IPvariables_t var_pred, Dvar_pred;
    var_pred.x.Resize(n, 1);      var_pred.y.Resize(m, 1);     var_pred.lam.Resize(m, 1);
    Dvar_pred.x.Resize(n, 1);    Dvar_pred.y.Resize(m, 1);    Dvar_pred.lam.Resize(m, 1);

    /*** find directions ***/
    // update y/lambda diagonal submatrix
    if (ADD_COMPLIANCE)
        for (auto diag_sel = 0; diag_sel < m; diag_sel++)
            BigMat.SetElement(
                n + diag_sel, n + diag_sel,
                var.y.GetElement(diag_sel, 0) / var.lam.GetElement(diag_sel, 0) + E.GetElement(diag_sel, diag_sel));
    else
        for (auto diag_sel = 0; diag_sel < m; diag_sel++)
            BigMat.SetElement( n + diag_sel, n + diag_sel, var.y.GetElement(diag_sel, 0) / var.lam.GetElement(diag_sel, 0) );

    BigMat.Compress();
    mumps_engine.SetProblem(BigMat, rhs_sol);
    if (mumps_engine.MumpsCall(ChMumpsEngine::ANALYZE_FACTORIZE))
        mumps_engine.PrintINFOG();
    // TODO: between different iterations only the bottom-right part of the matrix changes. This part is diagonal! Can't we avoid a FULL re-evalutation of all the factors?

    // fill 'rhs_sol' with rhs [-res.rd;-res.rp-y]
    for (auto row_sel = 0; row_sel < n; row_sel++)
        rhs_sol.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));
    for (auto row_sel = 0; row_sel < m; row_sel++)
        rhs_sol.SetElement(row_sel + n, 0, -res.rp(row_sel, 0) - var.y(row_sel, 0));

    // Solve the KKT system
    if (mumps_engine.MumpsCall(ChMumpsEngine::SOLVE))
        mumps_engine.PrintINFOG();
    if (mumps_engine.GetRINFOG(6) > 1e-6)
        std::cout << "Scaled residual norm of MUMPS call: " << mumps_engine.GetRINFOG(6) << std::endl;

    // Extract 'Dvar.x' and 'Dvar.lam' from 'sol'
    for (auto row_sel = 0; row_sel < n; row_sel++)
        Dvar_pred.x.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel, 0));
    for (auto row_sel = 0; row_sel < m; row_sel++)
        Dvar_pred.lam.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n, 0));

    // Calc 'Dvar.y' (it is also possible to evaluate Dvar.y as Dvar.y=(-lam°y+sigma*res.mu*e-y°Dvar.lam)./lam )
    multiplyA(Dvar_pred.x, Dvar_pred.y);  // Dvar.y = A*Dvar.x
    Dvar_pred.y += res.rp;
    if (ADD_COMPLIANCE) {
        E.MatrMultiply(Dvar_pred.lam, vectm);
        Dvar_pred.y += vectm;
    }

    /*** compute step lengths ***/
    // from 16.60 pag.482 from 14.32 pag.408 (remember that y>=0!)
    double alfa_pred_prim = find_Newton_step_length(var.y, Dvar_pred.y);
    double alfa_pred_dual = find_Newton_step_length(var.lam, Dvar_pred.lam);

    if (EQUAL_STEP_LENGTH) {
        double alfa_pred = std::min(alfa_pred_prim, alfa_pred_dual);
        alfa_pred_prim = alfa_pred;
        alfa_pred_dual = alfa_pred;
    }


    /*** make the prediction step ***/
    var_pred.y = Dvar_pred.y;        var_pred.y.MatrScale(alfa_pred_prim);        var_pred.y += var.y;
    var_pred.lam = Dvar_pred.lam;    var_pred.lam.MatrScale(alfa_pred_dual);    var_pred.lam += var.lam;

    /*** compute complementarity measure ***/
    double mu_pred = var_pred.y.MatrDot(var_pred.y, var_pred.lam) / m;  // from 16.56 pag.481

    if (ONLY_PREDICT) {
        var_pred.x = Dvar.x;        var_pred.x.MatrScale(alfa_pred_prim);        var_pred.x += var.x;        var.x = var_pred.x;
        var.y = var_pred.y;
        var.lam = var_pred.lam;

        res.rp.MatrScale(1 - alfa_pred_prim);

        multiplyG(Dvar.x, vectn);                          // vectn = G * Dvar.x
        vectn.MatrScale(alfa_pred_prim - alfa_pred_dual);  // vectn = (alfa_pred_prim - alfa_pred_dual) * (G * Dvar.x)
        res.rd.MatrScale(1 - alfa_pred_dual);
        res.rd += vectn;

        res.mu = mu_pred;

        res_nnorm.update_residual_status(res);

        return res_nnorm;
    }

    /*********************************************************************************/
    /******************************* Correction phase ********************************/
    /*********************************************************************************/

    /*** evaluate centering parameter ***/
    double sigma = std::pow(mu_pred / res.mu, 3.0);  // from 14.34 pag.408

    /*** find directions ***/
    for (auto row_sel = 0; row_sel < m; row_sel++)
        rhs_sol.SetElement(
            row_sel + n, 0,
            -res.rp(row_sel, 0) - var.y(row_sel, 0) +
                (sigma * res.mu - Dvar_pred.lam(row_sel, 0) * Dvar_pred.y(row_sel, 0)) / var.lam(row_sel, 0));

    // Solve the KKT system
    if (mumps_engine.MumpsCall(ChMumpsEngine::SOLVE))
        mumps_engine.PrintINFOG();
    if (verbose && mumps_engine.GetRINFOG(6) > 1e-6)
        std::cout << "Scaled residual norm of MUMPS call: " << mumps_engine.GetRINFOG(6) << std::endl;

    // Extract 'Dvar.x' and 'Dvar.lam' from 'sol'
    for (auto row_sel = 0; row_sel < n; row_sel++)
        Dvar.x.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel, 0));
    for (auto row_sel = 0; row_sel < m; row_sel++)
        Dvar.lam.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n, 0));

    // Calc 'Dvar.y' (it is also possible to evaluate Dvar.y as Dvar.y=(-lam°y+sigma*res.mu*e-y°Dvar.lam)./lam )
    multiplyA(Dvar.x, Dvar.y);  // Dvar.y = A*Dvar.x
    Dvar.y += res.rp;
    if (ADD_COMPLIANCE) {
        E.MatrMultiply(Dvar.lam, vectm);
        Dvar.y += vectm;
    }

    /*** step length correction ***/
    double tau = ADAPTIVE_ETA ? exp(-res.mu * m) * 0.1 + 0.9 : 0.95;  // exponential descent of tau

    /*** compute step lengths ***/
    double alfa_corr_prim = find_Newton_step_length(var.y, Dvar.y, tau);
    double alfa_corr_dual = find_Newton_step_length(var.lam, Dvar.lam, tau);

    if (EQUAL_STEP_LENGTH) {
        double alfa_corr = std::min(alfa_corr_prim, alfa_corr_dual);
        alfa_corr_prim = alfa_corr;
        alfa_corr_dual = alfa_corr;
    }

    IPvariables_t var_corr;
    var_corr.x.Resize(n, 1);
    var_corr.y.Resize(m, 1);
    var_corr.lam.Resize(m, 1);

    /*** make the correction step ***/
      var_corr.x = Dvar.x;      var_corr.x.MatrScale(alfa_corr_prim);        var_corr.x += var.x;        var.x = var_corr.x;
      var_corr.y = Dvar.y;      var_corr.y.MatrScale(alfa_corr_prim);        var_corr.y += var.y;        var.y = var_corr.y;
    var_corr.lam = Dvar.lam;    var_corr.lam.MatrScale(alfa_corr_dual);    var_corr.lam += var.lam;    var.lam = var_corr.lam;

    /********** Residuals update **********/
    res.rp.MatrScale(1 - alfa_corr_prim);
    res.rd.MatrScale(1 - alfa_corr_dual);
    res.mu = var.y.MatrDot(var.y, var.lam) / m;  // from 14.6 pag.395

    if (!EQUAL_STEP_LENGTH) {
        multiplyG(Dvar.x, vectn);                          // vectn = G*Dvar.x
        vectn.MatrScale(alfa_corr_prim - alfa_corr_dual);  // vectn = (alfa_pred_prim - alfa_pred_dual) * (G * Dvar.x)
        res.rd += vectn;
    }

    /*********************************************************************************/
    /******************************** Exit conditions ********************************/
    /*********************************************************************************/

    res_nnorm.update_residual_status(res);

    bool neg_y = false;
    bool neg_lam = false;
    for (int cont = 0; cont < m; cont++) {
        if (var.y(cont, 0) < 0)
            neg_y = true;
        if (var.lam(cont, 0) < 0)
            neg_lam = true;
    }

    if (neg_y)
        std::cout << "'y' has negative elements" << std::endl;

    if (neg_lam)
        std::cout << "'lam' has negative elements" << std::endl;

    if (verbose)
    {
        GetLog() << "InteriorPoint Results Pred+Corr\n";
        GetLog() << "Complementarity measure: " << res.mu << "\n";
        GetLog() << "|rd|/n (stationarity): " << res.rd.NormTwo() / n << "\n";
        GetLog() << "|rp|/m (constraint violation): " << res.rp.NormTwo() / m << "\n";
        GetLog() << "Objective Function: " << evaluate_objective_function() << "\n";
        GetLog() << "\n";
    }


#ifdef VTK_PLOT

    arr_call->InsertNextValue(solver_call * 100 + iteration_count);
    arr_rpnnorm->InsertNextValue(res_nnorm.rp_nnorm);
    arr_rdnnorm->InsertNextValue(res_nnorm.rd_nnorm);
    arr_mu->InsertNextValue(res.mu);

    if (arr_call->GetNumberOfValues() > 2) {
        chart->GetAxis(vtkAxis::BOTTOM)->SetMinimum(arr_call->GetValue(0));
        chart->GetAxis(vtkAxis::BOTTOM)->SetMaximum(arr_call->GetValue(arr_call->GetNumberOfValues() - 1));
        chart->GetAxis(vtkAxis::BOTTOM)->AutoScale();
        view->Render();
    }

    postprocess::ChVTKMatrixPlot vtk_matplot;
    vtk_matplot.PrintMatrix(BigMat);

#endif

    if (history_file.is_open())
        history_file << std::endl
                     << solver_call << ", " << iteration_count << ", " << res_nnorm.rp_nnorm << ", "
                     << res_nnorm.rd_nnorm << ", " << res.mu << ", " << evaluate_objective_function();

    DumpProblem("_end");

    return res_nnorm;
}

// Solve the KKT system in different modes: 'res.rp', 'res.rd', 'res.mu', 'x', 'y', 'lam' must be updated before calling
// this function 'sigma' is the centering parameter;
void ChInteriorPoint::KKTsolve(double sigma, bool apply_correction) {
    switch (KKT_solve_method) {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
            assert(false && "Standard mode to be fixed");
            // update lambda and y diagonal submatrices
            for (auto diag_sel = 0; diag_sel < m; diag_sel++) {
                BigMat.SetElement(n + m + diag_sel, n + diag_sel,
                                  var.lam.GetElement(diag_sel, 0));  // write lambda diagonal submatrix
                BigMat.SetElement(n + m + diag_sel, n + m + diag_sel,
                                  var.y.GetElement(diag_sel, 0));   // write y diagonal submatrix
                BigMat.SetElement(n + diag_sel, n + diag_sel, -1);  // write -identy_matrix diagonal submatrix
            }

            if (sigma != 0)  // rpd_corr
            {
                // I'm supposing that in 'res.rpd', since the previous call should have been without perturbation,
                // there is already y°lam
                vectm = Dvar.lam;  // I could use Dvar.lam directly, but it is not really clear
                vectm.MatrScale(Dvar.y);
                vectm.MatrInc(-sigma * res.mu);
                res.rpd += vectm;
            } else  // rpd_pred as (16.57 pag.481 suggests)
            {
                res.rpd = var.y;
                res.rpd.MatrScale(var.lam);
            }

            // Fill 'rhs_sol' with [-res.rd;-res.rp;-res.rpd]
            for (auto row_sel = 0; row_sel < n; row_sel++)
                rhs_sol.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));

            for (auto row_sel = 0; row_sel < m; row_sel++) {
                rhs_sol.SetElement(row_sel + n, 0, -res.rp.GetElement(row_sel, 0));
                rhs_sol.SetElement(row_sel + n + m, 0, -res.rpd.GetElement(row_sel, 0));
            }

            // Solve the KKT system
            BigMat.Compress();
            mumps_engine.SetProblem(BigMat, rhs_sol);
            printf("Mumps says: %d\n", mumps_engine.MumpsCall(ChMumpsEngine::SOLVE));

            // Extract 'Dvar.x', 'Dvar.y' and 'Dvar.lam' from 'sol'
            for (auto row_sel = 0; row_sel < n; row_sel++)
                Dvar.x.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel, 0));

            for (auto row_sel = 0; row_sel < m; row_sel++) {
                Dvar.y.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n, 0));
                Dvar.lam.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n + m, 0));
            }

            break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
            if (!apply_correction) {
                // update y/lambda diagonal submatrix
                if (ADD_COMPLIANCE)
                    for (auto diag_sel = 0; diag_sel < m; diag_sel++) {
                        BigMat.SetElement(n + diag_sel, n + diag_sel,
                                          var.y.GetElement(diag_sel, 0) / var.lam.GetElement(diag_sel, 0) +
                                              E.GetElement(diag_sel, diag_sel));
                    }
                else
                    for (auto diag_sel = 0; diag_sel < m; diag_sel++) {
                        BigMat.SetElement(n + diag_sel, n + diag_sel,
                                          var.y.GetElement(diag_sel, 0) / var.lam.GetElement(diag_sel, 0));
                    }

                BigMat.Compress();
                mumps_engine.SetProblem(BigMat, rhs_sol);
                mumps_engine.MumpsCall(ChMumpsEngine::ANALYZE_FACTORIZE);

                // fill 'rhs_sol' with rhs [-res.rd;-res.rp-y+sigma*res.mu/lam]
                for (auto row_sel = 0; row_sel < n; row_sel++)
                    rhs_sol.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));
            }

            for (auto row_sel = 0; row_sel < m; row_sel++)
                rhs_sol.SetElement(row_sel + n, 0,
                                   -res.rp(row_sel, 0) - var.y(row_sel, 0) + (apply_correction ? sigma * res.mu - Dvar.lam(row_sel, 0) * Dvar.y(row_sel, 0) : sigma * res.mu) / var.lam(row_sel, 0));

            // ExportArrayToFile(rhs_sol, "dump/rhs.txt");

            // Solve the KKT system
            mumps_engine.MumpsCall(ChMumpsEngine::SOLVE);
            mumps_engine.PrintINFOG();
            if (verbose) {
                double res_norm = mumps_engine.GetRINFOG(6);
                if (res_norm > 1e-6)
                    std::cout << "Scaled residual norm of MUMPS call: " << res_norm << std::endl;
            }

            // BigMat.ExportToDatFile("dump/COO.txt", true);
            // BigMat.ImportFromDatFile("COO.txt", true);
            // BigMat.ExportToDatFile("COO.txt", true);

            // ExportArrayToFile(rhs_sol, "dump/sol.txt");

            // Extract 'Dvar.x' and 'Dvar.lam' from 'sol'
            for (auto row_sel = 0; row_sel < n; row_sel++)
                Dvar.x.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel, 0));
            for (auto row_sel = 0; row_sel < m; row_sel++)
                Dvar.lam.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n, 0));

            // Calc 'Dvar.y' (it is also possible to evaluate Dvar.y as Dvar.y=(-lam°y+sigma*res.mu*e-y°Dvar.lam)./lam )
            multiplyA(Dvar.x, Dvar.y);  // Dvar.y = A*Dvar.x
            Dvar.y += res.rp;
            if (ADD_COMPLIANCE) {
                E.MatrMultiply(Dvar.lam, vectm);
                Dvar.y += vectm;
            }

            break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
            assert(false && "Normal mode to be fixed");
            for (auto row_sel = 0; row_sel < n; row_sel++) {
                for (auto col_sel = 0; col_sel < n; col_sel++) {
                    double temp = 0;
                    for (auto el_sel = 0; el_sel < m; el_sel++) {
                        temp += var.lam(el_sel, 0) / var.y(el_sel, 0) * SmallMat.GetElement(el_sel, row_sel) *
                                SmallMat.GetElement(el_sel, col_sel);
                        if (temp != 0)
                            BigMat.SetElement(row_sel, col_sel, temp, false);
                    }
                }
            }

            break;
    }
}

void ChInteriorPoint::set_starting_point(IP_STARTING_POINT_METHOD start_point_method, int n_old, int m_old) {
    switch (start_point_method) {
        case IP_STARTING_POINT_METHOD::STP1: {
            double infeas_dual_ratio = 0.1;  // TODO: dependant on n

            var.x.FillElem(1);
            var.y.FillElem(1);
            var.lam.FillElem(1);

            double duality_gap_calc = var.y.MatrDot(var.y, var.lam);  // [2] pag. 132
            double duality_gap = m;                                   // [2] pag. 132
            assert(duality_gap_calc == duality_gap);

            // norm of all residuals; [2] pag. 132
            residual_fullupdate();
            double res_norm = res.rp.MatrDot(res.rp, res.rp);
            res_norm += res.rp.MatrDot(res.rd, res.rd);
            res_norm = sqrt(res_norm);

            if (res_norm / duality_gap > infeas_dual_ratio) {
                double coeff = res_norm / (duality_gap * infeas_dual_ratio);
                var.x.MatrScale(coeff);
                var.y.MatrScale(coeff);
                var.lam.MatrScale(coeff);

                residual_fullupdate();
            }
        } break;

        case IP_STARTING_POINT_METHOD::STP2: {
            double threshold = 1;  // 'epsilon' in [2]

            if (!REUSE_OLD_SOLUTIONS || n != n_old) {
                // initialize x
                var.x.FillElem(1);
            }

            // initialize y and then lam
            multiplyA(var.x, vectm);
            vectm -= rhs.b;
            for (auto cont = 0; cont < m; cont++) {
                var.y(cont, 0) = vectm(cont, 0) > threshold ? vectm(cont, 0) : threshold;
                var.lam(cont, 0) = 1 / var.y(cont, 0);
            }

            residual_fullupdate();

        } break;

        case IP_STARTING_POINT_METHOD::NOCEDAL: {
            /********** Initialize IP algorithm **********/
            // Initial guess
            if (n_old != n || solver_call == 0 || !REUSE_OLD_SOLUTIONS)
                var.x.FillElem(1);  // TIP: every ChMatrix is initialized with zeros by default
            if (m_old != m || solver_call == 0 || !REUSE_OLD_SOLUTIONS)
                var.lam.FillElem(1);  // each element of lam will be at the denominator; avoid zeros!

            // since A is generally changed between calls, also with warm_start,
            // all the residuals and feasibility check must be redone
            multiplyA(var.x, var.y);  // y = A*x
            var.y -= rhs.b;

            // Calculate the residual
            residual_fullupdate();

            // Feasible starting Point (pag.484-485)
            KKTsolve(0, false);

            // x is accepted as it is
            var.y += Dvar.y;      // calculate y0
            var.lam += Dvar.lam;  // calculate lam0

            for (auto row_sel = 0; row_sel < m; row_sel++)
                var.y(row_sel) = abs(var.y(row_sel)) < 1 ? 1 : abs(var.y(row_sel));

            for (auto row_sel = 0; row_sel < m; row_sel++)
                var.lam(row_sel) = abs(var.lam(row_sel)) < 1 ? 1 : abs(var.lam(row_sel));

            // Update the residual considering the new values of 'y' and 'lam'
            residual_fullupdate();

        } break;

        case IP_STARTING_POINT_METHOD::NOCEDAL_WS: { /*Backup vectors*/
            ChMatrixDynamic<double> x_bkp(var.x);
            ChMatrixDynamic<double> y_bkp(var.y);
            ChMatrixDynamic<double> lam_bkp(var.lam);
            residual_fullupdate();
            double residual_value_bkp = res.rp.NormTwo() + res.rd.NormTwo() + res.mu * m;

            /********** Initialize IP algorithm **********/
            // Initial guess
            if (n_old != n || solver_call == 0)
                var.x.FillElem(1);  // TIP: every ChMatrix is initialized with zeros by default
            if (m_old != m || solver_call == 0)
                var.lam.FillElem(1);  // each element of lam will be at the denominator; avoid zeros!

            // since A is generally changed between calls, also with warm_start,
            // all the residuals and feasibility check must be redone
            multiplyA(var.x, var.y);  // y = A*x
            var.y -= rhs.b;

            // Calculate the residual
            residual_fullupdate();

            // Feasible starting Point (pag.484-485)
            KKTsolve(0, false);  // to obtain Dvar.x, Dvar.y, Dvar.lam called "affine"

            // x is accepted as it is
            var.y += Dvar.y;      // calculate y0
            var.lam += Dvar.lam;  // calculate lam0

            for (auto row_sel = 0; row_sel < m; row_sel++)
                var.y(row_sel) = abs(var.y(row_sel)) < 1 ? 1 : abs(var.y(row_sel));

            for (auto row_sel = 0; row_sel < m; row_sel++)
                var.lam(row_sel) = abs(var.lam(row_sel)) < 1 ? 1 : abs(var.lam(row_sel));

            // Update the residual considering the new values of 'y' and 'lam'
            residual_fullupdate();

            /* Check if restoring previous values would be better */
            double residual_value_new = res.rp.NormTwo() + res.rd.NormTwo() + res.mu * m;

            if (residual_value_bkp < residual_value_new) {
                var.x = x_bkp;
                var.y = y_bkp;
                var.lam = lam_bkp;
                residual_fullupdate();
            } else {
                std::cout << "Not WS\n";
            }

        } break;
        default:;
    }
}

/// Find the maximum step length, along the direction defined by \p Dvect, so that \p vect has no negative components;
/// It is applied to #lam and #y.
double ChInteriorPoint::find_Newton_step_length(const ChMatrix<double>& vect,
                                                const ChMatrix<double>& Dvect,
                                                double tau) {
    double alpha = 1;  // in this way alpha is clamped to a maximum of 1
    bool alpha_found = false;
    for (auto row_sel = 0; row_sel < vect.GetRows(); row_sel++) {
        if (Dvect(row_sel, 0) < 0) {
            double alfa_temp = - tau * vect(row_sel, 0) / Dvect(row_sel, 0);
            if (alfa_temp < alpha && alfa_temp > 0) {
                alpha_found = true;
                alpha = alfa_temp;
            }
        }
    }

    //if (!alpha_found || alpha <= 0)
    //    alpha = 0.001;

    return alpha;
}

double ChInteriorPoint::evaluate_objective_function() {
    multiplyG(var.x, vectn);
    auto obj_value = vectn.MatrDot(var.x, vectn);
    obj_value += rhs.c.MatrDot(var.x, rhs.c);

    return obj_value;
}


    /// Take care of adapting the size of matrices to \p n_new and \p m_new
void ChInteriorPoint::reset_dimensions(int n_new, int m_new) {
    if (n != n_new) {
        var.x.Resize(n_new, 1);
        Dvar.x.Resize(n_new, 1);
        rhs.c.Resize(n_new, 1);
        res.rd.Resize(n_new, 1);
        vectn.Resize(n_new, 1);
    }

    if (m != m_new) {
        var.y.Resize(m_new, 1);
        var.lam.Resize(m_new, 1);
        Dvar.y.Resize(m_new, 1);
        Dvar.lam.Resize(m_new, 1);
        rhs.b.Resize(m_new, 1);
        res.rp.Resize(m_new, 1);
        res.rpd.Resize(m_new, 1);
        vectm.Resize(m_new, 1);
    }

    SKIP_CONTACTS_UV ? sol_chrono.Resize(n_new + 3 * m_new, 1) : sol_chrono.Resize(n_new + m_new, 1);

    // BigMat and sol
    switch (KKT_solve_method) {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
            BigMat.Reset(2 * m_new + n_new, 2 * m_new + n_new, static_cast<int>(n_new * n_new * SPM_DEF_FULLNESS));
            rhs_sol.Resize(2 * m_new + n_new, 1);
            break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
            BigMat.Reset(n_new + m_new, n_new + m_new, static_cast<int>(n_new * n_new * SPM_DEF_FULLNESS));
            rhs_sol.Resize(n_new + m_new, 1);
            break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
            std::cout << std::endl << "Perturbed KKT system cannot be stored with 'NORMAL' method yet.";
            break;
    }

    n = n_new;
    m = m_new;
}

void ChInteriorPoint::DumpProblem(std::string suffix) {
    ExportArrayToFile(var.y, "dump/y" + suffix + ".txt");
    ExportArrayToFile(var.x, "dump/x" + suffix + ".txt");
    ExportArrayToFile(var.lam, "dump/lam" + suffix + ".txt");

    ExportArrayToFile(rhs.b, "dump/b" + suffix + ".txt");
    ExportArrayToFile(rhs.c, "dump/c" + suffix + ".txt");

    BigMat.Compress();
    BigMat.ExportToDatFile("dump/", 8);
}

void ChInteriorPoint::LoadProblem() {
    // ImportArrayFromFile(y, "dump/y.txt");
    // ImportArrayFromFile(x, "dump/x.txt");
    // ImportArrayFromFile(lam, "dump/lam.txt");

    ImportArrayFromFile(rhs.b, "dump/b.txt");
    ImportArrayFromFile(rhs.c, "dump/c.txt");

    // BigMat.ImportFromDatFile("dump/");
}

void ChInteriorPoint::DumpIPStatus(std::string suffix) const {
    ExportArrayToFile(var.y, "dump/y" + suffix + ".txt");
    ExportArrayToFile(var.x, "dump/x" + suffix + ".txt");
    ExportArrayToFile(var.lam, "dump/lam" + suffix + ".txt");

    ExportArrayToFile(Dvar.x, "dump/Dvar.x" + suffix + ".txt");
    ExportArrayToFile(Dvar.y, "dump/Dvar.y" + suffix + ".txt");
    ExportArrayToFile(Dvar.lam, "dump/Dvar.lam" + suffix + ".txt");

    ExportArrayToFile(rhs_sol, "dump/rhs_sol" + suffix + ".txt");
    // ExportArrayToFile(sol, "dump/sol" + suffix + ".txt");
}

void ChInteriorPoint::make_positive_definite() {
    if (m == 0)
        return;

    int offset_AT_col = n;
    if (KKT_solve_method == IP_KKT_SOLUTION_METHOD::STANDARD)
        offset_AT_col = n + m;

    BigMat.ForEachExistentValueInRange([](double* val) { *val *= -1; }, 0, n, offset_AT_col,
                                       BigMat.GetNumColumns() - 1);
}

// Perform moltiplication of A with vect_in: vect_out = A*vect_in
void ChInteriorPoint::multiplyA(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
    switch (KKT_solve_method) {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
            BigMat.MatrMultiplyClipped(vect_in, vect_out, n, n + m - 1, 0, n - 1, 0, 0);
            break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
            BigMat.MatrMultiplyClipped(vect_in, vect_out, n, n + m - 1, 0, n - 1, 0, 0);
            break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
            std::cout << std::endl << "A multiplication is not implemented in 'NORMAL' method yet.";
            break;
    }
}

// Perform moltiplication of -AT with vect_in: vect_out = -AT*vect_in (considers that in the top-right part there is
// already -A^T)
void ChInteriorPoint::multiplyNegAT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
    switch (KKT_solve_method) {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
            BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, n + m, n + 2 * m - 1, 0, 0);
            break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
            BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, n, n + m - 1, 0, 0);
            break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
            std::cout << std::endl << "AT multiplication is not implemented in 'NORMAL' method yet.";
            break;
    }
}

void ChInteriorPoint::multiplyG(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
    switch (KKT_solve_method) {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
            BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, 0, n - 1, 0, 0);
            break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
            BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, 0, n - 1, 0, 0);
            break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
            std::cout << std::endl << "G multiplication is not implemented in 'NORMAL' method yet.";
            break;
    }
}

void ChInteriorPoint::residual_fullupdate() {
    // Residual initialization (16.59 pag.482)

    // Dual Residual
    // res.rd = G*x + c - A^T*lam
    multiplyG(var.x, res.rd);  // res.rd = G*x
    res.rd += rhs.c;           // res.rd = G*x + c

    if (m > 0) {
        multiplyNegAT(var.lam, vectn);  // vectn = (-A^T)*lam
        res.rd += vectn;                // res.rd = (G*x + c) + (-A^T*lam)

        // Primal residual
        // res.rp = A*x - y - b
        multiplyA(var.x, res.rp);  // res.rp = A*x
        res.rp -= var.y;
        res.rp -= rhs.b;
        if (ADD_COMPLIANCE) {
            E.MatrMultiply(Dvar.lam, vectm);
            res.rp += vectm;
        }

        res.mu = var.y.MatrDot(var.y, var.lam) / m;
    } else {
        res.rp.FillElem(0);
        res.mu = 0;
    }
}

/// Export the IP variables in Chrono format
// TODO: FromVectorToUnknowns should accept const reference, but it doesn't. When it will be, we could return const
// ChMatrix<>&
ChMatrix<>& ChInteriorPoint::adapt_to_Chrono(ChMatrix<>& solution_vect) const {

    // copy 'x'
    for (auto row_sel = 0; row_sel < n; row_sel++)
        solution_vect(row_sel, 0) = var.x(row_sel, 0);

    // copy Lagrangian multipliers; skip tangential forces if needed
    if (SKIP_CONTACTS_UV) {
        // TODO: is 'lam' to be inverted?
        for (auto row_sel = 0; row_sel < m; row_sel++) {
            solution_vect(n + row_sel * 3, 0) = -var.lam(row_sel, 0);  // there will be an inversion inside FromVectorToUnknowns()
            solution_vect(n + row_sel * 3 + 1, 0) = 0;
            solution_vect(n + row_sel * 3 + 2, 0) = 0;
        }
    } else {
        for (auto row_sel = 0; row_sel < m; row_sel++)
            solution_vect(row_sel + n, 0) =
                -var.lam(row_sel);  // there will be an inversion inside FromVectorToUnknowns()
    }

    return solution_vect;
}

void ChInteriorPoint::RecordHistory(bool on_off, std::string filepath) {
    if (!history_file.is_open()) {
        history_file.open(filepath);
    }

    history_file << std::scientific << std::setprecision(3);
    history_file << "SolverCall, Iteration, res_nnorm.rp_nnorm, rd_nnormm, res.mu, obj_fun";

    print_history = true;
}

ChInteriorPoint::ChInteriorPoint() {
    mumps_engine.SetICNTL(11, 2);
    RecordHistory(true);

#ifdef VTK_PLOT
    // Set up the data table
    table = vtkSmartPointer<vtkTable>::New();

    arr_call = vtkSmartPointer<vtkIntArray>::New();
    arr_rpnnorm = vtkSmartPointer<vtkDoubleArray>::New();
    arr_rdnnorm = vtkSmartPointer<vtkDoubleArray>::New();
    arr_mu = vtkSmartPointer<vtkDoubleArray>::New();

    arr_call->SetName("call");
    arr_rpnnorm->SetName("res_nnorm.rp_nnorm");
    arr_rdnnorm->SetName("res_nnorm.rd_nnorm");
    arr_mu->SetName("res.mu");

    table->AddColumn(arr_call);
    table->AddColumn(arr_rpnnorm);
    table->AddColumn(arr_rdnnorm);
    table->AddColumn(arr_mu);

    // Set up the view
    view = vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    // Add multiple line plots, setting the colors etc
    chart = vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);
    chart->GetAxis(vtkAxis::LEFT)->SetLogScale(true);
    chart->GetAxis(vtkAxis::LEFT)->SetNotation(vtkAxis::SCIENTIFIC_NOTATION);
    chart->SetRenderEmpty(true);
    chart->SetShowLegend(true);

    vtkPlot* line;
    line = chart->AddPlot(vtkChart::LINE);
    vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CROSS);
    line->SetInputData(table, 0, 1);
    line->SetColor(255, 0, 0, 255);
    line->SetWidth(1.0);

    line = chart->AddPlot(vtkChart::LINE);
    vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CROSS);
    line->SetInputData(table, 0, 2);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(1.0);

    line = chart->AddPlot(vtkChart::LINE);
    vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CROSS);
    line->SetInputData(table, 0, 3);
    line->SetColor(0, 0, 255, 255);
    line->SetWidth(1.0);

// For dotted line, the line type can be from 2 to 5 for different dash/dot
// patterns (see enum in vtkPen containing DASH_LINE, value 2):
#ifndef WIN32
    line->GetPen()->SetLineType(vtkPen::DASH_LINE);
#endif
// (ifdef-ed out on Windows because DASH_LINE does not work on Windows
//  machines with built-in Intel HD graphics card...)

// view->GetRenderWindow()->SetMultiSamples(0);
// view->Render();

#endif
}

ChInteriorPoint::~ChInteriorPoint() {
    if (history_file.is_open())
        history_file.close();
}
}
