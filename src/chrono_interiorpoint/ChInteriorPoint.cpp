#include "ChInteriorPoint.h"
#include <algorithm>

#ifdef CHRONO_POSTPROCESS
#include "chrono_postprocess/ChGnuPlot.h"
#endif

//#define DEBUG_MODE
#define SKIP_CONTACTS_UV true
#define ADD_COMPLIANCE false
#define REUSE_OLD_SOLUTIONS true
#define STEPLENGTH_METHOD 0
//#define BYPASS_RESTITUTION_TERM

namespace chrono {

void PrintMatrix(ChMatrix<>& matrice) {
    for( auto i = 0; i < matrice.GetRows(); i++ )
    {
        for( auto j = 0; j < matrice.GetColumns(); j++ )
        {
            printf("%.1f ", matrice.GetElement(i, j));
        }
        printf("\n");
    }
}

void PrintCSR3Matrix(ChSparseMatrix& matrice) {
    for( auto i = 0; i < matrice.GetNumRows(); i++ )
    {
        for( auto j = 0; j < matrice.GetNumColumns(); j++ )
        {
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

    for( auto row_sel = 0; row_sel < mat.GetRows(); row_sel++ )
    {
        for( auto col_sel = 0; col_sel < mat.GetColumns(); col_sel++ )
        {
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
    for( row_sel = 0; row_sel < output_mat.GetRows(); row_sel++ )
    {
        my_file >> temp;
        output_mat.SetElement(row_sel, 0, temp);
    }
    my_file.close();
}

double ChInteriorPoint::Solve(ChSystemDescriptor& sysd) {
    solver_call++;

    // The problem to be solved is loaded into the main matrix that will be used to solve the various step of the IP
    // method The initial guess is modified in order to be feasible The residuals are computed

    /********** Load system **********/
    // TODO: dimensions generally change at each call; 'm' for sure, but maybe 'n' does not?
    auto n_old = n;
    auto m_old = m;
    reset_internal_dimensions(sysd.CountActiveVariables(), sysd.CountActiveConstraints(false, SKIP_CONTACTS_UV));

    // Load system matrix in 'BigMat', 'mumps_rhs', 'b' and 'c'
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

    switch( KKT_solve_method )
    {
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

#ifdef BYPASS_RESTITUTION_TERM
    rhs.b.FillElem(0); // TODO: WARNING: removing phi/h!!!
    //rhs.b*=1e-3; // TODO: WARNING: removing phi/h!!!
#endif

    /********** Check if system has constraints **********/
    if( m == 0 )  // if no constraints
    {
        ChMatrixDynamic<double> mumps_rhs(n, 1);

        // Fill 'mumps_rhs' with just Chrono's 'f' i.e. IP's '-c'
        for( auto row_sel = 0; row_sel < n; row_sel++ )
            mumps_rhs.SetElement(row_sel, 0, -rhs.c.GetElement(row_sel, 0));

        // Solve the KKT system
        BigMat.Compress();
        mumps_engine.SetProblem(BigMat, mumps_rhs);
        if( mumps_engine.MumpsCall(ChMumpsEngine::COMPLETE) )
            mumps_engine.PrintINFOG();

        if( verbose && mumps_engine.GetRINFOG(6) > 1e-6 )
            std::cout << "MUMPS scaled residual: " << mumps_engine.GetRINFOG(6) << std::endl;

        sysd.FromVectorToUnknowns(mumps_rhs);

        // Export variable so that can be used in the next iteration as starting point
        var.x.Resize(n, 1);
        var.x = mumps_rhs;

        if( verbose )
            std::cout << "IP call: " << solver_call << "; No constraints." << std::endl;

        return 0.0;
    }

    /********* The system DOES have constraints! Start Interior Point ********/
    ip_solver_call++;
    ip_timer.start();

    if( ADD_COMPLIANCE && m > 0 )
    {
        sysd.ConvertToMatrixForm(nullptr, nullptr, &E, nullptr, nullptr, nullptr, false, SKIP_CONTACTS_UV);
        E *= -1;
    }

    DumpProblem();

    set_starting_point(starting_point_method, n_old, m_old);

    for( iteration_count = 1; iteration_count < iteration_count_max; iteration_count++ )
    {

        iterate();
        iteration_count_tot++;

        /*********************************************************************************/
        /******************************** Exit conditions ********************************/
        /*********************************************************************************/


        if( verbose )
        {
            for (auto cont = 0; cont < m; cont++)
            {
                if (var.y(cont, 0) < 0 || var.lam(cont, 0) < 0)
                {
                    std::cout << "'y' or 'lambda' have negative elements" << std::endl;
                    break;
                }
            }
            GetLog() << "InteriorPoint | Call: " << solver_call << " Iter: " << iteration_count << "/" << iteration_count_max << "\n";
            GetLog() << "Complementarity Measure: " << res.mu << "\n";
            GetLog() << "|rd|/n (stationarity): " << res.rd.NormTwo() / n << "\n";
            GetLog() << "|rp|/m (constraint violation): " << res.rp.NormTwo() / m << "\n";
            GetLog() << "Objective Function: " << evaluate_objective_function() << "\n";
            GetLog() << "\n";
        }

        //DumpProblem("_end");

        if( res < res_nnorm_tol )
        {
            if( verbose )
                std::cout << "IP call: " << solver_call << "; iter: " << iteration_count << "/" << iteration_count_max
                << std::endl;


            break;
        }
    }

    ip_timer.stop();

    // Scatter the solution into the Chrono environment
    sysd.FromVectorToUnknowns(adapt_to_Chrono(sol_chrono));

    return 0.0;
}

// Iterating function
// output: (x, y, lam) are computed
// (res.rp, res.rd, res.mu) are updated based on most recent (x, y, lam)
// (res.rp, res.rd, res.mu, x, y, lam) are taken as they are
void ChInteriorPoint::iterate() {
    /*********************************************************************************/
    /***************************** Prediction Phase **********************************/
    /*********************************************************************************/

    setup_system_matrix(var);

    ChMatrixDynamic<double> mumps_rhs(n + m, 1);
    IPvariables_t Dvar_pred;
    Dvar_pred.x.Resize(n, 1);    Dvar_pred.y.Resize(m, 1);    Dvar_pred.lam.Resize(m, 1);

    // WARNING: the residual structure 'res' must be already updated at this point!

    // fill 'mumps_rhs' with rhs [-res.rd;-res.rp-y]
    for( auto row_sel = 0; row_sel < n; row_sel++ )
        mumps_rhs.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));
    for( auto row_sel = 0; row_sel < m; row_sel++ )
        mumps_rhs.SetElement(row_sel + n, 0, -res.rp(row_sel, 0) - var.y(row_sel, 0));

    makeNewtonStep(Dvar_pred, mumps_rhs, res);

#if USE_MESZEROS_STEPLENGTH
    double alfa_pred_prim, alfa_pred_dual;
    find_Newton_step_length(var, Dvar_pred, 1.0, alfa_pred_prim, alfa_pred_dual);
#else
    /*** compute step lengths ***/
    // from 16.60 pag.482 from 14.32 pag.408 (remember that y>=0!)
    auto alfa_pred_prim = find_Newton_step_length(var.y, Dvar_pred.y);
    auto alfa_pred_dual = find_Newton_step_length(var.lam, Dvar_pred.lam);

    if( EQUAL_STEP_LENGTH )
    {
        auto alfa_pred = std::min(alfa_pred_prim, alfa_pred_dual);
        alfa_pred_prim = alfa_pred;
        alfa_pred_dual = alfa_pred;
    }
#endif


    /*** make the prediction step ***/
    IPvariables_t var_pred;
    var_pred.x.Resize(n, 1);      var_pred.y.Resize(m, 1);     var_pred.lam.Resize(m, 1);

    var_pred.y = Dvar_pred.y;        var_pred.y.MatrScale(alfa_pred_prim);        var_pred.y += var.y;
    var_pred.lam = Dvar_pred.lam;    var_pred.lam.MatrScale(alfa_pred_dual);    var_pred.lam += var.lam;

    /*** compute complementarity measure ***/
    auto mu_pred = var_pred.y.MatrDot(var_pred.y, var_pred.lam) / m;  // from 16.56 pag.481

    if( ONLY_PREDICT )
    {
        var_pred.x = Dvar_pred.x;        var_pred.x.MatrScale(alfa_pred_prim);        var_pred.x += var.x;        var.x = var_pred.x;
        var.y = var_pred.y;
        var.lam = var_pred.lam;

        // update residuals
        res.rp.MatrScale(1 - alfa_pred_prim);

        multiplyG(Dvar_pred.x, vectn);                          // vectn = G * Dvar.x
        vectn.MatrScale(alfa_pred_prim - alfa_pred_dual);  // vectn = (alfa_pred_prim - alfa_pred_dual) * (G * Dvar.x)
        res.rd.MatrScale(1 - alfa_pred_dual);
        res.rd += vectn;

        res.mu = mu_pred;

        return;
    }

    /*********************************************************************************/
    /******************************* Correction phase ********************************/
    /*********************************************************************************/
    IPvariables_t Dvar;
    Dvar.x.Resize(n, 1);    Dvar.y.Resize(m, 1);    Dvar.lam.Resize(m, 1);

    /*** evaluate centering parameter ***/
    auto sigma = std::pow(mu_pred / res.mu, 3);  // from 14.34 pag.408

    /*** find directions ***/
    for( auto row_sel = 0; row_sel < n; row_sel++ )
        mumps_rhs.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));
    for( auto row_sel = 0; row_sel < m; row_sel++ )
        mumps_rhs.SetElement(row_sel + n, 0, -res.rp(row_sel, 0) - var.y(row_sel, 0) + (sigma * res.mu - Dvar_pred.lam(row_sel, 0) * Dvar_pred.y(row_sel, 0)) / var.lam(row_sel, 0));

    makeNewtonStep(Dvar, mumps_rhs, res);

    /*** step length correction ***/
    auto tau = ADAPTIVE_ETA ? exp(-res.mu * m) * 0.1 + 0.9 : 0.95;  // exponential descent of tau

#if USE_MESZEROS_STEPLENGTH
    double alfa_corr_prim, alfa_corr_dual;
    find_Newton_step_length(var, Dvar, tau, alfa_corr_prim, alfa_corr_dual);
#else
    /*** compute step lengths ***/
    auto alfa_corr_prim = find_Newton_step_length(var.y, Dvar.y, tau);
    auto alfa_corr_dual = find_Newton_step_length(var.lam, Dvar.lam, tau);

    if( EQUAL_STEP_LENGTH )
    {
        auto alfa_corr = std::min(alfa_corr_prim, alfa_corr_dual);
        alfa_corr_prim = alfa_corr;
        alfa_corr_dual = alfa_corr;
    }
#endif

    IPvariables_t var_corr;
    var_corr.x.Resize(n, 1);    var_corr.y.Resize(m, 1);    var_corr.lam.Resize(m, 1);

    /*** make the correction step ***/
    var_corr.x = Dvar.x;        var_corr.x.MatrScale(alfa_corr_prim);        var_corr.x += var.x;        var.x = var_corr.x;
    var_corr.y = Dvar.y;        var_corr.y.MatrScale(alfa_corr_prim);        var_corr.y += var.y;        var.y = var_corr.y;
    var_corr.lam = Dvar.lam;    var_corr.lam.MatrScale(alfa_corr_dual);    var_corr.lam += var.lam;    var.lam = var_corr.lam;

    /********** Residuals update **********/
    res.rp.MatrScale(1 - alfa_corr_prim);
    res.rd.MatrScale(1 - alfa_corr_dual);
    res.mu = var.y.MatrDot(var.y, var.lam) / m;  // from 14.6 pag.395

    if( !EQUAL_STEP_LENGTH )
    {
        multiplyG(Dvar.x, vectn);                          // vectn = G*Dvar.x
        vectn.MatrScale(alfa_corr_prim - alfa_corr_dual);  // vectn = (alfa_pred_prim - alfa_pred_dual) * (G * Dvar.x)
        res.rd += vectn;
    }


    /*** Data logging ***/
    if( print_history )
    {
        if( logfile_stream.is_open() )
        {
            logfile_stream << std::endl << res.rp.NormTwo() / m << ", " << res.rd.NormTwo() / n << ", " << res.mu << ", " << evaluate_objective_function();
        }
        else
        {
            logfile_stream.open(logfile_name + ".txt");

            if( logfile_stream.is_open() )
            {
                logfile_stream << std::scientific << std::setprecision(6);
                logfile_stream << "res_P res_D mu objfun";

                logfile_stream << std::endl << res.rp.NormTwo() / m << ", " << res.rd.NormTwo() / n << ", " << res.mu << ", " << evaluate_objective_function();

#ifdef CHRONO_POSTPROCESS
                postprocess::ChGnuPlot mplot(("__" + logfile_name + ".gpl").c_str());

                auto column_number = 3;
                mplot.SetGrid();

                mplot.SetCommand("bind s 'plot_history = 0'");
                mplot.SetCommand("bind 'Close' 'plot_history = 0'");
                mplot.SetCommand("plot_history = 1");
                mplot.SetCommand("set key autotitle columnhead");
                mplot.SetCommand("set terminal wxt size 800,600 enhanced font 'Verdana,10' title 'Residuals' persist");

                mplot.SetLabelX("iterations");
                mplot.SetLabelY("residuals");
                mplot.SetCommand("set logscale y");
                mplot.SetCommand("set format y \"%T\"");
                std::stringstream comm_template;
                comm_template << "plot ";
                comm_template << "'" << logfile_name << ".txt' using " << 1 << " with linespoints" << ",";
                comm_template << "'" << logfile_name << ".txt' using " << 2 << " with linespoints" << ",";
                comm_template << "'" << logfile_name << ".txt' using " << 3 << " with linespoints";
                mplot.SetCommand(comm_template.str().c_str());
                mplot.SetCommand("pause 0.1");
                mplot.SetCommand("if(plot_history==1) reread");
#endif
            }
            else
                throw std::exception(("Log file " + logfile_name + " cannot be opened\n").c_str());
        }
    }



}

void ChInteriorPoint::setup_system_matrix(const IPvariables_t& vars) {
    // update y/lambda diagonal submatrix
    if( ADD_COMPLIANCE )
        for( auto diag_sel = 0; diag_sel < m; diag_sel++ )
        {
            BigMat.SetElement(n + diag_sel, n + diag_sel,
                              vars.y.GetElement(diag_sel, 0) / vars.lam.GetElement(diag_sel, 0) + E.GetElement(diag_sel, diag_sel));
        }
    else
        for( auto diag_sel = 0; diag_sel < m; diag_sel++ )
        {
            BigMat.SetElement(n + diag_sel, n + diag_sel,
                              vars.y.GetElement(diag_sel, 0) / vars.lam.GetElement(diag_sel, 0));
        }

    // factorize the matrix
    BigMat.Compress();
    mumps_engine.SetMatrix(BigMat);
    for (auto loop_expand_workspace_size = 0; loop_expand_workspace_size < 5; ++loop_expand_workspace_size)
    {
        if (!mumps_engine.MumpsCall(ChMumpsEngine::ANALYZE_FACTORIZE))
            break;
        mumps_engine.SetICNTL(14, static_cast<int>(round(mumps_engine.GetICNTL(14)*1.5)));
        std::cout << "MUMPS work space will be allocated overestimating the estimate by " << mumps_engine.GetICNTL(14) << "%" << std::endl;
    }


}

void ChInteriorPoint::makeNewtonStep(IPvariables_t& Dvar_unknown, ChMatrix<>& rhs, const IPresidual_t& residuals) {
    // Solve the KKT system
    mumps_engine.SetRhsVector(rhs);
    if( mumps_engine.MumpsCall(ChMumpsEngine::SOLVE) )
        mumps_engine.PrintINFOG();
    if( mumps_engine.GetRINFOG(6) > 1e-6 )
        std::cout << "Scaled residual norm of MUMPS call: " << mumps_engine.GetRINFOG(6) << std::endl;

    // MUMPS uses its 'rhs' vector to store the solution. In order to clarify that:
    const ChMatrixDynamic<double>& sol{ rhs };

    // Extract 'x' and 'lam' from 'sol'
    for( auto row_sel = 0; row_sel < n; row_sel++ )
        Dvar_unknown.x.SetElement(row_sel, 0, sol.GetElement(row_sel, 0));
    for( auto row_sel = 0; row_sel < m; row_sel++ )
        Dvar_unknown.lam.SetElement(row_sel, 0, sol.GetElement(row_sel + n, 0));

    // Calc 'y' (it is also possible to evaluate y as (-lam°y+sigma*res.mu*e-y°lam)./lam )
    multiplyA(Dvar_unknown.x, Dvar_unknown.y);  // Dy = A*Dx
    Dvar_unknown.y += residuals.rp;
    if( ADD_COMPLIANCE )
    {
        E.MatrMultiply(Dvar_unknown.lam, vectm);
        Dvar_unknown.y += vectm;
    }
}

void ChInteriorPoint::set_starting_point(IP_STARTING_POINT_METHOD start_point_method, int n_old, int m_old) {
    IPvariables_t Dvar;
    Dvar.x.Resize(n, 1);    Dvar.y.Resize(m, 1);    Dvar.lam.Resize(m, 1);

    ChMatrixDynamic<double> mumps_rhs(n + m, 1);

    switch( start_point_method )
    {
        case IP_STARTING_POINT_METHOD::STP1:
        {
            auto infeas_dual_ratio = 0.1;  // TODO: dependant on n

            var.x.FillElem(1);
            var.y.FillElem(1);
            var.lam.FillElem(1);

            auto duality_gap_calc = var.y.MatrDot(var.y, var.lam);  // [2] pag. 132
            double duality_gap = m;                                   // [2] pag. 132
            assert(duality_gap_calc == duality_gap);

            // norm of all residuals; [2] pag. 132
            residual_fullupdate(res, var);
            auto res_norm = res.rp.MatrDot(res.rp, res.rp);
            res_norm += res.rp.MatrDot(res.rd, res.rd);
            res_norm = sqrt(res_norm);

            if( res_norm / duality_gap > infeas_dual_ratio )
            {
                auto coeff = res_norm / (duality_gap * infeas_dual_ratio);
                var.x.MatrScale(coeff);
                var.y.MatrScale(coeff);
                var.lam.MatrScale(coeff);

                residual_fullupdate(res, var);
            }
        } break;

        case IP_STARTING_POINT_METHOD::STP2:
        {
            double threshold = 1;  // 'epsilon' in [2]

            if( !REUSE_OLD_SOLUTIONS || n != n_old )
            {
                // initialize x
                var.x.FillElem(1);
            }

            // initialize y and then lam
            multiplyA(var.x, vectm);
            vectm -= rhs.b;
            for( auto cont = 0; cont < m; cont++ )
            {
                var.y(cont, 0) = vectm(cont, 0) > threshold ? vectm(cont, 0) : threshold;
                var.lam(cont, 0) = 1 / var.y(cont, 0);
            }

            residual_fullupdate(res, var);

        } break;

        case IP_STARTING_POINT_METHOD::NOCEDAL:
        {
            /********** Initialize IP algorithm **********/
            // Initial guess
            if( n_old != n || solver_call == 0 || !REUSE_OLD_SOLUTIONS )
                var.x.FillElem(1);  // TIP: every ChMatrix is initialized with zeros by default
            if( m_old != m || solver_call == 0 || !REUSE_OLD_SOLUTIONS )
                var.lam.FillElem(1);  // each element of lam will be at the denominator; avoid zeros!

            // since A is generally changed between calls, also with warm_start,
            // all the residuals and feasibility check must be redone
            multiplyA(var.x, var.y);  // y = A*x
            var.y -= rhs.b;

            // Calculate the residual
            residual_fullupdate(res, var);

            // Feasible starting Point (pag.484-485)
            setup_system_matrix(var);

            // fill 'mumps_rhs' with rhs [-res.rd;-res.rp-y]
            for( auto row_sel = 0; row_sel < n; row_sel++ )
                mumps_rhs.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));
            for( auto row_sel = 0; row_sel < m; row_sel++ )
                mumps_rhs.SetElement(row_sel + n, 0, -res.rp(row_sel, 0) - var.y(row_sel, 0));

            makeNewtonStep(Dvar, mumps_rhs, res);

            // x is accepted as it is
            var.y += Dvar.y;      // calculate y0
            var.lam += Dvar.lam;  // calculate lam0

            for( auto row_sel = 0; row_sel < m; row_sel++ )
                var.y(row_sel) = abs(var.y(row_sel)) < 1 ? 1 : abs(var.y(row_sel));

            for( auto row_sel = 0; row_sel < m; row_sel++ )
                var.lam(row_sel) = abs(var.lam(row_sel)) < 1 ? 1 : abs(var.lam(row_sel));

            // Update the residual considering the new values of 'y' and 'lam'
            residual_fullupdate(res, var);

        } break;

        case IP_STARTING_POINT_METHOD::NOCEDAL_WS:
        { /*Backup vectors*/
            ChMatrixDynamic<double> x_bkp(var.x);
            ChMatrixDynamic<double> y_bkp(var.y);
            ChMatrixDynamic<double> lam_bkp(var.lam);
            residual_fullupdate(res, var);
            auto residual_value_bkp = res.rp.NormTwo() + res.rd.NormTwo() + res.mu;

            /********** Initialize IP algorithm **********/
            // Initial guess
            if( n_old != n || solver_call == 0 )
                var.x.FillElem(1);  // TIP: every ChMatrix is initialized with zeros by default
            if( m_old != m || solver_call == 0 )
                var.lam.FillElem(1);  // each element of lam will be at the denominator; avoid zeros!

            // since A is generally changed between calls, also with warm_start,
            // all the residuals and feasibility check must be redone
            multiplyA(var.x, var.y);  // y = A*x
            var.y -= rhs.b;

            // Calculate the residual
            residual_fullupdate(res, var);

            setup_system_matrix(var);

            // Feasible starting Point (pag.484-485)

            // fill 'mumps_rhs' with rhs [-res.rd;-res.rp-y]
            for( auto row_sel = 0; row_sel < n; row_sel++ )
                mumps_rhs.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));
            for( auto row_sel = 0; row_sel < m; row_sel++ )
                mumps_rhs.SetElement(row_sel + n, 0, -res.rp(row_sel, 0) - var.y(row_sel, 0));

            makeNewtonStep(Dvar, mumps_rhs, res);


            // x is accepted as it is
            var.y += Dvar.y;      // calculate y0
            var.lam += Dvar.lam;  // calculate lam0

            for( auto row_sel = 0; row_sel < m; row_sel++ )
                var.y(row_sel) = abs(var.y(row_sel)) < 1 ? 1 : abs(var.y(row_sel));

            for( auto row_sel = 0; row_sel < m; row_sel++ )
                var.lam(row_sel) = abs(var.lam(row_sel)) < 1 ? 1 : abs(var.lam(row_sel));

            // Update the residual considering the new values of 'y' and 'lam'
            residual_fullupdate(res, var);

            /* Check if restoring previous values would be better */
            auto residual_value_new = res.rp.NormTwo() + res.rd.NormTwo() + res.mu;

            if( residual_value_bkp < residual_value_new )
            {
                var.x = x_bkp;
                var.y = y_bkp;
                var.lam = lam_bkp;
                residual_fullupdate(res, var);
            }
            else
            {
                std::cout << "Not WS\n";
            }

        } break;
        default:;
    }
}

/// Find the maximum step length, along the direction defined by \p Dvect, so that \p vect has no negative components;
/// Usually is called passing #var.lam and #var.y and their associated Delta.
double ChInteriorPoint::find_Newton_step_length(const ChMatrix<double>& vect,
                                                const ChMatrix<double>& Dvect,
                                                double tau) { //TODO: is 'tau' to be multiplied before or after checking for best alpha? if the latter can be done after calling the function
#if 0
    double alpha = 1;  // in this way alpha is clamped to a maximum of 1
    for( auto row_sel = 0; row_sel < vect.GetRows(); row_sel++ )
    {
        if( Dvect(row_sel, 0) < 0 )
        {
            auto alfa_temp = -tau * vect(row_sel, 0) / Dvect(row_sel, 0);
            if( alfa_temp < alpha && alfa_temp > 0 )
            {
                alpha = alfa_temp;
            }
        }
    }
    return alpha;

#else
    double alpha_inv = 1;
    for( auto row_sel = 0; row_sel < vect.GetRows(); ++row_sel )
    {
        auto alpha_inv_temp = -Dvect(row_sel, 0) / vect(row_sel, 0);
        if( alpha_inv_temp > alpha_inv )
        {
            alpha_inv = alpha_inv_temp;
        }
    }

    return tau / alpha_inv;
#endif


}

void ChInteriorPoint::find_Newton_step_length(const IPvariables_t& vars, const IPvariables_t& Dvars, double tau, double& alfa_prim, double& alfa_dual) const {

    //ExportArrayToFile(vars.x, "dump/x.txt");
    //ExportArrayToFile(vars.y, "dump/y.txt");
    //ExportArrayToFile(vars.lam, "dump/lam.txt");

    //ExportArrayToFile(Dvars.x, "dump/Dx.txt");
    //ExportArrayToFile(Dvars.y, "dump/Dy.txt");
    //ExportArrayToFile(Dvars.lam, "dump/Dlam.txt");


    // maximize the steplength the classical way (just for debug)
    auto alfa_prim_standard = find_Newton_step_length(vars.y, Dvars.y, 1);
    auto alfa_dual_standard = find_Newton_step_length(vars.lam, Dvars.lam, 1);

    // minimize primal residual (just for fun)
    ChMatrixDynamic<double> qp1(m, 1);
    multiplyA(Dvars.x, qp1);
    qp1 -= Dvars.y;
    auto alfa_prim_test = -qp1.MatrDot(res.rp, qp1) / qp1.MatrDot(qp1, qp1);

    // minimize dual residual
    ChMatrixDynamic<double> qd1(n, 1);
    ChMatrixDynamic<double> qd2(n, 1);

    multiplyG(Dvars.x, qd1);
    multiplyNegAT(Dvars.lam, qd2);

    ChMatrixNM<double, 2, 1> known_vect;
    known_vect[0][0] = -qd1.MatrDot(res.rd, qd1);
    known_vect[1][0] = -qd2.MatrDot(res.rd, qd2);

    ChMatrixNM<double, 2, 2> half_hessian_inverted;
    half_hessian_inverted[0][0] = qd2.MatrDot(qd2, qd2);
    half_hessian_inverted[0][1] = -qd1.MatrDot(qd1, qd2);
    half_hessian_inverted[1][0] = half_hessian_inverted[0][1];
    half_hessian_inverted[1][1] = qd1.MatrDot(qd1, qd1);
    half_hessian_inverted *= 1 / (half_hessian_inverted[1][1] * half_hessian_inverted[0][0] - half_hessian_inverted[0][1] * half_hessian_inverted[1][0]);

    ChMatrixNM<double, 2, 1> alfa_sol;
    alfa_sol.MatrMultiply(half_hessian_inverted, known_vect);

    switch( STEPLENGTH_METHOD )
    {
        case 0: // simple "damping"
        {
            alfa_prim = alfa_prim_standard;
            alfa_dual = alfa_dual_standard;
        }
        break;
        case 1: // "method 2" in [4]
        {
            alfa_prim = std::max(std::min(alfa_prim_standard, alfa_dual_standard), alfa_sol[0][0]);
            alfa_dual = std::max(0.0, std::min(alfa_sol[1][0], 1.0)); // simple clamping from 0 to 1
        }
        break;
        case 2:
        {
            alfa_prim = std::max(std::min(alfa_prim_standard, alfa_dual_standard), alfa_sol[0][0]);
            alfa_dual = std::max(alfa_dual_standard, std::min(alfa_sol[1][0], 1.0));
        }
        break;
        default:
        assert(0);
        break;
    }
    alfa_prim *= tau;
    alfa_dual *= tau;

    // obtain eigenvalues (just for fun)
    auto eig_max = half_hessian_inverted[1][1] + half_hessian_inverted[0][0] + sqrt(std::pow(half_hessian_inverted[1][1] - half_hessian_inverted[0][0], 2) + 4 * std::pow(half_hessian_inverted[0][1], 2));
    auto eig_min = half_hessian_inverted[1][1] + half_hessian_inverted[0][0] - sqrt(std::pow(half_hessian_inverted[1][1] - half_hessian_inverted[0][0], 2) + 4 * std::pow(half_hessian_inverted[0][1], 2));

    if( eig_max >= 0 && eig_min >= 0 )
    {
        // alfa_sol is a minimum for the dual residual norm
    }

}

double ChInteriorPoint::evaluate_objective_function() const {
    multiplyG(var.x, vectn);
    auto obj_value = vectn.MatrDot(var.x, vectn);
    obj_value += rhs.c.MatrDot(var.x, rhs.c);

    return obj_value;
}


/// Take care of adapting the size of matrices to \p n_new and \p m_new
void ChInteriorPoint::reset_internal_dimensions(int n_new, int m_new) {
    if( n != n_new )
    {
        var.x.Resize(n_new, 1);
        rhs.c.Resize(n_new, 1);
        res.rd.Resize(n_new, 1);
        vectn.Resize(n_new, 1);
    }

    if( m != m_new )
    {
        var.y.Resize(m_new, 1);
        var.lam.Resize(m_new, 1);
        rhs.b.Resize(m_new, 1);
        res.rp.Resize(m_new, 1);
        res.rpd.Resize(m_new, 1);
        vectm.Resize(m_new, 1);
    }

    SKIP_CONTACTS_UV ? sol_chrono.Resize(n_new + 3 * m_new, 1) : sol_chrono.Resize(n_new + m_new, 1);

    // BigMat and sol
    switch( KKT_solve_method )
    {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
        BigMat.Reset(2 * m_new + n_new, 2 * m_new + n_new, static_cast<int>(n_new * n_new * SPM_DEF_FULLNESS));
        break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
        BigMat.Reset(n_new + m_new, n_new + m_new, static_cast<int>(n_new * n_new * SPM_DEF_FULLNESS));
        break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
        std::cout << std::endl << "Perturbed KKT system cannot be stored with 'NORMAL' method yet.";
        break;
    }

    n = n_new;
    m = m_new;
}

void ChInteriorPoint::DumpProblem(std::string suffix) {
    CreateDirectory("dump", nullptr);
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

}

void ChInteriorPoint::make_positive_definite() {
    if( m == 0 )
        return;

    int offset_AT_col = n;
    if( KKT_solve_method == IP_KKT_SOLUTION_METHOD::STANDARD )
        offset_AT_col = n + m;

    BigMat.ForEachExistentValueInRange([](double* val) { *val *= -1; }, 0, n, offset_AT_col,
                                       BigMat.GetNumColumns() - 1);
}

// Perform moltiplication of A with vect_in: vect_out = A*vect_in
void ChInteriorPoint::multiplyA(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
    switch( KKT_solve_method )
    {
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
    switch( KKT_solve_method )
    {
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
    switch( KKT_solve_method )
    {
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


void ChInteriorPoint::residual_fullupdate(IPresidual_t& residuals, const IPvariables_t& variables) const {
    // Dual Residual
    // res.rd = G*x + c - A^T*lam
    multiplyG(variables.x, residuals.rd);  // res.rd = G*x
    residuals.rd += rhs.c;           // res.rd = G*x + c

    if( m > 0 )
    {
        multiplyNegAT(variables.lam, vectn);  // vectn = (-A^T)*lam
        residuals.rd += vectn;                // res.rd = (G*x + c) + (-A^T*lam)

        // Primal residual
        // res.rp = A*x - y - b
        multiplyA(variables.x, residuals.rp);  // res.rp = A*x
        residuals.rp -= variables.y;
        residuals.rp -= rhs.b;
        if( ADD_COMPLIANCE )
        {
            E.MatrMultiply(variables.lam, vectm);
            residuals.rp += vectm;
        }

        residuals.mu = variables.y.MatrDot(variables.y, variables.lam) / m;
    }
    else
    {
        residuals.rp.FillElem(0);
        residuals.mu = 0;
    }
}

/// Export the IP variables in Chrono format
// TODO: FromVectorToUnknowns should accept const reference, but it doesn't. When it will be, we could pass const ChMatrix<>& in the argument and as return
ChMatrix<>& ChInteriorPoint::adapt_to_Chrono(ChMatrix<>& solution_vect) const {

    // copy 'x'
    for( auto row_sel = 0; row_sel < n; row_sel++ )
        solution_vect(row_sel, 0) = var.x(row_sel, 0);

    // copy Lagrangian multipliers; skip tangential forces if needed
    if( SKIP_CONTACTS_UV )
    {
        for( auto row_sel = 0; row_sel < m; row_sel++ )
        {
            solution_vect(n + row_sel * 3, 0) = -var.lam(row_sel, 0);  // there will be an inversion inside FromVectorToUnknowns()
            solution_vect(n + row_sel * 3 + 1, 0) = 0;
            solution_vect(n + row_sel * 3 + 2, 0) = 0;
        }
    }
    else
    {
        for( auto row_sel = 0; row_sel < m; row_sel++ )
            solution_vect(row_sel + n, 0) =
            -var.lam(row_sel);  // there will be an inversion inside FromVectorToUnknowns()
    }

    return solution_vect;
}

void ChInteriorPoint::RecordHistory(bool on_off, std::string file_name) {
    print_history = on_off;

    if( print_history )
        logfile_name = file_name;
    else
    {
        logfile_stream.close();
        return;
    }

}

void ChInteriorPoint::Solve(const ChCOOMatrix& augmented_mat, const ChMatrix<double>& rhs_b, const ChMatrix<double>& rhs_c, ChMatrixDynamic<double>& var_x, ChMatrixDynamic<double>& var_y, ChMatrixDynamic<double>& var_lam) {
    // Guys, I know that it's the worst way to to it, but it's just for debug!
    // I promise that, if you want this feature to be left, I will change the settings in order to make moves instead of copies.
    BigMat = augmented_mat;
    rhs.b = rhs_b;
    rhs.c = rhs_c;

    n = rhs.c.GetRows();
    m = rhs.b.GetRows();

    var_x.Resize(n, 1);
    var_y.Resize(m, 1);
    var_lam.Resize(m, 1);

    reset_internal_dimensions(n, m);
    auto bkp_method = KKT_solve_method;

    KKT_solve_method = IP_KKT_SOLUTION_METHOD::AUGMENTED;


    /********** Check if system has constraints **********/
    if( m == 0 )  // if no constraints
    {
        ChMatrixDynamic<double> mumps_rhs_sol(n, 1);

        for( auto row_sel = 0; row_sel < n; row_sel++ )
            mumps_rhs_sol.SetElement(row_sel, 0, rhs.c.GetElement(row_sel, 0));

        // Solve the KKT system
        BigMat.Compress();
        mumps_engine.SetProblem(BigMat, mumps_rhs_sol);
        if( mumps_engine.MumpsCall(ChMumpsEngine::COMPLETE) )
            mumps_engine.PrintINFOG();

        if( verbose && mumps_engine.GetRINFOG(6) > 1e-6 )
            std::cout << "MUMPS scaled residual: " << mumps_engine.GetRINFOG(6) << std::endl;        

        if( verbose )
            std::cout << "IP called with no constraints. Switched to linear solver." << std::endl;

        return;
    }

    /********* The system DOES have constraints! Start Interior Point ********/

    set_starting_point(starting_point_method);

    for( iteration_count = 1; iteration_count < iteration_count_max; iteration_count++ )
    {

        iterate();

        /*********************************************************************************/
        /******************************** Exit conditions ********************************/
        /*********************************************************************************/

        for( auto cont = 0; cont < m; cont++ )
        {
            if( var.y(cont, 0) < 0 || var.lam(cont, 0) < 0 )
            {
                std::cout << "'y' or 'lambda' have negative elements" << std::endl;
                break;
            }
        }


        if( verbose )
        {
            GetLog() << "InteriorPoint | Call: " << solver_call << " Iter: " << iteration_count << "/" << iteration_count_max << "\n";
            GetLog() << "Complementarity Measure: " << res.mu << "\n";
            GetLog() << "|rd|/n (stationarity): " << res.rd.NormTwo() / n << "\n";
            GetLog() << "|rp|/m (constraint violation): " << res.rp.NormTwo() / m << "\n";
            GetLog() << "Objective Function: " << evaluate_objective_function() << "\n";
            GetLog() << "\n";
        }

        //DumpProblem("_end");

        if( res < res_nnorm_tol )
        {
            if( verbose )
                std::cout << "IP  iter: " << iteration_count << "/" << iteration_count_max
                << std::endl;


            break;
        }
    }


    // Scatter the solution
    var_x = var.x;
    var_x = var.y;
    var_x = var.lam;

    KKT_solve_method = bkp_method;

}

ChInteriorPoint::ChInteriorPoint() {
    mumps_engine.SetICNTL(11, 2);
    RecordHistory(true);
}

ChInteriorPoint::~ChInteriorPoint() {
    if( logfile_stream.is_open() )
        logfile_stream.close();
}
}
