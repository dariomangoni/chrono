#include "ChInteriorPoint.h"
#include <algorithm>
#include "solver/ChConstraintTwoTuplesFrictionT.h"

#ifdef CHRONO_POSTPROCESS
#include "chrono_postprocess/ChGnuPlot.h"
#endif

//#define DEBUG_MODE
#define SKIP_CONTACTS_UV true
#define ADD_COMPLIANCE false
#define REUSE_OLD_SOLUTIONS false
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
            printf("%3.1f ", matrice.GetElement(i, j));
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
	// for optimization of the starting point
    auto n_old = n;
    auto m_ineq_old = m_ineq;
    auto m_eq_old = m_eq;

	sysd.SortActiveConstraints(); // tell ChSystemDescriptor to write the Jacobian [equalities;inequalities]

	auto m_tuple = sysd.CountActiveConstraintsSorted(SKIP_CONTACTS_UV);
    reset_internal_dimensions(sysd.CountActiveVariables(), std::get<0>(m_tuple), std::get<2>(m_tuple), std::get<1>(m_tuple));
    sysd.ConvertToMatrixForm(&BigMat, nullptr, SKIP_CONTACTS_UV, false);

    make_positive_definite(); //TODO: make it an option


    sysd.ConvertToMatrixForm(nullptr, nullptr, nullptr, &rhs.c, &rhs.b, nullptr, false, SKIP_CONTACTS_UV);  // load f->c and b->b
    rhs.c.MatrScale(-1); // adapt to InteriorPoint convention
    rhs.b.MatrScale(-1); // adapt to InteriorPoint convention


#ifdef BYPASS_RESTITUTION_TERM
    rhs.b.FillElem(0); // TODO: WARNING: removing phi/h!!!
    //rhs.b*=1e-3; // TODO: WARNING: removing phi/h!!!
#endif

    /********** Check if system has inequality constraints **********/
    if( m_ineq == 0 )  // if no inequality constraints
    {
        ChMatrixDynamic<double> mumps_rhs(n + m_eq, 1);

        // Fill 'mumps_rhs' with just Chrono's 'f' i.e. IP's '-c'
        for( auto row_sel = 0; row_sel < n; row_sel++ )
            mumps_rhs.SetElement(row_sel, 0, -rhs.c.GetElement(row_sel, 0));
		for (auto row_sel = 0; row_sel < m_eq; row_sel++)
			mumps_rhs.SetElement( n + row_sel, 0, rhs.b.GetElement(row_sel, 0));

        // Solve the KKT system
        BigMat.Compress();
        mumps_engine.SetProblem(BigMat, mumps_rhs);
        if( mumps_engine.MumpsCall(ChMumpsEngine::COMPLETE) )
            mumps_engine.PrintINFOG();

        if( verbose && mumps_engine.GetRINFOG(6) > 1e-6 )
            std::cout << "MUMPS scaled residual: " << mumps_engine.GetRINFOG(6) << std::endl;

        sysd.FromVectorToUnknowns(mumps_rhs);

        // Export variable so that can be used in the next iteration as starting point
        var.v.Resize(n, 1);
        var.gamma.Resize(m_eq, 1);
		for (auto row_sel = 0; row_sel < n; row_sel++)
			var.v.SetElement(row_sel, 0, mumps_rhs.GetElement(row_sel, 0));
		for (auto row_sel = 0; row_sel < m_eq; row_sel++)
			var.gamma.SetElement(row_sel, 0, mumps_rhs.GetElement( n + row_sel, 0));

        if( verbose )
            std::cout << "IP call: " << solver_call << "; No inequality constraints." << std::endl;

        return 0.0;
    }

    /********* The system DOES have constraints! Start Interior Point ********/
	if (m_eq > 0)
	{
		rhs.b_eq.PasteClippedMatrix(rhs.b, 0, 0, m_eq, 1, 0, 0);
		rhs.b_ineq.PasteClippedMatrix(rhs.b, m_eq, 0, m_ineq, 1, 0, 0);
	}
	else
		rhs.b_ineq = rhs.b;


    ip_solver_call++;
    ip_timer.start();


    if( ADD_COMPLIANCE && m_ineq > 0 )
    {
        sysd.ConvertToMatrixForm(nullptr, nullptr, &E, nullptr, nullptr, nullptr, false, SKIP_CONTACTS_UV);
        E *= -1;
    }

    set_starting_point(starting_point_method, n_old, m_eq_old, m_ineq_old);

    for( iteration_count = 1; iteration_count < iteration_count_max; iteration_count++ )
    {

        iterate();
        iteration_count_tot++;

        /*********************************************************************************/
        /******************************** Exit conditions ********************************/
        /*********************************************************************************/


        if( verbose )
        {
			PrintIPStatus();
			residual_fullupdate(res, var);
			PrintIPStatus();
        }

        //DumpProblem("_end");

        if( res <= res_nnorm_tol )
        {
            if( verbose )
                std::cout << "IP call: " << solver_call << "; iter: " << iteration_count << "/" << iteration_count_max
                << std::endl;


            break;
        }
    }

    ip_timer.stop();

    // Scatter the solution into the Chrono environment
    sysd.FromVectorToUnknowns(adapt_to_Chrono(sysd, sol_chrono));

    return 0.0;
}

// Iterating function
void ChInteriorPoint::iterate() {
    /*********************************************************************************/
    /***************************** Prediction Phase **********************************/
    /*********************************************************************************/

    setup_system_matrix(var);

    ChMatrixDynamic<double> mumps_rhs(n + m_eq + m_ineq, 1);
    IPvariables_t Dvar_pred(n, m_eq, m_ineq);

    // WARNING: the residual structure 'res' must be already updated at this point!

    // fill 'mumps_rhs' with rhs [-res.rd;-res.rp_gamma;-res.rp_lambda-y]
    for( auto row_sel = 0; row_sel < n;      row_sel++ ) mumps_rhs.SetElement(row_sel,            0, -res.rd(row_sel, 0));
	for (auto row_sel = 0; row_sel < m_eq;   row_sel++ ) mumps_rhs.SetElement(row_sel + n,        0, -res.rp_gamma(row_sel, 0));
    for( auto row_sel = 0; row_sel < m_ineq; row_sel++ ) mumps_rhs.SetElement(row_sel + n + m_eq, 0, -res.rp_lambda(row_sel, 0) - var.y(row_sel, 0));

    makeNewtonStep(Dvar_pred, mumps_rhs, res);

    /*** compute step lengths ***/
    // from 16.60 pag.482 from 14.32 pag.408 (remember that y>=0!)
    auto alfa_pred_prim = find_Newton_step_length(var.y, Dvar_pred.y);
    auto alfa_pred_dual = find_Newton_step_length(var.lambda, Dvar_pred.lambda);

    if( EQUAL_STEP_LENGTH )
    {
        auto alfa_pred = std::min(alfa_pred_prim, alfa_pred_dual);
        alfa_pred_prim = alfa_pred;
        alfa_pred_dual = alfa_pred;
    }


    /*** make the prediction step ***/
	// update only variables needed for complementarity measure
    IPvariables_t var_pred(n, m_eq, m_ineq);

    var_pred.y = Dvar_pred.y;
	var_pred.lambda = Dvar_pred.lambda;

	var_pred.y.MatrScale(alfa_pred_prim);
	var_pred.lambda.MatrScale(alfa_pred_dual);

	var_pred.y += var.y;
	var_pred.lambda += var.lambda;

    /*** compute complementarity measure ***/
    auto mu_pred = var_pred.y.MatrDot(var_pred.y, var_pred.lambda) / m_ineq;  // from 16.56 pag.481

	if (mu_pred == 0)
		std::cout << "IP warning: complementarity measure is 0, but probably we are not at the optimal solution." << std::endl;

    if( ONLY_PREDICT )
    {
		// update remaining variables (v, gamma)
        var_pred.v = Dvar_pred.v;
    	var_pred.v.MatrScale(alfa_pred_prim);
    	var_pred.v += var.v;

		var_pred.gamma = Dvar_pred.gamma;
		var_pred.gamma.MatrScale(alfa_pred_dual); //TODO: check if dual or primal
		var_pred.gamma += var.gamma;

		// store in global variables
    	var.v = var_pred.v;
        var.y = var_pred.y;
        var.gamma = var_pred.gamma;
    	var.lambda = var_pred.lambda;

        // update residuals
        res.rp_lambda.MatrScale(1 - alfa_pred_prim);
        res.rp_gamma.MatrScale(1 - alfa_pred_prim);

		res.rd.MatrScale(1 - alfa_pred_dual);
		if (!EQUAL_STEP_LENGTH)
		{
			multiplyH(Dvar_pred.v, vectn);                          // vectn = H * Dvar.v
			vectn.MatrScale(alfa_pred_prim - alfa_pred_dual);  // vectn = (alfa_pred_prim - alfa_pred_dual) * (H * Dvar.v)
			res.rd += vectn;
		}

        res.mu = mu_pred;

        return;
    }

    /*********************************************************************************/
    /******************************* Correction phase ********************************/
    /*********************************************************************************/
    IPvariables_t Dvar(n, m_eq, m_ineq);

    /*** evaluate centering parameter ***/
    auto sigma = std::pow( mu_pred / res.mu, 3);  // from 14.34 pag.408

    /*** find directions ***/
	for (auto row_sel = 0; row_sel < n;      row_sel++) mumps_rhs.SetElement(row_sel,            0, -res.rd(row_sel, 0));
	for (auto row_sel = 0; row_sel < m_eq;   row_sel++) mumps_rhs.SetElement(row_sel + n,        0, -res.rp_gamma(row_sel, 0));
	for (auto row_sel = 0; row_sel < m_ineq; row_sel++) mumps_rhs.SetElement(row_sel + n + m_eq, 0, -res.rp_lambda(row_sel, 0) - var.y(row_sel, 0) + (sigma * res.mu - Dvar_pred.lambda(row_sel, 0) * Dvar_pred.y(row_sel, 0)) / var.lambda(row_sel, 0));

    makeNewtonStep(Dvar, mumps_rhs, res);

    /*** step length correction ***/
    auto tau = ADAPTIVE_ETA ? exp(-res.mu * m_ineq) * 0.1 + 0.9 : 0.95;  // exponential descent of tau

    /*** compute step lengths ***/
    auto alfa_corr_prim = find_Newton_step_length(var.y, Dvar.y, tau);
    auto alfa_corr_dual = find_Newton_step_length(var.lambda, Dvar.lambda, tau);

    if( EQUAL_STEP_LENGTH )
    {
        auto alfa_corr = std::min(alfa_corr_prim, alfa_corr_dual);
        alfa_corr_prim = alfa_corr;
        alfa_corr_dual = alfa_corr;
    }

    IPvariables_t var_corr(n, m_eq, m_ineq);

    /*** make the correction step ***/
    var_corr.v = Dvar.v;        var_corr.v.MatrScale(alfa_corr_prim);        var_corr.v += var.v;        var.v = var_corr.v;
    var_corr.y = Dvar.y;        var_corr.y.MatrScale(alfa_corr_prim);        var_corr.y += var.y;        var.y = var_corr.y;
    var_corr.gamma = Dvar.gamma;      var_corr.gamma.MatrScale(alfa_corr_dual);     var_corr.gamma += var.gamma;        var.gamma = var_corr.gamma;
    var_corr.lambda = Dvar.lambda;    var_corr.lambda.MatrScale(alfa_corr_dual);    var_corr.lambda += var.lambda;    var.lambda = var_corr.lambda;

    /********** Residuals update **********/
	res.mu = var.y.MatrDot(var.y, var.lambda) / m_ineq;  // from 14.6 pag.395
    res.rp_gamma.MatrScale(1 - alfa_corr_prim);
    res.rp_lambda.MatrScale(1 - alfa_corr_prim);

	res.rd.MatrScale(1 - alfa_corr_dual);
    if( !EQUAL_STEP_LENGTH )
    {
        multiplyH(Dvar.v, vectn);                          // vectn = H*Dvar.v
        vectn.MatrScale(alfa_corr_prim - alfa_corr_dual);  // vectn = (alfa_pred_prim - alfa_pred_dual) * (H * Dvar.v)
        res.rd += vectn;
    }


	/*********************************************************************************/
	/********************************* Data Logging **********************************/
	/*********************************************************************************/
    if( print_history )
    {
        if( logfile_stream.is_open() )
        {
            logfile_stream << std::endl << res.rp_lambda.NormTwo() / m_ineq << ", " << res.rd.NormTwo() / n << ", " << res.mu << ", " << evaluate_objective_function();
        }
        else
        {
            logfile_stream.open(logfile_name + ".txt");

            if( logfile_stream.is_open() )
            {
                logfile_stream << std::scientific << std::setprecision(6);
                logfile_stream << "res_P res_D mu objfun";

                logfile_stream << std::endl << res.rp_lambda.NormTwo() / m_ineq << ", " << res.rd.NormTwo() / n << ", " << res.mu << ", " << evaluate_objective_function();

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
		for (auto diag_sel = 0; diag_sel < m_ineq; ++diag_sel)
		{
			BigMat.SetElement(n + m_eq + diag_sel, n + m_eq + diag_sel, vars.y.GetElement(diag_sel, 0) / vars.lambda.GetElement(diag_sel, 0) + E.GetElement(diag_sel - n + m_eq, diag_sel - n + m_eq));
		}
    else
        for( auto diag_sel = 0; diag_sel < m_ineq; ++diag_sel )
        {
            BigMat.SetElement(n + m_eq + diag_sel, n + m_eq + diag_sel, vars.y.GetElement(diag_sel, 0) / vars.lambda.GetElement(diag_sel, 0));
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
    const ChMatrixDynamic<double>& sol = rhs;

    // Extract 'v' and 'lambda' from 'sol'
    for( auto row_sel = 0; row_sel < n;      row_sel++ ) Dvar_unknown.v.SetElement(row_sel, 0, sol.GetElement(row_sel, 0));
	for (auto row_sel = 0; row_sel < m_eq;   row_sel++ ) Dvar_unknown.gamma.SetElement(row_sel, 0, sol.GetElement(row_sel + n, 0));
    for( auto row_sel = 0; row_sel < m_ineq; row_sel++ ) Dvar_unknown.lambda.SetElement(row_sel, 0, sol.GetElement(row_sel + n + m_eq, 0));

    // Calc 'y' (it is also possible to evaluate y as (-lambda°y+sigma*res.mu*e-y°lambda)./lambda )
    multiplyIneqCon(Dvar_unknown.v, Dvar_unknown.y);  // Dy = IneqCon*Dv
    Dvar_unknown.y += residuals.rp_lambda; // Dy = (IneqCon*Dv) + rp_lambda
    if( ADD_COMPLIANCE )
    {
        E.MatrMultiply(Dvar_unknown.lambda, vectm_ineq);
        Dvar_unknown.y += vectm_ineq;
    }
}

void ChInteriorPoint::set_starting_point(IP_STARTING_POINT_METHOD start_point_method, int n_old, int m_eq_old, int m_ineq_old) {
    IPvariables_t Dvar(n, m_eq ,m_ineq);

    ChMatrixDynamic<double> mumps_rhs(n + m_eq + m_ineq, 1);

    switch( start_point_method )
    {
        case IP_STARTING_POINT_METHOD::STP1:
        {
            auto infeas_dual_ratio = 0.1;  // TODO: dependant on n

            var.v.FillElem(1);
            var.y.FillElem(1);
            var.gamma.FillElem(1);
            var.lambda.FillElem(1);

            auto duality_gap_calc = var.y.MatrDot(var.y, var.lambda);  // [2] pag. 132
            double duality_gap = m_ineq;                                   // [2] pag. 132
            assert(duality_gap_calc == duality_gap);

            // norm of all residuals; [2] pag. 132
            residual_fullupdate(res, var);
            auto res_norm = res.rp_lambda.MatrDot(res.rp_lambda, res.rp_lambda);
            res_norm += res.rp_lambda.MatrDot(res.rd, res.rd);
            res_norm = sqrt(res_norm);

            if( res_norm / duality_gap > infeas_dual_ratio )
            {
                auto coeff = res_norm / (duality_gap * infeas_dual_ratio);
                var.v.MatrScale(coeff);
                var.y.MatrScale(coeff);
                var.lambda.MatrScale(coeff);

                residual_fullupdate(res, var);
            }
        } break;

        case IP_STARTING_POINT_METHOD::STP2:
        {
            double threshold = 1;  // 'epsilon' in [2]

            if( !REUSE_OLD_SOLUTIONS || n != n_old )
            {
                // initialize v
                var.v.FillElem(1);
            }

            // initialize y and then lambda
            multiplyIneqCon(var.v, vectm_ineq);
            vectm_ineq -= rhs.b;
            for( auto cont = 0; cont < m_ineq; cont++ )
            {
                var.y(cont, 0) = vectm_ineq(cont, 0) > threshold ? vectm_ineq(cont, 0) : threshold;
                var.lambda(cont, 0) = 1 / var.y(cont, 0);
            }

            residual_fullupdate(res, var);

        } break;

        case IP_STARTING_POINT_METHOD::NOCEDAL:
        {
            /********** Initialize IP algorithm **********/
            // Initial guess
            if( n_old != n || solver_call == 0 || !REUSE_OLD_SOLUTIONS )
                var.v.FillElem(1);  // TIP: every ChMatrix is initialized with zeros by default
        	if( m_eq_old != m_eq || solver_call == 0 || !REUSE_OLD_SOLUTIONS )
                var.gamma.FillElem(1);  // TIP: every ChMatrix is initialized with zeros by default
            if( m_ineq_old != m_ineq || solver_call == 0 || !REUSE_OLD_SOLUTIONS )
                var.lambda.FillElem(1);  // each element of lambda will be at the denominator; avoid zeros!

            // since IneqCon is generally changed between calls, also with warm_start,
            // all the residuals and feasibility check must be redone
            multiplyIneqCon(var.v, var.y);
            var.y -= rhs.b_ineq;

            // Calculate the residual
            residual_fullupdate(res, var);

            // Feasible starting Point (pag.484-485)
            setup_system_matrix(var);

			// fill 'mumps_rhs' with rhs [-res.rd;-res.rp_gamma;-res.rp_lambda-y] like in prediction step
			for (auto row_sel = 0; row_sel < n; row_sel++) mumps_rhs.SetElement(row_sel, 0, -res.rd(row_sel, 0));
			for (auto row_sel = 0; row_sel < m_eq; row_sel++) mumps_rhs.SetElement(row_sel + n, 0, -res.rp_gamma(row_sel, 0));
			for (auto row_sel = 0; row_sel < m_ineq; row_sel++) mumps_rhs.SetElement(row_sel + n + m_eq, 0, -res.rp_lambda(row_sel, 0) - var.y(row_sel, 0));

            makeNewtonStep(Dvar, mumps_rhs, res);

            // v is accepted as it is
            var.y += Dvar.y;      // calculate y0
            var.gamma += Dvar.gamma;      // calculate gamma0
            var.lambda += Dvar.lambda;  // calculate lam0

            for( auto row_sel = 0; row_sel < m_ineq; row_sel++ )
                var.y(row_sel) = abs(var.y(row_sel)) < 1 ? 1 : abs(var.y(row_sel));

            for( auto row_sel = 0; row_sel < m_ineq; row_sel++ )
                var.lambda(row_sel) = abs(var.lambda(row_sel)) < 1 ? 1 : abs(var.lambda(row_sel));

            // Update the residual considering the new values of 'y', 'gamma' and 'lambda'
            residual_fullupdate(res, var);

        } break;

        case IP_STARTING_POINT_METHOD::NOCEDAL_WS:
        { /*Backup vectors*/
            ChMatrixDynamic<double> x_bkp(var.v);
            ChMatrixDynamic<double> y_bkp(var.y);
            ChMatrixDynamic<double> lam_bkp(var.lambda);
            residual_fullupdate(res, var);
            auto residual_value_bkp = res.rp_lambda.NormTwo() + res.rd.NormTwo() + res.mu;

            /********** Initialize IP algorithm **********/
            // Initial guess
            if( n_old != n || solver_call == 0 )
                var.v.FillElem(1);  // TIP: every ChMatrix is initialized with zeros by default
            if( m_ineq_old != m_ineq || solver_call == 0 )
                var.lambda.FillElem(1);  // each element of lambda will be at the denominator; avoid zeros!

            // since IneqCon is generally changed between calls, also with warm_start,
            // all the residuals and feasibility check must be redone
            multiplyIneqCon(var.v, var.y);  // y = IneqCon*v
            var.y -= rhs.b;

            // Calculate the residual
            residual_fullupdate(res, var);

            setup_system_matrix(var);

            // Feasible starting Point (pag.484-485)

            // fill 'mumps_rhs' with rhs [-res.rd;-res.rp_lambda-y]
            for( auto row_sel = 0; row_sel < n; row_sel++ )
                mumps_rhs.SetElement(row_sel, 0, -res.rd.GetElement(row_sel, 0));
            for( auto row_sel = 0; row_sel < m_ineq; row_sel++ )
                mumps_rhs.SetElement(row_sel + n, 0, -res.rp_lambda(row_sel, 0) - var.y(row_sel, 0));

            makeNewtonStep(Dvar, mumps_rhs, res);


            // v is accepted as it is
            var.y += Dvar.y;      // calculate y0
            var.lambda += Dvar.lambda;  // calculate lam0

            for( auto row_sel = 0; row_sel < m_ineq; row_sel++ )
                var.y(row_sel) = abs(var.y(row_sel)) < 1 ? 1 : abs(var.y(row_sel));

            for( auto row_sel = 0; row_sel < m_ineq; row_sel++ )
                var.lambda(row_sel) = abs(var.lambda(row_sel)) < 1 ? 1 : abs(var.lambda(row_sel));

            // Update the residual considering the new values of 'y' and 'lambda'
            residual_fullupdate(res, var);

            /* Check if restoring previous values would be better */
            auto residual_value_new = res.rp_lambda.NormTwo() + res.rd.NormTwo() + res.mu;

            if( residual_value_bkp < residual_value_new )
            {
                var.v = x_bkp;
                var.y = y_bkp;
                var.lambda = lam_bkp;
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
/// Usually is called passing #var.lambda and #var.y and their associated Delta.
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

    //ExportArrayToFile(vars.v, "dump/v.txt");
    //ExportArrayToFile(vars.y, "dump/y.txt");
    //ExportArrayToFile(vars.lambda, "dump/lambda.txt");

    //ExportArrayToFile(Dvars.v, "dump/Dv.txt");
    //ExportArrayToFile(Dvars.y, "dump/Dy.txt");
    //ExportArrayToFile(Dvars.lambda, "dump/Dlam.txt");


    // maximize the steplength the classical way (just for debug)
    auto alfa_prim_standard = find_Newton_step_length(vars.y, Dvars.y, 1);
    auto alfa_dual_standard = find_Newton_step_length(vars.lambda, Dvars.lambda, 1);

    // minimize primal residual (just for fun)
    ChMatrixDynamic<double> qp1(m_ineq, 1);
    multiplyIneqCon(Dvars.v, qp1);
    qp1 -= Dvars.y;
    auto alfa_prim_test = -qp1.MatrDot(res.rp_lambda, qp1) / qp1.MatrDot(qp1, qp1);

    // minimize dual residual
    ChMatrixDynamic<double> qd1(n, 1);
    ChMatrixDynamic<double> qd2(n, 1);

    multiplyH(Dvars.v, qd1);
    multiplyNegIneqConT(Dvars.lambda, qd2);

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
    multiplyH(var.v, vectn);
    auto obj_value = vectn.MatrDot(var.v, vectn);
    obj_value += rhs.c.MatrDot(var.v, rhs.c);

    return obj_value;
}


/// Take care of adapting the size of matrices to \p n_new and \p m_new
void ChInteriorPoint::reset_internal_dimensions(int n_new, int m_eq_new, int m_ineq_new, int m_ineq_full) {

	var.Resize(n_new, m_eq_new, m_ineq_new);
	rhs.Resize(n_new, m_eq_new, m_ineq_new);
	res.Resize(n_new, m_eq_new, m_ineq_new);

	// update mutables
    vectn.Resize(n_new, 1);
	vectm_eq.Resize(m_eq_new, 1);
    vectm_ineq.Resize(m_ineq_new, 1);
    sol_chrono.Resize(n_new + m_eq_new + m_ineq_full, 1);

    // BigMat and sol
	switch (KKT_solve_method)
	{
	case IP_KKT_SOLUTION_METHOD::STANDARD:
		std::cout << std::endl << "Perturbed KKT system cannot be stored with 'STANDARD	' method yet.";
		break;
	case IP_KKT_SOLUTION_METHOD::AUGMENTED:
		BigMat.Reset(n_new + m_eq_new + m_ineq_new, n_new + m_eq_new + m_ineq_new, static_cast<int>((n_new + m_eq_new + m_ineq_new) * (n_new + m_eq_new + m_ineq_new) * SPM_DEF_FULLNESS));
		break;
	case IP_KKT_SOLUTION_METHOD::NORMAL:
		std::cout << std::endl << "Perturbed KKT system cannot be stored with 'NORMAL' method yet.";
		break;
	}

	n = n_new;
	m_eq = m_eq_new; // constraints equalities lines
	m_ineq = m_ineq_new; // constraints inequalities lines (eventually excluding tangential contact forces)
}

void ChInteriorPoint::DumpProblem(std::string suffix) {
	CreateDirectory("dump", nullptr);
    ExportArrayToFile(var.y, "dump/var_y" + suffix + ".txt");
    ExportArrayToFile(var.v, "dump/var_v" + suffix + ".txt");
    ExportArrayToFile(var.gamma, "dump/var_gamma" + suffix + ".txt");
    ExportArrayToFile(var.lambda, "dump/var_lambda" + suffix + ".txt");

    ExportArrayToFile(rhs.b_eq, "dump/rhs_b_eq" + suffix + ".txt");
    ExportArrayToFile(rhs.b_ineq, "dump/rhs_b_ineq" + suffix + ".txt");
    ExportArrayToFile(rhs.c, "dump/rhs_c" + suffix + ".txt");

    BigMat.Compress();
    BigMat.ExportToDatFile("dump/", 8);
}

void ChInteriorPoint::LoadProblem() {
    // ImportArrayFromFile(y, "dump/y.txt");
    // ImportArrayFromFile(v, "dump/v.txt");
    // ImportArrayFromFile(lambda, "dump/lambda.txt");

    ImportArrayFromFile(rhs.b, "dump/b.txt");
    ImportArrayFromFile(rhs.c, "dump/c.txt");

    // BigMat.ImportFromDatFile("dump/");
}

void ChInteriorPoint::DumpIPStatus(std::string suffix) const {
    ExportArrayToFile(var.y, "dump/y" + suffix + ".txt");
    ExportArrayToFile(var.v, "dump/v" + suffix + ".txt");
    ExportArrayToFile(var.lambda, "dump/lambda" + suffix + ".txt");

}

void ChInteriorPoint::make_positive_definite() {
    if( m_ineq == 0 && m_eq == 0)
        return;

    BigMat.ForEachExistentValueInRange([](double* val) { *val *= -1; }, 0, n, n, n + m_eq + m_ineq - 1);
}

// Perform moltiplication of EqCon with vect_in: vect_out = EqCon*vect_in
void ChInteriorPoint::multiplyEqCon(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
	switch (KKT_solve_method)
	{
	case IP_KKT_SOLUTION_METHOD::STANDARD:
		std::cout << std::endl << "EqCon multiplication is not implemented in 'STANDARD' method yet.";
		break;
	case IP_KKT_SOLUTION_METHOD::AUGMENTED:
		BigMat.MatrMultiplyClipped(vect_in, vect_out, n, n + m_eq - 1,
													  0, n - 1,
													  0, 0);
		break;
	case IP_KKT_SOLUTION_METHOD::NORMAL:
		std::cout << std::endl << "EqCon multiplication is not implemented in 'NORMAL' method yet.";
		break;
	}
}

// Perform moltiplication of IneqCon with vect_in: vect_out = IneqCon*vect_in
void ChInteriorPoint::multiplyIneqCon(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
    switch( KKT_solve_method )
    {
        case IP_KKT_SOLUTION_METHOD::STANDARD:
		std::cout << std::endl << "IneqCon multiplication is not implemented in 'STANDARD' method yet.";
			break;
        case IP_KKT_SOLUTION_METHOD::AUGMENTED:
        BigMat.MatrMultiplyClipped(vect_in, vect_out, n + m_eq, n + m_eq + m_ineq - 1, // start | end row
			                                          0,        n - 1,                 // start | end column
			                                          0,        0);
			break;
        case IP_KKT_SOLUTION_METHOD::NORMAL:
        std::cout << std::endl << "IneqCon multiplication is not implemented in 'NORMAL' method yet.";
			break;
    }
}

// Perform moltiplication of -EqCon^T with vect_in: vect_out = -EqCon^T *vect_in (considers that in the top-right part there is
// already -EqCon^T)
void ChInteriorPoint::multiplyNegEqConT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
	switch (KKT_solve_method)
	{
	case IP_KKT_SOLUTION_METHOD::STANDARD:
		std::cout << std::endl << "NegEqCon multiplication is not implemented in 'STANDARD' method yet.";
		break;
	case IP_KKT_SOLUTION_METHOD::AUGMENTED:
		BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1,
			                                          n, n + m_eq - 1,
			                                          0, 0);
		break;
	case IP_KKT_SOLUTION_METHOD::NORMAL:
		std::cout << std::endl << "NegEqCon multiplication is not implemented in 'NORMAL' method yet.";
		break;
	}
}

// Perform moltiplication of -IneqCon^T with vect_in: vect_out = -IneqCon^T *vect_in (considers that in the top-right part there is
// already -IneqCon^T)
void ChInteriorPoint::multiplyNegIneqConT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
    switch( KKT_solve_method )
    {
	case IP_KKT_SOLUTION_METHOD::STANDARD:
		std::cout << std::endl << "NegIneqCon multiplication is not implemented in 'STANDARD' method yet.";
		break;
	case IP_KKT_SOLUTION_METHOD::AUGMENTED:
		BigMat.MatrMultiplyClipped(vect_in, vect_out, 0,        n - 1,
													  n + m_eq, n + m_eq + m_ineq - 1,
													  0, 0);
		break;
	case IP_KKT_SOLUTION_METHOD::NORMAL:
		std::cout << std::endl << "NegIneqCon multiplication is not implemented in 'NORMAL' method yet.";
		break;
    }
}



void ChInteriorPoint::multiplyH(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const {
	switch (KKT_solve_method)
	{
	case IP_KKT_SOLUTION_METHOD::STANDARD:
		std::cout << std::endl << "H multiplication is not implemented in 'STANDARD' method yet.";
		break;
	case IP_KKT_SOLUTION_METHOD::AUGMENTED:
		BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1,
			0, n - 1,
			0, 0);
		break;
	case IP_KKT_SOLUTION_METHOD::NORMAL:
		std::cout << std::endl << "H multiplication is not implemented in 'NORMAL' method yet.";
		break;
	}
}


void ChInteriorPoint::residual_fullupdate(IPresidual_t& residuals, const IPvariables_t& variables) const {
    // Dual Residual
    // res.rd = H*v + c - EqCon^T*gamma - IneqCon^T*lambda
    multiplyH(variables.v, residuals.rd);  // res.rd = H*v
	residuals.rd += rhs.c; // res.rd = (H*v) + c
	if (m_eq > 0)
	{
		multiplyNegEqConT(variables.gamma, vectn);
		residuals.rd += vectn; // res.rd = (H*v + c) + (-EqCon^T*gamma)

		// res.rp_gamma = EqCon*v - b_eq
		multiplyEqCon(variables.v, residuals.rp_gamma);  // res.rp_gamma = EqCon*v
		residuals.rp_gamma -= rhs.b_eq; // res.rp_gamma = (EqCon*v) - b_eq

	}
	else
		residuals.rp_gamma.FillElem(0);


    if( m_ineq > 0 )
    {
        multiplyNegIneqConT(variables.lambda, vectn);  // vectn = (-IneqCon^T)*lambda
        residuals.rd += vectn;                // res.rd = (H*v + c - EqCon^T*gamma) + (-IneqCon^T*lambda)

        // Primal residual
        // res.rp_lambda = IneqCon*v - y - b
        multiplyIneqCon(variables.v, residuals.rp_lambda);  // res.rp_lambda = IneqCon*v
		residuals.rp_lambda -= variables.y; // res.rp_lambda = (IneqCon*v) - y
        residuals.rp_lambda -= rhs.b_ineq; // res.rp_lambda = (IneqCon*v - y) - b_ineq;
        if( ADD_COMPLIANCE )
        {
            E.MatrMultiply(variables.lambda, vectm_ineq);
            residuals.rp_lambda += vectm_ineq;
        }

        residuals.mu = variables.y.MatrDot(variables.y, variables.lambda) / m_ineq;
    }
    else
    {
        residuals.rp_lambda.FillElem(0);
        residuals.mu = 0;
    }
}

/// Export the IP variables in Chrono format
// TODO: FromVectorToUnknowns should accept const reference, but it doesn't. When it will be, we could pass const ChMatrix<>& in the argument and as return
ChMatrix<>& ChInteriorPoint::adapt_to_Chrono(ChSystemDescriptor& sysd, ChMatrix<>& solution_vect) const {

    // copy 'v' as is
    for( auto row_sel = 0; row_sel < n; row_sel++ )
        solution_vect(row_sel, 0) = var.v(row_sel, 0);

	// copy 'gamma' as is
	for (auto row_sel = 0; row_sel < m_eq; row_sel++)
	{
		solution_vect(n + row_sel, 0) = -var.gamma(row_sel, 0);  // there will be an inversion inside FromVectorToUnknowns()
	}

	// copy 'lambda' as is
    if( SKIP_CONTACTS_UV )
    {
		auto ineq_sel = 0;
		auto row_sel = 0;
		auto con_list = sysd.GetConstraintsList();
		for (auto con_sel = 0;  con_sel < sysd.GetConstraintsList().size() && ineq_sel<m_ineq; ++con_sel)
        {
			if ((con_list[con_sel]->GetMode() == CONSTRAINT_UNILATERAL || con_list[con_sel]->GetMode() == CONSTRAINT_FRIC) && con_list[con_sel]->IsActive() )
			{
				if (dynamic_cast<ChConstraintTwoTuplesFrictionTall*>(con_list[con_sel]))
				{
					solution_vect(n + m_eq + row_sel, 0) = 0;
				}
				else
					solution_vect(n + m_eq + row_sel, 0) = -var.lambda(ineq_sel, 0);  // there will be an inversion inside FromVectorToUnknowns()

				++row_sel;
			}
        }
    }
    else
    {
        for( auto row_sel = 0; row_sel < m_ineq; row_sel++ )
            solution_vect(n + m_eq + row_sel, 0) = -var.lambda(row_sel);  // there will be an inversion inside FromVectorToUnknowns()
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

void ChInteriorPoint::Solve(const ChCOOMatrix& augmented_mat, const ChMatrix<double>& rhs_b, const ChMatrix<double>& rhs_c, ChMatrixDynamic<double>& var_x, ChMatrixDynamic<double>& var_y, ChMatrixDynamic<double>& var_gamma, ChMatrixDynamic<double>& var_lambda) {
    // Guys, I know that it's the worst way to to it, but it's just for debug!
    // I promise that, if you want this feature to be left, I will change the settings in order to make moves instead of copies.
    BigMat = augmented_mat;
    rhs.b = rhs_b;
    rhs.c = rhs_c;

    n = rhs.c.GetRows();
    m_ineq = rhs.b.GetRows();

    var_x.Resize(n, 1);
    var_gamma.Resize(m_eq, 1);
    var_y.Resize(m_ineq, 1);
    var_lambda.Resize(m_ineq, 1);

    reset_internal_dimensions(n, m_eq, m_ineq, m_ineq);
    auto bkp_method = KKT_solve_method;

    KKT_solve_method = IP_KKT_SOLUTION_METHOD::AUGMENTED;


    /********** Check if system has constraints **********/
    if( m_ineq == 0 )  // if no constraints
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

		if (verbose)
			PrintIPStatus();

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
    var_x = var.v;
    var_x = var.y;
    var_x = var.lambda;

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


void ChInteriorPoint::PrintIPStatus() const
{
	GetLog() << "InteriorPoint | Call: " << solver_call << " Iter: " << iteration_count << "/" << iteration_count_max << "\n";
	GetLog() << "Complementarity Measure: " << res.mu << "\n";
	GetLog() << "|rd|/n (stationarity): " << res.rd.NormTwo() / n << "\n";
	GetLog() << "|rp_gamma|/m_eq (BILAT violation): " << (m_eq == 0 ? 0 : res.rp_gamma.NormTwo() / m_eq) << "\n";
	GetLog() << "|rp_lambda|/m_ineq (UNILAT violation): " << res.rp_lambda.NormTwo() / m_ineq << "\n";
	GetLog() << "Objective Function: " << evaluate_objective_function() << "\n";
	GetLog() << "\n";
}
}
