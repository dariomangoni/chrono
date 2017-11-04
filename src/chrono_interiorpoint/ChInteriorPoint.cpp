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

#include "ChInteriorPoint.h"
#include "ChInteriorPointUtils.h"
#include "solver/ChConstraintTwoTuplesFrictionT.h"
#include <algorithm>
#include <iomanip>

#ifdef CHRONO_POSTPROCESS
#include "chrono_postprocess/ChGnuPlot.h"
#endif

#define DEBUG_TRAPS
#define DEBUG_CONSOLE
#define CORRECT_STARTING_POINTS
#define ADD_COMPLIANCE false
#define REUSE_OLD_SOLUTIONS false
//#define USE_MY_MAX_STEPLENGTHS
//#define BYPASS_RESTITUTION_TERM

//TODO: drive iterations toward the term that has not reached the tolerated level yet

namespace chrono {


double ChInteriorPoint::Solve(ChSystemDescriptor& sysd) {
    solver_call++;

	ip_timer_solve_assembly.start();
    /********** Load system **********/
	// for optimization of the starting point

    // set offsets so that equality constraints come first
    sysd.SortActiveConstraints();

    // get matrices size and identify the constraint mode for inequality constraints
    auto constr_list = sysd.GetConstraintsList();
    n = sysd.CountActiveVariables();
    m_eq = 0;
    auto m_ineq_tot = 0;
    auto nnz_scaling = 0;
    ineq_mode.clear();
    for (auto it = 0; it<constr_list.size(); ++it)
    {
        switch(constr_list[it]->GetMode())
        {
        case eChConstraintMode::CONSTRAINT_FREE:
            ++m_eq; break;
        case eChConstraintMode::CONSTRAINT_LOCK:
            ++m_eq; break;
        case eChConstraintMode::CONSTRAINT_UNILATERAL:
            ineq_mode.push_back(eChConstraintModeMOD::CONSTRAINT_UNILATERAL);
            nnz_scaling += 1;
            break;
        case eChConstraintMode::CONSTRAINT_FRIC:
            ++m_ineq_tot;
            if (!dynamic_cast<ChConstraintTwoTuplesFrictionTall*>(constr_list[it]))
            {
                ineq_mode.push_back(eChConstraintModeMOD::CONSTRAINT_FRIC_N);
                nnz_scaling += 1;
            }
            else
            // we assume that CONSTRAINT_FRIC has only N and UV type of constraints, otherwise if (enable_friction_cones && dynamic_cast<ChConstraintTwoTuplesFrictionTall*>(constr_list[it]))
            if (!skip_contacts_uv)
            {
                ineq_mode.push_back(eChConstraintModeMOD::CONSTRAINT_FRIC_UV);
                nnz_scaling += 4;
            }
            break;
        default:
            assert(0);
        }
    }
    m_ineq = static_cast<int>(ineq_mode.size());

    // scale variables
    var.Reset(n, m_eq, m_ineq);
    res.Reset(n, m_eq, m_ineq);
    yl_scaled.Reset(m_ineq, 1);

    rhs.Resize(n, m_eq, m_ineq);


    // update mutables
    vectn.Resize(n, 1);
    vectm_eq.Resize(m_eq, 1);
    vectm_ineq.Resize(m_ineq, 1);
    sol_chrono.Resize(n + m_eq + m_ineq_tot, 1);

    // Let the matrix acquire the information about ChSystem
    if (m_force_sparsity_pattern_update)
    {
        m_force_sparsity_pattern_update = false;

        ChSparsityPatternLearner sparsity_learner(n + m_eq + m_ineq, n + m_eq + m_ineq, true);
        sysd.ConvertToMatrixForm(&sparsity_learner, nullptr);
        for (auto cs = 0; cs < ineq_mode.size(); ++cs)
        {
            if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
            {
                sparsity_learner.SetElement(cs, cs, 0.0);
                continue;
            }


            if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
            {
                sparsity_learner.SetElement(cs+0, cs+0, 0.0);
                sparsity_learner.SetElement(cs+1, cs+0, 0.0);
                sparsity_learner.SetElement(cs+2, cs+0, 0.0);
                sparsity_learner.SetElement(cs+0, cs+1, 0.0);
                sparsity_learner.SetElement(cs+1, cs+1, 0.0);
                sparsity_learner.SetElement(cs+2, cs+1, 0.0);
                sparsity_learner.SetElement(cs+0, cs+2, 0.0);
                sparsity_learner.SetElement(cs+1, cs+2, 0.0);
                sparsity_learner.SetElement(cs+2, cs+2, 0.0);
            }
        }

        BigMat.LoadSparsityPattern(sparsity_learner);
    }
    else
    {
        // If an NNZ value for the underlying matrix was specified, perform an initial resizing, *before*
        // a call to ChSystemDescriptor::ConvertToMatrixForm(), to allow for possible size optimizations.
        // Otherwise, do this only at the first call, using the default sparsity fill-in.
        if (ip_solver_call == 0) {
            BigMat.Reset(n + m_eq + m_ineq, n + m_eq + m_ineq, static_cast<int>((n + m_eq + m_ineq) * ((n + m_eq + m_ineq) * SPM_DEF_FULLNESS)));
        }
    }

    sysd.ConvertToMatrixForm(&BigMat, nullptr, false, skip_contacts_uv, ADD_COMPLIANCE);

	if (!leverage_symmetry)
		make_positive_definite();

    sysd.ConvertToMatrixForm(nullptr, nullptr, nullptr, &rhs.c, &rhs.b, nullptr, false, skip_contacts_uv);  // load f->c and b->b
    rhs.c.MatrScale(-1); // adapt to InteriorPoint convention
    rhs.b.MatrScale(-1); // adapt to InteriorPoint convention

	ip_timer_solve_assembly.stop();
#ifdef BYPASS_RESTITUTION_TERM
    rhs.b.FillElem(0);
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

		if (leverage_symmetry) // then the gamma have flipped sign
			for (auto row_sel = 0; row_sel < m_eq; row_sel++)
				mumps_rhs(n + row_sel,0) *= -1;

        sysd.FromVectorToUnknowns(mumps_rhs);

        // Export variable so that can be used in the next iteration as starting point
        var.v.Resize(n, 1);
        var.gamma.Resize(m_eq, 1);
		for (auto row_sel = 0; row_sel < n; row_sel++)
			var.v.SetElement(row_sel, 0, mumps_rhs.GetElement(row_sel, 0));
		for (auto row_sel = 0; row_sel < m_eq; row_sel++)
			var.gamma.SetElement(row_sel, 0, leverage_symmetry ? -mumps_rhs.GetElement(n + row_sel, 0) : mumps_rhs.GetElement( n + row_sel, 0));

        if( verbose )
            std::cout << "IP call: " << solver_call << "; No inequality constraints." << std::endl;

        return 0.0;
    }

    /********* The system DOES have constraints! Start Interior Point ********/
	ip_solver_call++;
	ip_timer_solver_solvercall.start();

	if (m_eq > 0)
	{
		rhs.b_eq.PasteClippedMatrix(rhs.b, 0, 0, m_eq, 1, 0, 0);
		rhs.b_ineq.PasteClippedMatrix(rhs.b, m_eq, 0, m_ineq, 1, 0, 0);
	}
	else
		rhs.b_ineq = rhs.b;

    PrintCVXOPTproblem(sysd, "cvxopt_problem");


    if( ADD_COMPLIANCE && m_ineq > 0 )
    {
        sysd.ConvertToMatrixForm(nullptr, nullptr, &E, nullptr, nullptr, nullptr, false, skip_contacts_uv);
        E *= -1;
    }

    set_feasible_starting_point();

    if (verbose)
    {
        std::cout << "IP " << solver_call << "| n: " << n << "; m_eq: " << m_eq << "; m_ineq: " << m_ineq << std::endl;
        std::cout << "       || *** prediction *** ||      || correction  ||" << std::endl;
        std::cout << "  | it || a_pr | a_du |  mu  || tau  || a_pr | a_du ||  mu  |    resP    |   res_eq   |  res_ineq  |" << std::endl;
    }

    for( iteration_count = 0; iteration_count < iteration_count_max; iteration_count++ )
    {

        iterate();
        iteration_count_tot++;

        
        if( res <= res_nnorm_tol ) // check for exit conditions (less OR EQUAL in order to catch zero residual (for example, no bilateral))
        {
            if (verbose)
                PrintIPStatus();
            break;
        }
    }

    ip_timer_solver_solvercall.stop();

    // Scatter the solution into the Chrono environment
    sysd.FromVectorToUnknowns(adapt_to_Chrono(sysd, sol_chrono));

    return 0.0;
}

// Iterating function
void ChInteriorPoint::iterate() {
    /*********************************************************************************/
    /***************************** Prediction Phase **********************************/
    /*********************************************************************************/

    //DumpProblem("pre");
    // Paste scaling matrix W^T * W
    computeNesterovToddScalingMatrix(var);
    scaling_matrix.MatrMultiplyClipped(scaling_matrix, BigMat, 0, m_ineq - 1, 0, m_ineq - 1, 0, n + m_eq, true, 0, 0, n + m_eq);
    if (ADD_COMPLIANCE)
        for (auto diag_sel = 0; diag_sel < m_ineq; ++diag_sel)
            BigMat.SetElement(n + m_eq + diag_sel, n + m_eq + diag_sel, leverage_symmetry ? -E.GetElement(diag_sel - n + m_eq, diag_sel - n + m_eq) : +E.GetElement(diag_sel - n + m_eq, diag_sel - n + m_eq), false);
    factorize_system_matrix();
    scaling_matrix.MatrMultiply(var.lambda, yl_scaled);

    DumpProblem("post");

#ifdef DEBUG_TRAPS
    assert(is_valid(yl_scaled, m_ineq));
    assert(is_valid(var.y, m_ineq));
    assert(is_valid(var.lambda, m_ineq));
#endif

    // WARNING: the residual structure 'res' must be already updated at this point!

    // fill 'mumps_rhs' with rhs [-res.rd;-res.rp_gamma;-res.rp_lambda-y]
    mumps_rhs.Resize(n + m_eq + m_ineq, 1);
    for( auto row_sel = 0; row_sel < n;      row_sel++ ) mumps_rhs(row_sel           ) = -res.rd(row_sel);
	for( auto row_sel = 0; row_sel < m_eq;   row_sel++ ) mumps_rhs(row_sel + n       ) = -res.rp_gamma(row_sel);
    for( auto row_sel = 0; row_sel < m_ineq; row_sel++ ) mumps_rhs(row_sel + n + m_eq) = -res.rp_lambda(row_sel) - var.y(row_sel);

    Dvar_pred.Resize(n, m_eq, m_ineq);
    get_Newton_direction(Dvar_pred, mumps_rhs, res);

#ifdef DEBUG_TRAPS
    assert(is_valid(Dvar_pred.y, m_ineq));
    assert(is_valid(Dvar_pred.lambda, m_ineq));
#endif

    /*** compute step lengths ***/
    // from 16.60 pag.482 from 14.32 pag.408 (remember that y>=0!)
#ifdef USE_MY_MAX_STEPLENGTHS
    auto alfa_pred_prim = get_Newton_steplength_MAX_mod(var.y, Dvar_pred.y);
    auto alfa_pred_dual = get_Newton_steplength_MAX_mod(var.lambda, Dvar_pred.lambda);
#else
    auto alfa_pred_prim = get_Newton_steplength_MAX(var.y, Dvar_pred.y);
    auto alfa_pred_dual = get_Newton_steplength_MAX(var.lambda, Dvar_pred.lambda);
#endif

#ifdef DEBUG_CONSOLE
    printf("  |%4d", iteration_count);
    printf("||%-6.3g|%-6.3g", alfa_pred_prim, alfa_pred_dual);
    //printf("|%-6d|%-6d|", find_out_of_cone(var.y, Dvar_pred.y, alfa_pred_prim), find_out_of_cone(var.lambda, Dvar_pred.lambda, alfa_pred_dual));
#endif

    if( equal_step_lengths )
    {
        auto alfa_pred = std::min(alfa_pred_prim, alfa_pred_dual);
        alfa_pred_prim = alfa_pred;
        alfa_pred_dual = alfa_pred;
    }


    /*** make the prediction step ***/
	// update only variables needed for complementarity measure
    var_pred.Resize(n, m_eq, m_ineq);

    var_pred.y = Dvar_pred.y;
	var_pred.lambda = Dvar_pred.lambda;

	var_pred.y.MatrScale(alfa_pred_prim);
	var_pred.lambda.MatrScale(alfa_pred_dual);

	var_pred.y += var.y;
	var_pred.lambda += var.lambda;

    /*** compute complementarity measure ***/
    auto mu_pred = var_pred.y.MatrDot(var_pred.y, var_pred.lambda) / m_ineq;  // from 16.56 pag.481


    if( only_predict)
    {
		// update remaining variables (v, gamma)
        var_pred.v = Dvar_pred.v;
    	var_pred.v.MatrScale(alfa_pred_prim);
    	var_pred.v += var.v;

		var_pred.gamma = Dvar_pred.gamma;
		var_pred.gamma.MatrScale(alfa_pred_dual);
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
		if (!equal_step_lengths)
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

    /*** evaluate centering parameter ***/
    auto sigma = std::pow(std::max(0.0, std::min(1.0, mu_pred / res.mu)), 3.0); // from [5] pag. 12
    //auto sigma = std::pow(mu_pred / res.mu, 3.0); // from [5] pag. 12

    /*** step length correction ***/
    auto tau = adaptive_eta ? exp(-res.mu * m_ineq) * 0.1 + 0.9 : 0.95;  // exponential descent of tau

    auto rhs_corr = 0;
    //auto rhs_corr = sigma;

    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            vectm_ineq(cs) = sigma*res.mu - Dvar_pred.lambda(cs)*Dvar_pred.y(cs);
            continue;
        }
        
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_UV)
        {
            vectm_ineq(cs) = - Dvar_pred.lambda(cs)*Dvar_pred.y(cs);
            continue;
        }
    }

    inverse_Hadamard(yl_scaled, vectm_ineq, vectm_ineq);
    for (auto row_sel = 0; row_sel < m_ineq; row_sel++) mumps_rhs.SetElement(row_sel + n + m_eq, 0, -yl_scaled(row_sel) + vectm_ineq(row_sel));
    scaling_matrix.MatrMultiplyClipped(mumps_rhs, vectm_ineq, 0, m_ineq - 1, 0, m_ineq - 1, n + m_eq, 0, true, 0, 0, 0, true);


    /*** find directions ***/
	for (auto row_sel = 0; row_sel < n;      row_sel++) mumps_rhs(row_sel           ) =  -(1 - rhs_corr)*res.rd(row_sel);
	for (auto row_sel = 0; row_sel < m_eq;   row_sel++) mumps_rhs(row_sel + n       ) =  -(1 - rhs_corr)*res.rp_gamma(row_sel);
	for (auto row_sel = 0; row_sel < m_ineq; row_sel++) mumps_rhs(row_sel + n + m_eq) =  -(1 - rhs_corr)*res.rp_lambda(row_sel) + vectm_ineq(row_sel);

    Dvar.Resize(n, m_eq, m_ineq);
    get_Newton_direction(Dvar, mumps_rhs, res);


    /*** compute step lengths ***/
#ifdef USE_MY_MAX_STEPLENGTHS
    auto alfa_corr_prim = get_Newton_steplength_MAX_mod(var.y, Dvar.y);
    auto alfa_corr_dual = get_Newton_steplength_MAX_mod(var.lambda, Dvar.lambda);
#else
    auto alfa_corr_prim = get_Newton_steplength_MAX(var.y, Dvar.y);
    auto alfa_corr_dual = get_Newton_steplength_MAX(var.lambda, Dvar.lambda);
#endif


#ifdef DEBUG_CONSOLE
    printf("|%-6.0g", mu_pred);
    printf("||%-6.0g", tau);
    printf("||%-6.0g|%-6.0g", alfa_corr_prim, alfa_corr_dual);
    //printf("|%-6d|%-6d|", find_out_of_cone(var.y, Dvar.y, alfa_corr_prim), find_out_of_cone(var.lambda, Dvar.lambda, alfa_corr_dual));
#endif

    if( equal_step_lengths )
    {
        auto alfa_corr = std::min(alfa_corr_prim, alfa_corr_dual);
        alfa_corr_prim = alfa_corr;
        alfa_corr_dual = alfa_corr;
    }

    // apply correction
    alfa_corr_prim *= tau;
    alfa_corr_dual *= tau;


    /*** make the correction step ***/
    var_corr.Resize(n, m_eq, m_ineq);
    var_corr.v = Dvar.v;        var_corr.v.MatrScale(alfa_corr_prim);        var_corr.v += var.v;        var.v = var_corr.v;
    var_corr.y = Dvar.y;        var_corr.y.MatrScale(alfa_corr_prim);        var_corr.y += var.y;        var.y = var_corr.y;
    var_corr.gamma = Dvar.gamma;      var_corr.gamma.MatrScale(alfa_corr_dual);     var_corr.gamma += var.gamma;        var.gamma = var_corr.gamma;
    var_corr.lambda = Dvar.lambda;    var_corr.lambda.MatrScale(alfa_corr_dual);    var_corr.lambda += var.lambda;    var.lambda = var_corr.lambda;

    /********** Residuals update **********/
	res.mu = var.y.MatrDot(var.y, var.lambda) / m_ineq;  // from 14.6 pag.395
    res.rp_gamma.MatrScale(1 - alfa_corr_prim);
    res.rp_lambda.MatrScale(1 - alfa_corr_prim);

	res.rd.MatrScale(1 - alfa_corr_dual);
    if( !equal_step_lengths )
    {
        multiplyH(Dvar.v, vectn);                          // vectn = H*Dvar.v
        vectn.MatrScale(alfa_corr_prim - alfa_corr_dual);  // vectn = (alfa_pred_prim - alfa_pred_dual) * (H * Dvar.v)
        res.rd += vectn;
    }


    if( print_history )
        output_log();

#ifdef DEBUG_CONSOLE
    printf("||%-6.0g", res.mu);
    printf("|%-12.0g|%-12.0g|%-12.0g|", res.rd.NormTwo(), res.rp_gamma.NormTwo(), res.rp_lambda.NormTwo());
    printf("\n");
#endif


}

void ChInteriorPoint::factorize_system_matrix() {

    BigMat.Compress();
    mumps_engine.SetMatrix(BigMat);
    for (auto loop_expand_workspace_size = 0; loop_expand_workspace_size < 5; ++loop_expand_workspace_size)
    {
        if (mumps_engine.MumpsCall(ChMumpsEngine::ANALYZE_FACTORIZE) == 9)
        {
			mumps_engine.SetICNTL(14, static_cast<int>(round(mumps_engine.GetICNTL(14)*1.5)));
			std::cout << "MUMPS work space will be allocated overestimating the estimate by " << mumps_engine.GetICNTL(14) << "%" << std::endl;
        }
		else
			break;

    }

}

void ChInteriorPoint::get_Newton_direction(IPvariables_t& Dvar_unknown, ChMatrix<>& rhs, const IPresidual_t& residuals) {
    // Solve the KKT system
    mumps_engine.SetRhsVector(rhs);
    if( mumps_engine.MumpsCall(ChMumpsEngine::SOLVE) )
        mumps_engine.PrintINFOG();
    if( mumps_engine.GetRINFOG(6) > 1e-6 )
        std::cout << "Scaled residual norm of MUMPS call: " << mumps_engine.GetRINFOG(6) << std::endl;

    // MUMPS uses its 'rhs' vector to store the solution. In order to clarify that:
    const ChMatrixDynamic<double>& sol = rhs;

    // Extract 'v' and 'lambda' from 'sol'
    for (auto row_sel = 0; row_sel < n;      row_sel++ ) Dvar_unknown.v.SetElement(row_sel, 0, sol.GetElement(row_sel, 0));
	if (leverage_symmetry)
	{
		for (auto row_sel = 0; row_sel < m_eq; row_sel++) Dvar_unknown.gamma.SetElement(row_sel, 0, -sol.GetElement(row_sel + n, 0));
		for (auto row_sel = 0; row_sel < m_ineq; row_sel++) Dvar_unknown.lambda.SetElement(row_sel, 0, -sol.GetElement(row_sel + n + m_eq, 0));
	}
	else
	{
		for (auto row_sel = 0; row_sel < m_eq; row_sel++) Dvar_unknown.gamma.SetElement(row_sel, 0, sol.GetElement(row_sel + n, 0));
		for (auto row_sel = 0; row_sel < m_ineq; row_sel++) Dvar_unknown.lambda.SetElement(row_sel, 0, sol.GetElement(row_sel + n + m_eq, 0));
	}


    // Calc 'y' (it is also possible to evaluate y as (-lambda°y+sigma*res.mu*e-y°lambda)./lambda )
    multiplyIneqCon(Dvar_unknown.v, Dvar_unknown.y);  // Dy = IneqCon*Dv
    Dvar_unknown.y += residuals.rp_lambda; // Dy = (IneqCon*Dv) + rp_lambda
    if( ADD_COMPLIANCE )
    {
        E.MatrMultiply(Dvar_unknown.lambda, vectm_ineq);
        Dvar_unknown.y += vectm_ineq;
    }
}

double ChInteriorPoint::get_Newton_steplength_MAX(const ChMatrix<double>& pos, const ChMatrix<double>& dir) const
{
    // we are going to solve max{alfa | pos + alfa*dir >= 0 }
    double alfa_inv = 0.0;
    ChMatrixNM<double, 3, 1> a_pr;
    ChMatrixNM<double, 3, 1> pos_n; // normalized scaled variable
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            if (-dir(cs) / pos(cs) > alfa_inv)
            {
                alfa_inv = -dir(cs) / pos(cs);
            }

            continue;
        }


        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            auto pos_p = projection_on_polar_cone(pos, cs);

            pos_n(0) = pos(cs + 0) / pos_p;
            pos_n(1) = pos(cs + 1) / pos_p;
            pos_n(2) = pos(cs + 2) / pos_p;

            auto pd_p = pos_n(0)*dir(cs + 0) - pos_n(1)*dir(cs + 1) - pos_n(2)*dir(cs + 2);
            a_pr(0) = pd_p / pos_p;
            a_pr(1) = dir(cs + 1) - (pd_p + dir(cs + 0)) / (pos_n(0) + 1)*pos_n(1);
            a_pr(2) = dir(cs + 2) - (pd_p + dir(cs + 0)) / (pos_n(0) + 1)*pos_n(2);

            auto alfa_inv_temp = -a_pr(0) + sqrt(a_pr(1)*a_pr(1) + a_pr(2)*a_pr(2));

            if (alfa_inv_temp > alfa_inv)
                alfa_inv = alfa_inv_temp;

            continue;
        }
    } // end FOR loop

    return std::min(1.0 / alfa_inv, 1.0);
}

double ChInteriorPoint::get_Newton_steplength_MAX_mod(const ChMatrix<double>& pos, const ChMatrix<double>& dir, bool bound_from_zero_to_one) const
{
    // we are going to solve max{alfa | pos + alfa*dir >= 0 }
    double alfa = +1000;
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        auto pd = -pos(cs) / dir(cs);
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            if (dir(cs) < 0) // otherwise it gives no bounds to MAX of alfa
                alfa = std::min(alfa, pd);

            continue;
        }


        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            auto alfa_a = uTJv(dir, dir, cs);

            if (alfa_a > 0 && dir(cs) > 0)
                continue; // no MAX limit imposed

            auto alfa_k = uTJv(pos, dir, cs);
            auto alfa_c = uTJv(pos, pos, cs);

            auto alfa_new = (-alfa_k - sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;

            if (alfa_a < 0 && dir(cs) < 0)
            {
                //assert(pd > alfa_min);
                //auto alfa_max = (-alfa_k - sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;
                alfa = std::min({ alfa, alfa_new, pd });
                continue;
            }

            if (alfa_a > 0 && dir(cs) < 0)
            {
                //auto alfa_min = (-alfa_k - sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;
                alfa = std::min({ alfa, alfa_new, pd });
                continue;
            }

            if (alfa_a < 0 && dir(cs) > 0)
            {
                // assert(alfa_max>pd)
                //auto alfa_max = (-alfa_k - sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;
                alfa = std::min(alfa, alfa_new);
                continue;
            }

        }
    } // end FOR loop

    return bound_from_zero_to_one ? std::min(std::max(0.0, alfa), 1.0) : alfa;
}

double ChInteriorPoint::get_Newton_steplength_MIN(const ChMatrix<double>& pos, const ChMatrix<double>& dir, bool bound_from_zero_to_one) const
{
    // we are going to solve min{alfa | pos + alfa*dir >= 0 }
    double alfa = -1000;
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        auto pd = -pos(cs) / dir(cs);
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            if (dir(cs) > 0) // otherwise it gives no bounds to MIN of alfa
                alfa = std::max(alfa, pd);

            continue;
        }


        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            auto alfa_a = uTJv(dir, dir, cs);

            if (alfa_a > 0 && dir(cs) < 0)
                continue; // no MIN limit imposed

            auto alfa_k = uTJv(pos, dir, cs);
            auto alfa_c = uTJv(pos, pos, cs);

            auto alfa_new = (-alfa_k + sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;

            if (alfa_a < 0 && dir(cs) < 0)
            {
                //assert(pd > alfa_min);
                //auto alfa_min = (-alfa_k + sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;
                alfa = std::max({ alfa, alfa_new });
                continue;
            }

            if (alfa_a > 0 && dir(cs) > 0)
            {
                //auto alfa_max = (-alfa_k + sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;
                alfa = std::max({ alfa, alfa_new, pd });
                continue;
            }

            if (alfa_a < 0 && dir(cs) > 0)
            {
                // assert(alfa_max>pd)
                //auto alfa_min = (-alfa_k + sqrt(alfa_k*alfa_k - alfa_a*alfa_c)) / alfa_a;
                alfa = std::max({ alfa, alfa_new, pd });
                continue;
            }
        }
    } // end FOR loop

    return bound_from_zero_to_one ? std::min(std::max(0.0, alfa), 1.0) : alfa;
}

double ChInteriorPoint::get_Newton_steplength_MIN_along_cone_axis(const ChMatrix<double>& pos) const
{
    // we are going to solve min{alfa | pos + alfa*[1;0;0] >= 0 i.e. pos + alfa*[1;0;0] is inside the cone with equation u(0)>=norm(u(1:end))}
    double alfa_lb = -1000;
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            alfa_lb = std::max(alfa_lb, -pos(cs));
            continue;
        }


        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            alfa_lb = std::max(alfa_lb, -pos(cs) + sqrt(pos(cs + 1)*pos(cs + 1) + pos(cs + 2)*pos(cs + 2)));
            continue;
        }
    } // end FOR loop

    return alfa_lb;
}



void ChInteriorPoint::set_feasible_starting_point()
{
    // fill SE corner with identity matrix
    for (auto diag_sel = 0; diag_sel < m_ineq; ++diag_sel)
        BigMat.SetElement(n + m_eq + diag_sel, n + m_eq + diag_sel, leverage_symmetry ? -1 : +1, true);

    mumps_rhs.Resize(n + m_eq + m_ineq, 1);
    for (auto row_sel = 0; row_sel < n; row_sel++) mumps_rhs(row_sel) = -rhs.c(row_sel, 0);
    for (auto row_sel = 0; row_sel < m_eq; row_sel++) mumps_rhs(row_sel + n) = +rhs.b_eq(row_sel, 0);
    for (auto row_sel = 0; row_sel < m_ineq; row_sel++) mumps_rhs(row_sel + n + m_eq) = +rhs.b_ineq(row_sel, 0);

    BigMat.Compress();
    mumps_engine.SetProblem(BigMat, mumps_rhs);
    if (mumps_engine.MumpsCall(ChMumpsEngine::COMPLETE))
        mumps_engine.PrintINFOG();

    var.v.Resize(n, 1);
    var.gamma.Resize(m_eq, 1);
    var.lambda.Resize(m_ineq, 1);
    for (auto row_sel = 0; row_sel < n; row_sel++) var.v(row_sel) = mumps_rhs(row_sel);
    for (auto row_sel = 0; row_sel < m_eq; row_sel++) var.gamma(row_sel) = leverage_symmetry ? -mumps_rhs(n + row_sel) : mumps_rhs(n + row_sel);
    for (auto row_sel = 0; row_sel < m_ineq; row_sel++) vectm_ineq(row_sel) = leverage_symmetry ? -mumps_rhs(n + m_eq + row_sel) : mumps_rhs(n + m_eq + row_sel);


    auto alfa_p = get_Newton_steplength_MIN_along_cone_axis(-vectm_ineq);
    auto alfa_d = get_Newton_steplength_MIN_along_cone_axis(vectm_ineq);

    if (alfa_p >= -1e-8)
        for (auto cs = 0; cs < ineq_mode.size(); ++cs)
            var.y(cs) = ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_UV ? -vectm_ineq(cs) : -vectm_ineq(cs) + (1 + alfa_p);
    else
        var.y = -vectm_ineq;

    if (alfa_d >= -1e-8)
        for (auto cs = 0; cs < ineq_mode.size(); ++cs)
            var.lambda(cs) = ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_UV ? vectm_ineq(cs) : vectm_ineq(cs) + (1 + alfa_d);
    else
        var.lambda = vectm_ineq;


    residual_fullupdate(res, var);


#ifdef DEBUG_TRAPS
    if (find_out_of_cone(var.y) != -1) 
        std::cout << "y out of cone: " << find_out_of_cone(var.y) << std::endl;
    if (find_out_of_cone(var.lambda) != -1)
        std::cout << "lambda out of cone: " << find_out_of_cone(var.lambda) << std::endl;
    assert(is_valid(var.y, m_ineq));
    assert(is_valid(var.lambda, m_ineq));
#endif

}


double ChInteriorPoint::evaluate_objective_function() const {
    multiplyH(var.v, vectn);
    auto obj_value = vectn.MatrDot(var.v, vectn);
    obj_value += rhs.c.MatrDot(var.v, rhs.c);

    return obj_value;
}

double ChInteriorPoint::projection_on_polar_cone(const ChMatrixDynamic<double>& z, int offset)
{
    return std::sqrt( z(offset)*z(offset) - z(offset + 1)*z(offset + 1) - z(offset + 2)*z(offset + 2) );
}

inline double ChInteriorPoint::uTJv(const ChMatrix<double>& u, const ChMatrix<double>& v, int offset) 
{
    return u(offset)*v(offset) - u(offset + 1)*v(offset + 1) - u(offset + 2)*v(offset + 2);
}


void ChInteriorPoint::inverse_Hadamard(const ChMatrix<>& v1, const ChMatrix<>& v2, ChMatrix<>& v_out) const
{

    double sf = 0;
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            v_out(cs) = v2(cs) / v1(cs);
            continue;
        }


        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            sf = 1 / std::pow(projection_on_polar_cone(v1, cs), 2.0);

            // backup v2 (needed if v_out == v2)
            auto v2_0 = v2(cs + 0);
            auto v2_1 = v2(cs + 1);
            auto v2_2 = v2(cs + 2);

            v_out(cs + 0) = ( v1(cs + 0) * v2_0 - v1(cs + 1) * v2_1 - v1(cs + 2) * v2_2) / sf;
            v_out(cs + 1) = (-v1(cs + 1) * v2_0 + (sf + std::pow(v1(cs + 1), 2.0)) / v1(cs + 0) * v2_1 + v1(cs + 1)*v1(cs + 2) / v1(cs + 0) * v2_2) / sf;
            v_out(cs + 2) = (-v1(cs + 2) * v2_0 + v1(cs + 1)*v1(cs + 2) / v1(cs + 0) * v2_1 + (sf + std::pow(v1(cs + 2), 2.0)) / v1(cs + 0) * v2_2) / sf;
            cs += 2;
            continue;
        }
    }
    
}

void ChInteriorPoint::computeNesterovToddScalingMatrix(const IPvariables_t& var)
{
    // resize the scaling matrix
    auto nnz = 0;
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL)
            nnz++;

        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
            skip_contacts_uv ? ++nnz : nnz += 9;

        // skip UV, they have already been considered in the previous 'if'
    }
    scaling_matrix.Reset(m_ineq, m_ineq, nnz);

    // evaluate the scaling matrix
    for (auto cs = 0; cs<ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            scaling_matrix.SetElement(cs, cs, sqrt(var.y(cs))/ sqrt(var.lambda(cs)));
            continue;
        }

        
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            // grant that we are handling UV contact components
            assert(!skip_contacts_uv);
            // grant that the friction force N is followed by two UV
            assert(ineq_mode[cs + 1] == eChConstraintModeMOD::CONSTRAINT_FRIC_UV &&
                ineq_mode[cs + 2] == eChConstraintModeMOD::CONSTRAINT_FRIC_UV
                && "A friction reaction force along N is not followed by a couple of constraints on UV");

            // WARNING: assembling the scaling matrix is possible only if the vectors 'y' and 'lambda'
            // are STRICTLY inside the cone otherwise you will get NaN (if on the border) or complex numbers (if outside)

#ifdef CORRECT_STARTING_POINTS // there are multiple calls to uTJv that can be avoided
            auto alfa_y = -var.y(cs) + sqrt(var.y(cs + 1)*var.y(cs + 1) + var.y(cs + 2)*var.y(cs + 2));
            auto alfa_lam = -var.lambda(cs) + sqrt(var.lambda(cs + 1)*var.lambda(cs + 1) + var.lambda(cs + 2)*var.lambda(cs + 2));
            // if alfa is positive it means that we have to 
            if (alfa_y >= -1e-8 || var.y(cs)<=0)
                const_cast<ChMatrixDynamic<>&>(var.y)(cs) += std::max(alfa_y*1.000001, 1e-6);
            if (alfa_lam >= -1e-8 || var.lambda(cs) <= 0)
                const_cast<ChMatrixDynamic<>&>(var.lambda)(cs) += std::max(alfa_lam*1.000001, 1e-6);

#endif
            auto y_proj = std::sqrt(uTJv(var.y, var.y, cs));
            auto lam_proj = std::sqrt(uTJv(var.lambda, var.lambda, cs));
            auto gamma = std::sqrt((1 + (var.lambda(cs)*var.y(cs) + var.lambda(cs + 1)*var.y(cs + 1) + var.lambda(cs + 2)*var.y(cs + 2)) / y_proj / lam_proj) / 2);


            ChVector<double> nsp; // Normalized Scaling Point  i.e. \bar{w}_k
            nsp[0] = (var.y(cs    ) / y_proj + var.lambda(cs    ) / lam_proj) / 2 / gamma;
            nsp[1] = (var.y(cs + 1) / y_proj - var.lambda(cs + 1) / lam_proj) / 2 / gamma;
            nsp[2] = (var.y(cs + 2) / y_proj - var.lambda(cs + 2) / lam_proj) / 2 / gamma;

            ChMatrix33<double> nsm; // Normalized Scaling SubMatrix i.e. \bar{W}_k
            nsm[0][0] = nsp[0];
            nsm[0][1] = nsp[1];
            nsm[0][2] = nsp[2];

            nsm[1][0] = nsp[1];
            nsm[2][0] = nsp[2];

            nsm[1][1] = 1 + nsp[1] * nsp[1] / (nsp[0] + 1);
            nsm[1][2] =     nsp[1] * nsp[2] / (nsp[0] + 1);
            nsm[2][1] =     nsp[2] * nsp[1] / (nsp[0] + 1);
            nsm[2][2] = 1 + nsp[2] * nsp[2] / (nsp[0] + 1);


#ifdef DEBUG_TRAPS
            // // from [5] pag.9
            // deb = 2*J*(nsp*nsp')*J - J
            ChMatrix33<double> deb; 
            deb(0,0) = + 2 * nsp[0]*nsp[0] - 1;
            deb(0,1) = - 2 * nsp[0]*nsp[1];
            deb(0,2) = - 2 * nsp[0]*nsp[2];
            deb(1,0) = - 2 * nsp[1]*nsp[0];
            deb(1,1) = + 2 * nsp[1]*nsp[1] + 1;
            deb(1,2) = + 2 * nsp[1]*nsp[2];
            deb(2,0) = - 2 * nsp[2]*nsp[0];
            deb(2,1) = + 2 * nsp[2]*nsp[1];
            deb(2,2) = + 2 * nsp[2]*nsp[2] + 1;

            // deb = W * deb * W
            ChMatrix33<double> temp;
            temp.MatrMultiply(nsm, deb);
            deb.MatrMultiply(temp, nsm);

            // deb = identity            
            for (auto i = 0; i < 3; ++i)
                for (auto j = 0; j < 3; ++j)
                    i == j ? assert(std::abs(deb(i,j) - 1) < 1e-4) : assert(std::abs(deb(i, j)) < 1e-4);
#endif

            // De-normalized Scaling SubMatrix i.e. W_k
            nsm *= std::sqrt(y_proj / lam_proj);

            scaling_matrix.PasteClippedMatrix(nsm, 0, 0, 3, 3, cs, cs, true);

            cs += 2;
            continue;
        }

        // skip UV, they should have already been considered in the previous 'if'

#ifdef DEBUG_TRAPS
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_UV)
            continue;
        
        assert(false && "An unhandled constraint has been found.");
#endif


    } // end FOR loop

    scaling_matrix.Compress(); //TODO: it might be not needed

}


void ChInteriorPoint::SetUseSymmetry(bool val)
{
	leverage_symmetry = val;
	BigMat.SetType(val ? ChSparseMatrix::SYMMETRIC_INDEF : ChSparseMatrix::GENERAL);
	scaling_matrix.SetType(val ? ChSparseMatrix::SYMMETRIC_INDEF : ChSparseMatrix::GENERAL);
	mumps_engine.SetMatrixSymmetry(val ? ChMumpsEngine::mumps_SYM::SYMMETRIC_GENERAL : ChMumpsEngine::mumps_SYM::UNSYMMETRIC);
}

void ChInteriorPoint::DumpProblem(std::string suffix) {
    
    auto folder_name("dump_" + suffix);
	CreateDirectory(folder_name.c_str(), nullptr);
    ExportArrayToFile(var.y, folder_name+"/var_y.txt");
    ExportArrayToFile(var.v, folder_name+"/var_v.txt");
    ExportArrayToFile(var.gamma, folder_name+"/var_gamma.txt");
    ExportArrayToFile(var.lambda, folder_name+"/var_lambda.txt");

    ExportArrayToFile(rhs.b_eq, folder_name+"/rhs_b_eq.txt");
    ExportArrayToFile(rhs.b_ineq, folder_name+"/rhs_b_ineq.txt");
    ExportArrayToFile(rhs.c, folder_name+"/rhs_c.txt");

    ExportArrayToFile(res.rd, folder_name + "/res_rd.txt");
    ExportArrayToFile(res.rp_gamma, folder_name + "/res_rp_gamma.txt");
    ExportArrayToFile(res.rp_lambda, folder_name + "/res_rp_lambda.txt");

    BigMat.Compress();
    CreateDirectory((folder_name + "/BigMat/").c_str(), nullptr);
    BigMat.ExportToDatFile((folder_name+"/BigMat/").c_str(), 12);

    scaling_matrix.Compress();
    CreateDirectory((folder_name + "/scaling_mat/").c_str(), nullptr);
    scaling_matrix.ExportToDatFile((folder_name + "/scaling_mat/").c_str(), 12);

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

//TODO: it could be usefule to implement multiply_and_add function like multA(vi, vo, va, ci, ca) -> vo = ci*A*vi + ca*va;
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
		if (leverage_symmetry)
			vect_out.MatrScale(-1);
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
		if (leverage_symmetry)
			vect_out.MatrScale(-1);
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
		solution_vect(n + row_sel, 0) = -var.gamma(row_sel, 0); // there will be an inversion inside FromVectorToUnknowns()
	}

	// copy 'lambda', taking care of padding the space in case we skipped the UV contact forces
    for (auto cs = 0; cs<ineq_mode.size(); ++cs)
    {
        solution_vect(n + m_eq + cs, 0) = -var.lambda(cs, 0); // there will be an inversion inside FromVectorToUnknowns()

        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            solution_vect(n + m_eq + cs + 1, 0) = 0;
            solution_vect(n + m_eq + cs + 2, 0) = 0;
            cs += 2;
        }
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

void ChInteriorPoint::PrintCVXOPTproblem(ChSystemDescriptor& sysd, std::string folder)
{
    CreateDirectory(folder.c_str(), nullptr);
    std::ofstream file_out;
    file_out.open(folder + "/cvxopt_chrono.py");

    file_out << "from cvxopt import matrix, solvers, spmatrix" << std::endl;
    file_out << "dims = { 'l': ";
    int CONSTRAINT_UNILATERAL_count = 0;
    
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            ++CONSTRAINT_UNILATERAL_count;
            continue;
        }
    }

    file_out << CONSTRAINT_UNILATERAL_count << " , 'q' : [";
    if (!skip_contacts_uv)
    {
        for (auto cs = 0; cs < ineq_mode.size(); ++cs)
        {
            if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
            {
                file_out << "3,";
                continue;
            }
        }
    }

    file_out << "], 's' : []}" << std::endl;

    
    auto func_c =  [&file_out](int r, int c, double v) {file_out << c << "," ; };
    auto func_v =  [&file_out](int r, int c, double v) {file_out << v << "," ; };
    auto func_mv = [&file_out](int r, int c, double v) {file_out << -v << ","; };

    // Matrices
    file_out << "P = spmatrix([";
    auto func_rP = [&file_out](int r, int c, double v) {file_out << r << ","; };
    BigMat.ForEachExistentValueInRange(func_v, 0, n - 1, 0, n - 1);
    file_out << "], [";
    BigMat.ForEachExistentValueInRange(func_rP, 0, n-1, 0, n-1);
    file_out << "], [";
    BigMat.ForEachExistentValueInRange(func_c, 0, n-1, 0, n-1);
    file_out << "], ("<< n <<","<< n <<"))" << std::endl;

    if (m_eq > 0)
    {
        file_out << "A = spmatrix([";
        auto func_rA = [&file_out, this](int r, int c, double v) {file_out << r - n << ","; };
        BigMat.ForEachExistentValueInRange(func_v, n, n + m_eq - 1, 0, n);
        file_out << "], [";
        BigMat.ForEachExistentValueInRange(func_rA, n, n + m_eq - 1, 0, n);
        file_out << "], [";
        BigMat.ForEachExistentValueInRange(func_c, n, n + m_eq - 1, 0, n);
        file_out << "], (" << m_eq << "," << n << "))" << std::endl;
    }

    file_out << "G = spmatrix([";
    auto func_rG = [&file_out, this](int r, int c, double v) {file_out << r - n - m_eq << ","; };
    BigMat.ForEachExistentValueInRange(func_mv, n + m_eq, n + m_eq + m_ineq - 1, 0, n);
    file_out << "], [";
    BigMat.ForEachExistentValueInRange(func_rG, n + m_eq, n + m_eq + m_ineq - 1, 0, n);
    file_out << "], [";
    BigMat.ForEachExistentValueInRange(func_c, n + m_eq, n + m_eq + m_ineq - 1, 0, n);
    file_out << "], (" << m_ineq << "," << n << "))" << std::endl;

    // Vectors
    file_out << "c = matrix([";
    for (auto i = 0; i<n; ++i)
        file_out << rhs.c(i) << ",";
    file_out << "])" << std::endl;

    if (m_eq > 0)
    {
        file_out << "b_eq = matrix([";
        for (auto i = 0; i<m_eq; ++i)
            file_out << rhs.b_eq(i) << ",";
        file_out << "])" << std::endl;
    }
    

    file_out << "b_ineq = matrix([";
    for (auto i = 0; i<m_ineq; ++i)
        file_out << -rhs.b_ineq(i) << ",";
    file_out << "])" << std::endl;

    if (m_eq==0)
        file_out << "sol = solvers.coneqp(P, c, G, b_ineq, dims)" << std::endl;
    else
        file_out << "sol = solvers.coneqp(P, c, G, b_ineq, dims, A, b_eq)" << std::endl;

    // file_out << "print(sol['x'])" << std::endl;


    file_out.close();


}


int ChInteriorPoint::belongs_to_cone(const ChMatrix<double>& pos, int offset)
{
    if (pos(offset) > 0 && uTJv(pos, pos, offset) > 0)
        return 1;

    if (pos(offset) >= 0 && uTJv(pos, pos, offset) >= 0)
        return 0;

    return -1;
}

int ChInteriorPoint::belongs_to_cone(const ChMatrix<double>& pos, const ChMatrix<double>& dir, double alfa, int offset)
{
    ChMatrixNM<double, 3, 1> pos_temp;
    pos_temp(0) = pos(offset) + alfa*dir(offset);
    pos_temp(1) = pos(offset+1) + alfa*dir(offset+1);
    pos_temp(2) = pos(offset+2) + alfa*dir(offset+2);

    if (pos_temp(0) > 0 && uTJv(pos_temp, pos_temp, 0) > 0)
        return 1;

    if (pos_temp(0) >= 0 && uTJv(pos_temp, pos_temp, 0) >= 0)
        return 0;

    return -1;
}

int ChInteriorPoint::find_out_of_cone(const ChMatrix<double>& pos, const ChMatrix<double>& dir, double alfa)
{
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            if (pos(cs) + alfa*dir(cs) < 0)
                return cs;
            continue;
        }


        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
        {
            ChMatrixNM<double, 3, 1> pos_temp;
            pos_temp(0) = pos(cs) + alfa*dir(cs);
            pos_temp(1) = pos(cs + 1) + alfa*dir(cs + 1);
            pos_temp(2) = pos(cs + 2) + alfa*dir(cs + 2);

            if (pos_temp(0) < 0 || uTJv(pos_temp, pos_temp, 0) < 0)
                return cs;

        }
    }

    return -1;

}

int ChInteriorPoint::find_out_of_cone(const ChMatrix<double>& pos)
{
    for (auto cs = 0; cs < ineq_mode.size(); ++cs)
    {
        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_UNILATERAL || ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N && skip_contacts_uv)
        {
            if (pos(cs) < 0)
                return cs;

            continue;
        }

        if (ineq_mode[cs] == eChConstraintModeMOD::CONSTRAINT_FRIC_N)
            if (pos(cs) < 0 || uTJv(pos, pos, cs) < 0)
                return cs;

    }

    return -1;
}

    ChInteriorPoint::ChInteriorPoint() {
    mumps_engine.SetICNTL(11, 2);
    RecordHistory(true);
    BigMat.SetSparsityPatternLock(m_lock);
}

ChInteriorPoint::~ChInteriorPoint() {
    if( logfile_stream.is_open() )
        logfile_stream.close();
}


void ChInteriorPoint::PrintIPStatus() const
{
	GetLog() << "\nInteriorPoint | Call: " << solver_call << " Iter: " << iteration_count << "/" << iteration_count_max << "\n";
	GetLog() << "Complementarity Measure: " << res.mu << "\n";
	GetLog() << "|rd|/n (stationarity): " << res.rd.NormTwo() / n << "\n";
	GetLog() << "|rp_gamma|/m_eq (BILAT violation): " << (m_eq == 0 ? 0 : res.rp_gamma.NormTwo() / m_eq) << "\n";
	GetLog() << "|rp_lambda|/m_ineq (UNILAT violation): " << res.rp_lambda.NormTwo() / m_ineq << "\n";
	GetLog() << "Objective Function: " << evaluate_objective_function() << "\n";
	GetLog() << "\n";
}


void ChInteriorPoint::output_log()
{
    if (logfile_stream.is_open())
    {
        logfile_stream << std::endl << res.rd.NormTwo() / n << ", " << res.rp_gamma.NormTwo() / m_eq << ", " << res.rp_lambda.NormTwo() / m_ineq << ", " << res.mu << ", " << evaluate_objective_function();
    }
    else
    {
        logfile_stream.open(logfile_name + ".txt");

        if (logfile_stream.is_open())
        {
            logfile_stream << std::scientific << std::setprecision(6);
            logfile_stream << "r_{D} r_{P_{\\gamma}} r_{P_{\\lambda}} \\mu objfun";

            logfile_stream << std::endl << res.rd.NormTwo() / n << ", " << res.rp_gamma.NormTwo() / m_eq << ", " << res.rp_lambda.NormTwo() / m_ineq << ", " << res.mu << ", " << evaluate_objective_function();

#ifdef CHRONO_POSTPROCESS
            postprocess::ChGnuPlot mplot(("__" + logfile_name + ".gpl").c_str());
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
} // end namespace chrono
