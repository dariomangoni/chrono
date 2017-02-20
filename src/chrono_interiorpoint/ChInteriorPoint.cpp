#include "ChInteriorPoint.h"
#include <algorithm>


#define DEBUG_MODE true
#define SKIP_CONTACTS_UV true
#define ADD_COMPLIANCE false
#define REUSE_OLD_SOLUTIONS false

namespace chrono
{

	template<class ChMatrixIN>
	void PrintMatrix(ChMatrixIN& matrice){
		for (size_t i = 0; i < matrice.GetRows(); i++){
			for (size_t j = 0; j < matrice.GetColumns(); j++){
				printf("%.1f ", matrice.GetElement(i, j));
			}
			printf("\n");
		}
	}

	template <class matrix>
	void ExportArrayToFile(matrix mat, std::string filepath, int precision = 12)
	{
		std::ofstream ofile;
		ofile.open(filepath);
		ofile << std::scientific << std::setprecision(precision);

		for (size_t row_sel = 0; row_sel < mat.GetRows(); row_sel++)
		{
			for (size_t col_sel = 0; col_sel < mat.GetColumns(); col_sel++)
			{
				ofile << mat.GetElement(row_sel, col_sel);
			}

			ofile << std::endl;
		}

		ofile.close();
	}

	void ImportArrayFromFile(ChMatrix<>& output_mat, std::string filename)
	{
		std::ifstream my_file;
		my_file.open(filename);

		double temp;
		int row_sel = -1;
		for (row_sel = 0; row_sel < output_mat.GetRows(); row_sel++)
		{
			my_file >> temp;
			output_mat.SetElement(row_sel, 0, temp);
		}
		my_file.close();
	}
    // Chrono global matrix //TODO: check
    // | M   Cq'|*|q|- |f|= |0|
    // | Cq  -E | |l|  |b|  |c| 
    // after conversion
	// | M  -Cq'|*|q|- | f|= |0|
	// | Cq  -E | |l|  |-b|  |c| 

	double ChInteriorPoint::Solve(ChSystemDescriptor& sysd)
	{
		solver_call++;
        verbose = true;

        // The problem to be solved is loaded into the main matrix that will be used to solve the various step of the IP method
        // The initial guess is modified in order to be feasible
        // The residuals are computed

        /********** Load system **********/
        // TODO: dimensions generally change at each call; 'm' for sure, but maybe 'n' does not?
        auto n_old = n;
        auto m_old = m;
        double n_new = sysd.CountActiveVariables();
        double m_new = sysd.CountActiveConstraints(false, SKIP_CONTACTS_UV);
        reset_dimensions(n_new, m_new);

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
        switch (KKT_solve_method)
        {
            case STANDARD:
                sysd.ConvertToMatrixForm(&BigMat, nullptr, false, SKIP_CONTACTS_UV, 2);
                make_positive_definite();
                break;
            case AUGMENTED:
                sysd.ConvertToMatrixForm(&BigMat, nullptr, false, SKIP_CONTACTS_UV, 1);
                make_positive_definite();
                break;
            case NORMAL:
                sysd.ConvertToMatrixForm(&SmallMat, &BigMat, nullptr, nullptr, nullptr, nullptr, false, true);
                break;
        }

        sysd.ConvertToMatrixForm(nullptr, nullptr, nullptr, &c, &b, nullptr, false, SKIP_CONTACTS_UV); // load f->c and b->b
        c.MatrScale(-1); // adapt to InteriorPoint convention
        b.MatrScale(-1); // adapt to InteriorPoint convention


        /********** Check if constraints are found **********/
        if (m == 0) // if no constraints
        {
            // Fill 'rhs_sol' with [-rd;-rp-y-sigma*mu/lam]
            for (size_t row_sel = 0; row_sel < n; row_sel++)
                rhs_sol.SetElement(row_sel, 0, -c.GetElement(row_sel, 0));

            // Solve the KKT system
            BigMat.Compress();
            mumps_engine.SetProblem(BigMat, rhs_sol);
            printf("Mumps says: %d\n", mumps_engine.MumpsCall(ChMumpsEngine::COMPLETE));

            if (verbose)
            {
                std::cout << "Scaled residual norm of MUMPS call: " << mumps_engine.GetRINFOG(6) << std::endl;
            }

            std::cout << "IP call: " << solver_call << "; Switched to MUMPS because no constraints found" << std::endl;
            sysd.FromVectorToUnknowns(rhs_sol);
            return 0;

        }

        /********* the system DOES have contacts! ********/

        if (ADD_COMPLIANCE && m > 0)
        {
            sysd.ConvertToMatrixForm(nullptr, nullptr, &E, nullptr, nullptr, nullptr, false, SKIP_CONTACTS_UV);
            E *= -1;
        }

        starting_point_Nocedal(n_old, m_old);
        //starting_point_Nocedal_WS();
        //starting_point_STP1();
        //starting_point_STP2(n_old);


		SetTolerances(1e-7, 1e-8, 1e-8);

		for (iteration_count = 1; iteration_count < iteration_count_max && !check_exit_conditions(); iteration_count++)
		{
			iterate();
			update_history();
		}


		if (verbose)
		{
			std::cout << "IP call: " << solver_call << "; Iter: " << iteration_count << "/" << iteration_count_max << std::endl;
			VerifyKKTconditions(true);
		}


		generate_solution();
		sysd.FromVectorToUnknowns(sol_chrono); // the function expects [x;-l]


		return 0;
	}


	void ChInteriorPoint::starting_point_STP1() // from [2]
	{
		double infeas_dual_ratio = 0.1; // TODO: dependant on n
		
		x.FillElem(1);
		y.FillElem(1);
		lam.FillElem(1);
		
		double duality_gap_calc = y.MatrDot(y, lam); // [2] pag. 132
		double duality_gap = m; // [2] pag. 132
		assert(duality_gap_calc == duality_gap);

		
		// norm of all residuals; [2] pag. 132
		fullupdate_residual();
		double res_norm = rp.MatrDot(rp, rp); 
		res_norm += rp.MatrDot(rd, rd);
		res_norm = sqrt(res_norm);

		if (res_norm / duality_gap > infeas_dual_ratio)
		{
			double coeff = res_norm / (duality_gap * infeas_dual_ratio);
			x.MatrScale(coeff);
			y.MatrScale(coeff);
			lam.MatrScale(coeff);

			fullupdate_residual();
		}

	}

	void ChInteriorPoint::starting_point_STP2(int n_old)
	{
		double threshold = 1; // 'epsilon' in [2]

		if (n != n_old && !REUSE_OLD_SOLUTIONS)
		{
			// initialize x
			x.FillElem(1);
		}
		

		// initialize y and then lam
		multiplyA(x, vectm);
		vectm -= b;
		for (size_t cont = 0; cont < m; cont++)
		{
			y(cont, 0) = (vectm(cont, 0) > threshold) ? vectm(cont, 0) : threshold;
			lam(cont, 0) = 1 / y(cont, 0);
		}
		
		fullupdate_residual();

	}

	void ChInteriorPoint::starting_point_Nocedal_WS(int n_old, int m_old)
	{
		/*Backup vectors*/
		ChMatrixDynamic<double> x_bkp(x);
		ChMatrixDynamic<double> y_bkp(y);
		ChMatrixDynamic<double> lam_bkp(lam);
		fullupdate_residual();
		VerifyKKTconditions();
		double residual_value_bkp = IPres.rp_nnorm*m + IPres.rd_nnorm*n + mu*m;
		

		/********** Initialize IP algorithm **********/
		// Initial guess
		if (n_old != n || solver_call == 0)
			x.FillElem(1); // TIP: every ChMatrix is initialized with zeros by default
		if (m_old != m || solver_call == 0)
			lam.FillElem(1); // each element of lam will be at the denominator; avoid zeros!

		// since A is generally changed between calls, also with warm_start,
		// all the residuals and feasibility check must be redone
		multiplyA(x, y);  // y = A*x
		y -= b;

		// Calculate the residual
		fullupdate_residual();

		// Feasible starting Point (pag.484-485)
		KKTsolve(0, false); // to obtain Dx, Dy, Dlam called "affine"

		// x is accepted as it is
		y += Dy; // calculate y0
		lam += Dlam; // calculate lam0

		for (size_t row_sel = 0; row_sel < m; row_sel++)
			y(row_sel) = abs(y(row_sel)) < 1 ? 1 : abs(y(row_sel));

		for (size_t row_sel = 0; row_sel < m; row_sel++)
			lam(row_sel) = abs(lam(row_sel)) < 1 ? 1 : abs(lam(row_sel));

		// Update the residual considering the new values of 'y' and 'lam'
		fullupdate_residual();


		/* Check if restoring previous values would be better */
		VerifyKKTconditions();
		double residual_value_new = IPres.rp_nnorm*m + IPres.rd_nnorm*n + mu*m;

		if (residual_value_bkp<residual_value_new)
		{
			x = x_bkp;
			y = y_bkp;
			lam = lam_bkp;
			fullupdate_residual();
			VerifyKKTconditions();
			double residual_value_final = IPres.rp_nnorm*m + IPres.rd_nnorm*n + mu*m;
		}
		else
		{
			std::cout << "Not WS\n";
		}

	}

	void ChInteriorPoint::starting_point_Nocedal(int n_old, int m_old)
	{

		/********** Initialize IP algorithm **********/
		// Initial guess
		if (n_old != n || solver_call == 0 || !REUSE_OLD_SOLUTIONS)
			x.FillElem(1); // TIP: every ChMatrix is initialized with zeros by default
		if (m_old != m || solver_call == 0 || !REUSE_OLD_SOLUTIONS)
			lam.FillElem(1); // each element of lam will be at the denominator; avoid zeros!

		// since A is generally changed between calls, also with warm_start,
		// all the residuals and feasibility check must be redone
		multiplyA(x, y);  // y = A*x
		y -= b;

		// Calculate the residual
 		fullupdate_residual();

		// Feasible starting Point (pag.484-485)
		KKTsolve(0, false); // to obtain Dx, Dy, Dlam called "affine"

		// x is accepted as it is
		y += Dy; // calculate y0
		lam += Dlam; // calculate lam0

		for (size_t row_sel = 0; row_sel < m; row_sel++)
			y(row_sel) = abs(y(row_sel)) < 1 ? 1 : abs(y(row_sel));

		for (size_t row_sel = 0; row_sel < m; row_sel++)
			lam(row_sel) = abs(lam(row_sel)) < 1 ? 1 : abs(lam(row_sel));

		// Update the residual considering the new values of 'y' and 'lam'
		fullupdate_residual();

	}

	// Iterating function
	// output: (x, y, lam) are computed
	// (rp, rd, mu) are updated based on most recent (x, y, lam)

	// rp, rd, mu, x, y, lam are taken as they are
	void ChInteriorPoint::iterate()
	{
        /*********************************************************************************/
        /***************************** Prediction Phase **********************************/
        /*********************************************************************************/
#ifdef DEBUG_MODE
        GetLog() << "\n\n\n\n\n/*********************************************************************************/\n";
        GetLog() << "/******************* IP call: " << solver_call << "; iter: " << iteration_count << " ************************/\n";
        GetLog() << "/*********************************************************************************/\n";
#endif

#ifdef DEBUG_MODE
        GetLog() << "n: "<< n << " | m: " << m << "\n";
        GetLog() << "Starting with values:\n";
        GetLog() << "x: " << x << "\n";
        GetLog() << "y: " << y << "\n";
        GetLog() << "lambda: " << lam << "\n";
        GetLog() << "ResidualNorm primal: " << rp.NormTwo() << "\n";
        GetLog() << "ResidualNorm dual: " << rd.NormTwo() << "\n";
        GetLog() << "\n";
#endif

        /*** find directions ***/
		KKTsolve(0, false); // to obtain Dx, Dy, Dlam called "affine" aka "prediction"

#ifdef DEBUG_MODE
        GetLog() << "Direction (prediction):\n";
        GetLog() << "Dx: " << Dx << "\n";
        GetLog() << "Dy: " << Dy << "\n";
        GetLog() << "Dlambda: " << Dlam << "\n";
        GetLog() << "\n";
#endif

        /*** compute step lengths ***/
		// from 16.60 pag.482 from 14.32 pag.408 (remember that y>=0!)
		double alfa_pred_prim = find_Newton_step_length(y, Dy);
		double alfa_pred_dual = find_Newton_step_length(lam, Dlam);

		if (EQUAL_STEP_LENGTH)
		{
			double alfa_pred = std::min(alfa_pred_prim, alfa_pred_dual);
			alfa_pred_prim = alfa_pred;
			alfa_pred_dual = alfa_pred;
		}

#ifdef DEBUG_MODE
        GetLog() << "Step length 'alfa' (prediction):\n";
        GetLog() << "| primal: " << alfa_pred_prim << "\n";
        GetLog() << "| dual: " << alfa_pred_dual << "\n";
        GetLog() << "\n";
#endif

        /*** make the prediction step ***/
		  y_pred = Dy;     y_pred.MatrScale(alfa_pred_prim);   y_pred += y;
		lam_pred = Dlam; lam_pred.MatrScale(alfa_pred_dual); lam_pred += lam;

#ifdef DEBUG_MODE
        GetLog() << "Values (prediction):\n";
        GetLog() << "x_pred: " << x_pred << "\n";
        GetLog() << "y_pred: " << y_pred << "\n"; //TODO: is unchanged in prediction step???
        GetLog() << "lambda_pred: " << lam_pred << "\n";
#endif

        /*** compute complementarity measure ***/
		double mu_pred = y_pred.MatrDot(y_pred, lam_pred) / m; // from 14.33 pag.408 //TODO: why MatrDot is a member?

#ifdef DEBUG_MODE
        GetLog() << "Complementarity measure (prediction): " << mu_pred << "\n";
        GetLog() << "\n";
#endif

		if (ONLY_PREDICT)
		{
			x_pred = Dx; x_pred.MatrScale(alfa_pred_prim); x_pred += x; x = x_pred;
			y = y_pred;
			lam = lam_pred;
			
			rp.MatrScale(1 - alfa_pred_prim);
			
			multiplyG(Dx, vectn); // vectn = G * Dx
			vectn.MatrScale(alfa_pred_prim - alfa_pred_dual); // vectn = (alfa_pred_prim - alfa_pred_dual) * (G * Dx)
			rd.MatrScale(1 - alfa_pred_dual);
			rd += vectn;

			mu = mu_pred;

#ifdef DEBUG_MODE
            GetLog() << "InteriorPoint Results only_pred\n";
            GetLog() << "Complementarity measure: " << mu << "\n";
            GetLog() << "ResidualNorm primal: " << rp.NormTwo() << "\n";
            GetLog() << "ResidualNorm dual: " << rd.NormTwo() << "\n";
            GetLog() << "\n";
#endif

            return;
		}

		
        /*********************************************************************************/
		/******************************* Correction phase ********************************/
		/*********************************************************************************/

        /*** evaluate centering parameter ***/
		double sigma = std::pow(mu_pred / mu, 3.0); // from 14.34 pag.408

        /*** find directions ***/
		KKTsolve(sigma, true); // in 'Dlam' and 'Dy' there must be 'Dlam_pred' and 'Dy_pred'

#ifdef DEBUG_MODE
        GetLog() << "Direction (correction):\n";
        GetLog() << "centering param 'sigma': " << sigma << "\n";
        GetLog() << "Dx: " << Dx << "\n";
        GetLog() << "Dy: " << Dy << "\n";
        GetLog() << "Dlambda: " << Dlam << "\n";
        GetLog() << "\n";
#endif

        /*** step length correction ***/
		double eta = ADAPTIVE_ETA ? exp(-mu*m)*0.1 + 0.9 : 0.95; // exponential descent of eta

        /*** compute step lengths ***/
		double alfa_corr_prim = find_Newton_step_length(y, Dy, eta);
		double alfa_corr_dual = find_Newton_step_length(lam, Dlam, eta);

		if (EQUAL_STEP_LENGTH)
		{
			double alfa_corr = std::min(alfa_corr_prim, alfa_corr_dual);
			alfa_corr_prim = alfa_corr;
			alfa_corr_dual = alfa_corr;
		}

#ifdef DEBUG_MODE
        GetLog() << "Step length 'alfa' (correction):\n";
        GetLog() << "| primal: " << alfa_corr_prim << "\n";
        GetLog() << "| dual: " << alfa_corr_dual << "\n";
        GetLog() << "\n";
#endif

        /*** make the correction step ***/
		  x_corr = Dx;		  x_corr.MatrScale(alfa_corr_prim);		  x_corr += x;		  x = x_corr;
		  y_corr = Dy;		  y_corr.MatrScale(alfa_corr_prim);		  y_corr += y;		  y = y_corr;
		lam_corr = Dlam;	lam_corr.MatrScale(alfa_corr_dual);		lam_corr += lam;	lam = lam_corr;

#ifdef DEBUG_MODE
        GetLog() << "Values (correction):\n";
        GetLog() << "x: " << x << "\n";
        GetLog() << "y: " << y << "\n";
        GetLog() << "lambda: " << lam << "\n";
        GetLog() << "\n";
#endif

		/********** Residuals update **********/
		rp.MatrScale(1 - alfa_corr_prim);
		rd.MatrScale(1 - alfa_corr_dual);
		mu = y.MatrDot(y, lam) / m; // from 14.6 pag.395
		
		if (!EQUAL_STEP_LENGTH)
		{
			multiplyG(Dx, vectn); // vectn = G*Dx
			vectn.MatrScale(alfa_corr_prim - alfa_corr_dual); // vectn = (alfa_pred_prim - alfa_pred_dual) * (G * Dx)
			rd += vectn;
		}



#ifdef DEBUG_MODE
        GetLog() << "InteriorPoint Results pred+corr\n";
        GetLog() << "Complementarity measure: " << mu << "\n";
        GetLog() << "ResidualNorm primal: " << rp.NormTwo() << "\n";
        GetLog() << "ResidualNorm dual: " << rd.NormTwo() << "\n";
        GetLog() << "\n";
#endif

        DumpProblem();

	}

	// Solve the KKT system in different modes: 'rp', 'rd', 'mu', 'x', 'y', 'lam' must be updated before calling this function
    // 'sigma' is the centering parameter;
	void ChInteriorPoint::KKTsolve(double sigma, bool apply_correction)
	{

		switch (KKT_solve_method)
		{
		case STANDARD:
            assert(false && "Standard mode to be fixed");
			// update lambda and y diagonal submatrices
			for (size_t diag_sel = 0; diag_sel < m; diag_sel++)
			{
				BigMat.SetElement(n + m + diag_sel, n + diag_sel, lam.GetElement(diag_sel, 0)); // write lambda diagonal submatrix
				BigMat.SetElement(n + m + diag_sel, n + m + diag_sel, y.GetElement(diag_sel, 0)); // write y diagonal submatrix
				BigMat.SetElement(n + diag_sel, n + diag_sel, -1); // write -identy_matrix diagonal submatrix
			}

			if (sigma != 0) // rpd_corr
			{
				// I'm supposing that in 'rpd', since the previous call should have been without perturbation,
				// there is already y°lam
				vectm = Dlam; // I could use Dlam directly, but it is not really clear
				vectm.MatrScale(Dy);
				vectm.MatrInc(-sigma*mu);
				rpd += vectm;
			}
			else // rpd_pred as (16.57 pag.481 suggests)
			{
				rpd = y;
				rpd.MatrScale(lam);
			}

			// Fill 'rhs_sol' with [-rd;-rp;-rpd]
			for (size_t row_sel = 0; row_sel < n; row_sel++)
				rhs_sol.SetElement(row_sel, 0, -rd.GetElement(row_sel, 0));

			for (size_t row_sel = 0; row_sel < m; row_sel++)
			{
				rhs_sol.SetElement(row_sel + n, 0, -rp.GetElement(row_sel, 0));
				rhs_sol.SetElement(row_sel + n + m, 0, -rpd.GetElement(row_sel, 0));
			}


			// Solve the KKT system
            BigMat.Compress();
			mumps_engine.SetProblem(BigMat, rhs_sol);
			printf("Mumps says: %d\n", mumps_engine.MumpsCall(ChMumpsEngine::SOLVE));

			// Extract 'Dx', 'Dy' and 'Dlam' from 'sol'
			for (size_t row_sel = 0; row_sel < n; row_sel++)
				Dx.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel, 0));

			for (size_t row_sel = 0; row_sel < m; row_sel++)
			{
				Dy.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n, 0));
				Dlam.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n + m, 0));
			}

			break;
		case AUGMENTED:
            if (!apply_correction)
            {
                // update y/lambda diagonal submatrix
                if (ADD_COMPLIANCE)
                    for (size_t diag_sel = 0; diag_sel < m; diag_sel++)
                    {
                        BigMat.SetElement(n + diag_sel, n + diag_sel, y.GetElement(diag_sel, 0) / lam.GetElement(diag_sel, 0) + E.GetElement(diag_sel, diag_sel));
                    }
                else
                    for (size_t diag_sel = 0; diag_sel < m; diag_sel++)
                    {
                        BigMat.SetElement(n + diag_sel, n + diag_sel, y.GetElement(diag_sel, 0) / lam.GetElement(diag_sel, 0));
                    }

                BigMat.Compress();
                mumps_engine.SetProblem(BigMat, rhs_sol);
                mumps_engine.MumpsCall(ChMumpsEngine::ANALYZE_FACTORIZE);

                // fill 'rhs_sol' with rhs [-rd;-rp-y+sigma*mu/lam]
                for (size_t row_sel = 0; row_sel < n; row_sel++)
                    rhs_sol.SetElement(row_sel, 0, -rd.GetElement(row_sel, 0));
            }


			for (size_t row_sel = 0; row_sel < m; row_sel++)
				rhs_sol.SetElement(row_sel + n, 0, -rp(row_sel, 0) - y(row_sel, 0) + (apply_correction ? sigma*mu - Dlam(row_sel,0)*Dy(row_sel, 0) : sigma*mu) / lam(row_sel, 0));

			//ExportArrayToFile(rhs_sol, "dump/rhs.txt");

			// Solve the KKT system
			mumps_engine.MumpsCall(ChMumpsEngine::SOLVE);
            mumps_engine.PrintINFOG();
            if (verbose)
            {
                double res_norm = mumps_engine.GetRINFOG(6);
                if (res_norm > 1e-6)
                    std::cout << "Scaled residual norm of MUMPS call: " << res_norm << std::endl;
            }

			//BigMat.ExportToDatFile("dump/COO.txt", true);
			//BigMat.ImportFromDatFile("COO.txt", true);
			//BigMat.ExportToDatFile("COO.txt", true);

			//ExportArrayToFile(rhs_sol, "dump/sol.txt");


			// Extract 'Dx' and 'Dlam' from 'sol'
			for (size_t row_sel = 0; row_sel < n; row_sel++)
				Dx.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel, 0));
			for (size_t row_sel = 0; row_sel < m; row_sel++)
				Dlam.SetElement(row_sel, 0, rhs_sol.GetElement(row_sel + n, 0));

			// Calc 'Dy' (it is also possible to evaluate Dy as Dy=(-lam°y+sigma*mu*e-y°Dlam)./lam )
			multiplyA(Dx, Dy);  // Dy = A*Dx
			Dy += rp;
			if (ADD_COMPLIANCE)
			{
				E.MatrMultiply(Dlam, vectm);
				Dy += vectm;
			}
				

			break;
		case NORMAL:
            assert(false && "Normal mode to be fixed");
			for (size_t row_sel = 0; row_sel < n; row_sel++)
			{
				for (size_t col_sel = 0; col_sel < n; col_sel++)
				{
					int temp = 0;
					for (size_t el_sel = 0; el_sel < m; el_sel++)
					{
						temp += lam(el_sel, 0) / y(el_sel, 0) * SmallMat.GetElement(el_sel, row_sel) * SmallMat.GetElement(el_sel, col_sel);
						if (temp != 0)
							BigMat.SetElement(row_sel, col_sel, temp, false);
					}
				}
			}

			break;
		}
	}

	// this function obtain the maximum length, along the direction defined by Dvect, so that vect has no negative components;
	// it will be applied to the two vectors that are forced to be non-negative i.e. 'lam' and 'y'
	double ChInteriorPoint::find_Newton_step_length(const ChMatrix<double>& vect, const ChMatrix<double>& Dvect, double eta )
	{
		double alpha = 1;
		for (size_t row_sel = 0; row_sel < vect.GetRows(); row_sel++)
		{
			if (Dvect(row_sel,0)<0)
			{
				double alfa_temp = -eta * vect(row_sel,0) / Dvect(row_sel,0);
				if (alfa_temp < alpha)
					alpha = alfa_temp;
			}
		}

		return alpha>0 ? alpha : 0;
	}

	double ChInteriorPoint::evaluate_objective_function()
	{
		multiplyG(x, vectn);
		double obj_value = vectn.MatrDot(x, vectn);
		obj_value += c.MatrDot(x, c);

		return obj_value;
	}


    /// Take care of adapting the size of matrices to \p n_new and \p m_new
    void ChInteriorPoint::reset_dimensions(int n_new, int m_new)
    {
        if (n != n_new)
        {
            x.Resize(n_new, 1);
            x_pred.Resize(n_new, 1);
            x_corr.Resize(n_new, 1);
            Dx.Resize(n_new, 1);
            c.Resize(n_new, 1);
            rd.Resize(n_new, 1);
            vectn.Resize(n_new, 1);
        }

        if (m != m_new)
        {
            y.Resize(m_new, 1);
            lam.Resize(m_new, 1);
            y_pred.Resize(m_new, 1);
            lam_pred.Resize(m_new, 1);
            y_corr.Resize(m_new, 1);
            lam_corr.Resize(m_new, 1);
            Dy.Resize(m_new, 1);
            Dlam.Resize(m_new, 1);
            b.Resize(m_new, 1);
            rp.Resize(m_new, 1);
            rpd.Resize(m_new, 1);
            vectm.Resize(m_new, 1);
        }

        SKIP_CONTACTS_UV ? sol_chrono.Resize(n_new + 3 * m_new, 1) : sol_chrono.Resize(n_new + m_new, 1);

        // BigMat and sol
        switch (KKT_solve_method)
        {
            case STANDARD:
                BigMat.Reset(2 * m_new + n_new, 2 * m_new + n_new, static_cast<int>(n_new*n_new*SPM_DEF_FULLNESS));
                rhs_sol.Resize(2 * m_new + n_new, 1);
                break;
            case AUGMENTED:
                BigMat.Reset(n_new + m_new, n_new + m_new, static_cast<int>(n_new*n_new*SPM_DEF_FULLNESS));
                rhs_sol.Resize(n_new + m_new, 1);
                break;
            case NORMAL:
                std::cout << std::endl << "Perturbed KKT system cannot be stored with 'NORMAL' method yet.";
                break;
        }

        n = n_new;
        m = m_new;

    }


	void ChInteriorPoint::DumpProblem(std::string suffix)
	{
		ExportArrayToFile(y, "dump/y" + suffix + ".txt");
		ExportArrayToFile(x, "dump/x" + suffix + ".txt");
		ExportArrayToFile(lam, "dump/lam" + suffix + ".txt");

		ExportArrayToFile(b, "dump/b" + suffix + ".txt");
		ExportArrayToFile(c, "dump/c" + suffix + ".txt");

		BigMat.Compress();
		BigMat.ExportToDatFile("dump/", 3);
	}

	void ChInteriorPoint::LoadProblem()
	{
		//ImportArrayFromFile(y, "dump/y.txt");
		//ImportArrayFromFile(x, "dump/x.txt");
		//ImportArrayFromFile(lam, "dump/lam.txt");

		ImportArrayFromFile(b, "dump/b.txt");
		ImportArrayFromFile(c, "dump/c.txt");

		//BigMat.ImportFromDatFile("dump/");
	}

	void ChInteriorPoint::DumpIPStatus(std::string suffix) const
	{
		ExportArrayToFile(y, "dump/y" + suffix + ".txt");
		ExportArrayToFile(x, "dump/x" + suffix + ".txt");
		ExportArrayToFile(lam, "dump/lam" + suffix + ".txt");

		ExportArrayToFile(Dx, "dump/Dx" + suffix + ".txt");
		ExportArrayToFile(Dy, "dump/Dy" + suffix + ".txt");
		ExportArrayToFile(Dlam, "dump/Dlam" + suffix + ".txt");

		ExportArrayToFile(rhs_sol, "dump/rhs_sol" + suffix + ".txt");
		//ExportArrayToFile(sol, "dump/sol" + suffix + ".txt");
	}

	void ChInteriorPoint::make_positive_definite()
	{
        if (m == 0)
            return;

		int offset_AT_col = n;
		if (KKT_solve_method == STANDARD)
			offset_AT_col = n + m;

        BigMat.ForEachExistentValueInRange([](double* val) { *val *= -1; }, 0, n, offset_AT_col, BigMat.GetNumColumns()-1);

	}

	// Perform moltiplication of A with vect_in: vect_out = A*vect_in
	void ChInteriorPoint::multiplyA(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const
	{
		switch (KKT_solve_method)
		{
		case STANDARD:
			BigMat.MatrMultiplyClipped(vect_in, vect_out, n, n + m - 1, 0, n - 1, 0, 0);
			break;
		case AUGMENTED:
			BigMat.MatrMultiplyClipped(vect_in, vect_out, n, n + m - 1, 0, n - 1, 0, 0);
			break;
		case NORMAL:
			std::cout << std::endl << "A multiplication is not implemented in 'NORMAL' method yet.";
			break;
		}
	}

	// Perform moltiplication of -AT with vect_in: vect_out = -AT*vect_in
	void ChInteriorPoint::multiplyNegAT(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const
	{
		switch (KKT_solve_method)
		{
		case STANDARD:
			BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, n + m, n + 2*m - 1, 0, 0);
			break;
		case AUGMENTED:
			BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, n, n + m - 1, 0, 0);
			break;
		case NORMAL:
			std::cout << std::endl << "AT multiplication is not implemented in 'NORMAL' method yet.";
			break;
		}
	}

	void ChInteriorPoint::multiplyG(const ChMatrix<double>& vect_in, ChMatrix<double>& vect_out) const
	{
		switch (KKT_solve_method)
		{
		case STANDARD:
			BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, 0, n - 1, 0, 0);
			break;
		case AUGMENTED:
			BigMat.MatrMultiplyClipped(vect_in, vect_out, 0, n - 1, 0, n - 1, 0, 0);
			break;
		case NORMAL:
			std::cout << std::endl << "G multiplication is not implemented in 'NORMAL' method yet.";
			break;
		}
	}

	void ChInteriorPoint::normalize_Arows()
	{
		
		int offset_AT_cols = 0;
		switch (KKT_solve_method)
		{
		case STANDARD:
			offset_AT_cols = n;
			break;
		case AUGMENTED:
			offset_AT_cols = n + m;
			break;
		case NORMAL:
			std::cout << std::endl << "G multiplication is not implemented in 'NORMAL' method yet.";
			break;
		}

		double row_norm;

		auto* a_mat = BigMat.GetCOO_ValuesAddress();
		auto* ia_mat = BigMat.GetCOO_RowIndexAddress();
		auto* ja_mat = BigMat.GetCOO_ColIndexAddress();

		for (size_t row_sel = n; row_sel < n + m; row_sel++)
		{
			// Calc row norm from A in BigMat
			row_norm = 0;
			for (size_t col_sel = ia_mat[row_sel]; ja_mat[col_sel] < n && col_sel < ia_mat[row_sel + 1]; col_sel++)
			{
				row_norm += a_mat[col_sel] * a_mat[col_sel];
			}
			row_norm = sqrt(row_norm);

			// Normalize A and rewrite -AT
			for (size_t col_sel = ia_mat[row_sel]; ja_mat[col_sel] < n && col_sel < ia_mat[row_sel + 1]; col_sel++)
			{
				a_mat[col_sel] = a_mat[col_sel] / row_norm;
				BigMat.SetElement( ja_mat[col_sel], row_sel, -a_mat[col_sel] );
			}
		}
	}

	bool ChInteriorPoint::check_exit_conditions(bool only_mu)
	{

		if (false)
		{
			if (mu < mu_tolerance)
			{
				if (only_mu)
					return true;

				VerifyKKTconditions();
				if (
					IPres.rp_nnorm < IPtolerances.rp_nnorm &&
					IPres.rd_nnorm < IPtolerances.rd_nnorm
					)
					return true;
			}
		}

		if (true)
		{
			double relative_duality_gap = y.MatrDot(y, lam); // TODO: you can recycle 'mu'
			
			VerifyKKTconditions();
			if (relative_duality_gap < mu_tolerance &&
				IPres.rp_nnorm < IPtolerances.rp_nnorm &&
				IPres.rd_nnorm < IPtolerances.rd_nnorm
				)
				return true;

		}
		

		return false;
	}

	//void ChInteriorPoint::save_lam_mean()
	//{
	//	lam_mean = 0;
	//	for (size_t row_sel = 0; row_sel < m; row_sel++)
	//	{
	//		lam_mean += lam.GetElement(row_sel, 0);
	//	}
	//	lam_mean = lam_mean / m;
	//}

	/// Update rp, rd, and mu from x,y,lam and the problem matrices G, A and E
	void ChInteriorPoint::fullupdate_residual()
	{
			// Residual initialization (16.59 pag.482)
			// rp initialization; rp = A*x - y - b
			multiplyA(x, rp);  // rp = A*x
			rp -= y;
			rp -= b;
			if (ADD_COMPLIANCE)
			{
				E.MatrMultiply(Dlam, vectm);
				rp += vectm;
			}
				

			// rd initialization; rd = G*x - A^T*lam + c
			multiplyG(x, rd); // rd = G*x
			rd += c; // rd = G*x + c
			multiplyNegAT(lam, vectn); // vectn = (-A^T)*lam
			rd += vectn; // rd = (G*x + c) + (-A^T*lam)

			mu = y.MatrDot(y, lam) / m;
	}

	void ChInteriorPoint::update_history()
	{
		if (history_file.is_open())
		{
			history_file << std::endl << solver_call << ", " << iteration_count << ", " << IPres.rp_nnorm << ", " << IPres.rd_nnorm << ", " << mu;
		}
			
	}

	// copy solution from the variables used for IP to Chrono variables
	void ChInteriorPoint::generate_solution()
	{
		// copy 'x'
		for (size_t row_sel = 0; row_sel < n; row_sel++)
			sol_chrono(row_sel, 0) = x(row_sel,0);

		// copy Lagrangian multipliers; skip tangential forces if needed
		if (SKIP_CONTACTS_UV)
		{
			for (size_t row_sel = 0; row_sel < m; row_sel++)
			{
				sol_chrono(n + row_sel*3, 0) = -lam(row_sel,0); // there will be an inversion inside FromVectorToUnknowns()
				sol_chrono(n + row_sel*3 + 1, 0) = 0;
				sol_chrono(n + row_sel*3 + 2, 0) = 0;
			}
		}
		else
		{
			for (size_t row_sel = 0; row_sel < m; row_sel++)
				sol_chrono(row_sel + n, 0) = -lam(row_sel); // there will be an inversion inside FromVectorToUnknowns()
		}
	}


	bool ChInteriorPoint::check_feasibility(double tolerance)
	{
		VerifyKKTconditions();
		if (IPres.rp_nnorm < tolerance &&
			IPres.rd_nnorm < tolerance)
			return true;

		return false;
	}

	void ChInteriorPoint::PrintHistory(bool on_off, std::string filepath)
	{
		if (!history_file.is_open())
		{
			history_file.open(filepath);
		}

		history_file << std::scientific << std::setprecision(3);
		history_file << std::endl << "SolverCall" << ", " << "Iteration" << ", " << "rp_nnorm" << ", " << "rd_nnorm" << ", " << "mu";

		print_history = true;

	}

	void ChInteriorPoint::VerifyKKTconditions(bool print)
	{
		IPres.rp_nnorm = rp.NormTwo() / m;
		IPres.rd_nnorm = rd.NormTwo() / n;

		if (print)
		{
			std::cout << std::scientific << std::setprecision(1);
			std::cout << "|rp|/n: " << IPres.rp_nnorm
				<< "; |rd|/m: " << IPres.rd_nnorm
				<< "; mu: " << mu << std::endl;

			bool neg_y = false;
			bool neg_lam = false;
			for (int cont = 0; cont < m; cont++)
			{
				if (y(cont, 0) < 0)
					neg_y = true;
				if (lam(cont, 0) < 0)
					neg_lam = true;
			}

			if (neg_y)
				std::cout << "y has negative elements" << std::endl;

			if (neg_lam)
				std::cout << "lam has negative elements" << std::endl;

		}
	}

	ChInteriorPoint::ChInteriorPoint():
		m(0),
		n(0),
		solver_call(0),
		iteration_count(0),
		iteration_count_max(50),
		EQUAL_STEP_LENGTH(true),
		ADAPTIVE_ETA(true),
		ONLY_PREDICT(false),
		warm_start(true),
		warm_start_broken(false),
		IPtolerances(1e-12),
		mu_tolerance(1e-12),
		KKT_solve_method(AUGMENTED),
        mumps_engine(),
		mu(0)
	{
		mumps_engine.SetICNTL(11, 2);
		PrintHistory(true);
	}

	ChInteriorPoint::~ChInteriorPoint()
	{
		if (history_file.is_open())
			history_file.close();
	}
}
