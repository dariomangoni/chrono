
#include "ChModalAnalysisArpack.h"



void diagonal_matrix_vector_product(double const* const x, double* const y) {
  for (int i = 0; i < 1000; ++i) {
    y[i] = static_cast<double>(i + 1) * x[i];
  }
}


namespace chrono
{

    int chrono::ChArpackSolver::compute(int requested_eigval, double sigmar, double sigmai) const
    {

        ChSparseMatrix matA(100, 100);
        ChSparseMatrix matB(100, 100);

      a_int const N = matA.rows();;
      a_int const nev = 9;

      a_int const ncv = 2 * nev + 1;
      a_int const ldv = N;

      a_int const ldz = N + 1;

      a_int const lworkl = 3 * (ncv * ncv) + 6 * ncv;

      double const tol = 0.000001; // small tol => more stable checks after EV computation.
      double const sigma = 0.0;

      a_int const rvec = 1;

      std::vector<double> resid(N);
      std::vector<double> V(ncv * N);
      std::vector<double> workd(3 * N, 0.0);
      std::vector<double> workl(lworkl, 0.0);
      std::vector<double> workev(3*ncv, 0.0);
      std::vector<double> dr((nev + 1));
      std::vector<double> di((nev + 1));
      std::vector<double> z((N + 1) * (nev + 1));

      a_int ipntr[14];
      a_int iparam[11];

      iparam[0] = 1;
      iparam[2] = 10 * N;
      iparam[3] = 1;
      iparam[4] = 0;  // number of ev found by arpack.
      iparam[6] = 1;

      a_int info = 0, ido = 0;

      while (ido != 99) {
        arpack::naupd(ido, arpack::bmat::generalized, N,
                  arpack::which::smallest_magnitude, nev, tol, resid.data(),
                  ncv, V.data(), ldv, iparam, ipntr, workd.data(),
                  workl.data(), lworkl, info);

        diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]),
                                       &(workd[ipntr[1] - 1]));
      }

      // check number of ev found by arpack.
      if (iparam[4] < nev) { /*arpack may succeed to compute more EV than expected*/
        std::cout << "ERROR: iparam[4] " << iparam[4] << ", nev " << nev
                  << ", info " << info << std::endl;
        throw std::domain_error("Error inside ARPACK routines");
      }

      std::vector<a_int> select(ncv);
      for (int i = 0; i < ncv; i++) select[i] = 1;


      arpack::neupd(rvec, arpack::howmny::ritz_vectors, select.data(),
                  dr.data(), di.data(), z.data(), ldz,
                  sigmar, sigmai, workev.data(),
                  arpack::bmat::generalized, N,
                  arpack::which::smallest_magnitude, nev, tol, resid.data(),
                  ncv, V.data(), ldv, iparam, ipntr, workd.data(),
                    workl.data(), lworkl, info);

      for (int i = 0; i < nev; ++i) {
        double val = dr[i];
        double ref = static_cast<double>(N - (nev - 1) + i);
        double eps = std::fabs(val - ref);
        std::cout << val << " - " << ref << " - " << eps << std::endl;

        /*eigen value order: smallest -> biggest*/
        if (eps > 1e-6) {
          throw std::domain_error("Correct eigenvalues not computed");
        }
      }
      std::cout << "------\n";

    }
}


