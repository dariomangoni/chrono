
#include <chrono_mumps/ChMumpsEngine.h>
#include <chrono/core/ChMatrixDynamic.h>

using namespace chrono;


int main()
{
	int n = 3;
	ChCOOMatrix mat(n,n, true);
	ChMatrixDynamic<double> rhs(n,1);
	ChMumpsEngine mumps_engine;
 
	rhs(0) = 2;
	rhs(1) = 1;
	rhs(2) = 3;

	mat.SetElement(0, 0, 1.3);
	mat.SetElement(1, 1, 2.7);
	mat.SetElement(2, 2, 3.9);
 
    mat.Compress();
	mumps_engine.SetProblem(mat, rhs);
    auto return_value = mumps_engine.MumpsCall(ChMumpsEngine::mumps_JOB::COMPLETE);
    mumps_engine.PrintINFOG();
    printf("Mumps says: %d\n", return_value);


	for (int i = 0; i < n; i++)
	{
		printf("%f ", rhs(i));
 	}

 
 	getchar();
	return 0;
 }