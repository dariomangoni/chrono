
#include <chrono_mumps/ChMumpsEngine.h>
#include <chrono/core/ChMatrixDynamic.h>

using namespace chrono;


int main()
{
	int n = 3;
	ChCOOMatrix mat(n,n);
	ChMatrixDynamic<double> rhs(n,1);
	ChMumpsEngine mumps_engine;
 
	rhs(0) = 2;
	rhs(1) = 1;
	rhs(2) = 3;

	mat.SetElement(0, 0, 1.3);
	mat.SetElement(1, 1, 2.7);
	mat.SetElement(2, 2, 3.9);
 
	mumps_engine.Initialize();
	mumps_engine.SetProblem(mat, rhs);
	printf("Mumps says: %d\n",mumps_engine.MumpsCall());


	for (int i = 0; i < n; i++)
	{
		printf("%f ", rhs(i));
 	}

 
 	getchar();
	return 0;
 }