#include <iostream>
//#include <random>
#include <time.h>
#include <cstdlib>



//define function to integrate
double f(double x)
{
	return 1.5 * (1 - x * x);
}


// dividing axes into subgrid
double *linspace(double a, double b, int n)
{
	double *space = new double[n + 1];

	for (int i = 0; i < n + 1; i++)
	{
		space[i] = i * (b - a) / n;
	}

	return space;
}


/*double random_number()
{	//Mersenne twister (rng) for 2 b)
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> unif(0, 1);

	return unif(gen);
}

*/

int main()
{
	clock_t tStart = clock();

	double I = 0.0;  //value of integral
	const int d = 1; //number of dimensions

	// 2a)
	const int n = 6; //number of grid intervals

	double *x = linspace(0, 1, n);
	double dx = x[1] - x[0]; //grid spacing



	//initialize vector
	double x_vect[d][n + 1];

	for (int j = 0; j < d; j++)
	{
		for (int i = 0; i < n+1; i++)
		{
			x_vect[j][i] = x[i];
		}
	}

	//calculate integral using midpoint method
	for (int j = 0; j < d; j++)
	{
		for (int i = 0; i < n; i++)
		{

			I += f((x_vect[j][i + 1] + x_vect[j][i]) / 2) * dx / d;
		}
	}

	delete[] x;


	/*
	// 2b)
	const int N = 20000;

	//set up random vector
	float rand_vector[d][N];

	for (int j = 0; j < d; j++)
	{
		for (int i = 0; i < N; i++)
		{
			rand_vector[j][i] = random_number();

			I += f(rand_vector[j][i]) / N / d;
		}
	}
	*/


	std::cout << I << std::endl;

	std::cout << "Time taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s" << std::endl;
	std::cin.get();
	return 0;
}
