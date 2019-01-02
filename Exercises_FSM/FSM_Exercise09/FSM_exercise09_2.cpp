#include <iostream>
#include <vector>
#include <stdlib.h>
#include <random>
#include <cstdlib>

class Functor {
    public:
        virtual double operator()(double x) = 0;
        virtual ~Functor(){}
};

class TestFunction : public Functor {
    public:
        TestFunction(){}

        double operator()(double x) override {
            /*dimenstion = x.size();
            f = 0;
            for (int i = 0; i < dimenstion; i++) {
                f *= (3 / 2) * (1 - (x[i] * x[i]));
            }
            return f;*/
            return (3 / 2) * (1 - (x * x));
        }
};

//2a - Midpoint
double *subintervals(double a, double b, int n) {
    double *ang = new double[n + 1];
    
    for(int i = 0; i < n + 1; i++) {
        ang[i] = i * (b - a) / n;
    }

    return ang;
}

double Midpoint(Functor& f, double a = 0, double b = 1, int d = 1, int n = 6) {
    double *h = subintervals(a, b, n);
    double dh = h[1] - h[0];

    double vector[d][n + 1];
    double result = 0;

    for(int i = 0; i < d; i++) {
        for(int j = 0; j < n + 1; j++) {
            vector[i][j] = h[i];
        }
    }
    
    for(int i = 0; i < d; i++) {
        for(int j = 0; j < n; j++) {
            result += f(0.5 * (vector[i][j + 1] + vector[i][j])) * (dh / d);
        }
    }

    return result;
}


// 2b - Monte Carlo
double random_number() {	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> unif(0, 1);

	return unif(gen);
}

double MonteCarlo(Functor& f, int d = 1, int N = 1000) {
    double result = 0;
    double random_vector[d][N];
    for(int i = 0; i < d; i++) {
        for(int j = 0; j < N; j++) {
            random_vector[i][j] = random_number();
            result += f(random_vector[i][j]) / N / d;
        }
    }

    return result;
}


int main() {
    //int d = 10;
    static const int n = 6;
    static const int N = 20000;
    TestFunction func1;
    /*
    double testMid = Midpoint(func1, 0, 1, d, n);
    double testMonte = MonteCarlo(func1, d, N);
    std::cout << "Midpoint: " << testMid << ", Monte Carlo: " << testMonte << std::endl;
    */
    for(int i = 1; i <= 10; i++) {
        clock_t time_start = clock();
        double valueMid = Midpoint(func1, 0, 1, i, n);
        std::cout << "Midpoint: \t d = " << i << ", value: " << valueMid << ",\t time taken: " << (double)(clock() - time_start) / CLOCKS_PER_SEC << " s" << std::endl;
        time_start = clock();
        double valueMonte = MonteCarlo(func1, i, N);
        std::cout << "Monte Carlo: \t d = " << i << ", value: " << valueMonte << ",\t time taken: " << (double)(clock() - time_start) / CLOCKS_PER_SEC << " s" << std::endl;
    }

    return 0;
}