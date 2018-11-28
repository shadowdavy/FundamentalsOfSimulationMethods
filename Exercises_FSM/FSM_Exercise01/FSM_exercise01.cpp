#include <iostream>
#include <stdio.h>
#include <math.h>

// exercise 2.2 functor
class Foo {
    public:
        // use default construcor
        double operator()(double x) {
            //double result = (x + std::exp(-x) - 1.) / (x * x);
            if(x > 0 && x < 6.10352e-6) {
                return(0.5);
            } else {
                return((x + std::exp(-x) - 1.) / (x * x));
            }
        }
};

int main() {
    // exercise 1.1
    int i = 7;
    float y = 2 * (i / 2);
    float z = 2 * (i / 2.);
    printf("%e %e \n", y, z);

    // exercise 1.2
    double a = 1.0e17;
    double b = -1.0e17;
    double c = 1.0;
    double v = (a + b) + c;
    double w = a + (b + c);
    printf("%e %e \n", v, w);

    // exercise 1.3
    float x = 1e20;
    float s;
    s = x * x;
    printf("%e %e \n", x, s / x);
  
    // exercise 2
    // exercise 2.2
    double userinput;
    Foo functionValue;
    std::cin >> userinput;
    double functionOutput = functionValue(userinput);
    std::cout << functionOutput << std::endl;

    // exercise 2.3
    double j = 0.1;
    //j = 6.10352e-6;
    while(j > 1e-9) {
        functionOutput = functionValue(j);
        std::cout << "input value: " << j << ", output: " << functionOutput << std::endl;
        j = j / 2;
    }

    j = -0.1;
    while(j < -1e-9) {
        functionOutput = functionValue(j);
        std::cout << "input value: " << j << ", output: " << functionOutput << std::endl;
        j = j / 2;
    }
}

