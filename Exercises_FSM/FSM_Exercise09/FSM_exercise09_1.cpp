#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>


std::vector<double> randu(uint64_t initialseed, int iterations) {
    if (initialseed % 2 == 0) {
        throw "initial value needs to be odd";
    } else {
        static const uint64_t a = 65539;
        static const uint64_t m = 2147483648;  // 2^31

        std::vector<double> values;
        //values.reserve(iterations);
        values.push_back(initialseed);
        for (int i = 1; i < iterations; i++) {
            initialseed = values[i - 1];
            values.push_back((a * initialseed) % m);
            //if(i < 10) {std::cout << values[i - 1] << std::endl;}
        }

        for (int i = 0; i < iterations; i++) {
            values[i] = values[i] / m;
        }

        return values;
    }
}

std::vector<double> randuRestricted(uint64_t initialseed, int iterations, std::vector<double> x_interval, std::vector<double> y_interval) {
    if (initialseed % 2 == 0) {
        throw "initial value needs to be odd";
    } else {
        static const uint64_t a = 65539;
        static const uint64_t m = 2147483648; 

        std::vector<double> values_x;
        std::vector<double> values_y;
        uint64_t initial_2 = initialseed;
        uint64_t initial_1 = 0;
        double x_interval_small = x_interval[0] * m;
        double x_interval_large = x_interval[1] * m;
        double y_interval_small = y_interval[0] * m;
        double y_interval_large = y_interval[1] * m;

        int i = 0;
        while(i < iterations) {
            initial_1 = (a * initial_2) % m;
            initial_2 = (a * initial_1) % m;

            if (((x_interval_small < initial_1) && (initial_1 < x_interval_large)) && ((y_interval_small < initial_2) && (initial_2 < y_interval_large))) {
                values_x.push_back(initial_1 / double(m));
                values_y.push_back(initial_2 / double(m));
                i++;                
                //std::cout << "number in intervall: " << i << std::endl;
            }
        }
        std::vector<double> values;
        for(int i = 0; i < iterations; i++) {
            values.push_back(values_x[i]);
        }
        for(int i = 0; i < iterations; i++) {
            values.push_back(values_y[i]);
        }
        return values;
    }
}


int main() {
    std::vector<double> random_values = randu(1, 1000);
    std::ofstream output_file("./randu.txt");
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    std::copy(random_values.begin(), random_values.end(), output_iterator);

    std::vector<double> randomBla = randuRestricted(1, 1000, {0.2, 0.201}, {0.3, 0.301});
    std::ofstream output_file2("./randuRestricted.txt");
    std::ostream_iterator<double> output_iterator2(output_file2, "\n");
    std::copy(randomBla.begin(), randomBla.end(), output_iterator2);
    return 0;
};
