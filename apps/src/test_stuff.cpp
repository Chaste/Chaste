#include <iostream>
#include <typeinfo>
#include <vector>
#include <cmath>
#include "UblasIncludes.hpp"
#include "ImmersedBoundaryArray.hpp"

void test_strings()
{
    std::string greeting = "Hello World!";

    std::cout << greeting << std::endl;

    if (greeting == "Hello World!")
    {
        std::cout << "Greeting == Hello World!" << std::endl;
    }
    else
    {
        std::cout << "It was something else!" << std::endl;
    }
}

void test_vectors()
{
    std::vector<double> array;
    for (int i = 0; i < 10; i++)
    {
        array.push_back(double(i));
        std::cout << "element " << i << " " << array[i] << std::endl;
    }

    std::cout << array.size() << std::endl;

    std::vector<double> vec_1(2);
    vec_1[0] = 1.0; vec_1[1] = 2.0;

    std::cout << "vec_1 = (" << vec_1[0] << ", " << vec_1[1] << ")" << std::endl;

    std::vector<double> vec_2(2);
    vec_2 = vec_1;

    std::cout << "vec_2 = (" << vec_2[0] << ", " << vec_2[1] << ")" << std::endl;

    vec_1[0] = 0; vec_1[1] = 0;

    std::cout << "vec_1 = (" << vec_1[0] << ", " << vec_1[1] << ")" << std::endl;
    std::cout << "vec_2 = (" << vec_2[0] << ", " << vec_2[1] << ")" << std::endl;

}

void test_c_vectors()
{
    c_vector<double, 2> vec;
    vec[0] = 1; vec[1] = 1;

    std::cout << "vec = (" << vec[0] << ", " << vec[1] << ")" << std::endl;

    double magnitude;

    magnitude = norm_2(vec);

    std::cout << magnitude << std::endl;

    multi_array<double, 2> an_array;

    an_array.resize(extents[2][2]);
    an_array[0][0] = 1.0;
    an_array[1][0] = 0.0;
    an_array[0][1] = 0.0;
    an_array[1][1] = 1.0;

    std::cout << "array = [[" << an_array[0][0] << ", " << an_array[0][1] << "]" <<
        "\n\t\t [" << an_array[1][0] << ", " << an_array[1][1] << "]]" << std::endl;
}

void test_passing_vectors(double* array, int size)
{

    for (int i = 0; i < size; i++)
    {
        std::cout << "array[" << i << "] = " << array[i] << std::endl;
    }

}

void test_complex_numbers()
{
    std::complex<double> i(0.0, 1.0);
    std::complex<double> result = exp(M_PI * i);

    std::cout << result << std::endl;
}

void test_arrays()
{
    multi_array<double, 2> array;
    array.resize(extents[2][2]);
    array[0][0] = 0.0; array[0][1] = 1.0;
    array[1][0] = 1.0; array[1][1] = 0.0;
    //  array = array + array; <- this is not allowed
}

void test_c_vectors_2()
{
//    c_vector<double, 2> vec(1.0, 2.0);
}

int main()
{
//    double vec[5] = {1, 4, 3, 4, 5};
//    test_passing_vectors(vec, 5);
//    test_complex_numbers();
//    test_arrays();

    return 0;
}
