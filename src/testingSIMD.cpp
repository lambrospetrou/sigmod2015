/*

   Some examplary loops for GCC Auto-Vectorization

   by Thomas Hauth ( Thomas.Hauth@cern.ch )

   Compile with ( use at least gcc 4.7 ):
   g++ -Ofast -ftree-vectorizer-verbose=7 -march=native -std=c++11 -o autovect autovect.cpp

 */

#include <math.h>

#include <string>
#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include "include/aligned_allocator.hpp"
#include "include/LPUtils.hpp"

// Sturcture-Of-Array to hold coordinates
struct Vector3
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    // final result of the distance calcualtion
    std::vector<double> distance;

    void add( double _x, double _y, double _z )
    {
        x.push_back( _x );
        y.push_back( _y );
        z.push_back( _z );
        distance.push_back( 0.0f );
    }

};


template<typename T>
using aligned_vector = std::vector<T, aligned_allocator<T, 16>>;

int main()
{

    // Fixed Size Arrays

    typedef std::array<double, 10> DataArray;

    DataArray vect_a = { 0,1,2,3,4,5,6,7,8,9 };
    DataArray vect_b = {0.5,1,2,3,4,5,6,7,8,9 };
    DataArray vect_res_plain = { 0,0,0,0,0,0,0,0,0,0};   
    DataArray vect_res_lambda = { 0,0,0,0,0,0,0,0,0,0};

    typedef uint64_t aint __attribute__ ((__aligned__(16)));
    typedef std::vector<uint64_t> lparray;
    lparray lpArr = { 1,2,3,4,5,6,7,8,9,10 };
    auto kernelOR = [] (lparray const& a) -> uint64_t {
        aint ored;
        const unsigned int sz = a.size();
        for (unsigned int i=0; i<sz; ++i) {
            ored |= a[i];
        }
        return ored;
    };
    uint64_t orA = kernelOR(lpArr);
    std::cout << " or vect_a: " << orA << std::endl;
    uint64_t orB = lp::simd::or_all(lpArr);
    std::cout << " or vect_a: " << orB << std::endl;

    auto kernelExists = [] (lparray const& a, const uint64_t val) -> bool {
        aint mask = 0;
        const unsigned int sz = a.size();
        for (unsigned int i=0; i<sz; ++i) {
            mask |= a[i] == val;
        }
        return mask;
    };
    uint64_t exA = kernelExists(lpArr, 6);
    std::cout << " exists: " << exA << std::endl;
    exA = lp::simd::exists<uint64_t>(lpArr, 20);
    std::cout << " exists (no): " << exA << std::endl;
    exA = lp::simd::exists<uint64_t>(lpArr.data()+2, lpArr.size()-1, 3);
    std::cout << " exists: " << exA << std::endl;

    constexpr double cFixedMultiply = 23.5f;

    // simple loop vectorized
    // -- auto-vectorized --
    for( unsigned int i = 0; i < vect_a.size(); ++ i)
    {
        vect_res_plain[i] = vect_a[i]  + vect_b[i];
    }

    // Defining a compute kernel to encapsulate a specific computation
    auto kernel_multiply = 
        [ &cFixedMultiply ] // capture the constant by reference of the scope of the lambda expression
        ( DataArray const& a, DataArray const& b, DataArray & res ) // take 3 parameters by reference
        ->void // lambda function returns void
        {
            // simple loop vectorized
            // -- auto-vectorized --
            for( unsigned int i = 0; i < a.size(); ++ i)
            {
                res[i] = a[i] * b[i] * cFixedMultiply;
            }
        };

    // call the lambda function
    // this call is autovectorized
    kernel_multiply ( vect_a, vect_b, vect_res_lambda );


    // This kernel will be called multiple times and performs the quadrature
    auto kernel_square = 
        [] // capture nothing
        ( double const& a) // take 1 parameter by reference
        ->double // lambda function returns the square
        {
            return ( a * a );
        };

    // create struct and fill with dummy values
    Vector3 v3;
    for ( unsigned int i = 0; i < 50 ; ++ i)
    {
        v3.add( i * 1.1, i * 2.2,  i* 3.3 );   
    }

    // store the size in a local variable. This is needed to GCG
    // can estimate the loop iterations
    const unsigned int size = v3.x.size();

    // -- auto-vectorized --
    for ( unsigned int i = 0; i < size; ++ i)
    {
        v3.distance[i] = sqrt( kernel_square( v3.x[i] ) + kernel_square( v3.y[i] ) + kernel_square( v3.z[i] )) ;
    }

    // output the result, so GCC won't optimize the calculations away
    std::cout << std::endl << "Computation on std::array" << std::endl;
    for( unsigned int i = 0; i < vect_a.size(); ++ i)
    {
        std::cout << vect_res_plain[i] << std::endl;
        std::cout << vect_res_lambda[i] << std::endl;
    }

    std::cout << std::endl << "Computation on Structure-of-Array with variable sized std::vectors" << std::endl;
    for( unsigned int i = 0; i < v3.x.size(); ++ i)
    {
        std::cout << "sqrt( " << v3.x[i] << "^2 + " << v3.y[i] << "^2 + " << v3.z[i] << "^2 ) = " 
            << v3.distance[i] << std::endl;
    }

    return 0;
}
