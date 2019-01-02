/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTUBLASCUSTOMFUNCTIONS_HPP_
#define TESTUBLASCUSTOMFUNCTIONS_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasMatrixInclude.hpp"
#include <boost/numeric/ublas/io.hpp>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"

class TestUblasCustomFunctions : public CxxTest::TestSuite
{
public:

    void TestDeterminant()
    {
        using namespace boost::numeric::ublas;
        c_matrix<double, 1, 1> C;
        double OneOneDeterminant;
        C(0,0) = 5.6;
        OneOneDeterminant = Determinant(C);
        TS_ASSERT_DELTA( OneOneDeterminant, 5.6, 0.0000000001);

        c_matrix<double, 3, 3> A;
        A(0,0) = 2.4;
        A(0,1) = 5;
        A(0,2) = 5;
        A(1,0) = 5;
        A(1,1) = 6;
        A(1,2) = 7;
        A(2,0) = 6;
        A(2,1) = 8;
        A(2,2) = 9;
        double ThreeThreeDeterminant = Determinant(A);
        TS_ASSERT_DELTA( ThreeThreeDeterminant, 0.2, 0.0000000001);

        c_matrix<double, 2, 2> B;
        B(0,0) = 2.4;
        B(0,1) = 5;
        B(1,0) = 5;
        B(1,1) = 6;
        double TwoTwoDeterminant = Determinant(B);
        TS_ASSERT_DELTA( TwoTwoDeterminant, -10.6, 0.0000000001);

        c_matrix<double, 2, 1> D;
        D(0, 0) = 3.0;
        D(1, 0) = 4.0;
        double TwoOneDeterminant = Determinant(D);
        TS_ASSERT_DELTA(TwoOneDeterminant, 5.0, 0.000000001);

        c_matrix<double, 3, 1> E;
        E(0,0) = 4.0;
        E(1,0) = 4.0;
        E(2,0) = 2.0;
        double ThreeOneDeterminant = Determinant(E);
        TS_ASSERT_DELTA(ThreeOneDeterminant, 6.0, 0.000000001);

        c_matrix<double, 3, 2> F;
        F(0,0) = 4;
        F(1,0) = 2;
        F(2,0) = -1;
        F(0,1) = 3;
        F(1,1) = -2;
        F(2,1) = -1;

        double ThreeTwoDeterminant = Determinant(F);
        TS_ASSERT_DELTA(ThreeTwoDeterminant, 14.595, 0.001);
    }

    void TestSubDeterminant()
    {
        using namespace boost::numeric::ublas;

        c_matrix<double, 1, 1> C;
        double OneOneDeterminant;
        C(0,0) = 5.6;
        OneOneDeterminant = Determinant(C);
        TS_ASSERT_DELTA( OneOneDeterminant,
                         C(0,0)*SubDeterminant(C,0,0),
                         1e-10);

        c_matrix<double, 2, 2> B;
        B(0,0) = 2.4;
        B(0,1) = 5;
        B(1,0) = 5;
        B(1,1) = 6;
        double TwoTwoDeterminant = Determinant(B);


        TS_ASSERT_DELTA( TwoTwoDeterminant,
                         B(0,0)*SubDeterminant(B,0,0)-B(0,1)*SubDeterminant(B,0,1), 1e-10);
        TS_ASSERT_DELTA( TwoTwoDeterminant,
                         B(0,0)*SubDeterminant(B,0,0)-B(1,0)*SubDeterminant(B,1,0), 1e-10);
        c_matrix<double, 3, 3> A;
        A(0,0) = 2.4;
        A(0,1) = 5;
        A(0,2) = 5;
        A(1,0) = 5;
        A(1,1) = 6;
        A(1,2) = 7;
        A(2,0) = 6;
        A(2,1) = 8;
        A(2,2) = 9;
        double ThreeThreeDeterminant = Determinant(A);
        TS_ASSERT_DELTA( ThreeThreeDeterminant, 0.2, 1e-10);
        TS_ASSERT_DELTA( ThreeThreeDeterminant,
                         A(0,0)*SubDeterminant(A,0,0)
                         - A(1,0)*SubDeterminant(A,1,0)
                         + A(2,0)*SubDeterminant(A,2,0),
                         1e-10);
        TS_ASSERT_DELTA( ThreeThreeDeterminant,
                         - A(0,1)*SubDeterminant(A,0,1)
                         + A(1,1)*SubDeterminant(A,1,1)
                         - A(2,1)*SubDeterminant(A,2,1),
                         1e-10);
        TS_ASSERT_DELTA( ThreeThreeDeterminant,
                         A(0,0)*SubDeterminant(A,0,0)
                         - A(0,1)*SubDeterminant(A,0,1)
                         + A(0,2)*SubDeterminant(A,0,2),
                         1e-10);
        TS_ASSERT_DELTA( ThreeThreeDeterminant,
                         - A(1,0)*SubDeterminant(A,1,0)
                         + A(1,1)*SubDeterminant(A,1,1)
                         - A(1,2)*SubDeterminant(A,1,2),
                         1e-10);
    }

    void TestInverse()
    {
        using namespace boost::numeric::ublas;
        c_matrix<double, 1, 1> C;
        C(0,0) = 8;
        c_matrix<double, 1, 1> invC;
        invC = Inverse(C);
        TS_ASSERT_DELTA(invC(0,0), 0.125, 0.0000000001);
        c_matrix<double, 3, 3> A;
        A(0,0) = 2.4;
        A(0,1) = 5;
        A(0,2) = 5;
        A(1,0) = 5;
        A(1,1) = 6;
        A(1,2) = 7;
        A(2,0) = 6;
        A(2,1) = 8;
        A(2,2) = 9;
        c_matrix<double, 3, 3> invA;
        c_matrix<double, 3, 3> invAMatlab;
        invAMatlab(0,0) = -10.00;
        invAMatlab(0,1) = -25.00;
        invAMatlab(0,2) = 25.00;
        invAMatlab(1,0) = -15.00;
        invAMatlab(1,1) = -42.00;
        invAMatlab(1,2) = 41.00;
        invAMatlab(2,0) = 20.00;
        invAMatlab(2,1) = 54.00;
        invAMatlab(2,2) = -53.00;
        invA = Inverse(A);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                TS_ASSERT_DELTA( invA(i,j), invAMatlab(i,j), 0.0000000001);
            }
        }
        c_matrix<double, 2, 2> B;
        B(0,0) = 2.4;
        B(0,1) = 5;
        B(1,0) = 5;
        B(1,1) = 6;
        c_matrix<double, 2, 2> invB;
        invB = Inverse(B);
        c_matrix<double, 2, 2> invBMatlab;
        invBMatlab(0,0) = -0.5660;
        invBMatlab(0,1) = 0.4717;
        invBMatlab(1,0) = 0.4717;
        invBMatlab(1,1) = -0.2264;
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                TS_ASSERT_DELTA( invB(i,j), invBMatlab(i,j), 0.0001);
            }
        }

        c_matrix<double, 2, 1> D;
        D(0,0) = 2;
        D(1,0) = 1;
        c_matrix<double, 1, 2> invD;
        invD = Inverse(D);
        TS_ASSERT_DELTA(invD(0, 0), 0.4, 0.00000001);
        TS_ASSERT_DELTA(invD(0, 1), 0.2, 0.00000001);

        c_matrix<double, 3, 1> E;
        E(0,0) = 3;
        E(1,0) = 4;
        E(2,0) = 5;
        c_matrix<double, 1, 3> invE;
        invE = Inverse(E);
        TS_ASSERT_DELTA(invE(0, 0), 0.06, 0.00000001);
        TS_ASSERT_DELTA(invE(0, 1), 0.08, 0.00000001);
        TS_ASSERT_DELTA(invE(0, 2), 0.1, 0.00000001);

        c_matrix<double, 3, 2> F;
        F(0,0) = 4;
        F(1,0) = 2;
        F(2,0) = -1;
        F(0,1) = 3;
        F(1,1) = -2;
        F(2,1) = -1;
        c_matrix<double, 2, 3> invF;
        invF = Inverse(F);
        c_matrix<double, 2, 2> identity;
        identity = prod(invF, F);
        TS_ASSERT_DELTA(identity(0,0), 1.000, 1e-3);
        TS_ASSERT_DELTA(identity(1,1), 1.000, 1e-3);
        TS_ASSERT_DELTA(identity(0,1), 0.000, 1e-3);
        TS_ASSERT_DELTA(identity(1,0), 0.000, 1e-3);
    }

    void TestTraceAndSecondInvariant()
    {
        c_matrix<double, 1,1> a;
        a(0,0) = 13.03;
        TS_ASSERT_DELTA(Trace(a),13.03,1e-10);

        c_matrix<double, 2,2> b;
        b(0,0) = 13.03;
        b(1,0) = 3.03;
        b(0,1) = 3.03;
        b(1,1) = 165;
        TS_ASSERT_DELTA(Trace(b),13.03+165,1e-10);

        TS_ASSERT_DELTA(SecondInvariant(b), Determinant(b), 1e-12);

        // Symmetric 3 by 3 matrix.
        c_matrix<double, 3,3> c;
        c(0,0) = 13.03;
        c(1,0) = 1.3;
        c(2,0) = 2.3;
        c(0,1) = 1.3;
        c(1,1) = 45;
        c(2,1) = 9.9;
        c(0,2) = 2.3;
        c(1,2) = 9.9;
        c(2,2) = 34;

        c_matrix<double,3,3> c_squared = prod(trans(c),c);

        TS_ASSERT_DELTA(Trace(c),13.03+45+34,1e-10);
        TS_ASSERT_DELTA(SecondInvariant(c),0.5*(Trace(c)*Trace(c)-Trace(c_squared)),1e-10);

        c_matrix<double, 4,4> d = zero_matrix<double>(4,4);
        d(0,0) = 1.1;
        d(1,1) = 2.2;
        d(2,2) = 3.3;
        d(3,3) = 4.4;
        d(3,2) = 10;
        TS_ASSERT_DELTA(Trace(d),11.0,1e-10);
    }

    // Get a row from a matrix
    void TestUblasMatrixRow()
    {
        c_matrix<double,2,3> a;
        a(0,0) = 1;
        a(0,1) = 2;
        a(0,2) = 3;
        a(1,0) = 4;
        a(1,1) = 5;
        a(1,2) = 6;

        matrix_row< c_matrix<double,2,3> > row0(a, 0);
        matrix_row< c_matrix<double,2,3> > row1(a, 1);

        TS_ASSERT_EQUALS(a(0,0), row0(0));
        TS_ASSERT_EQUALS(a(0,1), row0(1));
        TS_ASSERT_EQUALS(a(0,2), row0(2));
        TS_ASSERT_EQUALS(a(1,0), row1(0));
        TS_ASSERT_EQUALS(a(1,1), row1(1));
        TS_ASSERT_EQUALS(a(1,2), row1(2));
    }

    void TestCreate_c_vector()
    {
        c_vector<double, 1> v1 = Create_c_vector(1);
        TS_ASSERT_EQUALS(v1[0], 1);

        c_vector<double, 2> v2 = Create_c_vector(1,2);
        TS_ASSERT_EQUALS(v2[0], 1);
        TS_ASSERT_EQUALS(v2[1], 2);

        c_vector<double, 3> v3 = Create_c_vector(1,2,3);
        TS_ASSERT_EQUALS(v3[0], 1);
        TS_ASSERT_EQUALS(v3[1], 2);
        TS_ASSERT_EQUALS(v3[2], 3);
    }

    void TestEigenVectorValueCalculation()
    {
        c_matrix<double, 3, 3> A;
        A(0,0) = 22.4;
        A(0,1) = 5;
        A(0,2) = 6;
        A(1,0) = 5;
        A(1,1) = 6;
        A(1,2) = 7;
        A(2,0) = 6;
        A(2,1) = 7;
        A(2,2) = 9;

        double smallest_eigenvalue = 0.3403300709235;

        c_vector<double, 3> eigenvector;

        eigenvector = CalculateEigenvectorForSmallestNonzeroEigenvalue(A);

        /*
         * It's important to check that eigenvector is non-zero. If eigenvector has been
         * zero-initialised and CalculateEigenvectorForSmallestNonzeroEigenvalue
         * did *nothing* then Au=lambda.u for any lambda.
         */
        double norm_of_eigenvector = norm_2(eigenvector);
        double delta = 1e-12;
        TS_ASSERT_DELTA(norm_of_eigenvector, 1.0, delta);

        // Check Au=lambda.u
        c_vector<double, 3> a_times_eigenvector = prod(A, eigenvector);

        TS_ASSERT_DELTA( a_times_eigenvector[0], smallest_eigenvalue*eigenvector[0], delta);
        TS_ASSERT_DELTA( a_times_eigenvector[1], smallest_eigenvalue*eigenvector[1], delta);
        TS_ASSERT_DELTA( a_times_eigenvector[2], smallest_eigenvalue*eigenvector[2], delta);

        //Diagonal matrix with eigenvalues 10.0, -0.5 and 0.0
        //Associated eigenvectors are x, y and z
        c_matrix<double, 3, 3> singular;
        singular(0,0) = 10.0;
        singular(0,1) = 0.0;
        singular(0,2) = 0.0;
        singular(1,0) = 0.0;
        singular(1,1) = 0.0;
        singular(1,2) = 0.0;
        singular(2,0) = 0.0;
        singular(2,1) = 0.0;
        singular(2,2) = 0.5;

        eigenvector = CalculateEigenvectorForSmallestNonzeroEigenvalue(singular);
        TS_ASSERT_DELTA( 0.0, eigenvector[0], delta);
        TS_ASSERT_DELTA( 0.0, eigenvector[1], delta);
        TS_ASSERT_DELTA( 1.0, eigenvector[2], delta);

        c_matrix<double, 3, 3> non_symmetric;
        non_symmetric(0,0) = 0.0;
        non_symmetric(0,1) = 0.0;
        non_symmetric(0,2) = 0.0;
        non_symmetric(1,0) = 0.0;
        non_symmetric(1,1) = 10.0;
        non_symmetric(1,2) = 0.0;
        non_symmetric(2,0) = 0.0;
        non_symmetric(2,1) = -100.0;
        non_symmetric(2,2) = 0.5;

        TS_ASSERT_THROWS_THIS(CalculateEigenvectorForSmallestNonzeroEigenvalue(non_symmetric),
                              "Matrix should be symmetric");

    }
};

#endif /*TESTUBLASCUSTOMFUNCTIONS_HPP_*/
