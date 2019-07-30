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

#ifndef TESTFOURTHORDERTENSOR_HPP_
#define TESTFOURTHORDERTENSOR_HPP_

#include <cxxtest/TestSuite.h>
#include "FourthOrderTensor.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestFourthOrderTensor : public CxxTest::TestSuite
{
public:

    void TestFourthOrderTensorAllSameDimensions()
    {
        FourthOrderTensor<2,2,2,2> x;

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(x(M,N,P,Q), 0.0, 1e-9);
                        x(M,N,P,Q) = M+N+P+Q;
                    }
                }
            }
        }

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(x(M,N,P,Q), M+N+P+Q, 1e-9);
                    }
                }
            }
        }

        FourthOrderTensor<3,3,3,3> y;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(y(M,N,P,Q), 0.0, 1e-9);
                        y(M,N,P,Q) = M+N+P+Q;
                    }
                }
            }
        }

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(y(M,N,P,Q), M+N+P+Q, 1e-9);
                    }
                }
            }
        }

        y.Zero();

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA(y(M,N,P,Q), 0.0, 1e-9);
                    }
                }
            }
        }
    }

    void TestSetAsContractionOnFirstDimensionSameDimensions()
    {
        FourthOrderTensor<3,3,3,3> X;
        c_matrix<double,3,3> A;
        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = 3*M+N+P+Q;
                    }
                }
            }
        }

        FourthOrderTensor<3,3,3,3> Y;
        Y.SetAsContractionOnFirstDimension<3>(A,X);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(Q+3*M+N+P)+(Q+N+P)*3*M, 1e-9);
                    }
                }
            }
        }
    }

    void TestSetAsContractionOnSecondDimensionSameDimensions()
    {
        FourthOrderTensor<3,3,3,3> X;
        c_matrix<double,3,3> A;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+3*N+P+Q;
                    }
                }
            }
        }

        FourthOrderTensor<3,3,3,3> Y;
        Y.SetAsContractionOnSecondDimension<3>(A,X);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(Q+M+3*N+P)+(Q+M+P)*3*N, 1e-9);
                    }
                }
            }
        }
    }

    void TestSetAsContractionOnThirdDimensionSameDimensions()
    {
        FourthOrderTensor<3,3,3,3> X;
        c_matrix<double,3,3> A;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+N+3*P+Q;
                    }
                }
            }
        }

        FourthOrderTensor<3,3,3,3> Y;
        Y.SetAsContractionOnThirdDimension<3>(A,X);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(Q+M+N+3*P)+(Q+M+N)*3*P, 1e-9);
                    }
                }
            }
        }
    }

    void TestSetAsContractionOnFourthDimensionSameDimensions()
    {
        FourthOrderTensor<3,3,3,3> X;
        c_matrix<double,3,3> A;

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                A(M,N) = M+N;
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        X(M,N,P,Q) = M+N+P+3*Q;
                    }
                }
            }
        }

        FourthOrderTensor<3,3,3,3> Y;
        Y.SetAsContractionOnFourthDimension<3>(A,X);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                for (unsigned P=0; P<3; P++)
                {
                    for (unsigned Q=0; Q<3; Q++)
                    {
                        TS_ASSERT_DELTA( Y(M,N,P,Q), 15+3*(3*Q+M+N+P)+(M+N+P)*3*Q, 1e-9);
                    }
                }
            }
        }
    }

    void TestFourthOrderTensorDifferentDimensions1()
    {
        FourthOrderTensor<2,3,1,1> X;

        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned k=0; k<1; k++)
                {
                    for (unsigned n=0; n<1; n++)
                    {
                        X(i,j,k,n) = (i==j);
                    }
                }
            }
        }

        c_matrix<double,1,2> A;
        A(0,0) = 2;
        A(0,1) = 3;

        FourthOrderTensor<1,3,1,1> Y;
        Y.SetAsContractionOnFirstDimension<2>(A,X);

        TS_ASSERT_DELTA( Y(0,0,0,0), A(0,0), 1e-8);
        TS_ASSERT_DELTA( Y(0,1,0,0), A(0,1), 1e-8);
        TS_ASSERT_DELTA( Y(0,2,0,0), 0.0,    1e-8);
    }

    void TestFourthOrderTensorDifferentDimensions2()
    {
        FourthOrderTensor<2,3,1,1> X;

        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                for (unsigned k=0; k<1; k++)
                {
                    for (unsigned n=0; n<1; n++)
                    {
                        X(i,j,k,n) = (i==j);
                    }
                }
            }
        }

        c_matrix<double,1,3> A;
        A(0,0) = 2;
        A(0,1) = 3;
        A(0,2) = 4;

        FourthOrderTensor<2,1,1,1> Y;
        Y.SetAsContractionOnSecondDimension<3>(A,X);

        TS_ASSERT_DELTA( Y(0,0,0,0), A(0,0), 1e-8);
        TS_ASSERT_DELTA( Y(1,0,0,0), A(0,1), 1e-8);
    }

    void TestFourthOrderTensorDifferentDimensions3()
    {
        FourthOrderTensor<1,1,2,3> X;

        for (unsigned i=0; i<1; i++)
        {
            for (unsigned j=0; j<1; j++)
            {
                for (unsigned k=0; k<2; k++)
                {
                    for (unsigned n=0; n<3; n++)
                    {
                        X(i,j,k,n) = (k==n);
                    }
                }
            }
        }

        c_matrix<double,1,2> A;
        A(0,0) = 2;
        A(0,1) = 3;

        FourthOrderTensor<1,1,1,3> Y;
        Y.SetAsContractionOnThirdDimension<2>(A,X);

        TS_ASSERT_DELTA( Y(0,0,0,0), A(0,0), 1e-8);
        TS_ASSERT_DELTA( Y(0,0,0,1), A(0,1), 1e-8);
        TS_ASSERT_DELTA( Y(0,0,0,2), 0.0,    1e-8);
    }

    void TestFourthOrderTensorDifferentDimensions4()
    {
        FourthOrderTensor<1,1,2,3> X;

        for (unsigned i=0; i<1; i++)
        {
            for (unsigned j=0; j<1; j++)
            {
                for (unsigned k=0; k<2; k++)
                {
                    for (unsigned n=0; n<3; n++)
                    {
                        X(i,j,k,n) = (k==n);
                    }
                }
            }
        }

        c_matrix<double,1,3> A;
        A(0,0) = 2;
        A(0,1) = 3;
        A(0,2) = 4;

        FourthOrderTensor<1,1,2,1> Y;
        Y.SetAsContractionOnFourthDimension<3>(A,X);

        TS_ASSERT_DELTA( Y(0,0,0,0), A(0,0), 1e-8);
        TS_ASSERT_DELTA( Y(0,0,1,0), A(0,1), 1e-8);
    }
};

#endif /*TESTFOURTHORDERTENSOR_HPP_*/
