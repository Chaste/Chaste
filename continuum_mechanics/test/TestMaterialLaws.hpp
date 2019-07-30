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

#ifndef TESTMATERIALLAWS_HPP_
#define TESTMATERIALLAWS_HPP_

#include <cxxtest/TestSuite.h>

#include <cassert>

#include "MooneyRivlinMaterialLaw.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "PolynomialMaterialLaw3d.hpp"
#include "PoleZeroMaterialLaw.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "SchmidCostaExponentialLaw2d.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "ToyCompressibleMaterialLaw.hpp"
#include "CompressibleExponentialLaw.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestMaterialLaws : public CxxTest::TestSuite
{
private:

    // Helper method for testing T=0 when C=I and p=p_zero_strain (incompressible case)
    template<unsigned DIM>
    void CheckZeroStressWhenNoDeformation(AbstractIncompressibleMaterialLaw<DIM>* pLaw)
    {
        c_matrix<double,DIM,DIM> C = identity_matrix<double>(DIM);
        c_matrix<double,DIM,DIM> invC = identity_matrix<double>(DIM);

        c_matrix<double,DIM,DIM> T;
        FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;

        double pressure = pLaw->GetZeroStrainPressure();

        pLaw->ComputeStressAndStressDerivative(C,invC,pressure,T,dTdE,false);

        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned N=0; N<DIM; N++)
            {
                TS_ASSERT_DELTA(T(M,N), 0.0, 1e-6);
            }
        }
    }

    // Helper method for testing T=0 when C=I (compressible case)
    template<unsigned DIM>
    void CheckZeroStressWhenNoDeformation(AbstractCompressibleMaterialLaw<DIM>* pLaw)
    {
        c_matrix<double,DIM,DIM> C = identity_matrix<double>(DIM);
        c_matrix<double,DIM,DIM> invC = identity_matrix<double>(DIM);

        c_matrix<double,DIM,DIM> T;
        FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;

        pLaw->ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,false);

        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned N=0; N<DIM; N++)
            {
                TS_ASSERT_DELTA(T(M,N), 0.0, 1e-6);
            }
        }
    }

    // Helper method for testing that the analytic dTdE code is correct, by
    // computing a numerical derivative of T
    //
    // Compressible version
    template<unsigned DIM>
    void CheckDTdEComputation(AbstractCompressibleMaterialLaw<DIM>* pLaw)
    {
        c_matrix<double,DIM,DIM> C;
        c_matrix<double,DIM,DIM> invC;
        C(0,0) = 1.1;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 0.9;
        if (DIM==3)
        {
            C(2,2) = 0.95;
            C(0,2) = C(2,0) = 0;
            C(1,2) = C(2,1) = 0;
        }

        invC = Inverse(C);

        c_matrix<double,DIM,DIM> T_base;
        FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;

        pLaw->ComputeStressAndStressDerivative(C,invC,0.0,T_base,dTdE,false);

        double h=0.00001;

        for (unsigned M=0; M<DIM; M++)
        {
            C(M,M) += h;     // just change C00 and C11 (and C22). Can't see how to compute numerical
                             // derivative of wrt C01,C10, given that C is assumed symmetric
            invC = Inverse(C);

            c_matrix<double,DIM,DIM> T;

            pLaw->ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

            for (unsigned P=0; P<DIM; P++)
            {
                for (unsigned Q=0; Q<DIM; Q++)
                {
                    double dtdc = (T(P,Q)-T_base(P,Q))/h;
                    //std::cout << P << Q << M << M << " " << dTdE(P,Q,M,M) << "\n";
                    double tol = std::max(fabs(dTdE(P,Q,M,M))*1e-3, 1e-6);
                    TS_ASSERT_DELTA(2*dtdc, dTdE(P,Q,M,M), tol);
                }
            }

            C(M,M) -= h;
        }
    }

    // Helper method for testing that the analytic dTdE code is correct, by
    // computing a numerical derivative of T
    //
    // Incompressible version
    template<unsigned DIM>
    void CheckDTdEComputation(AbstractIncompressibleMaterialLaw<DIM>* pLaw)
    {
        c_matrix<double,DIM,DIM> C;
        c_matrix<double,DIM,DIM> invC;
        C(0,0) = 1.06;
        C(0,1) = C(1,0) = 0.106;
        C(1,1) = 0.954;            // overall C satifies det(C) = 1
        if (DIM==3)
        {
            C(2,2) = 1.0;
            C(0,2) = C(2,0) = 0;
            C(1,2) = C(2,1) = 0;
        }

        invC = Inverse(C);

        c_matrix<double,DIM,DIM> T_base;
        FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE;

        double pressure = 1.0;

        pLaw->ComputeStressAndStressDerivative(C,invC,pressure,T_base,dTdE,false);

        double h=0.00001;

        for (unsigned M=0; M<DIM; M++)
        {
            C(M,M) += h;     // just change C00 and C11 (and C22). Can't see how to compute numerical
                             // derivative of wrt C01,C10, given that C is assumed symmetric
            invC = Inverse(C);

            c_matrix<double,DIM,DIM> T;

            pLaw->ComputeStressAndStressDerivative(C,invC,1.0,T,dTdE,true);

            for (unsigned P=0; P<DIM; P++)
            {
                for (unsigned Q=0; Q<DIM; Q++)
                {
                    double dtdc = (T(P,Q)-T_base(P,Q))/h;
                    //std::cout << P << Q << M << M << " " << dTdE(P,Q,M,M) << "\n";
                    double tol = std::max(fabs(dTdE(P,Q,M,M))*1e-3, 1e-6);
                    TS_ASSERT_DELTA(2*dtdc, dTdE(P,Q,M,M), tol);
                }
            }

            C(M,M) -= h;
        }
    }

    // Helper method for testing change of basis (implemented for 2d only)
    void CheckChangeOfBasis(AbstractMaterialLaw<2>* pLaw)
    {
        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.2;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 1.1;
        invC = Inverse(C);

        c_matrix<double,2,2> T_Xfibres;
        c_matrix<double,2,2> T_Yfibres;
        FourthOrderTensor<2,2,2,2> dTdE_Xfibres;
        FourthOrderTensor<2,2,2,2> dTdE_Yfibres;

        double p = 1.0;
        if (dynamic_cast<AbstractCompressibleMaterialLaw<2>*>(pLaw) != NULL)
        {
            p = 0.0; // ie if compressible, then should give p=0
        }

        pLaw->ComputeStressAndStressDerivative(C,invC,p,T_Xfibres,dTdE_Xfibres,true); // no change of basis no fibres in X-dir

        // Now assume fibres in Y-dir. first set up equivalent C
        C(0,0) = 1.1;
        C(1,1) = 1.2;
        invC = Inverse(C);

        // Change of basis matrix
        c_matrix<double,2,2> P_basis;
        P_basis(0,0) = P_basis(1,1) = 0.0;
        P_basis(1,0) = P_basis(0,1) = 1.0;

        pLaw->SetChangeOfBasisMatrix(P_basis);
        pLaw->ComputeStressAndStressDerivative(C,invC,p,T_Yfibres,dTdE_Yfibres,true);

        TS_ASSERT_DELTA(T_Xfibres(0,0), T_Yfibres(1,1), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(1,1), T_Yfibres(0,0), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(0,1), T_Yfibres(1,0), 1e-8);
        TS_ASSERT_DELTA(T_Xfibres(1,0), T_Yfibres(0,1), 1e-8);

        // dTdE_Xfibres(0,1,1,0) should be equal to dTdE_Yfibres(1,0,0,1), etc
        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(dTdE_Xfibres(M,N,P,Q), dTdE_Yfibres((M+1)%2,(N+1)%2,(P+1)%2,(Q+1)%2), 1e-8);
                    }
                }
            }
        }

        pLaw->ResetToNoChangeOfBasisMatrix();
    }

public:

    void TestMooneyRivlinLaw()
    {
        TS_ASSERT_THROWS_THIS(MooneyRivlinMaterialLaw<2> bad_mr_law(-3.0),"c1 must be positive in mooney-rivlin");
        TS_ASSERT_THROWS_THIS(MooneyRivlinMaterialLaw<3> bad_mr_law2(3.0),"Two parameters needed for 3d Mooney-Rivlin");

        double c1 = 2.0;

        MooneyRivlinMaterialLaw<2> law_2d(c1);

        CheckZeroStressWhenNoDeformation<2>(&law_2d);
        CheckDTdEComputation<2>(&law_2d);

        CheckChangeOfBasis(&law_2d);

        TS_ASSERT_DELTA(law_2d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(1.0,0.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1(1.0,0.0), 0.0, 1e-12);

        double c2 = 3.0;

        MooneyRivlinMaterialLaw<3> law_3d(c1, c2);

        CheckZeroStressWhenNoDeformation<3>(&law_3d);
        CheckDTdEComputation<3>(&law_3d);

        TS_ASSERT_DELTA(law_3d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.GetC2(), c2, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI1(1.0,0.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI2(1.0,0.0), c2, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_d2W_dI1(1.0,0.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_d2W_dI2(1.0,0.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_d2W_dI1I2(1.0,0.0), 0.0, 1e-12);

        TS_ASSERT_DELTA(law_3d.GetZeroStrainPressure(), 2*c1+4*c2, 1e-12);

        // Compute stress given a non-zero deformation
        c_matrix<double,3,3> F;
        F(0,0) = 3.0;
        F(0,1) = 1.0;
        F(1,0) = -1.0;
        F(0,2) = 2.0;
        F(2,0) = 1.0;
        F(1,1) = 6.0;
        F(1,2) = -1.0;
        F(2,1) = 1.5;
        F(2,2) = 0.5;

        c_matrix<double,3,3> C = prod(trans(F),F);

        double I1 =  Trace(C);

        c_matrix<double,3,3> invC = Inverse(C);

        double pressure = 5.0;

        c_matrix<double,3,3> T;
        c_matrix<double,3,3> T2;
        c_matrix<double,3,3> S;
        c_matrix<double,3,3> sigma;

        FourthOrderTensor<3,3,3,3> dTdE;

        law_3d.ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, true);
        law_3d.Compute1stPiolaKirchoffStress(F,pressure,S);
        law_3d.Compute2ndPiolaKirchoffStress(C,pressure,T2);
        law_3d.ComputeCauchyStress(F,pressure,sigma);

        c_matrix<double,3,3> FT = prod(F,T);
        c_matrix<double,3,3> F_T_tranF_over_detF = (1.0/Determinant(F))*prod(FT,trans(F));//F*T_as_unsym_tensor*transpose(F);

        c_matrix<double,3,3> T_transposeF = prod(T,trans(F));//T_as_unsym_tensor*transpose(F);

        // Check sigma is correct - sigma should be (1/detF) F * T * trans(F)
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(sigma(i,j), F_T_tranF_over_detF(i,j), 1e-12);
            }
        }

        // Check S is correct
        for (unsigned M=0; M<3; M++)
        {
            for (unsigned i=0; i<3; i++)
            {
                TS_ASSERT_DELTA(S(M,i), T_transposeF(M,i), 1e-12);
            }
        }

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                // Check we gave a symmetric C
                assert(C(M,N)==C(N,M));

                // Check the stress
                TS_ASSERT_DELTA(T(M,N), (2*c1+2*c2*I1)*(M==N) - 2*c2*C(M,N) - pressure*invC(M,N), 1e-12);

                // Check alternative computation of the stress
                TS_ASSERT_DELTA(T(M,N), T2(M,N), 1e-12);

                for (unsigned P=0;P<3;P++)
                {
                    for (unsigned Q=0;Q<3;Q++)
                    {
                        double true_val =   4*c2*((M==N)*(P==Q)-(M==P)*(N==Q))
                                            + 2*pressure*invC(M,P)*invC(Q,N);

                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), true_val, 1e-12);
                    }
                }
            }
        }

        law_3d.ScaleMaterialParameters(10);
        TS_ASSERT_DELTA(law_3d.GetC1(), c1/10, 1e-12);
        TS_ASSERT_DELTA(law_3d.GetC2(), c2/10, 1e-12);
    }

    void TestExponentialLaw()
    {
        TS_ASSERT_THROWS_THIS(ExponentialMaterialLaw<2> bad_exp_law(-3.0,1),"a must be positive");

        double a = 2.0;
        double b = 3.0;
        double I1 = 2.0;
        double I2 = 1.0;

        ExponentialMaterialLaw<2> exp_law_2d(a,b);

        CheckZeroStressWhenNoDeformation<2>(&exp_law_2d);
        CheckDTdEComputation<2>(&exp_law_2d);

        TS_ASSERT_DELTA(exp_law_2d.GetA(), a, 1e-12);
        TS_ASSERT_DELTA(exp_law_2d.GetB(), b, 1e-12);

        TS_ASSERT_DELTA(exp_law_2d.Get_dW_dI1(I1,I2),  a*b*exp(b*(I1-2)),   1e-12);
        TS_ASSERT_DELTA(exp_law_2d.Get_d2W_dI1(I1,I2), b*exp_law_2d.Get_dW_dI1(I1,I2), 1e-12);

        ExponentialMaterialLaw<3> exp_law_3d(a,b);

        CheckZeroStressWhenNoDeformation<3>(&exp_law_3d);
        CheckDTdEComputation<3>(&exp_law_3d);

        TS_ASSERT_DELTA(exp_law_3d.GetA(), a, 1e-12);
        TS_ASSERT_DELTA(exp_law_3d.GetB(), b, 1e-12);

        TS_ASSERT_DELTA(exp_law_3d.Get_dW_dI1(I1,I2),   a*b*exp(b*(I1-3)),   1e-12);
        TS_ASSERT_DELTA(exp_law_3d.Get_dW_dI2(I1,I2),   0.0,                 1e-12);

        TS_ASSERT_DELTA(exp_law_3d.Get_d2W_dI1(I1,I2),  b*exp_law_3d.Get_dW_dI1(I1,I2), 1e-12);
        TS_ASSERT_DELTA(exp_law_3d.Get_d2W_dI2(I1,I2),  0.0,                 1e-12);
        TS_ASSERT_DELTA(exp_law_3d.Get_d2W_dI1I2(I1,I2),0.0,                 1e-12);
    }

    // The polynomial material law with N=1 is exactly a Mooney-Rivlin law and shouldagree
    void TestPolynomialMaterialLawAgainstMooneyRivlin()
    {
        unsigned N = 1;
        std::vector< std::vector<double> > alpha = PolynomialMaterialLaw3d::GetZeroedAlpha(N);

        // test GetZeroedAlpha
        TS_ASSERT_EQUALS(alpha.size(),2u);
        TS_ASSERT_EQUALS(alpha[0].size(),2u);
        TS_ASSERT_EQUALS(alpha[1].size(),2u);

        TS_ASSERT_DELTA(alpha[0][0],0.0,1e-12);
        TS_ASSERT_DELTA(alpha[0][1],0.0,1e-12);
        TS_ASSERT_DELTA(alpha[1][0],0.0,1e-12);
        TS_ASSERT_DELTA(alpha[1][1],0.0,1e-12);

        double c1 = 3.0;
        double c2 = 2.0;

        alpha[1][0] = c1;
        alpha[0][1] = c2;

        PolynomialMaterialLaw3d poly_mr_law(N,alpha);
        MooneyRivlinMaterialLaw<3> mooney_rivlin_law(c1,c2);

        CheckZeroStressWhenNoDeformation<3>(&mooney_rivlin_law);
        CheckDTdEComputation<3>(&mooney_rivlin_law);

        double I1 = 4;
        double I2 = 2.4;

        TS_ASSERT_DELTA(mooney_rivlin_law.Get_dW_dI1(I1,I2),    poly_mr_law.Get_dW_dI1(I1,I2),    1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_dW_dI2(I1,I2),    poly_mr_law.Get_dW_dI2(I1,I2),    1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_d2W_dI1(I1,I2),   poly_mr_law.Get_d2W_dI1(I1,I2),   1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_d2W_dI2(I1,I2),   poly_mr_law.Get_d2W_dI2(I1,I2),   1e-12);
        TS_ASSERT_DELTA(mooney_rivlin_law.Get_d2W_dI1I2(I1,I2), poly_mr_law.Get_d2W_dI1I2(I1,I2), 1e-12);
    }

    // Test the Polynomial Material Law with a quadratic law
    //   W = c20 (I1-3)^2  +  c11 (I1-3)(I2-3)  +  c02 (I2-3)^2   -   p C^{-1}/2
    //
    // We test ComputeStressAndStressDerivative() with this law because it uses all the
    // bits in that method.
    void TestQuadraticPolynomialLaw()
    {
        unsigned param_n = 2;
        std::vector< std::vector<double> > alpha = PolynomialMaterialLaw3d::GetZeroedAlpha(param_n);

        double c20 = 3.0;
        double c11 = 4.0;
        double c02 = 5.0;

        alpha[2][0] = c20;
        alpha[1][1] = c11;
        alpha[0][2] = c02;

        PolynomialMaterialLaw3d poly_law(param_n, alpha);

        c_matrix<double,3,3> C;
        C(0,0) = 3.0;
        C(0,1) = 1.0;
        C(1,0) = 1.0;
        C(0,2) = 2.0;
        C(2,0) = 2.0;
        C(1,1) = 6.0;
        C(1,2) = -1.0;
        C(2,1) = -1.0;
        C(2,2) = 0.5;

        double I1 =  Trace(C);
        double I2 =  SecondInvariant(C);

        double true_dWdI1    = 2*c20*(I1-3) +   c11*(I2-3);
        double true_dWdI2    =   c11*(I1-3) + 2*c02*(I2-3);
        double true_d2WdI1   = 2*c20;
        double true_d2WdI1I2 =   c11;
        double true_d2WdI2   = 2*c02;

        TS_ASSERT_DELTA(true_dWdI1,    poly_law.Get_dW_dI1(I1,I2),    1e-12);
        TS_ASSERT_DELTA(true_dWdI2,    poly_law.Get_dW_dI2(I1,I2),    1e-12);
        TS_ASSERT_DELTA(true_d2WdI1,   poly_law.Get_d2W_dI1(I1,I2),   1e-12);
        TS_ASSERT_DELTA(true_d2WdI2,   poly_law.Get_d2W_dI2(I1,I2),   1e-12);
        TS_ASSERT_DELTA(true_d2WdI1I2, poly_law.Get_d2W_dI1I2(I1,I2), 1e-12);

        c_matrix<double,3,3> invC = Inverse(C);

        double pressure = 5.0;

        c_matrix<double,3,3> T;
        FourthOrderTensor<3,3,3,3> dTdE;

        poly_law.ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, true);

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                // Check we gave a symmetric C
                assert(C(M,N)==C(N,M));

                double dI1_dC_MN = (M==N);
                double dI2_dC_MN = I1*(M==N)-C(M,N);

                double true_val =   2 * true_dWdI1 * dI1_dC_MN
                                    + 2 * true_dWdI2 * dI2_dC_MN
                                    - pressure*invC(M,N);

                TS_ASSERT_DELTA(  T(M,N), true_val, 1e-12 );

                for (unsigned P=0;P<3;P++)
                {
                    for (unsigned Q=0;Q<3;Q++)
                    {
                        //see above
                        dI1_dC_MN = (M==N);
                        double dI1_dC_PQ = (P==Q);

                        double d2I1_dC2  = 0;

                        //see above
                        dI2_dC_MN = I1*(M==N)-C(M,N);
                        double dI2_dC_PQ = I1*(P==Q)-C(P,Q);

                        double d2I2_dC2  = (M==N)*(P==Q)-(M==P)*(N==Q);

                        double d_invC_dC = -invC(M,P)*invC(Q,N);

                        true_val =    4 * true_d2WdI1 * dI1_dC_MN * dI1_dC_PQ
                                             + 4 * true_dWdI1  * d2I1_dC2
                                             + 4 * true_d2WdI2 * dI2_dC_MN * dI2_dC_PQ
                                             + 4 * true_dWdI2  * d2I2_dC2
                                             + 4 * true_d2WdI1I2 * (dI1_dC_MN*dI2_dC_PQ + dI1_dC_PQ*dI2_dC_MN)
                                             - 2 * pressure * d_invC_dC;

                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), true_val, 1e-12);
                    }
                }
            }
        }

        // Cover GetAlpha
        TS_ASSERT_DELTA(poly_law.GetAlpha(0,1),0.0,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(0,2),c02,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(1,1),c11,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(1,0),0.0,1e-12);
        TS_ASSERT_DELTA(poly_law.GetAlpha(2,0),c20,1e-12);

        // Check exception thrown if N=0
        TS_ASSERT_THROWS_THIS(PolynomialMaterialLaw3d bad_poly1(0,alpha),"n must be positive");

        // Check exception thrown if alpha is not correctly sized
        std::vector< std::vector<double> > bad_alpha(2);
        bad_alpha[0].resize(1);
        bad_alpha[1].resize(1);
        TS_ASSERT_THROWS_THIS(PolynomialMaterialLaw3d bad_poly2(2,bad_alpha),"alpha not big enough");

        // Check exception thrown if alpha[p][q]!=0, when p+q>N
        alpha[2][2] = 1.0;
        TS_ASSERT_THROWS_THIS(PolynomialMaterialLaw3d bad_poly3(2,alpha),"alpha[2][2] should be zero, as p+q > 2");

        /*
         * Compute the stress given C=delta_{MN} and p=zero_strain_pressure,
         * obviously it should be zero.
         */
        c_matrix<double,3,3> identity_strain_3d = identity_matrix<double>(3);
        c_matrix<double,3,3> T_3d;

        poly_law.Compute2ndPiolaKirchoffStress(identity_strain_3d,
                                               poly_law.GetZeroStrainPressure(),
                                               T_3d);
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(T_3d(i,j),0.0,1e-12);
            }
        }
    }

    void TestPoleZeroMaterialLaw()
    {
        std::vector<std::vector<double> > k(2),a(2),b(2);
        for (unsigned i=0; i<2; i++)
        {
            k[i].resize(2);
            a[i].resize(2);
            b[i].resize(2);
        }

        k[0][0] = 1;
        k[1][0] = k[0][1] = 2;
        k[1][1] = 3;

        a[0][0] = 4;
        a[1][0] = a[0][1] = 5;
        a[1][1] = 6;

        b[0][0] = 7;
        b[1][0] = b[0][1] = 6;
        b[1][1] = 5;

        PoleZeroMaterialLaw<2> pole_zero_law(k,a,b);

        CheckZeroStressWhenNoDeformation<2>(&pole_zero_law);
        CheckDTdEComputation<2>(&pole_zero_law);

        CheckChangeOfBasis(&pole_zero_law);

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;

        c_matrix<double,2,2> T;
        FourthOrderTensor<2,2,2,2> dTdE;

//// currently been changed on that pole-zero law DOESN'T return T=0 if E<0

//        C(0,0) = 0.5;
//        C(0,1) = -0.1;
//        C(1,0) = -0.1;
//        C(1,1) = 0.5;
//        invC = Inverse(C);
//
//        // C such that E_MN < 0, p=0 => T=0, dTdE=0;
//        pole_zero_law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);
//
//        for (unsigned M=0; M<2; M++)
//        {
//            for (unsigned N=0; N<2; N++)
//            {
//                TS_ASSERT_DELTA(T(M,N), 0.0, 1e-9);
//                for (unsigned P=0; P<2; P++)
//                {
//                    for (unsigned Q=0; Q<2; Q++)
//                    {
//                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), 0.0, 1e-9);
//                    }
//                }
//            }
//        }

        // non-trivial deformation, (checking all components have such that E_MN < a_MN)
        C(0,0) = 2;
        C(0,1) = C(1,0) = 2;
        C(1,1) = 5;
        invC = Inverse(C);

        pole_zero_law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        // T_MN = ke(2+be/(a-e))/(a-e)^b
        double t00 = 1*0.5*3/6433.92969;
        double t10 = 2.0*1.0*(2+6.0/4.0)/4096.0;
        double t11 = 3.0*2.0*(2+5*2.0/4.0)/1024.0;

        TS_ASSERT_DELTA( T(0,0), t00, 1e-9 );
        TS_ASSERT_DELTA( T(1,0), t10, 1e-9 );
        TS_ASSERT_DELTA( T(0,1), t10, 1e-9 );
        TS_ASSERT_DELTA( T(1,1), t11, 1e-9 );

        // Test dTdE
        double dtde00 = ( 2*3.5*3.5 + 14*3.5 + 0.25*7*8 )/78815.6387;
        double dtde10 = 2*( 2*4*4 + 4*1*6*4 + 42)/65536.0;
        double dtde11 = 3*( 2*4*4 + 4*2*5*4 + 120)/16384.0;

        TS_ASSERT_DELTA(dTdE(0,0,0,0), dtde00, 1e-9);
        TS_ASSERT_DELTA(dTdE(0,1,0,1), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,0,1,0), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,1,1,1), dtde11, 1e-9);

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        if ((P!=M) || (Q!=N))
                        {
                            TS_ASSERT_DELTA(dTdE(M,N,P,Q), 0.0, 1e-9);
                        }
                    }
                }
            }
        }

        /*
         * Test the pressure terms in the stress and stress-deriv, by calling with
         * p=0 and p=1 and verifying the difference is what it should be.
         */
        c_matrix<double,2,2> T2;
        FourthOrderTensor<2,2,2,2> dTdE2;
        pole_zero_law.ComputeStressAndStressDerivative(C, invC, 0.0, T,  dTdE,  true);
        pole_zero_law.ComputeStressAndStressDerivative(C, invC, 1.0, T2, dTdE2, true);

        for (unsigned M=0; M<2; M++)
        {
            for (unsigned N=0; N<2; N++)
            {
                TS_ASSERT_DELTA(T(M,N) - T2(M,N), invC(M,N), 1e-6);
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        TS_ASSERT_DELTA(dTdE(M,N,P,Q)-dTdE2(M,N,P,Q), -2*invC(M,P)*invC(Q,N), 1e-6);
                    }
                }
            }
        }
    }

    void TestPoleZeroMaterialLaw3d()
    {
        std::vector<std::vector<double> > k(3),a(3),b(3);
        for (unsigned i=0; i<3; i++)
        {
            k[i].resize(3);
            a[i].resize(3);
            b[i].resize(3);
        }

        k[0][0] = 1;
        k[1][0] = k[0][1] = 2;
        k[0][2] = k[2][0] = 4;
        k[1][1] = 3;
        k[1][2] = k[2][1] = 5;
        k[2][2] = 3;

        a[0][0] = 4;
        a[1][0] = a[0][1] = 5;
        a[2][0] = a[0][2] = 6;
        a[1][1] = 6;
        a[2][1] = a[1][2] = 4;
        a[2][2] = 9;

        b[0][0] = 7;
        b[1][0] = b[0][1] = 6;
        b[2][0] = b[0][2] = 2;
        b[1][1] = 5;
        b[2][1] = b[1][2] = 4;
        b[2][2] = 2;

        PoleZeroMaterialLaw<3> pole_zero_law(k,a,b);

        CheckZeroStressWhenNoDeformation<3>(&pole_zero_law);
        CheckDTdEComputation<3>(&pole_zero_law);

        c_matrix<double,3,3> C;
        C(0,0) = 2;
        C(0,1) = C(1,0) = 2;
        C(0,2) = C(2,0) = 3;
        C(1,1) = 5;
        C(1,2) = C(2,1) = 4;
        C(2,2) = 3;

        c_matrix<double,3,3> invC = Inverse(C);

        c_matrix<double,3,3> T;
        FourthOrderTensor<3,3,3,3> dTdE;

        pole_zero_law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        // Same as previous test, except the t22 and dtde22 bits are new
        double t00 = 1*0.5*3/6433.92969;
        double t10 = 2.0*1.0*(2+6.0/4.0)/4096.0;
        double t11 = 3.0*2.0*(2+5*2.0/4.0)/1024.0;
        double t22 = 3.0*(2+2.0/8.0)/64.0;

        TS_ASSERT_DELTA( T(0,0), t00, 1e-9 );
        TS_ASSERT_DELTA( T(1,0), t10, 1e-9 );
        TS_ASSERT_DELTA( T(0,1), t10, 1e-9 );
        TS_ASSERT_DELTA( T(1,1), t11, 1e-9 );
        TS_ASSERT_DELTA( T(2,2), t22, 1e-9 );

        // Test dTdE
        double dtde00 = ( 2*3.5*3.5 + 14*3.5 + 0.25*7*8 )/78815.6387;
        double dtde10 = 2*( 2*4*4 + 4*1*6*4 + 42)/65536.0;
        double dtde11 = 3*( 2*4*4 + 4*2*5*4 + 120)/16384.0;
        double dtde22 = 3*( 2*64 + 4*1*2*8 + 2*3)/4096.0;

        TS_ASSERT_DELTA(dTdE(0,0,0,0), dtde00, 1e-9);
        TS_ASSERT_DELTA(dTdE(0,1,0,1), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,0,1,0), dtde10, 1e-9);
        TS_ASSERT_DELTA(dTdE(1,1,1,1), dtde11, 1e-9);
        TS_ASSERT_DELTA(dTdE(2,2,2,2), dtde22, 1e-9);

        pole_zero_law.ScaleMaterialParameters(10);
        TS_ASSERT_DELTA(pole_zero_law.mK[0][0], 0.1, 1e-12);
        TS_ASSERT_DELTA(pole_zero_law.mA[0][0], 4, 1e-12);
        TS_ASSERT_DELTA(pole_zero_law.mB[0][0], 7, 1e-12);
    }

    void TestNashHunterPoleZeroLaw3d()
    {
        NashHunterPoleZeroLaw<3> law;

        CheckZeroStressWhenNoDeformation<3>(&law);
        CheckDTdEComputation<3>(&law);

        c_matrix<double,3,3> C;
        c_matrix<double,3,3> invC;
        C(0,0) = 1.2;
        C(0,1) = C(1,0) = 0.1;
        C(0,2) = C(2,0) = 0.3;
        C(1,1) = 1.1;
        C(1,2) = C(2,1) = -0.1;
        C(2,2) = 1.3;
        invC = Inverse(C);

        c_matrix<double,3,3> T;
        FourthOrderTensor<3,3,3,3> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        /*
         * Hard-coded test to test nothing changes (and stresses are of
         * the correct magnitude, which is dependent on whether the params
         * have been entered a Pa or KPa).
         */
        TS_ASSERT_DELTA(T(0,0),2.0902,1e-3);

        // Pick a P such that P =/= P^T
        c_matrix<double,3,3> basis = identity_matrix<double>(3);
        basis(0,0) = 1/sqrt(2.0);
        basis(1,0) = 1/sqrt(2.0);
        basis(0,1) = -1/sqrt(2.0);
        basis(1,1) = 1/sqrt(2.0);

        law.SetChangeOfBasisMatrix(basis);
        law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        /*
         * Carefully checked that the P's and P^T's are the right way
         * round in the transformations, this is just a hardcoded check
         * that nothing has changed.
         */
        TS_ASSERT_DELTA(T(0,0),1.7052,1e-3);
    }

    void TestNashHunterPoleZeroLaw2d()
    {
        NashHunterPoleZeroLaw<2> law;

        CheckZeroStressWhenNoDeformation<2>(&law);
        CheckDTdEComputation<2>(&law);

        CheckChangeOfBasis(&law);

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.2;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 1.1;
        invC = Inverse(C);

        c_matrix<double,2,2> T;
        FourthOrderTensor<2,2,2,2> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true);

        /*
         * Hard-coded test to test nothing changes (and stresses are of
         * the correct magnitude, which is dependent on whether the params
         * have been entered a Pa or KPa).
         */
        TS_ASSERT_DELTA(T(0,0),2.0902,1e-3);

        C(0,0) = 10;
        TS_ASSERT_THROWS_CONTAINS(law.ComputeStressAndStressDerivative(C,invC,0.0,T,dTdE,true), "strain unacceptably large");
    }

    void TestSchmidCostaExponentialLaw()
    {
        SchmidCostaExponentialLaw2d law;

        CheckZeroStressWhenNoDeformation<2>(&law);
        CheckDTdEComputation<2>(&law);

        CheckChangeOfBasis(&law);

        double a = law.GetA();
        assert(a > 0);
        double bff = law.GetB()[0][0];
        double bfs = law.GetB()[0][1];
        double bsf = law.GetB()[1][0];
        double bss = law.GetB()[1][1];

        TS_ASSERT(bsf == bfs);

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.1;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 0.9;
        invC = Inverse(C);

        c_matrix<double,2,2> T_base;
        FourthOrderTensor<2,2,2,2> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T_base,dTdE,false);

        double e00 = 0.5*(C(0,0)-1);
        double e01 = 0.5*C(0,1);
        double e11 = 0.5*(C(1,1)-1);
        double Q = bff*e00*e00 + 2*bfs*e01*e01 + bss*e11*e11;
        TS_ASSERT_DELTA(T_base(0,0), a*exp(Q)*bff*e00/2, 1e-9);
        TS_ASSERT_DELTA(T_base(0,1), a*exp(Q)*bfs*e01/2, 1e-9);
        TS_ASSERT_DELTA(T_base(1,0), a*exp(Q)*bsf*e01/2, 1e-9);
        TS_ASSERT_DELTA(T_base(1,1), a*exp(Q)*bss*e11/2, 1e-9);
    }

    /*
     * Uses the simple material law W(I1,I2,I3) = c1(I1-3) + c2(I2-3) + c3(I3-1),
     * which may not correspond to a physically acceptable law but can still be
     * used to test the code.
     */
    void TestCompressibleLawsUsingToyCompressibleMaterialLaw()
    {
        double c1 = 2.0;

        ToyCompressibleMaterialLaw<2> law_2d(c1, 0.0, -c1);

        CheckZeroStressWhenNoDeformation<2>(&law_2d);
        CheckDTdEComputation<2>(&law_2d);

        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(2.0,1.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI3(2.0,1.0,1.0), -c1, 1e-12);

        double c2 = 3.0;

        ToyCompressibleMaterialLaw<3> law_3d(c1, c2, -c1-2*c2); // this choice of c3 gives zero stress for zero strain
        CheckZeroStressWhenNoDeformation<3>(&law_3d);
        CheckDTdEComputation<3>(&law_3d);

        TS_ASSERT_DELTA(law_3d.Get_dW_dI1(2.0,2.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI2(2.0,2.0,1.0), c2, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI3(2.0,2.0,1.0), -c1-2*c2, 1e-12);

        // Compute stress given a non-zero deformation
        c_matrix<double,3,3> F;
        F(0,0) = 3.0;
        F(0,1) = 1.0;
        F(1,0) = -1.0;
        F(0,2) = 2.0;
        F(2,0) = 1.0;
        F(1,1) = 6.0;
        F(1,2) = -1.0;
        F(2,1) = 1.5;
        F(2,2) = 0.5;

        c_matrix<double,3,3> C = prod(trans(F),F);

        double I1 =  Trace(C);

        c_matrix<double,3,3> invC = Inverse(C);

        double I3 = Determinant(C);

        c_matrix<double,3,3> T;
        c_matrix<double,3,3> T2;
        c_matrix<double,3,3> S;
        c_matrix<double,3,3> sigma;

        FourthOrderTensor<3,3,3,3> dTdE;

        law_3d.ComputeStressAndStressDerivative(C, invC, 0.0, T, dTdE, true);
        law_3d.Compute1stPiolaKirchoffStress(F,0.0,S);
        law_3d.Compute2ndPiolaKirchoffStress(C,0.0,T2);
        law_3d.ComputeCauchyStress(F,0.0,sigma);

        c_matrix<double,3,3> FT = prod(F,T);
        c_matrix<double,3,3> F_T_tranF_over_detF = (1.0/Determinant(F))*prod(FT,trans(F));//F*T_as_unsym_tensor*transpose(F);

        c_matrix<double,3,3> T_transposeF = prod(T,trans(F));//T_as_unsym_tensor*transpose(F);

        // Check sigma is correct - sigma should be (1/detF) F * T * trans(F)
        for (unsigned i=0; i<3; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_DELTA(sigma(i,j), F_T_tranF_over_detF(i,j), 1e-12);
            }
        }

        // Check S is correct
        for (unsigned M=0; M<3; M++)
        {
            for (unsigned i=0; i<3; i++)
            {
                TS_ASSERT_DELTA(S(M,i), T_transposeF(M,i), 1e-12);
            }
        }

        for (unsigned M=0; M<3; M++)
        {
            for (unsigned N=0; N<3; N++)
            {
                // Check we gave a symmetric C
                assert(C(M,N)==C(N,M));

                // Check the stress
                TS_ASSERT_DELTA(T(M,N), (2*c1+2*c2*I1)*(M==N) - 2*c2*C(M,N) - 2*(c1+2*c2)*I3*invC(M,N), 1e-12);

                // Check alternative computation of the stress
                TS_ASSERT_DELTA(T(M,N), T2(M,N), 1e-12);

                for (unsigned P=0;P<3;P++)
                {
                    for (unsigned Q=0;Q<3;Q++)
                    {
                        double true_val =   4*c2*((M==N)*(P==Q)-(M==P)*(N==Q))
                                          - 4*(c1+2*c2)*I3*(invC(M,N)*invC(P,Q) - invC(M,P)*invC(Q,N) );

                        TS_ASSERT_DELTA(dTdE(M,N,P,Q), true_val, 1e-12);
                    }
                }
            }
        }
    }

    void TestCompressibleMooneyRivlinLaw()
    {
        double c1 = 3.0;
        double c3 = 2.0;

        CompressibleMooneyRivlinMaterialLaw<2> law_2d(c1, c3);

        CheckZeroStressWhenNoDeformation<2>(&law_2d);
        CheckDTdEComputation<2>(&law_2d);

        CheckChangeOfBasis(&law_2d);

        TS_ASSERT_DELTA(law_2d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.GetC3(), c3, 1e-12);

        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(2.0,2.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI2(2.0,2.0,1.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI3(2.0,2.0,1.0), -c1, 1e-12);

        double i1 = 2.4;
        double i3 = 0.6;
        TS_ASSERT_DELTA(law_2d.Get_dW_dI1(i1,2.0,i3), c1/sqrt(i3), 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI2(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_dW_dI3(i1,2.0,i3), -c1*i1*pow(i3,-1.5)/2 + c3 - c3*pow(i3,-0.5), 1e-12);

        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1I3(i1,2.0,i3), -0.5*c1*pow(i3,-1.5), 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI3(i1,2.0,i3),0.75*c1*i1*pow(i3,-2.5) + 0.5*c3*pow(i3,-1.5), 1e-12);

        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI2(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI1I2(i1,2.0,i3), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_2d.Get_d2W_dI2I3(i1,2.0,i3), 0.0, 1e-12);

        law_2d.ScaleMaterialParameters(10);
        TS_ASSERT_DELTA(law_2d.GetC1(), c1/10, 1e-12);
        TS_ASSERT_DELTA(law_2d.GetC3(), c3/10, 1e-12);

        CompressibleMooneyRivlinMaterialLaw<3> law_3d(c1, c3);

        CheckZeroStressWhenNoDeformation<3>(&law_3d);
        CheckDTdEComputation<3>(&law_3d);

        TS_ASSERT_DELTA(law_3d.GetC1(), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.GetC3(), c3, 1e-12);

        TS_ASSERT_DELTA(law_3d.Get_dW_dI1(3.0,3.0,1.0), c1, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI2(3.0,3.0,1.0), 0.0, 1e-12);
        TS_ASSERT_DELTA(law_3d.Get_dW_dI3(3.0,3.0,1.0), -c1, 1e-12);
    }

    void TestCompressibleExponentialLaw()
    {
        CompressibleExponentialLaw<2> law;

        CheckZeroStressWhenNoDeformation<2>(&law);
        CheckDTdEComputation<2>(&law);

        CheckChangeOfBasis(&law);

        double a = law.GetA();
        assert(a > 0);
        double bff = law.GetB()[0][0]; // i.e. b_{fibre,fibre}
        double bfs = law.GetB()[0][1];
        double bsf = law.GetB()[1][0];
        double bss = law.GetB()[1][1];

        double c = law.GetCompressibilityParam();

        TS_ASSERT(bsf == bfs);

        c_matrix<double,2,2> C;
        c_matrix<double,2,2> invC;
        C(0,0) = 1.1;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 0.9;
        invC = Inverse(C);

        double I3 = Determinant(C);
        double J = sqrt(I3);
        double w3  = c*log(J)/(2*J);

        c_matrix<double,2,2> T_base;
        FourthOrderTensor<2,2,2,2> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T_base,dTdE,false);

        double e00 = 0.5*(C(0,0)-1);
        double e01 = 0.5*C(0,1);
        double e11 = 0.5*(C(1,1)-1);
        double Q = bff*e00*e00 + 2*bfs*e01*e01 + bss*e11*e11;
        TS_ASSERT_DELTA(T_base(0,0), a*exp(Q)*bff*e00 + 2*w3*I3*invC(0,0), 1e-9);
        TS_ASSERT_DELTA(T_base(0,1), a*exp(Q)*bfs*e01 + 2*w3*I3*invC(0,1), 1e-9);
        TS_ASSERT_DELTA(T_base(1,0), a*exp(Q)*bsf*e01 + 2*w3*I3*invC(1,0), 1e-9);
        TS_ASSERT_DELTA(T_base(1,1), a*exp(Q)*bss*e11 + 2*w3*I3*invC(1,1), 1e-9);

        CheckDTdEComputation<2>(&law);
    }

    void TestCompressibleExponentialLaw3d()
    {
        CompressibleExponentialLaw<3> law;

        CheckZeroStressWhenNoDeformation<3>(&law);
        CheckDTdEComputation<3>(&law);

        double a = law.GetA();
        assert(a > 0);
        double bff = law.GetB()[0][0]; // i.e. b_{fibre,fibre}
        double bfs = law.GetB()[0][1];
        double bsf = law.GetB()[1][0];
        double bss = law.GetB()[1][1];

        double c = law.GetCompressibilityParam();

        double bnn = law.GetB()[2][2];
        double bsn = law.GetB()[1][2];
        double bns = law.GetB()[2][1];
        double bfn = law.GetB()[0][2];
        double bnf = law.GetB()[2][0];

        TS_ASSERT(bsf == bfs);
        TS_ASSERT(bfn == bnf);
        TS_ASSERT(bsn == bns);

        c_matrix<double,3,3> C;
        c_matrix<double,3,3> invC;
        C(0,0) = 1.1;
        C(0,1) = C(1,0) = 0.1;
        C(1,1) = 0.9;
        C(0,2) = C(2,0) = 0.05;
        C(1,2) = C(2,1) = 0.01;
        C(2,2) = 0.95;
        invC = Inverse(C);

        double I3 = Determinant(C);
        double J = sqrt(I3);
        double w3  = c*log(J)/(2*J);

        c_matrix<double,3,3> T_base;
        FourthOrderTensor<3,3,3,3> dTdE;

        law.ComputeStressAndStressDerivative(C,invC,0.0,T_base,dTdE,false);

        double e00 = 0.5*(C(0,0)-1);
        double e11 = 0.5*(C(1,1)-1);
        double e22 = 0.5*(C(2,2)-1);
        double e01 = 0.5*C(0,1);
        double e12 = 0.5*C(1,2);
        double e02 = 0.5*C(0,2);
        double Q = bff*e00*e00 + bss*e11*e11 + bnn*e22*e22 + 2*bfs*e01*e01 + 2*bfn*e02*e02 + 2*bsn*e12*e12;
        TS_ASSERT_DELTA(T_base(0,0), a*exp(Q)*bff*e00 + 2*w3*I3*invC(0,0), 1e-9);
        TS_ASSERT_DELTA(T_base(0,1), a*exp(Q)*bfs*e01 + 2*w3*I3*invC(0,1), 1e-9);
        TS_ASSERT_DELTA(T_base(1,0), a*exp(Q)*bsf*e01 + 2*w3*I3*invC(1,0), 1e-9);
        TS_ASSERT_DELTA(T_base(1,1), a*exp(Q)*bss*e11 + 2*w3*I3*invC(1,1), 1e-9);

        TS_ASSERT_DELTA(T_base(0,2), a*exp(Q)*bfs*e02 + 2*w3*I3*invC(0,2), 1e-9);
        TS_ASSERT_DELTA(T_base(1,2), a*exp(Q)*bsf*e12 + 2*w3*I3*invC(1,2), 1e-9);
        TS_ASSERT_DELTA(T_base(2,2), a*exp(Q)*bss*e22 + 2*w3*I3*invC(2,2), 1e-9);
    }
};

#endif /*TESTMATERIALLAWS_HPP_*/
