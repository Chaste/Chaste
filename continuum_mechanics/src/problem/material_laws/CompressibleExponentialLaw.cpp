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

#include "CompressibleExponentialLaw.hpp"

template<unsigned DIM>
CompressibleExponentialLaw<DIM>::CompressibleExponentialLaw()
{
    mA = 0.88;  // kPa

    double bff = 18.5; // dimensionless
    double bss = 3.58; // dimensionless
    double bnn = 3.58; // dimensionless
    double bfn = 2.8;  // etc
    double bfs = 2.8;
    double bsn = 2.8;

    mCompressibilityParam = 100.0;

    mB.resize(DIM);
    for (unsigned i=0; i<DIM; i++)
    {
        mB[i].resize(DIM);
    }

    mB[0][0] = bff;
    mB[0][1] = mB[1][0] = bfs;
    mB[1][1] = bss;

    if (DIM > 2)
    {
        mB[2][2] = bnn;
        mB[0][2] = mB[2][0] = bfn;
        mB[2][1] = mB[1][2] = bsn;
    }

    for (unsigned M=0; M<DIM; M++)
    {
        for (unsigned N=0; N<DIM; N++)
        {
            mIdentity(M,N) = M==N ? 1.0 : 0.0;
        }
    }
}

template<unsigned DIM>
void CompressibleExponentialLaw<DIM>::ComputeStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                                                       c_matrix<double,DIM,DIM>& rInvC,
                                                                       double                pressure /* not used */,
                                                                       c_matrix<double,DIM,DIM>& rT,
                                                                       FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                                                       bool                  computeDTdE)
{
    static c_matrix<double,DIM,DIM> C_transformed;
    static c_matrix<double,DIM,DIM> invC_transformed;

    // The material law parameters are set up assuming the fibre direction is (1,0,0)
    // and sheet direction is (0,1,0), so we have to transform C,inv(C),and T.
    // Let P be the change-of-basis matrix P = (\mathbf{m}_f, \mathbf{m}_s, \mathbf{m}_n).
    // The transformed C for the fibre/sheet basis is C* = P^T C P.
    // We then compute T* = T*(C*), and then compute T = P T* P^T.

    this->ComputeTransformedDeformationTensor(rC, rInvC, C_transformed, invC_transformed);

    // Compute T*

    c_matrix<double,DIM,DIM> E = 0.5*(C_transformed - mIdentity);

    double QQ = 0;
    for (unsigned M=0; M<DIM; M++)
    {
        for (unsigned N=0; N<DIM; N++)
        {
            QQ += mB[M][N]*E(M,N)*E(M,N);
        }
    }
    assert(QQ < 10.0);///\todo #2193 This line is to trap for large deformations which lead to blow up in the exponential Uysk model
    double multiplier = mA*exp(QQ);
    rDTdE.Zero();

    double J = sqrt(Determinant(rC));

    for (unsigned M=0; M<DIM; M++)
    {
        for (unsigned N=0; N<DIM; N++)
        {
            rT(M,N) = multiplier*mB[M][N]*E(M,N) + mCompressibilityParam * J*log(J)*invC_transformed(M,N);

            if (computeDTdE)
            {
                for (unsigned P=0; P<DIM; P++)
                {
                    for (unsigned Q=0; Q<DIM; Q++)
                    {
                        rDTdE(M,N,P,Q) =    multiplier * mB[M][N] * (M==P)*(N==Q)
                                         +  2*multiplier*mB[M][N]*mB[P][Q]*E(M,N)*E(P,Q)
                                         +  mCompressibilityParam * (J*log(J) + J) * invC_transformed(M,N) * invC_transformed(P,Q)
                                         -  mCompressibilityParam * 2*J*log(J) * invC_transformed(M,P) * invC_transformed(Q,N);
                    }
                }
            }
        }
    }

    // Now do:   T = P T* P^T   and   dTdE_{MNPQ}  =  P_{Mm}P_{Nn}P_{Pp}P_{Qq} dT*dE*_{mnpq}
    this->TransformStressAndStressDerivative(rT, rDTdE, computeDTdE);
}

// Explicit instantiation
template class CompressibleExponentialLaw<2>;
template class CompressibleExponentialLaw<3>;
