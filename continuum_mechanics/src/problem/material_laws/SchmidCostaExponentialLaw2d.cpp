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

#include "SchmidCostaExponentialLaw2d.hpp"

SchmidCostaExponentialLaw2d::SchmidCostaExponentialLaw2d()
{
    mA = 0.221;    // kiloPascals, presumably, although the paper doesn't say.
                   // gives results matching Pole-zero anyway.
                   // Obtained from Table 1 of Schmid reference (see class doxygen), the mu (mean) value.

    double bff = 42.5; // dimensionless
    double bfs = 11.0; // dimensionless
    double bss = 18.6; // dimensionless

    mB.resize(2);
    mB[0].resize(2);
    mB[1].resize(2);

    mB[0][0] = bff;
    mB[0][1] = bfs;
    mB[1][0] = bfs;
    mB[1][1] = bss;

    for (unsigned M=0; M<2; M++)
    {
        for (unsigned N=0; N<2; N++)
        {
            mIdentity(M,N) = M==N ? 1.0 : 0.0;
        }
    }
}

void SchmidCostaExponentialLaw2d::ComputeStressAndStressDerivative(c_matrix<double,2,2>& rC,
                                                                   c_matrix<double,2,2>& rInvC,
                                                                   double                pressure,
                                                                   c_matrix<double,2,2>& rT,
                                                                   FourthOrderTensor<2,2,2,2>& rDTdE,
                                                                   bool                  computeDTdE)
{
    static c_matrix<double,2,2> C_transformed;
    static c_matrix<double,2,2> invC_transformed;

    // The material law parameters are set up assuming the fibre direction is (1,0,0)
    // and sheet direction is (0,1,0), so we have to transform C,inv(C),and T.
    // Let P be the change-of-basis matrix P = (\mathbf{m}_f, \mathbf{m}_s, \mathbf{m}_n).
    // The transformed C for the fibre/sheet basis is C* = P^T C P.
    // We then compute T* = T*(C*), and then compute T = P T* P^T.

    ComputeTransformedDeformationTensor(rC, rInvC, C_transformed, invC_transformed);

    // Compute T*

    c_matrix<double,2,2> E = 0.5*(C_transformed - mIdentity);

    double QQ = 0;
    for (unsigned M=0; M<2; M++)
    {
        for (unsigned N=0; N<2; N++)
        {
            QQ += mB[M][N]*E(M,N)*E(M,N);
        }
    }

    double multiplier = mA*exp(QQ)/2;
    rDTdE.Zero();

    for (unsigned M=0; M<2; M++)
    {
        for (unsigned N=0; N<2; N++)
        {
            rT(M,N) = multiplier*mB[M][N]*E(M,N) - pressure*invC_transformed(M,N);

            if (computeDTdE)
            {
                for (unsigned P=0; P<2; P++)
                {
                    for (unsigned Q=0; Q<2; Q++)
                    {
                        rDTdE(M,N,P,Q) =   multiplier * mB[M][N] * (M==P)*(N==Q)
                                        +  2*multiplier*mB[M][N]*mB[P][Q]*E(M,N)*E(P,Q)
                                        +  2*pressure*invC_transformed(M,P)*invC_transformed(Q,N);
                    }
                }
            }
        }
    }

    // Now do:   T = P T* P^T   and   dTdE_{MNPQ}  =  P_{Mm}P_{Nn}P_{Pp}P_{Qq} dT*dE*_{mnpq}
    this->TransformStressAndStressDerivative(rT, rDTdE, computeDTdE);
}

double SchmidCostaExponentialLaw2d::GetA()
{
    return mA;
}

std::vector<std::vector<double> > SchmidCostaExponentialLaw2d::GetB()
{
    return mB;
}

double SchmidCostaExponentialLaw2d::GetZeroStrainPressure()
{
    return 0.0;
}
