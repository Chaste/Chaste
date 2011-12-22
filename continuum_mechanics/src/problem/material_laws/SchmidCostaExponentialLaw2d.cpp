/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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

    double Q = 0;
    for (unsigned M=0; M<2; M++)
    {
        for (unsigned N=0; N<2; N++)
        {
            Q += mB[M][N]*E(M,N)*E(M,N);
        }
    }

    double multiplier = mA*exp(Q)/2;
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
