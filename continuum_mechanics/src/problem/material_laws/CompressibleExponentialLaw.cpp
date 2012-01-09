/*

Copyright (C) University of Oxford, 2005-2012

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

    ComputeTransformedDeformationTensor(rC, rInvC, C_transformed, invC_transformed);

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

    double multiplier = mA*exp(QQ)/2;
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

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class CompressibleExponentialLaw<2>;
template class CompressibleExponentialLaw<3>;
