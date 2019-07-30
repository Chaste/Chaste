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

#include "PoleZeroMaterialLaw.hpp"

template<unsigned DIM>
PoleZeroMaterialLaw<DIM>::PoleZeroMaterialLaw()
{
}

template<unsigned DIM>
void PoleZeroMaterialLaw<DIM>::SetParameters(std::vector<std::vector<double> > k,
                                             std::vector<std::vector<double> > a,
                                             std::vector<std::vector<double> > b)
{
    if (DIM!=2 && DIM !=3)
    {
        EXCEPTION("Can only have a 2 or 3d incompressible pole-zero law");
    }

    assert(k.size()==DIM);
    assert(a.size()==DIM);
    assert(b.size()==DIM);

    for (unsigned i=0; i<DIM; i++)
    {
        assert(k[i].size()==DIM);
        assert(a[i].size()==DIM);
        assert(b[i].size()==DIM);

        for (unsigned j=0; j<DIM; j++)
        {
            assert( k[i][j] = k[j][i] );
            assert( a[i][j] = a[j][i] );
            assert( b[i][j] = b[j][i] );
        }
    }

    mK = k;
    mA = a;
    mB = b;

    for (unsigned M=0; M<DIM; M++)
    {
        for (unsigned N=0; N<DIM; N++)
        {
            mIdentity(M,N) = M==N ? 1.0 : 0.0;
        }
    }
}

template<unsigned DIM>
PoleZeroMaterialLaw<DIM>::PoleZeroMaterialLaw(std::vector<std::vector<double> > k,
                                              std::vector<std::vector<double> > a,
                                              std::vector<std::vector<double> > b)
{
    SetParameters(k,a,b);
}

template<unsigned DIM>
void PoleZeroMaterialLaw<DIM>::ComputeStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                                                c_matrix<double,DIM,DIM>& rInvC,
                                                                double                    pressure,
                                                                c_matrix<double,DIM,DIM>& rT,
                                                                FourthOrderTensor<DIM,DIM,DIM,DIM>&   rDTdE,
                                                                bool                      computeDTdE)
{
    static c_matrix<double,DIM,DIM> C_transformed;
    static c_matrix<double,DIM,DIM> invC_transformed;

    // The material law parameters are set up assuming the fibre direction is (1,0,0)
    // and sheet direction is (0,1,0), so we have to transform C,inv(C),and T.
    // Let P be the change-of-basis matrix P = (\mathbf{m}_f, \mathbf{m}_s, \mathbf{m}_n).
    // The transformed C for the fibre/sheet basis is C* = P^T C P.
    // We then compute T* = T*(C*), and then compute T = P T* P^T.

    this->ComputeTransformedDeformationTensor(rC, rInvC, C_transformed, invC_transformed);

    // compute T*

    c_matrix<double,DIM,DIM> E = 0.5*(C_transformed - mIdentity);

    for (unsigned M=0; M<DIM; M++)
    {
        for (unsigned N=0; N<DIM; N++)
        {
            double e = E(M,N);
            {
                double b = mB[M][N];
                double a = mA[M][N];
                double k = mK[M][N];

                //if this fails one of the strain values got too large for the law
                if (e>=a)
                {
                    EXCEPTION("E_{MN} >= a_{MN} - strain unacceptably large for model");
                }

                rT(M,N) =   k
                          * e
                          * (2*(a-e) + b*e)
                          * pow(a-e,-b-1)
                          - pressure*invC_transformed(M,N);
            }
        }
    }

    if (computeDTdE)
    {
        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned N=0; N<DIM; N++)
            {
                for (unsigned P=0; P<DIM; P++)
                {
                    for (unsigned Q=0; Q<DIM; Q++)
                    {
                        rDTdE(M,N,P,Q) = 2 * pressure * invC_transformed(M,P) * invC_transformed(Q,N);
                    }
                }

                double e = E(M,N);

                double b = mB[M][N];
                double a = mA[M][N];
                double k = mK[M][N];

                rDTdE(M,N,M,N) +=   k
                                  * pow(a-e, -b-2)
                                  * (
                                       2*(a-e)*(a-e)
                                     + 4*b*e*(a-e)
                                     + b*(b+1)*e*e
                                    );
            }
        }
    }

    // Now do:   T = P T* P^T   and   dTdE_{MNPQ}  =  P_{Mm}P_{Nn}P_{Pp}P_{Qq} dT*dE*_{mnpq}
    this->TransformStressAndStressDerivative(rT, rDTdE, computeDTdE);
}

template<unsigned DIM>
double PoleZeroMaterialLaw<DIM>::GetZeroStrainPressure()
{
    return 0.0;
}

template<unsigned DIM>
void PoleZeroMaterialLaw<DIM>::ScaleMaterialParameters(double scaleFactor)
{
    assert(scaleFactor > 0.0);
    for (unsigned i=0; i<mK.size(); i++)
    {
        for (unsigned j=0; j<mK[i].size(); j++)
        {
            mK[i][j] /= scaleFactor;
        }
    }
}

// Explicit instantiation
template class PoleZeroMaterialLaw<2>;
template class PoleZeroMaterialLaw<3>;
