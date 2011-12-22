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

    ComputeTransformedDeformationTensor(rC, rInvC, C_transformed, invC_transformed);

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

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class PoleZeroMaterialLaw<2>;
template class PoleZeroMaterialLaw<3>;
