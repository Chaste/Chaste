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

#include "AbstractIsotropicIncompressibleMaterialLaw.hpp"

template<unsigned DIM>
AbstractIsotropicIncompressibleMaterialLaw<DIM>::~AbstractIsotropicIncompressibleMaterialLaw()
{
}

template<unsigned DIM>
void AbstractIsotropicIncompressibleMaterialLaw<DIM>::ComputeStressAndStressDerivative(
        c_matrix<double,DIM,DIM>& rC,
        c_matrix<double,DIM,DIM>& rInvC,
        double                    pressure,
        c_matrix<double,DIM,DIM>& rT,
        FourthOrderTensor<DIM,DIM,DIM,DIM>&   rDTdE,
        bool                      computeDTdE)
{
    /*
     * This is covered, but gcov doesn't see this as being covered
     * for some reason, maybe because of optimisations.
     */
    #define COVERAGE_IGNORE
    assert((DIM==2) || (DIM==3));
    #undef COVERAGE_IGNORE

    static c_matrix<double,DIM,DIM> identity = identity_matrix<double>(DIM);

    double I1 = Trace(rC);
    double I2 = SecondInvariant(rC);

    double  w1 = Get_dW_dI1(I1, I2);
    double  w2; // only computed if DIM==3

    // Compute stress:  **** See FiniteElementImplementations document. ****
    //
    //  T = dW_dE
    //    = 2 * w1 * dI1_dC_MN   +   2 * w2 * dI1_dC_MN  -  p * invC
    //    = 2 * w1 * delta_MN    +   2 * w2 * (I1 delta_MN - C_MN)  -  p * invC
    //
    //  (where w1 = dW/dI1, etc).

    rT = 2*w1*identity - pressure*rInvC;
    if (DIM==3)
    {
        w2 = Get_dW_dI2(I1, I2);
        rT += 2*w2*(I1*identity - rC);
    }

    // Compute stress derivative if required:   **** See FiniteElementImplementations document. ****
    //
    // The stress derivative dT_{MN}/dE_{PQ} is
    //
    //  dT_dE =    4 * w11 * dI1_dC_MN * dI1_dC_PQ
    //           + 4 * w1  * d2I1_dC2
    //           + 4 * w22 * dI2_dC_MN * dI2_dC_PQ
    //           + 4 * w2  * d2I2_dC2
    //           + 4 * w12 * (dI1_dC_MN*dI2_dC_PQ + dI1_dC_PQ*dI2_dC_MN)
    //           - 2 * pressure * d_invC_dC;
    //
    // where
    //   dI1_dC_MN = (M==N); // ie delta_{MN}
    //   dI1_dC_PQ = (P==Q);
    //   d2I1_dC2  = 0;
    //
    //   dI2_dC_MN = I1*(M==N)-C[M][N];
    //   dI2_dC_PQ = I1*(P==Q)-C[P][Q];
    //   d2I2_dC2  = (M==N)*(P==Q)-(M==P)*(N==Q);
    //
    //   d_invC_dC = -invC[M][P]*invC[Q][N];
    //
    if (computeDTdE)
    {
        double w11 = Get_d2W_dI1(I1,I2);

        double w12;
        double w22;

        if (DIM==3)
        {
            w22 = Get_d2W_dI2(I1, I2);
            w12 = Get_d2W_dI1I2(I1, I2);
        }

        for (unsigned M=0; M<DIM; M++)
        {
            for (unsigned N=0; N<DIM; N++)
            {
                for (unsigned P=0; P<DIM; P++)
                {
                    for (unsigned Q=0; Q<DIM; Q++)
                    {
                        rDTdE(M,N,P,Q) =   4 * w11  * (M==N) * (P==Q)
                                         + 2 * pressure * rInvC(M,P) * rInvC(Q,N);

                        if (DIM==3)
                        {
                            rDTdE(M,N,P,Q) +=   4 * w22   * (I1*(M==N) - rC(M,N)) * (I1*(P==Q) - rC(P,Q))
                                              + 4 * w2    * ((M==N)*(P==Q) - (M==P)*(N==Q))
                                              + 4 * w12 * ((M==N)*(I1*(P==Q) - rC(P,Q)) + (P==Q)*(I1*(M==N) - rC(M,N)));
                        }
                    }
                }
            }
        }
    }
}

template<>
double AbstractIsotropicIncompressibleMaterialLaw<2>::GetZeroStrainPressure()
{
    return 2*Get_dW_dI1(2,0);
}

template<>
double AbstractIsotropicIncompressibleMaterialLaw<3>::GetZeroStrainPressure()
{
    return 2*Get_dW_dI1(3,3) + 4*Get_dW_dI2(3,3);
}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class AbstractIsotropicIncompressibleMaterialLaw<2>;
template class AbstractIsotropicIncompressibleMaterialLaw<3>;
