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

#include "AbstractMaterialLaw.hpp"

template<unsigned DIM>
AbstractMaterialLaw<DIM>::AbstractMaterialLaw()
    : mpChangeOfBasisMatrix(nullptr)
{
}

template<unsigned DIM>
void AbstractMaterialLaw<DIM>::ComputeCauchyStress(c_matrix<double,DIM,DIM>& rF,
                                                   double pressure,
                                                   c_matrix<double,DIM,DIM>& rSigma)
{
    double detF = Determinant(rF);

    c_matrix<double,DIM,DIM> C = prod(trans(rF), rF);
    c_matrix<double,DIM,DIM> invC = Inverse(C);

    c_matrix<double,DIM,DIM> T;

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE; // not filled in, made static for efficiency

    ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, false);

    /*
     * Looping is probably more eficient then doing rSigma = (1/detF)*rF*T*transpose(rF),
     * which doesn't seem to compile anyway, as rF is a Tensor<2,DIM> and T is a
     * SymmetricTensor<2,DIM>.
     */
    for (unsigned i=0; i<DIM; i++)
    {
        for (unsigned j=0; j<DIM; j++)
        {
            rSigma(i,j) = 0.0;
            for (unsigned M=0; M<DIM; M++)
            {
                for (unsigned N=0; N<DIM; N++)
                {
                    rSigma(i,j) += rF(i,M)*T(M,N)*rF(j,N);
                }
            }
            rSigma(i,j) /= detF;
        }
    }
}

template<unsigned DIM>
void AbstractMaterialLaw<DIM>::Compute1stPiolaKirchoffStress(c_matrix<double,DIM,DIM>& rF,
                                                             double pressure,
                                                             c_matrix<double,DIM,DIM>& rS)
{
    c_matrix<double,DIM,DIM> C = prod(trans(rF), rF);
    c_matrix<double,DIM,DIM> invC = Inverse(C);

    c_matrix<double,DIM,DIM> T;

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE; // not filled in, made static for efficiency

    ComputeStressAndStressDerivative(C, invC, pressure, T, dTdE, false);

    rS = prod(T, trans(rF));
}

template<unsigned DIM>
void AbstractMaterialLaw<DIM>::Compute2ndPiolaKirchoffStress(c_matrix<double,DIM,DIM>& rC,
                                                             double pressure,
                                                             c_matrix<double,DIM,DIM>& rT)
{
    c_matrix<double,DIM,DIM> invC = Inverse(rC);

    static FourthOrderTensor<DIM,DIM,DIM,DIM> dTdE; // not filled in, made static for efficiency

    ComputeStressAndStressDerivative(rC, invC, pressure, rT, dTdE, false);
}

// LCOV_EXCL_START
template<unsigned DIM>
void AbstractMaterialLaw<DIM>::ScaleMaterialParameters(double scaleFactor)
{
    EXCEPTION("[the material law you are using]::ScaleMaterialParameters() has not been implemented\n");
}
// LCOV_EXCL_STOP

template<unsigned DIM>
void AbstractMaterialLaw<DIM>::SetChangeOfBasisMatrix(c_matrix<double,DIM,DIM>& rChangeOfBasisMatrix)
{
    mpChangeOfBasisMatrix = &rChangeOfBasisMatrix;
}

template<unsigned DIM>
void AbstractMaterialLaw<DIM>::ResetToNoChangeOfBasisMatrix()
{
    mpChangeOfBasisMatrix = nullptr;
}

template<unsigned DIM>
void AbstractMaterialLaw<DIM>::ComputeTransformedDeformationTensor(c_matrix<double,DIM,DIM>& rC, c_matrix<double,DIM,DIM>& rInvC,
                                                                   c_matrix<double,DIM,DIM>& rCTransformed, c_matrix<double,DIM,DIM>& rInvCTransformed)
{
    // Writing the local coordinate system as fibre/sheet/normal, as in cardiac problems..

    // Let P be the change-of-basis matrix P = (\mathbf{m}_f, \mathbf{m}_s, \mathbf{m}_n).
    // The transformed C for the fibre/sheet basis is C* = P^T C P.
    if (mpChangeOfBasisMatrix)
    {
        // C* = P^T C P, and ditto inv(C)
        rCTransformed = prod(trans(*mpChangeOfBasisMatrix),(c_matrix<double,DIM,DIM>)prod(rC,*mpChangeOfBasisMatrix));         // C*    = P^T C    P
        rInvCTransformed = prod(trans(*mpChangeOfBasisMatrix),(c_matrix<double,DIM,DIM>)prod(rInvC,*mpChangeOfBasisMatrix));   // invC* = P^T invC P
    }
    else
    {
        rCTransformed = rC;
        rInvCTransformed = rInvC;
    }
}

template<unsigned DIM>
void AbstractMaterialLaw<DIM>::TransformStressAndStressDerivative(c_matrix<double,DIM,DIM>& rT,
                                                                  FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                                                  bool transformDTdE)
{
    //  T = P T* P^T   and   dTdE_{MNPQ}  =  P_{Mm}P_{Nn}P_{Pp}P_{Qq} dT*dE*_{mnpq}
    if (mpChangeOfBasisMatrix)
    {
        static c_matrix<double,DIM,DIM> T_transformed_times_Ptrans;
        T_transformed_times_Ptrans = prod(rT, trans(*mpChangeOfBasisMatrix));

        rT = prod(*mpChangeOfBasisMatrix, T_transformed_times_Ptrans);  // T = P T* P^T

        // dTdE_{MNPQ}  =  P_{Mm}P_{Nn}P_{Pp}P_{Qq} dT*dE*_{mnpq}
        if (transformDTdE)
        {
            static FourthOrderTensor<DIM,DIM,DIM,DIM> temp;
            temp.template SetAsContractionOnFirstDimension<DIM>(*mpChangeOfBasisMatrix, rDTdE);
            rDTdE.template SetAsContractionOnSecondDimension<DIM>(*mpChangeOfBasisMatrix, temp);
            temp.template SetAsContractionOnThirdDimension<DIM>(*mpChangeOfBasisMatrix, rDTdE);
            rDTdE.template SetAsContractionOnFourthDimension<DIM>(*mpChangeOfBasisMatrix, temp);
        }
    }
}

// Explicit instantiation
template class AbstractMaterialLaw<2>;
template class AbstractMaterialLaw<3>;
