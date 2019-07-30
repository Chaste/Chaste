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


#ifndef ABSTRACTMATERIALLAW_HPP_
#define ABSTRACTMATERIALLAW_HPP_

#include "UblasCustomFunctions.hpp"
#include <cassert>
#include <vector>
#include "Exception.hpp"
#include "FourthOrderTensor.hpp"

/**
 *  AbstractMaterialLaw
 *
 *  A hyper-elastic material law for finite elasticity
 *
 *  The law is given by a strain energy function W(E), where E is the strain, such
 *  that the (2nd Piola-Kirchhoff) stress T = dW/dE
 */
template<unsigned DIM>
class AbstractMaterialLaw
{
protected:
    /**
     *  Some material laws are based on a particular local set of preferred directions, eg anisotropic
     *  cardiac laws, which use the fibre sheet and normal directions. This matrix can defines the
     *  orientation and should set before T and dTdE are computed.
     *
     *  The change of matrix should have the form P = [a_f a_s a_n], where each a_i is a vector.
     */
    c_matrix<double,DIM,DIM>* mpChangeOfBasisMatrix;

    /**
     *  Transform the input (C and inv(C), where C in the deformation tensor) to the local coordinate system
     *  @param rC deformation tensor C (input)
     *  @param rInvC inverse of C (input)
     *  @param rCTransformed P^T C P (output)
     *  @param rInvCTransformed P^T inv(C) P (output)
     */
    void ComputeTransformedDeformationTensor(c_matrix<double,DIM,DIM>& rC, c_matrix<double,DIM,DIM>& rInvC,
                                             c_matrix<double,DIM,DIM>& rCTransformed, c_matrix<double,DIM,DIM>& rInvCTransformed);

    /**
     *  Transform the output (T and dTdE) back to the original coordinate system
     *  @param rT stress being computed
     *  @param rDTdE the stress derivative to be transformed (assuming
     *    the final parameter is true)
     *  @param transformDTdE a boolean flag saying whether the stress derivative is
     *    to be transformed or not
     */
    void TransformStressAndStressDerivative(c_matrix<double,DIM,DIM>& rT,
                                            FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                            bool transformDTdE);

public:

    /** Constuctor */
    AbstractMaterialLaw();

    /** Destructor */
    virtual ~AbstractMaterialLaw()
    {
    }

    /**
     *  Compute the (2nd Piola Kirchoff) stress T and the stress derivative dT/dE for
     *  a given strain.
     *
     *  NOTE: the strain E is not expected to be passed in, instead the Lagrangian
     *  deformation tensor C is required (recall, E = 0.5(C-I))
     *
     *  dTdE is a fourth-order tensor, where dTdE(M,N,P,Q) = dT^{MN}/dE_{PQ}
     *
     *  @param rC The Lagrangian deformation tensor (F^T F)
     *  @param rInvC The inverse of C. Should be computed by the user.
     *  @param pressure the current pressure
     *  @param rT the stress will be returned in this parameter
     *  @param rDTdE the stress derivative will be returned in this parameter, assuming
     *    the final parameter is true
     *  @param computeDTdE a boolean flag saying whether the stress derivative is
     *    required or not.
     */
    virtual void ComputeStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                                  c_matrix<double,DIM,DIM>& rInvC,
                                                  double                    pressure,
                                                  c_matrix<double,DIM,DIM>& rT,
                                                  FourthOrderTensor<DIM,DIM,DIM,DIM>&   rDTdE,
                                                  bool                      computeDTdE)=0;

    /**
     *  Compute the Cauchy stress (the true stress), given the deformation gradient
     *  F and the pressure. The Cauchy stress is given by
     *
     *  sigma^{ij} = (1/detF) F^i_M T^{MN} F^j_N
     *
     *  where T is the 2nd Piola Kirchoff stress, dW/dE
     *
     *  @param rF the deformation gradient
     *  @param pressure the pressure
     *  @param rSigma an empty matrix, which will be filled in with the Cauchy stress
     *
     *  Note: the compute the material part of the stress (the pressure-independent
     *  part), just pass in pressure=0.0
     */
    void ComputeCauchyStress(c_matrix<double,DIM,DIM>& rF, double pressure, c_matrix<double,DIM,DIM>& rSigma);

    /**
     *  Compute the 1st Piola Kirchoff stress, given the deformation gradient F
     *  and the pressure. The 1st Piola Kirchoff stress given by
     *
     *  S^{Mi} = T^{MN} F^i_M,
     *
     *  where T is the 2nd PK stress, dW/dE.
     *
     *  Note that this stress is not symmetric and the least useful of the three
     *  stresses.
     *
     *  @param rF the deformation gradient
     *  @param pressure the pressure
     *  @param rS an empty matrix, which will be filled in with the stress
     *
     *  Note: the compute the material part of the stress (the pressure-independent
     *  part), just pass in pressure=0.0
     */
    void Compute1stPiolaKirchoffStress(c_matrix<double,DIM,DIM>& rF, double pressure, c_matrix<double,DIM,DIM>& rS);

    /**
     *  Compute the 2nd Piola Kirchoff stress, given the deformation tensor C
     *  and the pressure. The 2nd Piola Kirchoff stress given by
     *
     *  T^{MN} = dW/dE_{MN} = 2dW/dC_{MN}
     *
     *  @param rC the Lagrange deformation tensor (C=F^T F), *not* F, and *not* E
     *  @param pressure the pressure
     *  @param rT an empty matrix, which will be filled in with the stress
     *
     *  Note: to compute the material part of the stress (the pressure-independent
     *  part), just pass in pressure=0.0
     */
    void Compute2ndPiolaKirchoffStress(c_matrix<double,DIM,DIM>& rC, double pressure, c_matrix<double,DIM,DIM>& rT);


    /**
     *  Set a scale factor by which (dimensional) material parameters are scaled. This method
     *  can be optionally implemented in the child class; if no implementation is made an
     *  exception is thrown. A scale factor may be used/needed to improve GMRES convergence.
     *  Note that is a material law is scaled like this any dimensionally equivalent terms
     *  (eg gravity, tractions, active tensions) must also be scaled. Also, computed pressure
     *  will come out scaled.
     *
     *  @param scaleFactor  the scale factor
     */
    virtual void ScaleMaterialParameters(double scaleFactor);

    /**
     *  Some material laws (eg pole-zero) may have preferred directions (eg fibre direction),
     *  but be implemented to assume the preferred directions are parallel to the X-axis etc.
     *  Call this with the change of basis matrix and C will be transformed from the Lagrangian
     *  coordinate system to the appropriate coordinate system before used to calculate T, which
     *  will then be transformed from the appropriate coordinate system back to the Lagrangian
     *  coordinate system before being returned, as will dTdE
     *
     *  Note that no copy of this matrix is taken, so the original matrix must persist whilst
     *  this class is used. Call ResetToNoChangeOfBasisMatrix() if necessary.
     *
     *  The change of matrix should have the form (writing the preferred directions as fibre,
     *  sheet and normal, as in heart simulations): P = [a_f a_s a_n], where each a_i is a vector.
     *
     *  @param rChangeOfBasisMatrix Change of basis matrix.
     */
    void SetChangeOfBasisMatrix(c_matrix<double,DIM,DIM>& rChangeOfBasisMatrix);

    /**
     *  Reset back to no change of basis matrix
     */
    void ResetToNoChangeOfBasisMatrix();
};

#endif /*ABSTRACTMATERIALLAW_HPP_*/
