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


#ifndef ABSTRACTISOTROPICINCOMPRESSIBLEMATERIALLAW_HPP_
#define ABSTRACTISOTROPICINCOMPRESSIBLEMATERIALLAW_HPP_

#include "AbstractIncompressibleMaterialLaw.hpp"

/**
 *  AbstractIsotropicIncompressibleMaterialLaw
 *
 *  An isotropic incompressible hyperelastic material law for finite elastiticy
 *
 *  The law is given by a strain energy function W(I1,I2,I3), where I_i are the principal
 *  invariants of C, the Lagrangian deformation tensor. (I1=trace(C), I2=trace(C)^2-trace(C^2),
 *  I3=det(C)). Since it is incompressible, the full strain energy has the form
 *  W^{full} = W(I_1,I_2) - p/2 C^{-1}
 *
 *  Note: only dimension equals 2 or 3 should be permitted.
 */
template<unsigned DIM>
class AbstractIsotropicIncompressibleMaterialLaw : public AbstractIncompressibleMaterialLaw<DIM>
{
protected:

    /**
     * @return the first derivative dW/dI1.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_dW_dI1(double I1, double I2)=0;

    /**
     * @return the first derivative dW/dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_dW_dI2(double I1, double I2)=0;

    /**
     * @return the second derivative d^2W/dI1^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_d2W_dI1(double I1, double I2)=0;

    /**
     * @return the second derivative d^2W/dI2^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_d2W_dI2(double I1, double I2)=0;

    /**
     * @return the second derivative d^2W/dI1dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_d2W_dI1I2(double I1, double I2)=0;

public:

    /**
     *  Compute the (2nd Piola Kirchoff) stress T and the stress derivative dT/dE
     *  for a given strain.
     *
     *  NOTE: the strain E is not expected to be passed in, instead the Lagrangian
     *  deformation tensor C is required (recall, E = 0.5(C-I)
     *
     *  dT/dE is a fourth-order tensor, where dT/dE(M,N,P,Q) = dT^{MN}/dE_{PQ}
     *
     *  @param rC The Lagrangian deformation tensor (F^T F)
     *  @param rInvC The inverse of C. Should be computed by the user. (Change this?)
     *  @param pressure the current pressure
     *  @param rT the stress will be returned in this parameter
     *  @param rDTdE the stress derivative will be returned in this parameter, assuming
     *    the final parameter is true
     *  @param computeDTdE a boolean flag saying whether the stress derivative is
     *    required or not.
     *
     *  This is the implemtation for an isotropic material law, so the stress etc is
     *  computed by calling methods returning dW/dI1, dW/dI2 etc.
     */
    void ComputeStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                          c_matrix<double,DIM,DIM>& rInvC,
                                          double                    pressure,
                                          c_matrix<double,DIM,DIM>& rT,
                                          FourthOrderTensor<DIM,DIM,DIM,DIM>&   rDTdE,
                                          bool                      computeDTdE);

    /**
     * Destructor.
     */
    virtual ~AbstractIsotropicIncompressibleMaterialLaw();

    /**
     *  @return the pressure corresponding to E=0, ie corresponding to C=identity
     *
     *  Since T = 2*Get_dW_dI1 identity + 4*Get_dW_dI2 (I1*identity - C) - p inverse(C),
     *  this is equal to 2*Get_dW_dI1(3,3) + 4*Get_dW_dI2(3,3) in 3D
     */
    double GetZeroStrainPressure();
};

#endif /*ABSTRACTISOTROPICINCOMPRESSIBLEMATERIALLAW_HPP_*/
