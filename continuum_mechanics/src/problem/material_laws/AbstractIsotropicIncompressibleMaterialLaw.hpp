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
     * Get the first derivative dW/dI1.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_dW_dI1(double I1, double I2)=0;

    /**
     * Get the first derivative dW/dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_dW_dI2(double I1, double I2)=0;

    /**
     * Get the second derivative d^2W/dI1^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_d2W_dI1(double I1, double I2)=0;

    /**
     * Get the second derivative d^2W/dI2^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    virtual double Get_d2W_dI2(double I1, double I2)=0;

    /**
     * Get the second derivative d^2W/dI1dI2.
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
     *  Get the pressure corresponding to E=0, ie corresponding to C=identity
     *
     *  Since T = 2*Get_dW_dI1 identity + 4*Get_dW_dI2 (I1*identity - C) - p inverse(C),
     *  this is equal to 2*Get_dW_dI1(3,3) + 4*Get_dW_dI2(3,3) in 3D
     */
    double GetZeroStrainPressure();
};

#endif /*ABSTRACTISOTROPICINCOMPRESSIBLEMATERIALLAW_HPP_*/
