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


#ifndef EXPONENTIALMATERIALLAW_HPP_
#define EXPONENTIALMATERIALLAW_HPP_

#include "AbstractIsotropicIncompressibleMaterialLaw.hpp"
#include "Exception.hpp"

/**
 *  ExponentialMaterialLaw
 *
 *  An exponential isotropic incompressible hyperelastic material law for finite
 *  elasticity
 *
 *  The law is given by a strain energy function
 *      W(I_1,I_2,I_3) = a exp( b(I_1-3) ) - p/2 C^{-1}
 *  in 3d, or
 *      W(I_1,I_2,I_3) = a exp( b(I_1-2) ) - p/2 C^{-1}
 *  in 2d.
 *
 *  Here I_i are the principal invariants of C, the Lagrangian deformation tensor.
 *  (I1=trace(C), I2=trace(C)^2-trace(C^2), I3=det(C)).

 *  Note: only dimension equals 2 or 3 is permitted.
 */

template<unsigned DIM>
class ExponentialMaterialLaw : public AbstractIsotropicIncompressibleMaterialLaw<DIM>
{
private:

    /** Parameter a. */
    double mA;

    /** Parameter b. */
    double mB;

public:

    /**
     * @return the first derivative dW/dI1.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI1(double I1, double I2);

    /**
     * @return the first derivative dW/dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_dW_dI2(double I1, double I2);

    /**
     * @return the second derivative d^2W/dI1^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1(double I1, double I2);

    /**
     * @return the second derivative d^2W/dI2^2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI2(double I1, double I2);

    /**
     * @return the second derivative d^2W/dI1dI2.
     *
     * \todo The name of this method should not include underscores.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     */
    double Get_d2W_dI1I2(double I1, double I2);

    /** @return  mA. */
    double GetA();

    /** @return  mB. */
    double GetB();

public:

    /**
     * Constructor, taking in the parameters a and b. a must be positive.
     *
     * @param a the parameter a
     * @param b the parameter b
     */
    ExponentialMaterialLaw(double a, double b);
};

#endif /*EXPONENTIALMATERIALLAW_HPP_*/
