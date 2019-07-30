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


#ifndef COMPRESSIBLEMOONEYRIVLINMATERIALLAW_HPP_
#define COMPRESSIBLEMOONEYRIVLINMATERIALLAW_HPP_

#include "AbstractIsotropicCompressibleMaterialLaw.hpp"
#include "Exception.hpp"

/**
 *  CompressibleMooneyRivlinMaterialLaw
 *
 *  A Mooney-Rivlin isotropic compressible hyperelastic material law for finite
 *  elasticity
 *
 *  The law is given by a strain energy function
 *      W(I1,I2,I3) = c1 ( dev(I1)-3 )  +  c3(J-1)^2
 *
 *  where (assuming Ii are the principal invariants of C, the Lagrangian deformation tensor,
 *  I1=trace(C), I2=0.5(trace(C)^2-trace(C^2)), I3=det(C)):
 *      J = det(F) = sqrt(I3)
 *      dev(I1) = I1 * J^(-2/DIM)  is the first invariant of the deviatoric part of C
 *
 *  Note T(E=0) = 0 regardless of choice of c1, c3.
 *
 *  NOTE: this is really just a NEO-HOOKEAN law at present - the c2 (dev(I2)-3) term hasn't been
 *  added yet...
 *
 */
template<unsigned DIM>
class CompressibleMooneyRivlinMaterialLaw : public AbstractIsotropicCompressibleMaterialLaw<DIM>
{
private:

    /** Parameter c1. */
    double mC1;

    /** Parameter c3 */
    double mC3;

    /** Initialised to -1.0/DIM*/
    static const double msMinusOneOverDimension;

public:

    /**
     * @return the first derivative dW/dI1.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_dW_dI1(double I1, double I2, double I3)
    {
        return mC1 * pow(I3, msMinusOneOverDimension);
    }

    /**
     * @return the first derivative dW/dI2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_dW_dI2(double I1, double I2, double I3)
    {
        return 0.0;
    }

    /**
     * @return the first derivative dW/dI3.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_dW_dI3(double I1, double I2, double I3)
    {
        return     mC1*I1*msMinusOneOverDimension*pow(I3,msMinusOneOverDimension - 1)
                +  mC3*(1 - pow(I3,-0.5));
    }

    /**
     * @return the second derivative d^2W/dI1^2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI1(double I1, double I2, double I3)
    {
        return 0.0;
    }


    /**
     * @return the second derivative d^2W/dI2^2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI2(double I1, double I2, double I3)
    {
        return 0.0;
    }


    /**
     * @return the second derivative d^2W/dI3^2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI3(double I1, double I2, double I3)
    {
        return    mC1*I1*msMinusOneOverDimension*(msMinusOneOverDimension - 1)*pow(I3,msMinusOneOverDimension - 2)
                + 0.5*mC3*pow(I3,-1.5);
    }

    /**
     * @return the second derivative d^2W/dI2dI3.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI2I3(double I1, double I2, double I3)
    {
        return 0.0;
    }


    /**
     * @return the second derivative d^2W/dI1dI3.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI1I3(double I1, double I2, double I3)
    {
        return mC1*msMinusOneOverDimension*pow(I3,msMinusOneOverDimension-1);
    }


    /**
     * @return the second derivative d^2W/dI1dI2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI1I2(double I1, double I2, double I3)
    {
        return 0.0;
    }


    /** @return  mC1. */
    double GetC1()
    {
        return mC1;
    }

    /** @return  mC3. */
    double GetC3()
    {
        return mC3;
    }

    /**
     * Constructor, taking in parameters c1 and c3.
     *
     * @param c1 parameter c1
     * @param c3 parameter c3
     */
    CompressibleMooneyRivlinMaterialLaw(double c1, double c3)
    {
        assert(c1 > 0.0);
        assert(c3 > 0.0);
        mC1 = c1;
        mC3 = c3;
    }

    /**
     * Scale the dimensional material parameters.
     *
     * @param scaleFactor
     */
    void ScaleMaterialParameters(double scaleFactor)
    {
        assert(scaleFactor > 0.0);
        mC1 /= scaleFactor;
        mC3 /= scaleFactor;
    }
};

/** Initialised to -1.0/DIM  */
template<unsigned DIM>
const double CompressibleMooneyRivlinMaterialLaw<DIM>::msMinusOneOverDimension = -1.0/DIM;

#endif /*COMPRESSIBLEMOONEYRIVLINMATERIALLAW_HPP_*/
