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

    /** -1.0/DIM */
    static const double mMinusOneOverDimension = -1.0/DIM;

public:

    /**
     * Get the first derivative dW/dI1.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_dW_dI1(double I1, double I2, double I3)
    {
        return mC1 * pow(I3, mMinusOneOverDimension);
    }

    /**
     * Get the first derivative dW/dI2.
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
     * Get the first derivative dW/dI3.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_dW_dI3(double I1, double I2, double I3)
    {
        return     mC1*I1*mMinusOneOverDimension*pow(I3,mMinusOneOverDimension - 1)
                +  mC3*(1 - pow(I3,-0.5));
    }

    /**
     * Get the second derivative d^2W/dI1^2.
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
     * Get the second derivative d^2W/dI2^2.
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
     * Get the second derivative d^2W/dI3^2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI3(double I1, double I2, double I3)
    {
        return    mC1*I1*mMinusOneOverDimension*(mMinusOneOverDimension - 1)*pow(I3,mMinusOneOverDimension - 2)
                + 0.5*mC3*pow(I3,-1.5);
    }



    /**
     * Get the second derivative d^2W/dI2dI3.
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
     * Get the second derivative d^2W/dI1dI3.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI1I3(double I1, double I2, double I3)
    {
        return mC1*mMinusOneOverDimension*pow(I3,mMinusOneOverDimension-1);
    }


    /**
     * Get the second derivative d^2W/dI1dI2.
     *
     * @param I1 first principal invariant of C
     * @param I2 second principal invariant of C
     * @param I3 third principal invariant of C
     */
    double Get_d2W_dI1I2(double I1, double I2, double I3)
    {
        return 0.0;
    }


    /** Get method for mC1. */
    double GetC1()
    {
        return mC1;
    }

    /** Get method for mC3. */
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


#endif /*COMPRESSIBLEMOONEYRIVLINMATERIALLAW_HPP_*/
