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

#ifndef TOYCOMPRESSIBLEMATERIALLAW_HPP_
#define TOYCOMPRESSIBLEMATERIALLAW_HPP_

#include "AbstractIsotropicCompressibleMaterialLaw.hpp"

/*
 * Simple material law W(I1,I2,I3) = c1(I1-3) + c2(I2-3) + c3(I3-1),
 * which may not correspond to a physically acceptable law but can
 * still be used to test the code.
 */
template<unsigned DIM>
class ToyCompressibleMaterialLaw : public AbstractIsotropicCompressibleMaterialLaw<DIM>
{
private:

    double mC1;
    double mC2;
    double mC3;

public:

    double Get_dW_dI1(double I1, double I2, double I3)
    {
        return mC1;
    }

    double Get_dW_dI2(double I1, double I2, double I3)
    {
        return mC2;
    }

    double Get_dW_dI3(double I1, double I2, double I3)
    {
        return mC3;
    }

    double Get_d2W_dI1(double I1, double I2, double I3)
    {
        return 0.0;
    }

    double Get_d2W_dI2(double I1, double I2, double I3)
    {
        return 0.0;
    }

    double Get_d2W_dI3(double I1, double I2, double I3)
    {
        return 0.0;
    }

    double Get_d2W_dI2I3(double I1, double I2, double I3)
    {
        return 0.0;
    }

    double Get_d2W_dI1I3(double I1, double I2, double I3)
    {
        return 0.0;
    }

    double Get_d2W_dI1I2(double I1, double I2, double I3)
    {
        return 0.0;
    }

    ToyCompressibleMaterialLaw(double c1, double c2, double c3)
    {
        assert(c1 > 0.0);
        assert(DIM!=2 || c2==0.0);
        mC1 = c1;
        mC2 = c2;
        mC3 = c3;
        if (DIM==3 && fabs(c1+2*c2+c3)>1e-8)
        {
            EXCEPTION("c1+2*c2+c3 should be equal to zero");
        }
    }
};

#endif /* TOYCOMPRESSIBLEMATERIALLAW_HPP_ */
