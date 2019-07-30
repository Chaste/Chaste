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
