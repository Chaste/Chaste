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

#ifndef _FAKEBATHCONTRACTIONMODEL_HPP_
#define _FAKEBATHCONTRACTIONMODEL_HPP_

#include "AbstractAlgebraicContractionModel.hpp"


/**
 * A fake contraction model that returns zero active tension always.
 * It is useful to be placed in bath nodes within EM simulations
 */
class FakeBathContractionModel : public AbstractAlgebraicContractionModel
{

public:
    /**
     * Constructor.
     */
    FakeBathContractionModel();

    /**
     *  Set the input parameters (none of which are used in this model)
     *  @param rInputParameters Structure containing voltage, [Ca]_i.
     */
    void SetInputParameters(ContractionModelInputParameters& rInputParameters);

    /**
     *  Set the fibre stretch and stretch-rate (none of these are used)
     *  @param stretch stretch in the fibre direction
     *  @param stretchRate rate of change of stretch in the fibre direction (not used in this model).
     */
    void SetStretchAndStretchRate(double stretch, double stretchRate);

    /**
     *  @return the active tension given the current stretch and time
     *  It will return zero always (this is a fake contraction modoel)
     */
    double GetActiveTension();

    /**
     *  @return model is stretch-rate-independent
     */
    bool IsStretchDependent()
    {
        return false;
    }

    /**
     *  @return model is stretch-rate-independent
     */
    bool IsStretchRateDependent()
    {
        return false;
    }
};


#endif /*_FAKEBATHCONTRACTIONMODEL_HPP_*/

