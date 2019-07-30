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


#include "NonPhysiologicalContractionModel.hpp"
#include <math.h>

#include <iostream>

NonPhysiologicalContractionModel::NonPhysiologicalContractionModel(unsigned option)
    : AbstractAlgebraicContractionModel()
{
    assert(option>=1);
    assert(option<=3);

    mOption = option;
    mStretch = 1.0;
}


void NonPhysiologicalContractionModel::SetInputParameters(ContractionModelInputParameters& rInputParameters)
{
}

void NonPhysiologicalContractionModel::SetStretchAndStretchRate(double stretch, double stretchRate)
{
    mStretch = stretch;
}

double NonPhysiologicalContractionModel::GetActiveTension()
{
    if (mOption==1)
    {
        // If solving mechanics problem, results using implicit and explicit solvers are identical, as they should be
        return fabs(sin(mTime));
    }
    else if (mOption==2)
    {
        // small error between implicit and explicit at lowest dt
        // next dt, small difference between explicit at lowest dt, mostly due to first timestep
        // largest dt (1ms) completely wrong after first timestep => solution translated across
        return fabs(5*mStretch*sin(mTime));
    }
    else
    {
        // same conclusions as for above
        return fabs(5*exp(1-mStretch)*sin(mTime));
    }
}

