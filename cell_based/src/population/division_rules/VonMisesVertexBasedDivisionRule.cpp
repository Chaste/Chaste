/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "VonMisesVertexBasedDivisionRule.hpp"
#include "MathsCustomFunctions.hpp"

template<unsigned DIM>
VonMisesVertexBasedDivisionRule<DIM>::VonMisesVertexBasedDivisionRule()
   : AbstractVertexBasedDivisionRule<DIM>(),
     mMeanParameter(0.0),
     mConcentrationParameter(1.0)
{
}

template<unsigned DIM>
VonMisesVertexBasedDivisionRule<DIM>::~VonMisesVertexBasedDivisionRule()
{
}

template<unsigned DIM>
double VonMisesVertexBasedDivisionRule<DIM>::GetMeanParameter()
{
    return mMeanParameter;
}

template<unsigned DIM>
double VonMisesVertexBasedDivisionRule<DIM>::GetConcentrationParameter()
{
    return mConcentrationParameter;
}

template<unsigned DIM>
void VonMisesVertexBasedDivisionRule<DIM>::SetMeanParameter(double meanParameter)
{
    mMeanParameter = meanParameter;
}

template<unsigned DIM>
void VonMisesVertexBasedDivisionRule<DIM>::SetConcentrationParameter(double concentrationParameter)
{
    assert(concentrationParameter > 0);
    mConcentrationParameter = concentrationParameter;
}

template <unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> VonMisesVertexBasedDivisionRule<SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation)
{
    // Algorithm taken from: Fisher, Statistical Analysis of Circular Data, CUP (1993)

    double a = 1.0 + sqrt(1.0 + 4.0*mConcentrationParameter*mConcentrationParameter);
    double b = (a - sqrt(2.0*a)) / (2.0*mConcentrationParameter);
    double r = (1.0 + b*b) / (2.0*b);

    unsigned counter = 0;
    double theta = 0;
    while (counter <= 0)
    {
        double u1 = RandomNumberGenerator::Instance()->ranf();
        double u2 = RandomNumberGenerator::Instance()->ranf();
        double u3 = RandomNumberGenerator::Instance()->ranf();

        double z = cos(M_PI*u1);
        double f = (1.0 + r*z) / (r + z);
        double c = mConcentrationParameter*(r - f);

        if (((c*(2.0 - c) - u2) > 0.0) || ((log(c/u2) + 1.0 - c) > 0.0))
        {
            theta = fmod((Signum(u3 - 0.5)*std::acos(f)) + mMeanParameter, 2*M_PI);
            counter += 1;
        }
    }

    c_vector<double, SPACE_DIM> vector;
    vector(0) = cos(theta);
    vector(1) = sin(theta);

    return vector;
}

// Explicit instantiation
template class VonMisesVertexBasedDivisionRule<1>;
template class VonMisesVertexBasedDivisionRule<2>;
template class VonMisesVertexBasedDivisionRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VonMisesVertexBasedDivisionRule)
