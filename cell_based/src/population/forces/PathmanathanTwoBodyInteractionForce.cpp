/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "PathmanathanTwoBodyInteractionForce.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::PathmanathanTwoBodyInteractionForce()
   : AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mAlpha(5.0)
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::~PathmanathanTwoBodyInteractionForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::CalculateLinkInteraction(double overlap,
                                                                                                                   double restLength,
                                                                                                                   const c_vector<double, SPACE_DIM>& rUnitDifference,
                                                                                                                   double multiplicationFactor)
{
    // Pathmanathan force: logarithmic repulsion / exponential attraction
    if (overlap <= 0)
    {
        assert(overlap > -restLength);
        return multiplicationFactor * this->mSpringStiffness * rUnitDifference * restLength * log(1.0 + overlap/restLength);
    }
    else
    {
        return multiplicationFactor * this->mSpringStiffness * rUnitDifference * overlap * exp(-mAlpha * overlap/restLength);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetAlpha()
{
    return mAlpha;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::SetAlpha(double alpha)
{
    assert(alpha > 0.0);
    mAlpha = alpha;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringStiffness>" << this->mSpringStiffness << "</SpringStiffness>\n";
    *rParamsFile << "\t\t\t<DivisionRestingSpringLength>" << this->mDivisionRestingSpringLength << "</DivisionRestingSpringLength>\n";
    *rParamsFile << "\t\t\t<SpringGrowthDuration>" << this->mSpringGrowthDuration << "</SpringGrowthDuration>\n";
    *rParamsFile << "\t\t\t<Alpha>" << mAlpha << "</Alpha>\n";

    // Call method on direct parent class
    AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PathmanathanTwoBodyInteractionForce<1,1>;
template class PathmanathanTwoBodyInteractionForce<1,2>;
template class PathmanathanTwoBodyInteractionForce<2,2>;
template class PathmanathanTwoBodyInteractionForce<1,3>;
template class PathmanathanTwoBodyInteractionForce<2,3>;
template class PathmanathanTwoBodyInteractionForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PathmanathanTwoBodyInteractionForce)
