/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "SemRegionalForce.hpp"

#include "SemEnumerations.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SemRegionalForce<ELEMENT_DIM,SPACE_DIM>::SemRegionalForce()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> SemRegionalForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    const double distance =  rCellPopulation.rGetMesh().GetDistanceBetweenNodes(nodeAGlobalIndex, nodeBGlobalIndex);

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    const c_vector<double, SPACE_DIM>& r_node_a_location = p_node_a->rGetLocation();
    const c_vector<double, SPACE_DIM>& r_node_b_location = p_node_b->rGetLocation();

    const SemNodeRegion node_a_region = static_cast<SemNodeRegion>(p_node_a->GetRegion());
    const SemNodeRegion node_b_region = static_cast<SemNodeRegion>(p_node_b->GetRegion());

    if (distance < 0.5)
    {
        const double rest_length = 0.5 * (mRestLengths[node_a_region] + mRestLengths[node_b_region]);
        const double spring_const = 0.5 * (mSpringConstants[node_a_region] + mSpringConstants[node_b_region]);

        const double overlap = distance - rest_length;
        c_vector<double, SPACE_DIM> calculated_force = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location) * (spring_const * overlap / distance);
        
        return calculated_force;
    }
    else
    {
        return zero_vector<double>(SPACE_DIM);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SemRegionalForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // *rParamsFile << "\t\t\t<MeinekeSpringStiffness>" << mMeinekeSpringStiffness << "</MeinekeSpringStiffness>\n";
    // *rParamsFile << "\t\t\t<MeinekeDivisionRestingSpringLength>" << mMeinekeDivisionRestingSpringLength << "</MeinekeDivisionRestingSpringLength>\n";
    // *rParamsFile << "\t\t\t<MeinekeSpringGrowthDuration>" << mMeinekeSpringGrowthDuration << "</MeinekeSpringGrowthDuration>\n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SemRegionalForce<1,1>;
template class SemRegionalForce<1,2>;
template class SemRegionalForce<2,2>;
template class SemRegionalForce<1,3>;
template class SemRegionalForce<2,3>;
template class SemRegionalForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SemRegionalForce)
