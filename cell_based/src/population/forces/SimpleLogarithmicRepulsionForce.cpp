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

#include "SimpleLogarithmicRepulsionForce.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SimpleLogarithmicRepulsionForce<ELEMENT_DIM, SPACE_DIM>::SimpleLogarithmicRepulsionForce()
    : AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>(),
      mRepulsionParameter(15.0)
{
    if constexpr (SPACE_DIM == 1)
    {
        mRepulsionParameter = 30.0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double SimpleLogarithmicRepulsionForce<ELEMENT_DIM, SPACE_DIM>::GetRepulsionParameter()
{
    return mRepulsionParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SimpleLogarithmicRepulsionForce<ELEMENT_DIM, SPACE_DIM>::SetRepulsionParameter(double repulsionParameter)
{
    assert(repulsionParameter > 0.0);
    mRepulsionParameter = repulsionParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> SimpleLogarithmicRepulsionForce<ELEMENT_DIM, SPACE_DIM>::CalculateForceBetweenNodes(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
{
    // Throw an exception if not using a NodeBasedCellPopulation
    if (dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SimpleLogarithmicRepulsionForce is to be used with a NodeBasedCellPopulation only");
    }

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    c_vector<double, SPACE_DIM> unit_difference =
        rCellPopulation.rGetMesh().GetVectorFromAtoB(p_node_a->rGetLocation(), p_node_b->rGetLocation());

    double distance = norm_2(unit_difference);
    double rest_length = p_node_a->GetRadius() + p_node_b->GetRadius();
    double overlap = distance - rest_length;

    if (overlap >= 0.0)
    {
        return zero_vector<double>(SPACE_DIM);
    }

    unit_difference /= distance;

    // Logarithmic repulsion (overlap is negative, cells are compressed)
    // log(x+1) is undefined for x<=-1; overlap/rest_length > -1 since distance > 0
    assert(overlap > -rest_length);
    return mRepulsionParameter * unit_difference * rest_length * log(1.0 + overlap/rest_length);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SimpleLogarithmicRepulsionForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RepulsionParameter>" << mRepulsionParameter << "</RepulsionParameter>\n";

    // Call direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SimpleLogarithmicRepulsionForce<1,1>;
template class SimpleLogarithmicRepulsionForce<1,2>;
template class SimpleLogarithmicRepulsionForce<2,2>;
template class SimpleLogarithmicRepulsionForce<1,3>;
template class SimpleLogarithmicRepulsionForce<2,3>;
template class SimpleLogarithmicRepulsionForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SimpleLogarithmicRepulsionForce)
