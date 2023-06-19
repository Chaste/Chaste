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

#include "SlidingBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
SlidingBoundaryCondition<DIM>::SlidingBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                    double threshold)
        : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
          mThreshold(threshold)
{
}

template<unsigned DIM>
double SlidingBoundaryCondition<DIM>::GetThreshold() const
{
    return mThreshold;
}

template<unsigned DIM>
void SlidingBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("SlidingBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    ChasteCuboid<DIM> bounds = this->mpCellPopulation->rGetMesh().CalculateBoundingBox();
    double x_min = bounds.rGetLowerCorner()[0];

    // Loop over every node
    for (unsigned node_index=0; node_index<this->mpCellPopulation->GetNumNodes(); node_index++)
    {
        Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

        if (p_node->IsBoundaryNode())
        {
            c_vector<double, DIM> old_node_location;
            old_node_location = rOldLocations.find(p_node)->second;

            // If the node lies on the left, then revert its x coordinate
            if (p_node->rGetLocation()[0] < x_min + mThreshold)
            {
                p_node->rGetModifiableLocation()[0] = old_node_location[0];
            }
        }
    }
}

template<unsigned DIM>
bool SlidingBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    // Without further information, it is not possible to directly verify the boundary condition
    return true;
}

template<unsigned DIM>
void SlidingBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Threshold>" << mThreshold << "</Threshold>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class SlidingBoundaryCondition<1>;
template class SlidingBoundaryCondition<2>;
template class SlidingBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SlidingBoundaryCondition)
