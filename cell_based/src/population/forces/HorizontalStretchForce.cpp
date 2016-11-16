/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "HorizontalStretchForce.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
HorizontalStretchForce<DIM>::HorizontalStretchForce()
    : AbstractForce<DIM>(),
      mForceMagnitude(1.0),
      mRelativeWidth(0.1)
{
}

template<unsigned DIM>
HorizontalStretchForce<DIM>::~HorizontalStretchForce()
{
}

template<unsigned DIM>
void HorizontalStretchForce<DIM>::SetForceMagnitude(double forceMagnitude)
{
    mForceMagnitude = forceMagnitude;
}

template<unsigned DIM>
void HorizontalStretchForce<DIM>::SetRelativeWidth(double relativeWidth)
{
    mRelativeWidth = relativeWidth;
}

template<unsigned DIM>
void HorizontalStretchForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("HorizontalStretchForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();

    const unsigned num_nodes = p_cell_population->GetNumNodes();
    c_vector<double, DIM> force_in_positive_y = zero_vector<double>(DIM);
    force_in_positive_y[0] = mForceMagnitude;
    c_vector<double, DIM> force_in_negative_y = -force_in_positive_y;

    ChasteCuboid<DIM> bounding_box_tmp = rMesh.CalculateBoundingBox();
    const double min_x = bounding_box_tmp.rGetLowerCorner()[0];
    const double max_x = bounding_box_tmp.rGetUpperCorner()[0];
    assert(max_x > min_x);
    const double width_x = max_x - min_x;
    const double width = width_x*mRelativeWidth;
    const double left_bound = min_x + width;
    const double right_bound = max_x - width;
    // could have built two more Cuboids to check if the point is inside,
    // but will have a lot of wrapping which reduce efficiency.
    // Iterate over nodes in the cell population
    for (unsigned node_global_index=0; node_global_index<num_nodes; ++node_global_index)
    {
        Node<3>* p_this_node = p_cell_population->GetNode(node_global_index);
//            if ( ! p_this_node->IsBoundaryNode())
//                continue;
        const double loc_x = p_this_node->rGetLocation()[0];
        if (loc_x < left_bound)
        {
            p_this_node->ClearAppliedForce();
            p_this_node->AddAppliedForceContribution(force_in_negative_y);
        }
        else if (loc_x > right_bound)
        {
            p_this_node->ClearAppliedForce();
            p_this_node->AddAppliedForceContribution(force_in_positive_y);
        }
    }
}

template<unsigned DIM>
void HorizontalStretchForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ForceMagnitude>" << mForceMagnitude << "</ForceMagnitude>\n";
    *rParamsFile << "\t\t\t<RelativeWidth>" << mRelativeWidth << "</RelativeWidth>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HorizontalStretchForce)
