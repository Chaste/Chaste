/*

Copyright (c) 2005-2017, University of Oxford.
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
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned DIM>
HorizontalStretchForce<DIM>::HorizontalStretchForce(const double ForceMagnitude, const double RelativeWidth)
        : AbstractForce<DIM>(),
          mForceMagnitude(ForceMagnitude),
          mRelativeWidth(RelativeWidth)
{
    SetUpForceVectorUsingMagnitude();
}

template <unsigned DIM>
HorizontalStretchForce<DIM>::~HorizontalStretchForce()
{
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::SetUpForceVectorUsingMagnitude()
{
    mForceVector = zero_vector<double>(DIM);
    mForceVector[0] = mForceMagnitude;
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::SetForceMagnitude(double forceMagnitude)
{
    mForceMagnitude = forceMagnitude;
    SetUpForceVectorUsingMagnitude();
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::SetRelativeWidth(double relativeWidth)
{
    mRelativeWidth = relativeWidth;
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    std::set<Node<DIM>*> moving_nodes;
    for (typename std::set<VertexElement<DIM, DIM>*>::iterator it = mMovingPinnedElements.begin(); it != mMovingPinnedElements.end(); ++it)
    {
        VertexElement<DIM, DIM>* p_elem = *it;
        for (unsigned i = 0; i < p_elem->GetNumNodes(); ++i)
        {
            moving_nodes.insert(p_elem->GetNode(i));
        }
    }

    for (typename std::set<Node<DIM>*>::iterator p_p_node = moving_nodes.begin(); p_p_node != moving_nodes.end(); ++p_p_node)
    {
        (*p_p_node)->ClearAppliedForce();
        (*p_p_node)->AddAppliedForceContribution(mForceVector);
    }

    std::set<Node<DIM>*> non_moving_nodes;
    for (typename std::set<VertexElement<DIM, DIM>*>::iterator it = mNonMovingPinnedElements.begin(); it != mNonMovingPinnedElements.end(); ++it)
    {
        VertexElement<DIM, DIM>* p_elem = *it;
        for (unsigned i = 0; i < p_elem->GetNumNodes(); ++i)
        {
            non_moving_nodes.insert(p_elem->GetNode(i));
        }
    }

    for (typename std::set<Node<DIM>*>::iterator p_p_node = non_moving_nodes.begin(); p_p_node != non_moving_nodes.end(); ++p_p_node)
    {
        (*p_p_node)->ClearAppliedForce();
    }
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ForceMagnitude>" << mForceMagnitude << "</ForceMagnitude>\n";
    *rParamsFile << "\t\t\t<RelativeWidth>" << mRelativeWidth << "</RelativeWidth>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::SetUpPinnedElements(AbstractCellPopulation<DIM>& rCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("HorizontalStretchForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    const MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();

    ChasteCuboid<DIM> bounding_box_tmp = rMesh.CalculateBoundingBox();
    const double min_x = bounding_box_tmp.rGetLowerCorner()[0];
    const double max_x = bounding_box_tmp.rGetUpperCorner()[0];
    assert(max_x > min_x);
    const double width_x = max_x - min_x;
    const double width = width_x * mRelativeWidth;
    const double left_bound = min_x + width;
    const double right_bound = max_x - width;

    // Iterate over nodes in the cell population
    for (unsigned elem_index = 0; elem_index < rMesh.GetNumElements(); ++elem_index)
    {
        VertexElement<DIM, DIM>* p_elem = rMesh.GetElement(elem_index);
        const double loc_x = p_elem->GetCentroid()[0];
        bool is_elem_pinned = false;
        if (loc_x < left_bound)
        {
            mNonMovingPinnedElements.insert(p_elem);
            is_elem_pinned = true;
        }
        else if (loc_x > right_bound)
        {
            mMovingPinnedElements.insert(p_elem);
            is_elem_pinned = true;
        }

        if (is_elem_pinned)
        {
            CellPtr p_cell_tmp = rCellPopulation.GetCellUsingLocationIndex(p_elem->GetIndex());
            NoCellCycleModel* p_no_cell_cycle_model = new NoCellCycleModel;
            p_cell_tmp->SetCellCycleModel(p_no_cell_cycle_model);
            p_cell_tmp->InitialiseCellCycleModel();
        }
    }
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::ClearPinnedElements()
{
    mNonMovingPinnedElements.clear();
    mMovingPinnedElements.clear();
}

template <unsigned DIM>
std::set<VertexElement<DIM, DIM>*>& HorizontalStretchForce<DIM>::GetPinnedElements(bool isMovingElements)
{
    return isMovingElements ? mMovingPinnedElements : mNonMovingPinnedElements;
}

template <unsigned DIM>
void HorizontalStretchForce<DIM>::SetForceVector(const c_vector<double, DIM>& ForceVector)
{
    mForceVector = ForceVector;
}

template <unsigned DIM>
c_vector<double, DIM>& HorizontalStretchForce<DIM>::rGetForceVector()
{
    return mForceVector;
}

template class HorizontalStretchForce<1>;
template class HorizontalStretchForce<2>;
template class HorizontalStretchForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HorizontalStretchForce)
