/*

Copyright (c) 2005-2024, University of Oxford.
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

#include "AbstractCentreBasedCellPopulation.hpp"
#include "AttractingPlaneBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "UblasIncludes.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::AttractingPlaneBoundaryCondition(
        AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation,
                                                    c_vector<double, SPACE_DIM> point,
                                                    c_vector<double, SPACE_DIM> normal)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation),
          mPointOnPlane(normal),
          mUseJiggledNodesOnPlane(false),
          mAttractionThreshold(0.1)
{
    assert(norm_2(normal) > 0.0);
    mNormalToPlane = normal/norm_2(normal);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetPointOnPlane() const
{
    return mPointOnPlane;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetNormalToPlane() const
{
    return mNormalToPlane;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::SetUseJiggledNodesOnPlane(bool useJiggledNodesOnPlane)
{
    mUseJiggledNodesOnPlane = useJiggledNodesOnPlane;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::GetUseJiggledNodesOnPlane()
{
    return mUseJiggledNodesOnPlane;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation)==NULL)
    {
        EXCEPTION("AttractingPlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation))
            || (SPACE_DIM==ELEMENT_DIM && (dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation))) );

    // This is a magic number
    double max_jiggle = 1e-4;

    if (SPACE_DIM != 1)
    {
        if (dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation))
        {
            for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                 cell_iter != this->mpCellPopulation->End();
                 ++cell_iter)
            {
                unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

                double signed_distance = inner_prod(node_location - mPointOnPlane, mNormalToPlane);
                if (signed_distance > 0.0 or std::abs(signed_distance) < mAttractionThreshold)
                {
                    // For the closest point on the plane we travel from node_location the signed_distance in the direction of -mNormalToPlane
                    c_vector<double, SPACE_DIM> nearest_point;
                    if (mUseJiggledNodesOnPlane)
                    {
                        nearest_point = node_location - (signed_distance+max_jiggle*RandomNumberGenerator::Instance()->ranf())*mNormalToPlane;
                    }
                    else
                    {
                        nearest_point = node_location - signed_distance*mNormalToPlane;
                    }
                    p_node->rGetModifiableLocation() = nearest_point;
                }
            }
        }
        else
        {
            assert(SPACE_DIM == ELEMENT_DIM); // LCOV_EXCL_LINE
            assert(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));

            // Iterate over all nodes and update their positions according to the boundary conditions
            unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

                double signed_distance = inner_prod(node_location - mPointOnPlane, mNormalToPlane);
                if (signed_distance > 0.0 or std::abs(signed_distance) < mAttractionThreshold)
                {
                    // For the closest point on the plane we travel from node_location the signed_distance in the direction of -mNormalToPlane
                    c_vector<double, SPACE_DIM> nearest_point;
                    if (mUseJiggledNodesOnPlane)
                    {
                        nearest_point = node_location - (signed_distance+max_jiggle*RandomNumberGenerator::Instance()->ranf())*mNormalToPlane;
                    }
                    else
                    {
                        nearest_point = node_location - signed_distance*mNormalToPlane;
                    }
                    p_node->rGetModifiableLocation() = nearest_point;
                }
            }
        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        //AttractingPlaneBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    if (SPACE_DIM == 1)
    {
        EXCEPTION("AttractingPlaneBoundaryCondition is not implemented in 1D");
    }
    else
    {
        for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, SPACE_DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

            if (inner_prod(cell_location - mPointOnPlane, mNormalToPlane) > 0.0)
            {
                condition_satisfied = false;
                break;
            }
        }
    }

    return condition_satisfied;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::SetPointOnPlane(const c_vector<double, SPACE_DIM>& rPoint)
{
    mPointOnPlane = rPoint;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::SetAttractionThreshold(double attractionThreshold)
{
    mAttractionThreshold = attractionThreshold;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AttractingPlaneBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PointOnPlane>";
    for (unsigned index=0; index != SPACE_DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnPlane[index] << ",";
    }
    *rParamsFile << mPointOnPlane[SPACE_DIM-1] << "</PointOnPlane>\n";

    *rParamsFile << "\t\t\t<NormalToPlane>";
    for (unsigned index=0; index != SPACE_DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mNormalToPlane[index] << ",";
    }
    *rParamsFile << mNormalToPlane[SPACE_DIM-1] << "</NormalToPlane>\n";
    *rParamsFile << "\t\t\t<UseJiggledNodesOnPlane>" << mUseJiggledNodesOnPlane << "</UseJiggledNodesOnPlane>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class AttractingPlaneBoundaryCondition<1,1>;
template class AttractingPlaneBoundaryCondition<1,2>;
template class AttractingPlaneBoundaryCondition<2,2>;
template class AttractingPlaneBoundaryCondition<1,3>;
template class AttractingPlaneBoundaryCondition<2,3>;
template class AttractingPlaneBoundaryCondition<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AttractingPlaneBoundaryCondition)
