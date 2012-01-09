/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "PlaneBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
PlaneBoundaryCondition<DIM>::PlaneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                    c_vector<double, DIM> point,
                                                    c_vector<double, DIM> normal)
        : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
          mPointOnPlane(point)
{
    assert(norm_2(normal) > 0.0);
    mNormalToPlane = normal/norm_2(normal);
}

template<unsigned DIM>
const c_vector<double, DIM>& PlaneBoundaryCondition<DIM>::rGetPointOnPlane() const
{
    return mPointOnPlane;
}

template<unsigned DIM>
const c_vector<double, DIM>& PlaneBoundaryCondition<DIM>::rGetNormalToPlane() const
{
    return mNormalToPlane;
}

template<unsigned DIM>
void PlaneBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::vector< c_vector<double, DIM> >& rOldLocations)
{
    //TODO Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<DIM>*>(this->mpCellPopulation)==NULL)
    {
        EXCEPTION("PlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(this->mpCellPopulation))
            || (dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation)) );

    if (DIM==2)
    {
        if(dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(this->mpCellPopulation))
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                 cell_iter != this->mpCellPopulation->End();
                 ++cell_iter)
            {
                c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

                unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
                Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

                if (inner_prod(cell_location - mPointOnPlane,mNormalToPlane) > 0.0)
                {
                    c_vector<double, 2> tangent;
                    tangent(0) = -mNormalToPlane(1);
                    tangent(1) = mNormalToPlane(0);

                    c_vector<double, 2> intersection = mPointOnPlane + inner_prod(tangent,cell_location- mPointOnPlane)*tangent;

                    p_node->rGetModifiableLocation() = intersection;
                }
            }
        }
        else
        {
            assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation));

            VertexBasedCellPopulation<DIM>* pStaticCastCellPopulation = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);

            // Iterate over all nodes and update their positions according to the boundary conditions
            unsigned num_nodes = pStaticCastCellPopulation->GetNumNodes();
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<DIM>* p_node = pStaticCastCellPopulation->GetNode(node_index);
                c_vector<double, DIM> node_location = p_node->rGetLocation();

                if (inner_prod(node_location - mPointOnPlane,mNormalToPlane) > 0.0)
                {
                    c_vector<double, 2> tangent;
                    tangent(0) = -mNormalToPlane(1);
                    tangent(1) = mNormalToPlane(0);

                    c_vector<double, 2> intersection = mPointOnPlane + inner_prod(tangent,node_location- mPointOnPlane)*tangent;

                    p_node->rGetModifiableLocation() = intersection;
                }
            }
        }

    }
    else
    {
        EXCEPTION("PlaneBoundaryCondition is not yet implemented in 1D or 3D");
    }
}

template<unsigned DIM>
bool PlaneBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    if (DIM==2)
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

            if (inner_prod(cell_location - mPointOnPlane,mNormalToPlane) > 0.0)
            {
                condition_satisfied = false;
                break;
            }
        }
    }
    else
    {
        EXCEPTION("PlaneBoundaryCondition is not yet implemented in 1D or 3D");
    }

    return condition_satisfied;
}

template<unsigned DIM>
void PlaneBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PointOnPlane>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnPlane[index] << ",";
    }
    *rParamsFile << mPointOnPlane[DIM-1] << "</PointOnPlane>\n";

    *rParamsFile << "\t\t\t<NormalToPlane>";
     for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
     {
         *rParamsFile << mNormalToPlane[index] << ",";
     }
     *rParamsFile << mNormalToPlane[DIM-1] << "</NormalToPlane>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class PlaneBoundaryCondition<1>;
template class PlaneBoundaryCondition<2>;
template class PlaneBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlaneBoundaryCondition)
