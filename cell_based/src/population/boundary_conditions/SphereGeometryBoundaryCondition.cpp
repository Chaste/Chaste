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

#include "SphereGeometryBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
SphereGeometryBoundaryCondition<DIM>::SphereGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                                      c_vector<double, DIM> centre,
                                                                      double radius,
                                                                      double distance)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mCentreOfSphere(centre),
      mRadiusOfSphere(radius),
      mMaximumDistance(distance)
{
    assert(mRadiusOfSphere > 0.0);
    assert(mMaximumDistance > 0.0);

    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == NULL)
    {
        EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
    if (DIM == 1)
    {
        EXCEPTION("This boundary condition is not implemented in 1D.");
    }
}

template<unsigned DIM>
const c_vector<double, DIM>& SphereGeometryBoundaryCondition<DIM>::rGetCentreOfSphere() const
{
    return mCentreOfSphere;
}

template<unsigned DIM>
double SphereGeometryBoundaryCondition<DIM>::GetRadiusOfSphere() const
{
    return mRadiusOfSphere;
}

template<unsigned DIM>
void SphereGeometryBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::vector<c_vector<double, DIM> >& rOldLocations)
{
    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Find the radial distance between this cell and the surface of the sphere
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        double radius = norm_2(cell_location - mCentreOfSphere);

        // If the cell is too far from the surface of the sphere...
        if (fabs(radius - mRadiusOfSphere) > mMaximumDistance)
        {
            // ...move the cell back onto the surface of the sphere
            c_vector<double, DIM> location_on_sphere =
                mCentreOfSphere + mRadiusOfSphere*(cell_location - mCentreOfSphere)/radius;

            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            p_node->rGetModifiableLocation() = location_on_sphere;
        }
    }
}

template<unsigned DIM>
bool SphereGeometryBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Find the radial distance between this cell and the surface of the sphere
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        double radius = norm_2(cell_location - mCentreOfSphere);

        // If the cell is too far from the surface of the sphere...
        if (fabs(radius - mRadiusOfSphere) > mMaximumDistance)
        {
            // ...then the boundary condition is not satisfied
            condition_satisfied = false;
            break;
        }
    }
    return condition_satisfied;
}

template<unsigned DIM>
void SphereGeometryBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CentreOfSphere>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mCentreOfSphere[index] << ",";
    }
    *rParamsFile << mCentreOfSphere[DIM-1] << "</CentreOfSphere>\n";

    *rParamsFile << "\t\t\t<RadiusOfSphere>" << mRadiusOfSphere << "</RadiusOfSphere>\n";
    *rParamsFile << "\t\t\t<MaximumDistance>" << mMaximumDistance << "</MaximumDistance>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SphereGeometryBoundaryCondition<1>;
template class SphereGeometryBoundaryCondition<2>;
template class SphereGeometryBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SphereGeometryBoundaryCondition)
