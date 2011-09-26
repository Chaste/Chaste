/*

Copyright (C) University of Oxford, 2005-2011

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

#include "GeneralisedPeriodicLinearSpringForce.hpp"

template<unsigned DIM>
GeneralisedPeriodicLinearSpringForce<DIM>::GeneralisedPeriodicLinearSpringForce()
   : AbstractPeriodicTwoBodyInteractionForce<DIM>(),
     mMeinekeSpringStiffness(15.0),        // denoted by mu in Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
     mMeinekeDivisionRestingSpringLength(0.5),
     mMeinekeSpringGrowthDuration(1.0)
{
    /**
     * \todo once AbstractPeriodicTwoBodyInteractionForce (and hence this class) is extended
     * to cope with 1D/3D, reset mMeinekeSpringStiffness to 30.0 if DIM==1 here (see also
     * GeneralisedLinearSpringForce)
     */
}

template<unsigned DIM>
double GeneralisedPeriodicLinearSpringForce<DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                                                     unsigned nodeBGlobalIndex,
                                                                                     AbstractCellPopulation<DIM>& rCellPopulation,
                                                                                     bool isCloserThanRestLength)
{
    return 1.0;
}

template<unsigned DIM>
GeneralisedPeriodicLinearSpringForce<DIM>::~GeneralisedPeriodicLinearSpringForce()
{
}

template<unsigned DIM>
c_vector<double, DIM> GeneralisedPeriodicLinearSpringForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<DIM>& rCellPopulation)
{
    ///\todo extend force class to cope with a NodeBasedCellPopulation (#1856)
    assert(rCellPopulation.IsMeshBasedCellPopulation());

    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = this->mpExtendedMesh->GetNode(nodeAGlobalIndex)->rGetLocation();
    c_vector<double, DIM> node_b_location = this->mpExtendedMesh->GetNode(nodeBGlobalIndex)->rGetLocation();

    /*
     * Get the unit vector parallel to the line joining the two nodes.
     *
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     *
     * todo: Though, if this is a periodic force class, we won't be overloading this method?
     */
    c_vector<double, DIM> unit_difference = this->mpExtendedMesh->GetVectorFromAtoB(node_a_location, node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutOffLength has been set, then there is zero force between
     * two nodes located a distance apart greater than mMechanicsCutOffLength in AbstractTwoBodyInteractionForce.
     */
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(DIM);
        }
    }

    // Calculate the rest length of the spring connecting the two nodes

    double rest_length = 1.0;

    ///\todo Extend force class to cope with newly divided cells (#1856)

    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    ///\todo Extend force class to cope with apoptotic cells (#1856)

    rest_length = a_rest_length + b_rest_length;
    assert(rest_length <= 1.0+1e-12);

    bool is_closer_than_rest_length = (distance_between_nodes - rest_length <= 0);

    // Although in this class the 'spring constant' is a constant parameter, in
    // subclasses it can depend on properties of each of the cells
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(this->mExtendedMeshNodeIndexMap[nodeAGlobalIndex], this->mExtendedMeshNodeIndexMap[nodeBGlobalIndex], rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = mMeinekeSpringStiffness;
    double overlap = distance_between_nodes - rest_length;

    return multiplication_factor * spring_stiffness * unit_difference * overlap;
}

template<unsigned DIM>
double GeneralisedPeriodicLinearSpringForce<DIM>::GetMeinekeSpringStiffness()
{
    return mMeinekeSpringStiffness;
}
template<unsigned DIM>
double GeneralisedPeriodicLinearSpringForce<DIM>::GetMeinekeDivisionRestingSpringLength()
{
    return mMeinekeDivisionRestingSpringLength;
}
template<unsigned DIM>
double GeneralisedPeriodicLinearSpringForce<DIM>::GetMeinekeSpringGrowthDuration()
{
    return mMeinekeSpringGrowthDuration;
}

template<unsigned DIM>
void GeneralisedPeriodicLinearSpringForce<DIM>::SetMeinekeSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mMeinekeSpringStiffness = springStiffness;
}

template<unsigned DIM>
void GeneralisedPeriodicLinearSpringForce<DIM>::SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength <= 1.0);
    assert(divisionRestingSpringLength >= 0.0);

    mMeinekeDivisionRestingSpringLength = divisionRestingSpringLength;
}

template<unsigned DIM>
void GeneralisedPeriodicLinearSpringForce<DIM>::SetMeinekeSpringGrowthDuration(double springGrowthDuration)
{
    assert(springGrowthDuration >= 0.0);

    mMeinekeSpringGrowthDuration = springGrowthDuration;
}

template<unsigned DIM>
void GeneralisedPeriodicLinearSpringForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MeinekeSpringStiffness>" << mMeinekeSpringStiffness << "</MeinekeSpringStiffness>\n";
    *rParamsFile << "\t\t\t<MeinekeDivisionRestingSpringLength>" << mMeinekeDivisionRestingSpringLength << "</MeinekeDivisionRestingSpringLength>\n";
    *rParamsFile << "\t\t\t<MeinekeSpringGrowthDuration>" << mMeinekeSpringGrowthDuration << "</MeinekeSpringGrowthDuration>\n";

    // Call method on direct parent class
    AbstractPeriodicTwoBodyInteractionForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class GeneralisedPeriodicLinearSpringForce<1>;
template class GeneralisedPeriodicLinearSpringForce<2>;
template class GeneralisedPeriodicLinearSpringForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeneralisedPeriodicLinearSpringForce)
