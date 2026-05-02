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

#include "AbstractVariableSizeTwoBodyInteractionForce.hpp"

#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::AbstractVariableSizeTwoBodyInteractionForce()
    : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
      mSpringStiffness(15.0),
      mDivisionRestingSpringLength(0.5),
      mSpringGrowthDuration(1.0)
{
    if (SPACE_DIM == 1)
    {
        mSpringStiffness = 30.0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::~AbstractVariableSizeTwoBodyInteractionForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{
    return 1.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    const c_vector<double, SPACE_DIM>& r_node_a_location = p_node_a->rGetLocation();
    const c_vector<double, SPACE_DIM>& r_node_b_location = p_node_b->rGetLocation();

    // Get the node radii for a NodeBasedCellPopulation
    double node_a_radius = 0.0;
    double node_b_radius = 0.0;

    if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        node_a_radius = p_node_a->GetRadius();
        node_b_radius = p_node_b->GetRadius();
    }

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, SPACE_DIM> unit_difference;
    /*
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);

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
            return zero_vector<double>(SPACE_DIM); // c_vector<double,SPACE_DIM>() is not guaranteed to be fresh memory
        }
    }

    /*
     * Calculate the rest length of the spring connecting the two nodes with a default
     * value of 1.0.
     */
    double rest_length_final = 1.0;

    if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        rest_length_final = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex, nodeBGlobalIndex);
    }
    else if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        rest_length_final = node_a_radius + node_b_radius;
    }

    double rest_length = rest_length_final;

    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();

    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));

    /*
     * If the cells are both newly divided, then the rest length of the spring
     * connecting them grows linearly with time, until 1 hour after division.
     */
    if (ageA < mSpringGrowthDuration && ageB < mSpringGrowthDuration)
    {
        AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

        if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        {
            // Spring rest length increases from a small value to the normal rest length over 1 hour
            double lambda = mDivisionRestingSpringLength;
            rest_length = lambda + (rest_length_final - lambda) * ageA/mSpringGrowthDuration;
        }
        if (ageA + SimulationTime::Instance()->GetTimeStep() >= mSpringGrowthDuration)
        {
            // This spring is about to go out of scope
            p_static_cast_cell_population->UnmarkSpring(cell_pair);
        }
    }

    /*
     * For apoptosis, progressively reduce the radius of the cell.
     */
    double a_rest_length = rest_length*0.5;
    double b_rest_length = a_rest_length;

    if (bool(dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        assert(node_a_radius > 0 && node_b_radius > 0);
        a_rest_length = (node_a_radius/(node_a_radius+node_b_radius))*rest_length;
        b_rest_length = (node_b_radius/(node_a_radius+node_b_radius))*rest_length;
    }

    /*
     * If either of the cells has begun apoptosis, then the length of the spring
     * connecting them decreases linearly with time.
     */
    if (p_cell_A->HasApoptosisBegun())
    {
        double time_until_death_a = p_cell_A->GetTimeUntilDeath();
        a_rest_length = a_rest_length * time_until_death_a / p_cell_A->GetApoptosisTime();
    }
    if (p_cell_B->HasApoptosisBegun())
    {
        double time_until_death_b = p_cell_B->GetTimeUntilDeath();
        b_rest_length = b_rest_length * time_until_death_b / p_cell_B->GetApoptosisTime();
    }

    rest_length = a_rest_length + b_rest_length;

    double overlap = distance_between_nodes - rest_length;
    bool is_closer_than_rest_length = (overlap <= 0);
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex,
                                                                              nodeBGlobalIndex,
                                                                              rCellPopulation,
                                                                              is_closer_than_rest_length);

    // TODO issue 1017 This should be rest_length not rest_length_final, keeping it the same for now to maintain backwards compatibility with existing force laws, but we should eventually change it and update the force laws accordingly.
    return CalculateLinkInteraction(overlap, rest_length_final, unit_difference, multiplication_factor);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetSpringStiffness()
{
    return mSpringStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetDivisionRestingSpringLength()
{
    return mDivisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetSpringGrowthDuration()
{
    return mSpringGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::SetSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mSpringStiffness = springStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::SetDivisionRestingSpringLength(double divisionRestingSpringLength)
{
    assert(divisionRestingSpringLength <= 1.0);
    assert(divisionRestingSpringLength >= 0.0);

    mDivisionRestingSpringLength = divisionRestingSpringLength;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::SetSpringGrowthDuration(double springGrowthDuration)
{
    assert(springGrowthDuration >= 0.0);

    mSpringGrowthDuration = springGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractVariableSizeTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractVariableSizeTwoBodyInteractionForce<1,1>;
template class AbstractVariableSizeTwoBodyInteractionForce<1,2>;
template class AbstractVariableSizeTwoBodyInteractionForce<2,2>;
template class AbstractVariableSizeTwoBodyInteractionForce<1,3>;
template class AbstractVariableSizeTwoBodyInteractionForce<2,3>;
template class AbstractVariableSizeTwoBodyInteractionForce<3,3>;
