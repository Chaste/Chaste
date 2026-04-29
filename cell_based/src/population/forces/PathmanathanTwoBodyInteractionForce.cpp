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

#include "PathmanathanTwoBodyInteractionForce.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::PathmanathanTwoBodyInteractionForce()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mSpringStiffness(15.0),
     mAlpha(5.0)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::~PathmanathanTwoBodyInteractionForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
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

    // Get the node radii
    double node_a_radius = p_node_a->GetRadius();
    double node_b_radius = p_node_b->GetRadius();
    assert(node_a_radius > 0 && node_b_radius > 0);

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

    // Rest length is the sum of the two cell radii
    double rest_length = node_a_radius + node_b_radius;

    double overlap = distance_between_nodes - rest_length;

    bool is_closer_than_rest_length = (overlap <= 0);
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);

    // Logarithmic repulsion (cells closer than rest length, overlap is negative)
    if (overlap <= 0)
    {
        //log(x+1) is undefined for x<=-1
        assert(overlap > -rest_length);
        return multiplication_factor * mSpringStiffness * unit_difference * rest_length * log(1.0 + overlap/rest_length);
    }
    else
    {
        // Exponential attraction (cells further than rest length, overlap is positive)
        return multiplication_factor * mSpringStiffness * unit_difference * overlap * exp(-mAlpha * overlap/rest_length);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                                                                              unsigned nodeBGlobalIndex,
                                                                                                              AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                                                                              bool isCloserThanRestLength)
{
    return 1.0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetSpringStiffness()
{
    return mSpringStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::SetSpringStiffness(double springStiffness)
{
    assert(springStiffness > 0.0);
    mSpringStiffness = springStiffness;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::GetAlpha()
{
    return mAlpha;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::SetAlpha(double alpha)
{
    assert(alpha > 0.0);
    mAlpha = alpha;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PathmanathanTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SpringStiffness>" << mSpringStiffness << "</SpringStiffness>\n";
    *rParamsFile << "\t\t\t<Alpha>" << mAlpha << "</Alpha>\n";

    // Call method on direct parent class
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PathmanathanTwoBodyInteractionForce<1,1>;
template class PathmanathanTwoBodyInteractionForce<1,2>;
template class PathmanathanTwoBodyInteractionForce<2,2>;
template class PathmanathanTwoBodyInteractionForce<1,3>;
template class PathmanathanTwoBodyInteractionForce<2,3>;
template class PathmanathanTwoBodyInteractionForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PathmanathanTwoBodyInteractionForce)
