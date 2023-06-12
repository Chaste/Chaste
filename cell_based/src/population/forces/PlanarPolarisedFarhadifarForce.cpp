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

#include "PlanarPolarisedFarhadifarForce.hpp"

template<unsigned DIM>
PlanarPolarisedFarhadifarForce<DIM>::PlanarPolarisedFarhadifarForce()
   : FarhadifarForce<DIM>(),
     mPlanarPolarisedLineTensionMultiplier(2.0)
{
}

template<unsigned DIM>
PlanarPolarisedFarhadifarForce<DIM>::~PlanarPolarisedFarhadifarForce()
{
}

template<unsigned DIM>
double PlanarPolarisedFarhadifarForce<DIM>::GetPlanarPolarisedLineTensionMultiplier()
{
    return mPlanarPolarisedLineTensionMultiplier;
}

template<unsigned DIM>
void PlanarPolarisedFarhadifarForce<DIM>::SetPlanarPolarisedLineTensionMultiplier(double planarPolarisedLineTensionMultiplier)
{
    mPlanarPolarisedLineTensionMultiplier = planarPolarisedLineTensionMultiplier;
}

template<unsigned DIM>
double PlanarPolarisedFarhadifarForce<DIM>::GetLineTensionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // Since each internal edge is visited twice in the loop above, we have to use half the line tension parameter
    // for each visit.
    double line_tension_parameter_in_calculation = this->mLineTensionParameter/2.0;

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        line_tension_parameter_in_calculation = this->mBoundaryLineTensionParameter;
    }

    // Get the vector between the two vertices
    c_vector<double, 2> vector = pNodeB->rGetLocation() - pNodeA->rGetLocation();
    double theta = atan2(fabs(vector(1)), fabs(vector(0)));

    if ((theta > 0.25*M_PI) && (theta < 0.75*M_PI))
    {
        line_tension_parameter_in_calculation *= mPlanarPolarisedLineTensionMultiplier;
    }

    return line_tension_parameter_in_calculation;
}

template<unsigned DIM>
double PlanarPolarisedFarhadifarForce<DIM>::GetLineTensionParameter()
{
    return this->mLineTensionParameter;
}

template<unsigned DIM>
double PlanarPolarisedFarhadifarForce<DIM>::GetBoundaryLineTensionParameter()
{
    return this->mBoundaryLineTensionParameter;
}


template<unsigned DIM>
void PlanarPolarisedFarhadifarForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PlanarPolarisedLineTensionMultiplier>" << mPlanarPolarisedLineTensionMultiplier << "</PlanarPolarisedLineTensionMultiplier>\n";

    // Call method on direct parent class
    FarhadifarForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PlanarPolarisedFarhadifarForce<1>;
template class PlanarPolarisedFarhadifarForce<2>;
template class PlanarPolarisedFarhadifarForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlanarPolarisedFarhadifarForce)
