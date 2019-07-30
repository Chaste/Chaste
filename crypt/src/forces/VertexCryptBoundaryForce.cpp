/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "VertexCryptBoundaryForce.hpp"
#include "MathsCustomFunctions.hpp"

template<unsigned DIM>
VertexCryptBoundaryForce<DIM>::VertexCryptBoundaryForce(double forceStrength)
   : AbstractForce<DIM>(),
     mForceStrength(forceStrength)
{
    // We don't want the force to act in the wrong direction
    assert(mForceStrength > 0.0);
}

template<unsigned DIM>
VertexCryptBoundaryForce<DIM>::~VertexCryptBoundaryForce()
{
}

template<unsigned DIM>
void VertexCryptBoundaryForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("VertexCryptBoundaryForce is to be used with VertexBasedCellPopulations only");
    }

    // Iterate over nodes
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = p_cell_population->rGetMesh().GetNodeIteratorBegin();
         node_iter != p_cell_population->rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        double y = node_iter->rGetLocation()[1]; // y-coordinate of node

        // If the node lies below the line y=0, then add the boundary force contribution to the node forces
        if (y < 0.0)
        {
            c_vector<double, DIM> boundary_force = zero_vector<double>(DIM);
            boundary_force[1] = mForceStrength*SmallPow(y, 2);

            node_iter->AddAppliedForceContribution(boundary_force);
        }
    }
}

template<unsigned DIM>
double VertexCryptBoundaryForce<DIM>::GetForceStrength() const
{
    return mForceStrength;
}

template<unsigned DIM>
void VertexCryptBoundaryForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ForceStrength>" << mForceStrength << "</ForceStrength>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class VertexCryptBoundaryForce<1>;
template class VertexCryptBoundaryForce<2>;
template class VertexCryptBoundaryForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexCryptBoundaryForce)
