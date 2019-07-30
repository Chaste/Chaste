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

#include "DiffusionCaUpdateRule.hpp"

template<unsigned DIM>
DiffusionCaUpdateRule<DIM>::DiffusionCaUpdateRule()
    : AbstractCaUpdateRule<DIM>(),
      mDiffusionParameter(0.5)
{
}

template<unsigned DIM>
DiffusionCaUpdateRule<DIM>::~DiffusionCaUpdateRule()
{
}

template<unsigned DIM>
double DiffusionCaUpdateRule<DIM>::EvaluateProbability(unsigned currentNodeIndex,
                                                               unsigned targetNodeIndex,
                                                               CaBasedCellPopulation<DIM>& rCellPopulation,
                                                               double dt,
                                                               double deltaX,
                                                               CellPtr cell)
{
   c_vector<double, DIM> node_index_location = rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
   c_vector<double, DIM> node_neighbour_location = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation();

   return (mDiffusionParameter*dt/(2* pow(norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(node_index_location, node_neighbour_location)), 2)));
}

template<unsigned DIM>
double DiffusionCaUpdateRule<DIM>::GetDiffusionParameter()
{
    return mDiffusionParameter;
}

template<unsigned DIM>
void DiffusionCaUpdateRule<DIM>::SetDiffusionParameter(double diffusionParameter)
{
    mDiffusionParameter = diffusionParameter;
}

template<unsigned DIM>
void DiffusionCaUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DiffusionParameter>" << mDiffusionParameter << "</DiffusionParameter>\n";

    // Call method on direct parent class
    AbstractCaUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class DiffusionCaUpdateRule<1u>;
template class DiffusionCaUpdateRule<2u>;
template class DiffusionCaUpdateRule<3u>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionCaUpdateRule)
