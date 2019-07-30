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

#include "AbstractCorrectionTermAssembler.hpp"
#include <typeinfo>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractCorrectionTermAssembler(
        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
        AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* pTissue)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,false,CARDIAC>(pMesh,pTissue)
{
    // Work out which elements can do SVI
    mElementsCanDoSvi.resize(pMesh->GetNumElements(), true);
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = pMesh->GetElementIteratorBegin();
         iter != pMesh->GetElementIteratorEnd();
         ++iter)
    {
        Element<ELEMENT_DIM, SPACE_DIM>& r_element = *iter;
        if (r_element.GetOwnership())
        {
            unsigned element_index = r_element.GetIndex();
            // If bath element, don't try to use SVI
            if (HeartRegionCode::IsRegionBath(r_element.GetUnsignedAttribute()))
            {
                mElementsCanDoSvi[element_index] = false;
                continue;
            }
            // See if the nodes in this element all use the same cell model
            unsigned node_zero = r_element.GetNodeGlobalIndex(0);
            AbstractCardiacCellInterface* p_cell_zero = this->mpCardiacTissue->GetCardiacCellOrHaloCell(node_zero);
            const std::type_info& r_zero_info = typeid(*p_cell_zero);
            // Check the other nodes match. If they don't, no SVI
            for (unsigned local_index=1; local_index<r_element.GetNumNodes(); local_index++)
            {
                unsigned global_index = r_element.GetNodeGlobalIndex(local_index);
                AbstractCardiacCellInterface* p_cell = this->mpCardiacTissue->GetCardiacCellOrHaloCell(global_index);
                const std::type_info& r_info = typeid(*p_cell);
                if (r_zero_info != r_info)
                {
                    mElementsCanDoSvi[element_index] = false;
                    break;
                }
            }
        }
    }
    // Note: the mStateVariables std::vector is resized if correction will be applied to a given element
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ResetInterpolatedQuantities()
{
    // reset ionic current, and state variables
    mIionicInterp = 0;
    for (unsigned i=0; i<mStateVariablesAtQuadPoint.size(); i++)
    {
        mStateVariablesAtQuadPoint[i] = 0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::IncrementInterpolatedQuantities(
            double phiI, const Node<SPACE_DIM>* pNode)
{
    // interpolate ionic current
    unsigned node_global_index = pNode->GetIndex();
    mIionicInterp  += phiI * this->mpCardiacTissue->rGetIionicCacheReplicated()[ node_global_index ];
    // and state variables
    std::vector<double> state_vars = this->mpCardiacTissue->GetCardiacCellOrHaloCell(node_global_index)->GetStdVecStateVariables();
    for (unsigned i=0; i<mStateVariablesAtQuadPoint.size(); i++)
    {
        mStateVariablesAtQuadPoint[i] += phiI * state_vars[i];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ElementAssemblyCriterion(Element<ELEMENT_DIM,SPACE_DIM>& rElement)
{
    // Check that SVI is allowed on this element
    if (!mElementsCanDoSvi[rElement.GetIndex()])
    {
        return false;
    }
    double DELTA_IIONIC = 1; // tolerance

    //The criterion and the correction both need the ionic cache, so we better make sure that it's up-to-date
    assert(this->mpCardiacTissue->GetDoCacheReplication());
    ReplicatableVector& r_cache = this->mpCardiacTissue->rGetIionicCacheReplicated();

    double diionic = fabs(r_cache[rElement.GetNodeGlobalIndex(0)] - r_cache[rElement.GetNodeGlobalIndex(1)]);

    if (ELEMENT_DIM > 1)
    {
        diionic = std::max(diionic, fabs(r_cache[rElement.GetNodeGlobalIndex(0)] - r_cache[rElement.GetNodeGlobalIndex(2)]) );
        diionic = std::max(diionic, fabs(r_cache[rElement.GetNodeGlobalIndex(1)] - r_cache[rElement.GetNodeGlobalIndex(2)]) );
    }

    if (ELEMENT_DIM > 2)
    {
        diionic = std::max(diionic, fabs(r_cache[rElement.GetNodeGlobalIndex(0)] - r_cache[rElement.GetNodeGlobalIndex(3)]) );
        diionic = std::max(diionic, fabs(r_cache[rElement.GetNodeGlobalIndex(1)] - r_cache[rElement.GetNodeGlobalIndex(3)]) );
        diionic = std::max(diionic, fabs(r_cache[rElement.GetNodeGlobalIndex(2)] - r_cache[rElement.GetNodeGlobalIndex(3)]) );
    }

    bool will_assemble = (diionic > DELTA_IIONIC);

    if (will_assemble)
    {
        unsigned any_node = rElement.GetNodeGlobalIndex(0);
        mStateVariablesAtQuadPoint.resize(this->mpCardiacTissue->GetCardiacCellOrHaloCell(any_node)->GetNumberOfStateVariables() );
    }

    return will_assemble;
}

// Explicit instantiation
template class AbstractCorrectionTermAssembler<1,1,1>;
template class AbstractCorrectionTermAssembler<1,2,1>;
template class AbstractCorrectionTermAssembler<1,3,1>;
template class AbstractCorrectionTermAssembler<2,2,1>;
template class AbstractCorrectionTermAssembler<3,3,1>;
template class AbstractCorrectionTermAssembler<1,1,2>;
template class AbstractCorrectionTermAssembler<2,2,2>;
template class AbstractCorrectionTermAssembler<3,3,2>;
