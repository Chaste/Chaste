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

#include "AbstractCorrectionTermAssembler.hpp"
#include <typeinfo>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractCorrectionTermAssembler(
        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
        AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
        unsigned numQuadPoints)
    : AbstractCardiacFeVolumeIntegralAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM,true,false,CARDIAC>(pMesh,pTissue,numQuadPoints)
{
    // Work out which elements have the same cell at every node, and hence can have SVI done
    mElementsHasIdenticalCellModels.resize(pMesh->GetNumElements(), true);
    for (typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator iter = pMesh->GetElementIteratorBegin();
         iter != pMesh->GetElementIteratorEnd();
         ++iter)
    {
        Element<ELEMENT_DIM, SPACE_DIM>& r_element = *iter;
        if (r_element.GetOwnership())
        {
            unsigned node_zero = r_element.GetNodeGlobalIndex(0);
            AbstractCardiacCell* p_cell_zero = this->mpCardiacTissue->GetCardiacCellOrHaloCell(node_zero);
            const std::type_info& r_zero_info = typeid(*p_cell_zero);
            // Check the other nodes match
            for (unsigned local_index=1; local_index<r_element.GetNumNodes(); local_index++)
            {
                unsigned global_index = r_element.GetNodeGlobalIndex(local_index);
                AbstractCardiacCell* p_cell = this->mpCardiacTissue->GetCardiacCellOrHaloCell(global_index);
                const std::type_info& r_info = typeid(*p_cell);
                if (r_zero_info != r_info)
                {
                    mElementsHasIdenticalCellModels[r_element.GetIndex()] = false;
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
    for(unsigned i=0; i<mStateVariablesAtQuadPoint.size(); i++)
    {
        mStateVariablesAtQuadPoint[i] = 0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::IncrementInterpolatedQuantities(
            double phiI, const Node<SPACE_DIM>* pNode)
{
    // interpolate ionic current, and state variables

    unsigned node_global_index = pNode->GetIndex();

    mIionicInterp  += phiI * this->mpCardiacTissue->rGetIionicCacheReplicated()[ node_global_index ];
    for (unsigned i=0; i<mStateVariablesAtQuadPoint.size(); i++)
    {
        mStateVariablesAtQuadPoint[i] += phiI * this->mpCardiacTissue->GetCardiacCellOrHaloCell(node_global_index)->rGetStateVariables()[i];
    }
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ElementAssemblyCriterion(Element<ELEMENT_DIM,SPACE_DIM>& rElement)
{
    // if element doesn't have identical cell models, can't do SVI.
    if (!mElementsHasIdenticalCellModels[rElement.GetIndex()])
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
        mStateVariablesAtQuadPoint.resize(this->mpCardiacTissue->GetCardiacCellOrHaloCell(any_node)->rGetStateVariables().size());
    }

    return will_assemble;
}

///////////////////////////////////////////////////////
// explicit instantiation
///////////////////////////////////////////////////////

template class AbstractCorrectionTermAssembler<1,1,1>;
template class AbstractCorrectionTermAssembler<1,2,1>;
template class AbstractCorrectionTermAssembler<1,3,1>;
template class AbstractCorrectionTermAssembler<2,2,1>;
template class AbstractCorrectionTermAssembler<3,3,1>;
template class AbstractCorrectionTermAssembler<1,1,2>;
template class AbstractCorrectionTermAssembler<2,2,2>;
template class AbstractCorrectionTermAssembler<3,3,2>;
