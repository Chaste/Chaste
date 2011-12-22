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


#include "AbstractPurkinjeCellFactory.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::AbstractPurkinjeCellFactory()
    : AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>(),
      mpMixedDimensionMesh(NULL)
{
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    mpMixedDimensionMesh = dynamic_cast<MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>*>(pMesh);
    if (mpMixedDimensionMesh ==NULL)
    {
        EXCEPTION("AbstractPurkinjeCellFactory must take a MixedDimensionMesh");
    }
    mLocalPurkinjeNodes.clear();
    for (typename MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>::CableElementIterator iter = mpMixedDimensionMesh->GetCableElementIteratorBegin();
          iter != mpMixedDimensionMesh->GetCableElementIteratorEnd();
          ++iter)
    {
        mLocalPurkinjeNodes.insert((*iter)->GetNodeGlobalIndex(0u));
        mLocalPurkinjeNodes.insert((*iter)->GetNodeGlobalIndex(1u));
    }
    AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>::SetMesh(pMesh);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCardiacCell*  AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::CreatePurkinjeCellForNode(
    unsigned nodeIndex)
{
    if(mLocalPurkinjeNodes.count(nodeIndex)>0)
    {
        return CreatePurkinjeCellForTissueNode(nodeIndex);
    }
    else
    {
        return new FakeBathCell(this->mpSolver, this->mpZeroStimulus);
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>::GetMixedDimensionMesh()
{
    if (mpMixedDimensionMesh == NULL)
    {
        EXCEPTION("The mixed dimension mesh object has not been set in the cell factory");
    }
    return mpMixedDimensionMesh;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractPurkinjeCellFactory<1,1>;
template class AbstractPurkinjeCellFactory<2,2>;
template class AbstractPurkinjeCellFactory<3,3>;
template class AbstractPurkinjeCellFactory<1,2>;
template class AbstractPurkinjeCellFactory<1,3>;
