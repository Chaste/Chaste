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

#include "ChasteNodesList.hpp"

template <unsigned SPACE_DIM>
ChasteNodesList<SPACE_DIM>::ChasteNodesList(const std::vector<Node<SPACE_DIM>*> rNodesList, bool ownNodes)
    : mListOfNodes(rNodesList),
      mOwnNodes(ownNodes)
{
}

template <unsigned SPACE_DIM>
ChasteNodesList<SPACE_DIM>::~ChasteNodesList()
{
    if (mOwnNodes)
    {
        for (unsigned i=0; i<mListOfNodes.size(); i++)
        {
            delete mListOfNodes[i];
        }
    }
}

template <unsigned SPACE_DIM>
bool ChasteNodesList<SPACE_DIM>::DoesContain(const ChastePoint<SPACE_DIM>& rPointToCheck) const
{
    bool returned_value = false;
    for (unsigned index = 0; index < mListOfNodes.size(); index++)
    {
        if (mListOfNodes[index]->GetPoint().IsSamePoint(rPointToCheck))
        {
            returned_value = true;
            break;
        }
    }

    return returned_value;
}

template <unsigned SPACE_DIM>
const std::vector< Node<SPACE_DIM>*>& ChasteNodesList<SPACE_DIM>::rGetNodesList() const
{
    return mListOfNodes;
}

template <unsigned SPACE_DIM>
unsigned ChasteNodesList<SPACE_DIM>::GetSize() const
{
    return mListOfNodes.size();
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ChasteNodesList<1>;
template class ChasteNodesList<2>;
template class ChasteNodesList<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChasteNodesList)
