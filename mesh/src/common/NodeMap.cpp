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


#include "NodeMap.hpp"
#include "Exception.hpp"


///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


NodeMap::NodeMap(unsigned size)
{
    // this used to be reserve, but this acts oddly:
    // eg: mMap.reserve(2); mMap[0]=1;
    // runs and mMap[0] returns 1, but mMap.size() returns 0
    mMap.resize(size);
}

void NodeMap::Resize(unsigned size)
{
    mMap.resize(size);
}

void NodeMap::ResetToIdentity()
{
    for (unsigned oldIndex=0; oldIndex<mMap.size(); oldIndex++)
    {
        mMap[oldIndex] = oldIndex;
    }
}

void NodeMap::SetNewIndex(unsigned oldIndex, unsigned newIndex)
{
    mMap[oldIndex] = newIndex;
}

void NodeMap::SetDeleted(unsigned index)
{
    mMap[index] = UINT_MAX;
}

bool NodeMap::IsDeleted(unsigned index)
{
    return (mMap[index] == UINT_MAX);
}

unsigned NodeMap::GetNewIndex(unsigned oldIndex) const
{
    if (mMap[oldIndex] == UINT_MAX)
    {
        EXCEPTION("Node has been deleted");
    }
    return (unsigned) mMap[oldIndex];
}

bool NodeMap::IsIdentityMap()
{
    for (unsigned i=0; i<mMap.size(); i++)
    {
        if (mMap[i] != i)
        {
            return false;
        }
    }
    return true;
}

unsigned NodeMap::Size()
{
    return mMap.size();
}
