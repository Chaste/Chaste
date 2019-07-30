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


#include "NodeMap.hpp"
#include "Exception.hpp"
#include <cassert>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////


NodeMap::NodeMap(unsigned size)
    : mIsIdentity(true)
{
    mMap.resize(size);
}

void NodeMap::Resize(unsigned size)
{
    mMap.resize(size);
}

void NodeMap::ResetToIdentity()
{
    for (unsigned index=0; index<mMap.size(); index++)
    {
        mMap[index] = index;  // This is needed if we later diverge from identity (e.g. delete nodes)
    }
    mIsIdentity = true;
}

void NodeMap::SetNewIndex(unsigned oldIndex, unsigned newIndex)
{
    assert(oldIndex<mMap.size());
    mMap[oldIndex] = newIndex;
    if (mIsIdentity && (oldIndex != newIndex))
    {
        mIsIdentity = false;
    }
}

void NodeMap::SetDeleted(unsigned index)
{
    mMap[index] = UINT_MAX;
    mIsIdentity = false;
}

bool NodeMap::IsDeleted(unsigned index)
{
    return (mMap[index] == UINT_MAX);
}

unsigned NodeMap::GetNewIndex(unsigned oldIndex) const
{
    if (mIsIdentity)
    {
        return oldIndex;
    }
    if (mMap[oldIndex] == UINT_MAX)
    {
        EXCEPTION("Node has been deleted");
    }
    return mMap[oldIndex];
}

bool NodeMap::IsIdentityMap()
{
    return mIsIdentity;
}

unsigned NodeMap::GetSize()
{
    return mMap.size();
}
