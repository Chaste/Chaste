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

#include "CellVecData.hpp"

CellVecData::CellVecData()
    : mFreeVecs(false)
{
}

CellVecData::CellVecData(const CellVecData& rAnotherCellVecData)
    : mFreeVecs(true)
{
    std::vector<std::string> keys = rAnotherCellVecData.GetKeys();

    for (std::vector<std::string>::iterator iter = keys.begin(); iter != keys.end(); ++iter)
    {
        std::string map_key = *iter;
        Vec map_value = rAnotherCellVecData.GetItem(map_key);

        Vec map_value_copy;
        VecDuplicate(map_value, &map_value_copy);
        VecCopy(map_value, map_value_copy);

        mCellVecData[map_key] = map_value_copy;
    }
}

CellVecData::CellVecData(const std::map<std::string, Vec>& rCellVecDataMap)
    : mCellVecData(rCellVecDataMap),
      mFreeVecs(true)
{
}

CellVecData::~CellVecData()
{
    // If the object was loaded from a checkpoint, the Vecs in the map need freeing. Otherwise is the user's responsibility.
    if (mFreeVecs)
    {
        for (std::map<std::string, Vec>::iterator iter = mCellVecData.begin(); iter != mCellVecData.end(); ++iter)
        {
            PetscTools::Destroy(iter->second);
        }
    }
}

void CellVecData::SetItem(const std::string& rVariableName, Vec data)
{
    mCellVecData[rVariableName] = data;
}

Vec CellVecData::GetItem(const std::string& rVariableName) const
{
    /*
     * Note that mCellVecData[rVariableName] is not const. If rVariableName is not
     * a key, then mCellVecData[rVariableName] will create a new item in the map
     * and increase the size by one.  Using a const_iterator ensures that the
     * map remains const.
     */
    std::map<std::string, Vec>::const_iterator it = mCellVecData.find(rVariableName);
    if (it == mCellVecData.end())
    {
        EXCEPTION("The item " << rVariableName << " is not stored");
    }
    return it->second;
}

unsigned CellVecData::GetNumItems() const
{
    return mCellVecData.size();
}

std::vector<std::string> CellVecData::GetKeys() const
{
    std::vector<std::string> keys;
    for (std::map<std::string, Vec>::const_iterator it = mCellVecData.begin(); it != mCellVecData.end(); ++it)
    {
        keys.push_back(it->first);
    }

    // From STL documentation we assume that the iterator is returning sorted keys sort(keys.begin(), keys.end());
    return keys;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellVecData)
