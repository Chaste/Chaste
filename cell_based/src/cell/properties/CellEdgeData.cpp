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

#include "CellEdgeData.hpp"

CellEdgeData::~CellEdgeData()
{
}

void CellEdgeData::SetItem(const std::string &rVariableName, const std::vector<double> &data)
{
    this->mCellEdgeData[rVariableName] = data;
}

std::vector<double> CellEdgeData::GetItem(const std::string &rVariableName) const
{
    /*
     * Note that mCellData[rVariableName] is not const. If rVariableName is not
     * a key, then mCellData[rVariableName] will create a new item in the map
     * and increase the size by one.  Using a const_iterator ensures that the
     * map remains const.
     */
    std::map<std::string, std::vector<double> >::const_iterator it = mCellEdgeData.find(rVariableName);
    if (it == mCellEdgeData.end())
    {
        EXCEPTION("The item " << rVariableName << " is not stored");
    }
    return(it->second);
}

double CellEdgeData::GetItemAtIndex(const std::string &rVariableName, const unsigned int index)
{
    /*
     * Note that mCellData[rVariableName] is not const. If rVariableName is not
     * a key, then mCellData[rVariableName] will create a new item in the map
     * and increase the size by one.  Using a const_iterator ensures that the
     * map remains const.
     */
    std::map<std::string, std::vector<double> >::const_iterator it = mCellEdgeData.find(rVariableName);
    if (it == mCellEdgeData.end())
    {
        EXCEPTION("The item " << rVariableName << " is not stored");
    }
    if (it->second.size() <= index)
    {
        EXCEPTION("The item " << rVariableName << " does not have index " << index);
    }

    return(it->second[index]);
}

unsigned CellEdgeData::GetNumItems() const {
    return mCellEdgeData.size();
}

std::vector<std::string> CellEdgeData::GetKeys() const {
    std::vector<std::string> keys;
    for (auto& kv: mCellEdgeData)
    {
        keys.push_back(kv.first);
    }
    return keys;
}
