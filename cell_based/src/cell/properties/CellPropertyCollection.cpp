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

#include "CellPropertyCollection.hpp"

CellPropertyCollection::CellPropertyCollection()
    : mpCellPropertyRegistry(nullptr)
{
}

CellPropertyRegistry* CellPropertyCollection::GetCellPropertyRegistry()
{
    if (!mpCellPropertyRegistry)
    {
        /*
         * If no cell population has taken ownership of the registry and assigned
         * this class a pointer, then store a pointer to the singleton.
         */
        mpCellPropertyRegistry = CellPropertyRegistry::Instance();
    }
    return mpCellPropertyRegistry;
}

void CellPropertyCollection::SetCellPropertyRegistry(CellPropertyRegistry* pRegistry)
{
    mpCellPropertyRegistry = pRegistry;
}

void CellPropertyCollection::AddProperty(const boost::shared_ptr<AbstractCellProperty>& rProp)
{
    if (HasProperty(rProp))
    {
        EXCEPTION("That property object is already in the collection.");
    }
    mProperties.insert(rProp);
}

bool CellPropertyCollection::HasProperty(const boost::shared_ptr<AbstractCellProperty>& rProp) const
{
    return (mProperties.find(rProp) != mProperties.end());
}

void CellPropertyCollection::RemoveProperty(const boost::shared_ptr<AbstractCellProperty>& rProp)
{
    IteratorType it = mProperties.find(rProp);
    if (it == mProperties.end())
    {
        EXCEPTION("Collection does not contain the given property.");
    }
    else
    {
        mProperties.erase(it);
    }
}

unsigned CellPropertyCollection::GetSize() const
{
    return mProperties.size();
}

CellPropertyCollection::Iterator CellPropertyCollection::Begin()
{
    return mProperties.begin();
}

CellPropertyCollection::Iterator CellPropertyCollection::End()
{
    return mProperties.end();
}

boost::shared_ptr<AbstractCellProperty> CellPropertyCollection::GetProperty() const
{
    if (GetSize() == 1)
    {
        return *mProperties.begin();
    }
    else
    {
        EXCEPTION("Can only call GetProperty on a collection of size 1.");
    }
}
