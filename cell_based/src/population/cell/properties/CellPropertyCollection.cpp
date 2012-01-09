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

#include "CellPropertyCollection.hpp"

CellPropertyCollection::CellPropertyCollection()
    : mpCellPropertyRegistry(NULL)
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
