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

#include <algorithm>

#include "CellPropertyRegistry.hpp"
#include "Exception.hpp"

CellPropertyRegistry* CellPropertyRegistry::mpInstance = NULL;

CellPropertyRegistry* CellPropertyRegistry::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CellPropertyRegistry;
    }
    return mpInstance;
}

const std::vector<boost::shared_ptr<AbstractCellProperty> >& CellPropertyRegistry::rGetAllCellProperties()
{
    return mCellProperties;
}

void CellPropertyRegistry::Clear()
{
    mCellProperties.clear();
    mOrderingHasBeenSpecified = false;
}

CellPropertyRegistry::CellPropertyRegistry()
    : mOrderingHasBeenSpecified(false)
{
}

CellPropertyRegistry* CellPropertyRegistry::TakeOwnership()
{
    mpInstance = NULL;
    return this;
}

void CellPropertyRegistry::SpecifyOrdering(const std::vector<boost::shared_ptr<AbstractCellProperty> >& rOrdering)
{
    if (mOrderingHasBeenSpecified)
    {
        EXCEPTION("An ordering has already been specified.");
    }

    std::vector<boost::shared_ptr<AbstractCellProperty> > temp_vector = rOrdering;
    for (unsigned i=0; i<mCellProperties.size(); i++)
    {
        std::vector<boost::shared_ptr<AbstractCellProperty> >::const_iterator it
            = find(rOrdering.begin(), rOrdering.end(), mCellProperties[i]);
        if (it == rOrdering.end())
        {
            temp_vector.push_back(mCellProperties[i]);
        }
    }
    mCellProperties = temp_vector;

    mOrderingHasBeenSpecified = true;
}

bool CellPropertyRegistry::HasOrderingBeenSpecified()
{
    return mOrderingHasBeenSpecified;
}
