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

#include <algorithm>

#include "CellPropertyRegistry.hpp"
#include "Exception.hpp"

CellPropertyRegistry* CellPropertyRegistry::mpInstance = nullptr;

CellPropertyRegistry* CellPropertyRegistry::Instance()
{
    if (mpInstance == nullptr)
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
    mpInstance = nullptr;
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
