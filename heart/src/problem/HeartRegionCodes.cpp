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
#include "HeartRegionCodes.hpp"
#include "HeartConfig.hpp"
#include <cassert>


HeartRegionType HeartRegionCode::GetValidTissueId()
{
    // Returns the identifier of the first region defined as tissue
    assert(!HeartConfig::Instance()->rGetTissueIdentifiers().empty());
    return *HeartConfig::Instance()->rGetTissueIdentifiers().begin();
}

HeartRegionType HeartRegionCode::GetValidBathId()
{
    // Returns the identifier of the first region defined as bath
    assert(!HeartConfig::Instance()->rGetBathIdentifiers().empty());
    return *HeartConfig::Instance()->rGetBathIdentifiers().begin();
}

bool HeartRegionCode::IsRegionTissue(HeartRegionType regionId)
{
    return (HeartConfig::Instance()->rGetTissueIdentifiers().find(regionId) != HeartConfig::Instance()->rGetTissueIdentifiers().end());
}

bool HeartRegionCode::IsRegionBath(HeartRegionType regionId)
{
    return (HeartConfig::Instance()->rGetBathIdentifiers().find(regionId) != HeartConfig::Instance()->rGetBathIdentifiers().end());
}
