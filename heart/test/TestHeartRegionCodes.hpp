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
#ifndef TESTHEARTREGIONCODES_HPP_
#define TESTHEARTREGIONCODES_HPP_

#include "HeartRegionCodes.hpp"
#include "HeartConfig.hpp"

class TestHeartRegionCodes : public CxxTest::TestSuite
{
    public:

    void TestRegions()
    {
        std::set<unsigned> tissue_ids;
        tissue_ids.insert(0);
        tissue_ids.insert(10);

        std::set<unsigned> bath_ids;
        bath_ids.insert(1);
        bath_ids.insert(2);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids);

        TS_ASSERT(HeartRegionCode::IsRegionTissue(0));
        TS_ASSERT(HeartRegionCode::IsRegionTissue(10));
        TS_ASSERT(!HeartRegionCode::IsRegionBath(0));
        TS_ASSERT(!HeartRegionCode::IsRegionBath(10));

        TS_ASSERT(HeartRegionCode::IsRegionBath(1));
        TS_ASSERT(HeartRegionCode::IsRegionBath(2));
        TS_ASSERT(!HeartRegionCode::IsRegionTissue(1));
        TS_ASSERT(!HeartRegionCode::IsRegionTissue(2));

////////////////////////////////

        int tissue_region_one;//, tissue_region_two;
        int bath_region_one;//, bath_region_two;

        tissue_region_one = HeartRegionCode::GetValidTissueId();

        bath_region_one = HeartRegionCode::GetValidBathId();

        TS_ASSERT(HeartRegionCode::IsRegionTissue(tissue_region_one));
        TS_ASSERT(!HeartRegionCode::IsRegionBath(tissue_region_one));


        TS_ASSERT(!HeartRegionCode::IsRegionTissue(bath_region_one));
        TS_ASSERT(HeartRegionCode::IsRegionBath(bath_region_one));

    }

    void TestExceptions() throw (Exception)
    {
        std::set<unsigned> tissue_ids;
        tissue_ids.insert(0);
        tissue_ids.insert(10);

        std::set<unsigned> bath_ids;
        //Empty set
        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids),
            "Identifying set must be non-empty");


        bath_ids.insert(1);
        bath_ids.insert(10);  //Overlaps with previous

        TS_ASSERT_THROWS_THIS(HeartConfig::Instance()->SetTissueAndBathIdentifiers(tissue_ids, bath_ids),
            "Tissue identifiers and bath identifiers overlap");

     }
};

#endif /*TESTHEARTREGIONCODES_HPP_*/
