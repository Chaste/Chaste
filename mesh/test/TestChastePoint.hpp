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


#ifndef _TESTCHASTEPOINT_HPP_
#define _TESTCHASTEPOINT_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "ChastePoint.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestChastePoint : public CxxTest::TestSuite
{
public:

    /**
     * Test that values set at a coordinate are the same when accessed.
     * Also check constructors.
     * We only test dimensions from 1 to 3 inclusive.
     */
    void TestSetAndGetCoordinate()
    {
        ChastePoint<1> point1;

        double value = 12.0;
        int index = 0;
        point1.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point1[index], 1e-12);
        TS_ASSERT_DELTA(value, point1.GetWithDefault(index), 1e-12);
        TS_ASSERT_DELTA(0.0, point1.GetWithDefault(1), 1e-12);
        TS_ASSERT_DELTA(0.0, point1.GetWithDefault(2), 1e-12);

        ChastePoint<2> point2;
        point2.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point2[index], 1e-12);

        index = 1;
        value = -13.56;
        point2.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point2[index], 1e-12);
        TS_ASSERT_DELTA(0.0, point2.GetWithDefault(2), 1e-12);

        ChastePoint<3> point3;
        index = 0;
        point3.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point3[index], 1e-12);

        index = 1;
        value = 1e5;
        point3.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point3[index], 1e-12);

        index = 2;
        value = 1e-5;
        point3.SetCoordinate(index, value);
        TS_ASSERT_DELTA(value, point3[index], 1e-12);
        TS_ASSERT_DELTA(0.0, point1.GetWithDefault(3), 1e-12);

        ChastePoint<1> point4(1);
        TS_ASSERT_DELTA(point4[0], 1, 1e-12);

        ChastePoint<2> point5(2,3);
        TS_ASSERT_DELTA(point5[0], 2, 1e-12);
        TS_ASSERT_DELTA(point5[1], 3, 1e-12);

        ChastePoint<3> point6(4,5,6);
        TS_ASSERT_DELTA(point6[0], 4, 1e-12);
        TS_ASSERT_DELTA(point6[1], 5, 1e-12);
        TS_ASSERT_DELTA(point6[2], 6, 1e-12);

        ChastePoint<1> point7;
        TS_ASSERT_DELTA(point7[0], 0, 1e-12);
    }

    void TestGetLocation()
    {
        ChastePoint<3> point(1.0, 2.0, 3.0);

        c_vector<double, 3>& point_location = point.rGetLocation();
        TS_ASSERT_EQUALS(point_location(1), 2.0);
        point_location(0) = 0;

        const c_vector<double, 3>& const_point_location = point.rGetLocation();

        for (int i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(const_point_location[i], point[i], 1e-7);
        }
    }

    void TestSameChastePoints()
    {
        ChastePoint<3> point1(4,5,6);
        ChastePoint<3> point2(4,5,6);
        ChastePoint<3> point3(12,5,6);

        TS_ASSERT(point1.IsSamePoint(point2));
        TS_ASSERT(!point1.IsSamePoint(point3));
    }
    void TestZeroDimPoint()
    {
        ChastePoint<0> zero_dim_point;
        TS_ASSERT_THROWS_THIS(zero_dim_point[0], "Zero-dimensional point has no data");
    }

    void TestCreateFromCvector()
    {
        c_vector<double, 1> location;
        location[0] = 34.0;
        ChastePoint<1> point(location);
        TS_ASSERT_EQUALS(point[0], 34.0);
    }

    void TestCreateFromStdVector()
    {
        std::vector<double> location;
        location.push_back(10.0);
        location.push_back(20.0);
        location.push_back(30.0);

        ChastePoint<3> point(location);
        TS_ASSERT_EQUALS(point[0], 10.0);
        TS_ASSERT_EQUALS(point[1], 20.0);
        TS_ASSERT_EQUALS(point[2], 30.0);
    }

    void TestArchivingPoint()
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "points.arch";

        // Create and archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ChastePoint<3>* const p_point_3d = new ChastePoint<3>(-3.0, -2.0, -1.0);
            ChastePoint<2>* const p_point_2d = new ChastePoint<2>(-33.0, -22.0);
            ChastePoint<1>* const p_point_1d = new ChastePoint<1>(-185.0);

            // Should always archive a pointer
            output_arch << p_point_3d;
            output_arch << p_point_2d;
            output_arch << p_point_1d;

            delete p_point_3d;
            delete p_point_2d;
            delete p_point_1d;

        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create pointer to regions
            ChastePoint<3>* p_point_3d;
            ChastePoint<2>* p_point_2d;
            ChastePoint<1>* p_point_1d;
            input_arch >> p_point_3d;
            input_arch >> p_point_2d;
            input_arch >> p_point_1d;

            TS_ASSERT_EQUALS((*p_point_3d)[0], -3.0);
            TS_ASSERT_EQUALS((*p_point_3d)[1], -2.0);
            TS_ASSERT_EQUALS((*p_point_3d)[2], -1.0);
            TS_ASSERT_EQUALS((*p_point_2d)[0], -33.0);
            TS_ASSERT_EQUALS((*p_point_2d)[1], -22.0);
            TS_ASSERT_EQUALS((*p_point_1d)[0], -185.0);

            delete p_point_3d;
            delete p_point_2d;
            delete p_point_1d;
        }
    }
};

#endif //_TESTCHASTEPOINT_HPP_
