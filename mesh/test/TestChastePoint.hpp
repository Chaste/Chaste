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


#ifndef _TESTCHASTEPOINT_HPP_
#define _TESTCHASTEPOINT_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "ChastePoint.hpp"

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
