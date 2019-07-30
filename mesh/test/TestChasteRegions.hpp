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

#ifndef TESTCHASTEREGIONS_HPP_
#define TESTCHASTEREGIONS_HPP_

#include "CheckpointArchiveTypes.hpp"

#include <cxxtest/TestSuite.h>
#include "AbstractChasteRegion.hpp"
#include "ChasteCuboid.hpp"
#include "ChasteEllipsoid.hpp"
#include "ChasteNodesList.hpp"
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "OutputFileHandler.hpp"
#include "PetscTools.hpp"

#include "PetscSetupAndFinalize.hpp"

typedef ChasteCuboid<3> POINT_3D;

class TestRegions : public CxxTest::TestSuite
{
public:
    void TestCuboidCreationAndContained()
    {
        //3D case
        ChastePoint<3> point_a(-3, -3, -3);
        ChastePoint<3> point_b(3, 3, 3);
        ChastePoint<3> point_inside(0, 0, 0);
        ChastePoint<3> point_outside(-4, -4, -4);

        ChasteCuboid<3> cuboid_a_b(point_a, point_b);

        TS_ASSERT_THROWS_THIS(POINT_3D cuboid_b_a(point_b, point_a), "Attempt to create a cuboid with MinCorner greater than MaxCorner in some dimension");

        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_inside), true);
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_a), true); //boundary is considered in
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_b), true); //boundary is considered in
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(point_outside), false);

        // A point that is just outside the cuboid counts as inside to deal with rounding errors
        double just = 3.00000000000000008882; // taken from error in cuboid mesh generation
        ChastePoint<3> just_outside(just, just, just);
        TS_ASSERT_EQUALS(cuboid_a_b.DoesContain(just_outside), true);

        //2D case
        ChastePoint<2> point_a_2d(-3, -3);
        ChastePoint<2> point_b_2d(3, 3);
        ChastePoint<2> point_inside_2d(0, 0);
        ChastePoint<2> point_outside_2d(-4, -4);
        ChastePoint<2> point_outside_2d_one_coord(0, -4);

        ChasteCuboid<2> cuboid_a_b_2d(point_a_2d, point_b_2d);

        TS_ASSERT_EQUALS(cuboid_a_b_2d.DoesContain(point_inside_2d), true);
        TS_ASSERT_EQUALS(cuboid_a_b_2d.DoesContain(point_a_2d), true); //boundary is considered in
        TS_ASSERT_EQUALS(cuboid_a_b_2d.DoesContain(point_b_2d), true); //boundary is considered in
        TS_ASSERT_EQUALS(cuboid_a_b_2d.DoesContain(point_outside_2d), false);
        TS_ASSERT_EQUALS(cuboid_a_b_2d.DoesContain(point_outside_2d_one_coord), false); //test the algroithm when one coord is in but the other is out

        //1D case
        ChastePoint<1> point_a_1d(-3);
        ChastePoint<1> point_b_1d(3);
        ChastePoint<1> point_inside_1d(0);
        ChastePoint<1> point_outside_1d(-4);

        ChasteCuboid<1> cuboid_a_b_1d(point_a_1d, point_b_1d);

        TS_ASSERT_EQUALS(cuboid_a_b_1d.DoesContain(point_inside_1d), true);
        TS_ASSERT_EQUALS(cuboid_a_b_1d.DoesContain(point_a_1d), true); //boundary is considered in
        TS_ASSERT_EQUALS(cuboid_a_b_1d.DoesContain(point_b_1d), true); //boundary is considered in
        TS_ASSERT_EQUALS(cuboid_a_b_1d.DoesContain(point_outside_1d), false);

        // Mixed dimensional cases

        ///\ To do The following test case has been deprecated as of r7769 because the rest of the code assumes that the point are in 3D in a 3D cuboid
        //        ChastePoint<2> two_d_point_in(0.0, 0.0);
        //        ChastePoint<1> one_d_point_in(0.0);
        //
        //        ChasteCuboid<3> cuboid_3_2(point_a, point_b);
        //        ChasteCuboid<3> cuboid_3_1(point_a, point_b);
        //
        //        TS_ASSERT(cuboid_3_2.DoesContain(two_d_point_in));
        //        TS_ASSERT(cuboid_3_1.DoesContain(one_d_point_in));
        //
        //        ChastePoint<2> two_d_point_out(-4.0, -4.0);
        //        ChastePoint<1> one_d_point_out(-4.0);
        //        TS_ASSERT_EQUALS(cuboid_3_2.DoesContain(two_d_point_out), false);
        //        TS_ASSERT_EQUALS(cuboid_3_1.DoesContain(one_d_point_out), false);

        ChastePoint<3> upper = cuboid_a_b.rGetUpperCorner();
        c_vector<double, 3> diff_upper = upper.rGetLocation() - point_b.rGetLocation();
        TS_ASSERT_DELTA(norm_2(diff_upper), 0.0, 1e-10);
        ChastePoint<3> lower = cuboid_a_b.rGetLowerCorner();
        c_vector<double, 3> diff_lower = lower.rGetLocation() - point_a.rGetLocation();
        TS_ASSERT_DELTA(norm_2(diff_lower), 0.0, 1e-10);

        // Coverage
        TS_ASSERT_THROWS_NOTHING(cuboid_a_b.Destroy());
    }

    void TestNodesList()
    {
        ChastePoint<3> point_a(-3, -3, -3);
        ChastePoint<3> point_b(3, 3, 3);
        ChastePoint<3> point_c(9, 4, 7);

        Node<3> first_node(0u, point_a);
        Node<3> second_node(1u, point_b);
        Node<3> third_node(2u, point_c);

        std::vector<Node<3u>*> array_of_nodes;

        array_of_nodes.push_back(&first_node);
        array_of_nodes.push_back(&second_node);
        array_of_nodes.push_back(&third_node);

        //declare the object
        ChasteNodesList<3> nodes_list(array_of_nodes);

        ChastePoint<3> test_point_contained(9, 4, 7);
        ChastePoint<3> test_point_non_contained(10, 4, 7);

        TS_ASSERT_EQUALS(nodes_list.DoesContain(test_point_contained), true);
        TS_ASSERT_EQUALS(nodes_list.DoesContain(test_point_non_contained), false);

        TS_ASSERT_EQUALS(nodes_list.GetSize(), 3u);
    }

    void TestEllipsoidCreationAndContained()
    {
        ChastePoint<3> centre_3D(0, 0, 0);
        ChastePoint<3> radii_3D(2, 4, 6);
        ChastePoint<3> bad_radii_3D(-2, 4, 6);
        ChastePoint<3> point_inside_x_3D(1, 0, 0);
        ChastePoint<3> point_inside_y_3D(0, 2, 0);
        ChastePoint<3> point_inside_z_3D(0, 0, 3);
        ChastePoint<3> point_outside_x_3D(3, 0, 0);
        ChastePoint<3> point_outside_y_3D(0, 5, 0);
        ChastePoint<3> point_outside_z_3D(0, 0, 7);

        ChasteEllipsoid<3> ellipsoid_3D(centre_3D, radii_3D);

        TS_ASSERT_THROWS_THIS(ChasteEllipsoid<3> bad_ellipsoid_3D(centre_3D, bad_radii_3D), "Attempted to create an ellipsoid with a negative radius");

        TS_ASSERT_EQUALS(ellipsoid_3D.DoesContain(point_inside_x_3D), true);
        TS_ASSERT_EQUALS(ellipsoid_3D.DoesContain(point_inside_y_3D), true);
        TS_ASSERT_EQUALS(ellipsoid_3D.DoesContain(point_inside_z_3D), true);

        TS_ASSERT_EQUALS(ellipsoid_3D.DoesContain(point_outside_x_3D), false);
        TS_ASSERT_EQUALS(ellipsoid_3D.DoesContain(point_outside_y_3D), false);
        TS_ASSERT_EQUALS(ellipsoid_3D.DoesContain(point_outside_z_3D), false);

        // A point that is just outside the ellipsoid counts as inside to deal with rounding errors
        ChastePoint<3> just_outside_3D(0, 0, 6.00000000000000008882); // taken from error in cuboid mesh generation
        TS_ASSERT_EQUALS(ellipsoid_3D.DoesContain(just_outside_3D), true);

        ChastePoint<3> returned_centre = ellipsoid_3D.rGetCentre();
        c_vector<double, 3> diff_centre = returned_centre.rGetLocation() - centre_3D.rGetLocation();
        TS_ASSERT_DELTA(norm_2(diff_centre), 0.0, 1e-10);
        ChastePoint<3> returned_radii = ellipsoid_3D.rGetRadii();
        c_vector<double, 3> diff_radii = returned_radii.rGetLocation() - radii_3D.rGetLocation();
        TS_ASSERT_DELTA(norm_2(diff_radii), 0.0, 1e-10);

        // Lower dimensional cases
        ChastePoint<2> centre_2D(0, 0);
        ChastePoint<2> radii_2D(2, 4);
        ChastePoint<2> bad_radii_2D(-2, 4);
        ChastePoint<2> point_inside_x_2D(1, 0);
        ChastePoint<2> point_inside_y_2D(0, 2);
        ChastePoint<2> point_outside_x_2D(3, 0);
        ChastePoint<2> point_outside_y_2D(0, 5);
        ChasteEllipsoid<2> ellipsoid_2D(centre_2D, radii_2D);
        TS_ASSERT_EQUALS(ellipsoid_2D.DoesContain(point_inside_x_2D), true);
        TS_ASSERT_EQUALS(ellipsoid_2D.DoesContain(point_inside_y_2D), true);
        TS_ASSERT_EQUALS(ellipsoid_2D.DoesContain(point_outside_x_2D), false);
        TS_ASSERT_EQUALS(ellipsoid_2D.DoesContain(point_outside_y_2D), false);

        ChastePoint<1> centre_1D(0);
        ChastePoint<1> radii_1D(2);
        ChastePoint<1> bad_radii_1D(-2);
        ChastePoint<1> point_inside_x_1D(1);
        ChastePoint<1> point_outside_x_1D(3);
        ChasteEllipsoid<1> ellipsoid_1D(centre_1D, radii_1D);
        TS_ASSERT_EQUALS(ellipsoid_1D.DoesContain(point_inside_x_1D), true);
        TS_ASSERT_EQUALS(ellipsoid_1D.DoesContain(point_outside_x_1D), false);
    }

    void TestArchivingRegions()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("archive", false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "regions.arch";

        // Create and archive
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            ChastePoint<3> point_a(-3, -3, -3);
            ChastePoint<3> point_b(3, 3, 3);

            AbstractChasteRegion<3>* const p_cuboid = new ChasteCuboid<3>(point_a, point_b);

            ChastePoint<3> centre(0, 0, 0);
            ChastePoint<3> radii(2, 4, 6);

            AbstractChasteRegion<3>* const p_ellipsoid = new ChasteEllipsoid<3>(centre, radii);

            Node<3>* p_node_1 = new Node<3>(0u, false, 13, 13, 13);
            Node<3>* p_node_2 = new Node<3>(1u, false, 4, 5, 6);
            Node<3>* p_node_3 = new Node<3>(2u, false, -13, -13, -13);
            std::vector<Node<3>*> node_list;
            node_list.push_back(p_node_1);
            node_list.push_back(p_node_2);
            node_list.push_back(p_node_3);

            AbstractChasteRegion<3>* const p_node_list = new ChasteNodesList<3>(node_list);

            // Should always archive a pointer
            output_arch << p_cuboid;
            output_arch << p_ellipsoid;
            output_arch << p_node_list;

            delete p_cuboid;
            delete p_ellipsoid;
            delete p_node_list;
            delete p_node_1;
            delete p_node_2;
            delete p_node_3;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create pointer to regions
            AbstractChasteRegion<3>* p_cuboid;
            AbstractChasteRegion<3>* p_ellipsoid;
            AbstractChasteRegion<3>* p_node_list;

            input_arch >> p_cuboid;
            input_arch >> p_ellipsoid;
            input_arch >> p_node_list;

            ChastePoint<3> point_inside_cuboid(0, 0, 0);
            ChastePoint<3> point_just_inside_cuboid(-2.99, -2.99, -2.99);
            ChastePoint<3> point_outside_cuboid(-4, -4, -4);
            ChastePoint<3> point_just_outside_cuboid(-3.01, -3.01, -3.01);

            ChastePoint<3> point_inside_x_ellipsoid(1, 0, 0);
            ChastePoint<3> point_inside_y_ellipsoid(0, 2, 0);
            ChastePoint<3> point_inside_z_ellipsoid(0, 0, 3);
            ChastePoint<3> point_outside_x_ellipsoid(3, 0, 0);
            ChastePoint<3> point_outside_y_ellipsoid(0, 5, 0);
            ChastePoint<3> point_outside_z_ellipsoid(0, 0, 7);

            ChastePoint<3> node_1_in_the_list(-13, -13, -13);
            ChastePoint<3> node_2_in_the_list(13, 13, 13);
            ChastePoint<3> node_not_in_the_list(-5, 4, 0);

            TS_ASSERT_EQUALS(p_cuboid->DoesContain(point_inside_cuboid), true);
            TS_ASSERT_EQUALS(p_cuboid->DoesContain(point_outside_cuboid), false);
            TS_ASSERT_EQUALS(p_cuboid->DoesContain(point_just_inside_cuboid), true);
            TS_ASSERT_EQUALS(p_cuboid->DoesContain(point_just_outside_cuboid), false);

            TS_ASSERT_EQUALS(p_ellipsoid->DoesContain(point_inside_x_ellipsoid), true);
            TS_ASSERT_EQUALS(p_ellipsoid->DoesContain(point_inside_y_ellipsoid), true);
            TS_ASSERT_EQUALS(p_ellipsoid->DoesContain(point_inside_z_ellipsoid), true);

            TS_ASSERT_EQUALS(p_ellipsoid->DoesContain(point_outside_x_ellipsoid), false);
            TS_ASSERT_EQUALS(p_ellipsoid->DoesContain(point_outside_y_ellipsoid), false);
            TS_ASSERT_EQUALS(p_ellipsoid->DoesContain(point_outside_z_ellipsoid), false);

            TS_ASSERT_EQUALS(p_node_list->DoesContain(node_1_in_the_list), true);
            TS_ASSERT_EQUALS(p_node_list->DoesContain(node_2_in_the_list), true);
            TS_ASSERT_EQUALS(p_node_list->DoesContain(node_not_in_the_list), false);

            delete p_cuboid;
            delete p_ellipsoid;
            delete p_node_list;
        }
    }
};

#endif /*TESTCHASTEREGIONS_HPP_*/
