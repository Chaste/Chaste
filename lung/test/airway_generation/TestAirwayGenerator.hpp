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

#ifndef TESTAIRWAYGENERATOR_HPP_
#define TESTAIRWAYGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "AirwayGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "TetrahedralMesh.hpp"
#include "FileFinder.hpp"
#include "UblasIncludes.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>

#ifdef CHASTE_VTK

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include "vtkVersion.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkSphereSource.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkSTLReader.h"

#endif //CHASTE_VTK

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestAirwayGenerator : public CxxTest::TestSuite
{
public:
    void TestCreatePointCloud()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere(50);

        AirwayGenerator generator(sphere);

        double test_point[3];
        test_point[0] = 1.0;
        test_point[1] = 1.0;
        test_point[2] = 0.0;
        TS_ASSERT(!generator.IsInsideLobeSurface(test_point));
        TS_ASSERT_DELTA(generator.DistanceFromLobeSurface(test_point), sqrt(2.0) - 1, 5e-3);

        vtkSmartPointer<vtkPolyData> point_data = generator.CreatePointCloud(std::pow(4*M_PI/3/50, 1.0/3.0)); //Gives spacing for ~50 points

        //Check that all points created are inside the sphere
        std::cout << point_data->GetNumberOfPoints() << std::endl;
        TS_ASSERT_EQUALS(point_data->GetNumberOfPoints(), 51);
        for (int i = 0; i < point_data->GetNumberOfPoints(); ++i)
        {
            double coords[3];
            point_data->GetPoint(i, coords);

            TS_ASSERT_LESS_THAN(std::sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]), 1.0);
        }

        //Check that the centre of mass is close to zero
        double centre[3];
        generator.GetCentreOfMass(point_data, centre);

        TS_ASSERT_DELTA(centre[0], 0.0, 5e-2);
        TS_ASSERT_DELTA(centre[1], 0.0, 5e-2);
        TS_ASSERT_DELTA(centre[2], 0.0, 5e-2);

        std::set<unsigned>& invalid_ids = generator.GetInvalidIds();
        TS_ASSERT_EQUALS(invalid_ids.size(), 0u);
#endif
    }

    void TestCreatePointCloudOtherMethods()
            {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

        EXIT_IF_PARALLEL;

        {
            vtkSmartPointer<vtkPolyData> sphere = CreateSphere(50);

            AirwayGenerator generator(sphere);
            vtkSmartPointer<vtkPolyData> point_data = generator.CreatePointCloudUsingTargetPoints(50);

            TS_ASSERT_EQUALS(point_data->GetNumberOfPoints(), 51);

        }

        {
            vtkSmartPointer<vtkPolyData> sphere = CreateSphere(50);

            AirwayGenerator generator(sphere);
            vtkSmartPointer<vtkPolyData> point_data = generator.CreatePointCloudUsingTargetVolume(4*M_PI/3/50);

            TS_ASSERT_EQUALS(point_data->GetNumberOfPoints(), 51);

        }
#endif
            }

    void TestSplitPointCloud()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();

        AirwayGenerator generator(sphere);
        vtkSmartPointer<vtkPolyData> point_cloud = generator.CreatePointCloudUsingTargetPoints(50);

        double normal[3] = {0.0, 1.0, 0.0 };
        double origin[3] = {0.0, 0.0, 0.0};
        vtkSmartPointer<vtkPolyData> top_half_sphere = generator.SplitPointCloud(point_cloud,
                                                                                 normal,
                                                                                 origin,
                                                                                 false);

        TS_ASSERT_EQUALS(top_half_sphere->GetNumberOfPoints(), 21);

        for (int i = 0; i < top_half_sphere->GetNumberOfPoints(); ++i)
        {
            double coords[3];
            top_half_sphere->GetPoint(i, coords);

            TS_ASSERT_LESS_THAN(0.0, coords[1]);
        }

        top_half_sphere = generator.SplitPointCloud(point_cloud,
                                                    normal,
                                                    origin,
                                                    true);

        TS_ASSERT_EQUALS(top_half_sphere->GetNumberOfPoints(), 28);

        for (int i = 0; i < top_half_sphere->GetNumberOfPoints(); ++i)
        {
            double coords[3];
            top_half_sphere->GetPoint(i, coords);

            TS_ASSERT_LESS_THAN(coords[1], 0.0);
        }

#endif
    }

    void TestAddInitialApex()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();
        AirwayGenerator generator(sphere);

        double origin[3] = {0.0, 1.0, 0.0};
        double direction[3] = {0.0, 1.0, 0.0};
        double parent_direction[3] = {1.0, 0.0, 0.0};

        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 0);
        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 1);
        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 1);
        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 3);

        TS_ASSERT_EQUALS(generator.GetGenerations().size(), 30u);
        TS_ASSERT_EQUALS(generator.GetGenerations()[0].GetApices().size(), 1u);
        TS_ASSERT_EQUALS(generator.GetGenerations()[1].GetApices().size(), 2u);
        TS_ASSERT_EQUALS(generator.GetGenerations()[2].GetApices().size(), 0u);
        TS_ASSERT_EQUALS(generator.GetGenerations()[3].GetApices().size(), 1u);

        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 4);

        //The generator only supports generations up to 30
        TS_ASSERT_THROWS_THIS(generator.AddInitialApex(origin, direction, parent_direction, 10.0, 35),
                              "Error: Airway generation can only generate up to 30 generations.");
#endif
    }


    void TestInsertBranch()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere(50);

        double min_branch_length = 1.2;
        unsigned point_limit = 5;
        double angle_limit = 90.0;
        double branching_fraction = 0.5;

        AirwayGenerator generator(sphere, min_branch_length, point_limit, angle_limit, branching_fraction);


        double origin[3] = {0.0, 1.0, 0.0};
        double direction[3] = {0.0, -1.0, 0.0};
        double parent_direction[3] = {1.0, 0.0, 0.0};

        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 0); //creates a single node for the start of the apex

        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 1);

        //Create a single branch, should end up 50% of the way towards the centre
        double end_location[3];
        vtkSmartPointer<vtkPolyData> p_point_cloud = generator.CreatePointCloudUsingTargetPoints(100);
        unsigned branch_end_id = generator.InsertBranch(p_point_cloud, 0, direction, end_location);

        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 2);
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfCells(), 1);
        TS_ASSERT_EQUALS(branch_end_id, 1u);

        TS_ASSERT_DELTA(end_location[0], 0.0, 1e-2);
        TS_ASSERT_DELTA(end_location[1], 0.5, 1e-2);
        TS_ASSERT_DELTA(end_location[2], 0.0, 1e-2);

        double inserted_end_location[3];
        generator.GetAirwayTree()->GetPoints()->GetPoint(1, inserted_end_location);

        TS_ASSERT_DELTA(end_location[0], inserted_end_location[0], 1e-8);
        TS_ASSERT_DELTA(end_location[1], inserted_end_location[1], 1e-8);
        TS_ASSERT_DELTA(end_location[2], inserted_end_location[2], 1e-8);
#endif
    }

    void TestGrowApex()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere(100);

        double min_branch_length = 0.01;
        unsigned point_limit = 5;
        double angle_limit = 180.0;
        double branching_fraction = 1.0;

        AirwayGenerator generator(sphere, min_branch_length, point_limit, angle_limit, branching_fraction);

        double origin[3] = {0.0, 1.0001, 0.0}; //to avoid bias in the splitting plane
        double direction[3] = {1.0, 0.0, 0.0};
        double parent_direction[3] = {1.0, 0.0, 0.0};

        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 0);
        Apex& apex = generator.GetGenerations()[0].GetApices()[0];
        apex.mPointCloud = generator.CreatePointCloudUsingTargetPoints(5000);

        generator.GrowApex(apex);

        //Test that two branches were generated
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 3);
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfCells(), 2);

        //Test the exact generated locations
        double inserted_end_location1[3];
        generator.GetAirwayTree()->GetPoints()->GetPoint(1, inserted_end_location1);

        double inserted_end_location2[3];
        generator.GetAirwayTree()->GetPoints()->GetPoint(2, inserted_end_location2);

        //Generated points should have x ~= 0, y ~= 0.0,  z ~= +/- 3/8
        TS_ASSERT_DELTA(inserted_end_location1[0], 0.0, 1e-2);
        TS_ASSERT_DELTA(inserted_end_location1[1], 0.0, 1e-2);
        TS_ASSERT_DELTA(inserted_end_location1[2], -3.0/8.0, 1e-2);
        TS_ASSERT_DELTA(inserted_end_location2[0], 0.0, 1e-2);
        TS_ASSERT_DELTA(inserted_end_location2[1], 0.0, 1e-2);
        TS_ASSERT_DELTA(inserted_end_location2[2], 3.0/8.0, 1e-2);

        //Test that two child apices were created in generation 1
        TS_ASSERT_EQUALS(generator.GetGenerations()[1].GetApices().size(), 2u);

        //Check the start ids, parent directions and current directions of the two apices
        Apex& new_apex_1 = generator.GetGenerations()[1].GetApices()[0];
        Apex& new_apex_2 = generator.GetGenerations()[1].GetApices()[1];

        TS_ASSERT_EQUALS(new_apex_1.mStartId, 1);
        TS_ASSERT_EQUALS(new_apex_2.mStartId, 2);

        double branch_length = std::sqrt(9.0/64.0 + 1.0);

        TS_ASSERT_DELTA(new_apex_1.mOriginalDirection[0], 0.0, 1e-2);
        TS_ASSERT_DELTA(new_apex_1.mOriginalDirection[1], -1/branch_length, 1e-2);
        TS_ASSERT_DELTA(new_apex_1.mOriginalDirection[2], -3.0/8.0/branch_length, 1e-2);
        TS_ASSERT_DELTA(new_apex_2.mOriginalDirection[0], 0.0, 1e-2);
        TS_ASSERT_DELTA(new_apex_2.mOriginalDirection[1], -1/branch_length, 1e-2);
        TS_ASSERT_DELTA(new_apex_2.mOriginalDirection[2], 3.0/8.0/branch_length, 1e-2);

#endif
    }

    void TestInvalidateClosestPoint()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere(50);

        AirwayGenerator generator(sphere);
        vtkSmartPointer<vtkPolyData> cloud = generator.CreatePointCloudUsingTargetPoints(100);

        double invalid1[3] = {-1.0, -1.0, -0.3}; //Expect this to be the first point
        double invalid2[3] = {1.0, 1.0, 0.7}; //Expect this to be the last point
        double invalid3[3] = {0.0, 0.0, 0.0}; //Expect this to be roughly the middle inserted point

        generator.InvalidateClosestPoint(invalid1, cloud);
        generator.InvalidateClosestPoint(invalid2, cloud);
        generator.InvalidateClosestPoint(invalid3, cloud);

        TS_ASSERT_EQUALS(cloud->GetNumberOfPoints(), 93);

        //Assert that the correct points were invalidated and that no others were
        std::set<unsigned>& invalid_ids = generator.GetInvalidIds();
        TS_ASSERT_EQUALS(invalid_ids.size(), 3u);

        TS_ASSERT(invalid_ids.count(0));
        TS_ASSERT(invalid_ids.count(92));
        TS_ASSERT(invalid_ids.count(49));

#endif
    }

    void TestGrowTerminalLengthApex()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();

        AirwayGenerator generator(sphere, 2.0); //Long limit stops growth after one generation

        double origin[3] = {0.0, 1.0001, 0.0}; //to avoid bias in the splitting plane
        double direction[3] = {1.0, 0.0, 0.0};
        double parent_direction[3] = {1.0, 0.0, 0.0};

        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 0);
        Apex& apex = generator.GetGenerations()[0].GetApices()[0];
        apex.mPointCloud = generator.CreatePointCloudUsingTargetPoints(50);

        generator.GrowApex(apex);

        //Test that two branches were generated
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 3);
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfCells(), 2);

        //Test that no apices were added and that points were invalidated
        TS_ASSERT_EQUALS(generator.GetGenerations()[1].GetApices().size(), 0u);
        TS_ASSERT_EQUALS(generator.GetInvalidIds().size(), 2u);
#endif
    }
//
//    void xxxTestGrowTerminalPointsApex()
//    {
//#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
//        EXIT_IF_PARALLEL;
//
//        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();
//
//        AirwayGenerator generator(sphere, 0.1, 51);
//
//        double origin[3] = {0.0, 1.0, 0.0};
//        double direction[3] = {0.0, -1.0, 0.0};
//
//        generator.AddInitialApex(origin, direction, 10.0, 1, generator.CreatePointCloud(50));
//        generator.GrowCurrentApex();
//
//        //Test that the branch was generated
//        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 2);
//        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfCells(), 1);
//
//        //Test that no child apices were created
//        TS_ASSERT_EQUALS(generator.GetApices().size(), 0u);
//#endif
//    }
//
//
//
    void TestCheckAngleAndLength()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();

        //Check the angle code works for right angles
        {
            AirwayGenerator generator(sphere, 0.1, 5, 45, 0.5);

            double origin[3] = {0.0, 1.0, 0.0};
            double direction[3] = {0.0, -1.0, 0.0};

            double centre[3] = {1.0, 1.0, 0.0}; //Requires rotating and scaling

            generator.AddInitialApex(origin, direction, direction, 10.0, 1);

            generator.CheckBranchAngleLengthAndAdjust(0, direction, centre);

            TS_ASSERT_DELTA(centre[0], 0.5*1.0/std::sqrt(2.0), 1e-6);
            TS_ASSERT_DELTA(centre[1], 1.0 - 0.5*1/std::sqrt(2.0), 1e-6);
            TS_ASSERT_DELTA(centre[2], 0.0, 1e-6);
        }

        //Check the angle code works for acute angles
        {
            AirwayGenerator generator(sphere, 0.1, 5, 0.0, 1.0); //No branching angle allowed

            double origin[3] = {0.0, 1.0, 0.0};
            double direction[3] = {0.0, -1.0, 0.0};

            double centre[3] = {1/std::sqrt(2.0), 1.0-1/std::sqrt(2.0), 0.0};

            generator.AddInitialApex(origin, direction, direction, 10.0, 1);

            generator.CheckBranchAngleLengthAndAdjust(0, direction, centre);

            TS_ASSERT_DELTA(centre[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(centre[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(centre[2], 0.0, 1e-6);
        }

        //check the code works for obtuse angles
        {
            AirwayGenerator generator(sphere, 0.1, 5, 90, 1.0);

            double origin[3] = {0.0, 1.0, 0.0};
            double direction[3] = {0.0, -1.0, 0.0};

            double centre[3] = {1.0, 2.0, 0.0}; //results in a vector of length sqrt(2.0)

            generator.AddInitialApex(origin, direction, direction, 10.0, 1);

            generator.CheckBranchAngleLengthAndAdjust(0, direction, centre);

            TS_ASSERT_DELTA(centre[0], std::sqrt(2.0), 1e-6);
            TS_ASSERT_DELTA(centre[1], 1.0, 1e-6);
            TS_ASSERT_DELTA(centre[2], 0.0, 1e-6);
        }
#endif
    }

    void TestOutsideHostVolume()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();

        AirwayGenerator generator(sphere, 2, 5, 0.5);

        double origin[3] = {0.0, 1.0, 0.0};
        double direction[3] = {1.0, 0.0, 0.0};
        generator.CreatePointCloudUsingTargetPoints(50);

        generator.AddInitialApex(origin, direction, direction, 10.0, 1);

        //This should fail the host volume test and not generate
        generator.Generate();

        //Test that the branch was not generated
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 1);
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfCells(), 0);
#endif
    }

    void TestGenerate()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;
        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();

        double min_branch_length = 0.1;
        unsigned point_limit = 1;
        double angle_limit = 180.0;
        double branching_fraction = 0.4;

        AirwayGenerator generator(sphere, min_branch_length, point_limit, angle_limit, branching_fraction);

        double origin[3] = {0.0, 1.0, 0.0};
        double direction[3] = {1.0, 0.0, 0.0};
        double parent_direction[3] = {1.0, 0.0, 0.0};

        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 0);
        generator.CreatePointCloudUsingTargetPoints(50);

        generator.Generate();

        //Test that the tree was generated
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfPoints(), 97);
        TS_ASSERT_EQUALS(generator.GetAirwayTree()->GetNumberOfCells(), 96);

        //to visualise
        /*OutputFileHandler output("TestAirwayGenerator");
        std::string output_file = output.GetOutputDirectoryFullPath() + "/test.vtp";

        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(output_file.c_str());
#if VTK_MAJOR_VERSION >= 6
        writer->SetInputData(generator.GetAirwayTree());
#else
        writer->SetInput(generator.GetAirwayTree());
#endif
        writer->Write();*/
#endif
    }

    void TestHorsfieldOrder()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;

        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();

        AirwayGenerator generator(sphere, 0.1, 10);
        generator.CreatePointCloudUsingTargetPoints(50);

        double origin[3] = {0.0, 1.0, 0.0};
        double direction[3] = {0.0, -1.0, 0.0};

        generator.AddInitialApex(origin, direction, direction, 10.0, 1);
        generator.Generate();

        generator.CalculateHorsfieldOrder();

        //Test resulting orders, these have been calculated by hand for this tree
        vtkSmartPointer<vtkDataArray> order = generator.GetAirwayTree()->GetPointData()->GetArray("horsfield_order");
        TS_ASSERT_DELTA(order->GetTuple1(0), 4, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(1), 3, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(2), 3, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(3), 2, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(4), 1, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(5), 2, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(6), 2, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(7), 1, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(8), 1, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(9), 1, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(10), 1, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(11), 1, 1e-6);
        TS_ASSERT_DELTA(order->GetTuple1(12), 1, 1e-6);

        generator.CalculateRadii(1.6);

        vtkSmartPointer<vtkDataArray> radii = generator.GetAirwayTree()->GetPointData()->GetArray("radius");

        double log_max_diameter = std::log10(20.0);
        double max_order = 4.0;
        double log_diameter_ratio = std::log10(1.6);

        for (int i = 0; i < radii->GetSize(); ++i)
        {
            TS_ASSERT_DELTA(radii->GetTuple1(i), std::pow(10.0, log_diameter_ratio*(order->GetTuple1(i) - max_order) + log_max_diameter)/2.0, 1e-4);
        }

        generator.MarkStartIds();
        vtkSmartPointer<vtkDataArray> start_ids = generator.GetAirwayTree()->GetPointData()->GetArray("start_id");

        TS_ASSERT_DELTA(start_ids->GetTuple1(0), 1.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(1), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(2), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(3), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(4), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(5), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(6), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(7), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(8), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(9), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(10), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(11), 0.0, 1e-4);
        TS_ASSERT_DELTA(start_ids->GetTuple1(12), 0.0, 1e-4);

        //To visualise
        /*OutputFileHandler output("TestAirwayGenerator");
        std::string output_file = output.GetOutputDirectoryFullPath() + "/test.vtp";

        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer->SetFileName(output_file.c_str());
#if VTK_MAJOR_VERSION >= 6
        writer->SetInputData(generator.GetAirwayTree());
#else
        writer->SetInput(generator.GetAirwayTree());
#endif
        writer->Write();*/
#endif
    }


    void TestGenerateDecomposedAirways()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;
        vtkSmartPointer<vtkPolyData> sphere = CreateSphere();

        double min_branch_length = 0.1;
        unsigned point_limit = 1;
        double angle_limit = 180.0;
        double branching_fraction = 0.4;

        AirwayGenerator generator(sphere, min_branch_length, point_limit, angle_limit, branching_fraction);

        double origin1[3] = {0.0, 1.0, 0.0};
        double origin2[3] = {1.0, 0.0, 0.0};
        double origin3[3] = {0.0, 0.0, 1.0};
        double direction[3] = {1.0, 0.0, 0.0};
        double parent_direction[3] = {1.0, 0.0, 0.0};

        generator.AddInitialApex(origin1, direction, parent_direction, 10.0, 0);
        generator.AddInitialApex(origin2, direction, parent_direction, 10.0, 0);
        generator.AddInitialApex(origin3, direction, parent_direction, 10.0, 0);
        generator.CreatePointCloudUsingTargetPoints(50);

        generator.Generate();

        generator.CalculateHorsfieldOrder();
        generator.CalculateRadii(1.15);

        generator.WriteDecomposedAirways("TestAirwayGenerator", "decomposed_test");

        //We expect three separate files to be written.
        FileFinder finder_0("TestAirwayGenerator/decomposed_test_0.vtu", RelativeTo::ChasteTestOutput);
        FileFinder finder_1("TestAirwayGenerator/decomposed_test_1.vtu", RelativeTo::ChasteTestOutput);
        FileFinder finder_2("TestAirwayGenerator/decomposed_test_2.vtu", RelativeTo::ChasteTestOutput);
        FileFinder finder_3("TestAirwayGenerator/decomposed_test_0.node", RelativeTo::ChasteTestOutput);
        FileFinder finder_4("TestAirwayGenerator/decomposed_test_1.node", RelativeTo::ChasteTestOutput);
        FileFinder finder_5("TestAirwayGenerator/decomposed_test_2.node", RelativeTo::ChasteTestOutput);
        FileFinder finder_6("TestAirwayGenerator/decomposed_test_0.edge", RelativeTo::ChasteTestOutput);
        FileFinder finder_7("TestAirwayGenerator/decomposed_test_1.edge", RelativeTo::ChasteTestOutput);
        FileFinder finder_8("TestAirwayGenerator/decomposed_test_2.edge", RelativeTo::ChasteTestOutput);

        TS_ASSERT(finder_0.Exists() && !finder_0.IsEmpty());
        TS_ASSERT(finder_1.Exists() && !finder_1.IsEmpty());
        TS_ASSERT(finder_2.Exists() && !finder_2.IsEmpty());
        TS_ASSERT(finder_3.Exists() && !finder_3.IsEmpty());
        TS_ASSERT(finder_4.Exists() && !finder_4.IsEmpty());
        TS_ASSERT(finder_5.Exists() && !finder_5.IsEmpty());
        TS_ASSERT(finder_6.Exists() && !finder_6.IsEmpty());
        TS_ASSERT(finder_7.Exists() && !finder_7.IsEmpty());
        TS_ASSERT(finder_8.Exists() && !finder_8.IsEmpty());
#endif
    }


    void TestCalculateLobeVolume()
    {
    #if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
       EXIT_IF_PARALLEL;

       vtkSmartPointer<vtkPolyData> sphere = CreateSphere(100);

       AirwayGenerator generator(sphere);

       //The sphere is coarsely meshed, hence relatively large tolerance
       TS_ASSERT_DELTA(generator.CalculateLobeVolume(), 4.0/3.0*M_PI, 1e-2);

    #endif
    }

    void TestEndBranchDistanceLimit()
    {
    #if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
        EXIT_IF_PARALLEL;
        vtkSmartPointer<vtkPolyData> sphere = CreateSphere(100, 100.0); //lung sized sphere to test distance limit

        double min_branch_length = 0.1;
        unsigned point_limit = 1;
        double angle_limit = 180.0;
        double branching_fraction = 0.4;
        bool applyPointReassignmentLimit = true;

        AirwayGenerator generator(sphere, min_branch_length, point_limit, angle_limit, branching_fraction, applyPointReassignmentLimit);

        std::deque<AirwayGeneration>& generations = generator.GetGenerations();

        //The bounding box of the sphere is (-100,-100,-100) and (100,100,100), the diagonal of the bounding box is of size
        //sqrt(dx^2 + dy^2 + dz^2) where dx = dy = dz = 200.0
        double lobe_bound_size = sqrt(120000);

        for (unsigned generation_number = 0; generation_number < generations.size(); ++generation_number)
        {
            double scale_distance_limit = lobe_bound_size/30; //30 is the expected total number of generations
            TS_ASSERT_DELTA(generations[generation_number].GetDistributionRadius(), std::max(lobe_bound_size - scale_distance_limit*generation_number, 5.0), 1e-1);
        }

        double origin[3] = {0.0, 1.0, 0.0};
        double direction[3] = {1.0, 0.0, 0.0};
        double parent_direction[3] = {1.0, 0.0, 0.0};

        //Check that all points are distributed to an apex in generation 0
        generator.AddInitialApex(origin, direction, parent_direction, 10.0, 0);

        vtkSmartPointer<vtkPolyData> cloud = generator.CreatePointCloudUsingTargetPoints(1000);
        std::set<unsigned> invalid_ids;
        generations[0].DistributeGrowthPoints(cloud, invalid_ids);
        unsigned assigned_points = generations[0].GetApices()[0].mPointCloud->GetNumberOfPoints();
        unsigned total_points = cloud->GetNumberOfPoints();
        TS_ASSERT_EQUALS(assigned_points, total_points);

        //Check that the number of distributed points becomes less with each generation
        unsigned previous_assigned_points = assigned_points;
        for (unsigned generation_number = 1; generation_number < generations.size(); ++generation_number)
        {
            generator.AddInitialApex(origin, direction, parent_direction, 10.0, generation_number);

            generations[generation_number].DistributeGrowthPoints(cloud, invalid_ids);
            assigned_points = generations[generation_number].GetApices()[0].mPointCloud->GetNumberOfPoints();
            TS_ASSERT_LESS_THAN_EQUALS(assigned_points, previous_assigned_points);
            previous_assigned_points = assigned_points;
        }
    #endif
    }

    void TestDummyClassCoverage()
    {
#if !(defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6))
       EXIT_IF_PARALLEL;

       AirwayGenerator generator;

    #endif
    }

private:
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
    vtkSmartPointer<vtkPolyData> CreateSphere(unsigned resolution = 18, double radius = 1.0)
    {
        vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetRadius(radius);
        sphere->SetThetaResolution(resolution);
        sphere->SetPhiResolution(resolution);
        sphere->Update();

        vtkSmartPointer<vtkPolyData> sphere_data = sphere->GetOutput();
        return sphere_data;
    }
#endif
};

#endif /* TESTAIRWAYGENERATOR_HPP_ */
