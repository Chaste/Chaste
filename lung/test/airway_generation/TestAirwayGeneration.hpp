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

#ifndef TESTAIRWAYGENERATION_HPP_
#define TESTAIRWAYGENERATION_HPP_

#include <cxxtest/TestSuite.h>
#include "AirwayGeneration.hpp"
#include "OutputFileHandler.hpp"
#include "TetrahedralMesh.hpp"

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

class TestAirwayGeneration : public CxxTest::TestSuite
{
public:
    void TestGeneration()
    {
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)

        EXIT_IF_PARALLEL;

        //We are only testing point assignment, so directions are arbitrary
        double current_direction[3] = {1.0, 0.0, 0.0};
        double parent_direction[3] = {0.0, 1.0, 0.0};

        //Initial apex locations
        double loc0[3] = {0.0, 0.0, 0.5};
        double loc1[3] = {1.0, 0.0, 0.5};
        double loc2[3] = {0.0, 1.0, 0.5};
        double loc3[3] = {1.0, 1.0, 0.5};

        //Check that we can distribute growth points
        {
            AirwayGeneration generation(5);

            generation.AddApex(0, loc0, current_direction, parent_direction);
            generation.AddApex(1, loc1, current_direction, parent_direction);
            generation.AddApex(2, loc2, current_direction, parent_direction);
            generation.AddApex(3, loc3, current_direction, parent_direction);

            std::set<unsigned> no_invalid;
            generation.DistributeGrowthPoints(CreatePointCube(), no_invalid);

            TS_ASSERT_EQUALS(generation.GetApices().size(), 4u);
            TS_ASSERT_EQUALS(generation.GetApices()[0].mPointCloud->GetNumberOfPoints(), 250);
            TS_ASSERT_EQUALS(generation.GetApices()[1].mPointCloud->GetNumberOfPoints(), 250);
            TS_ASSERT_EQUALS(generation.GetApices()[2].mPointCloud->GetNumberOfPoints(), 250);
            TS_ASSERT_EQUALS(generation.GetApices()[3].mPointCloud->GetNumberOfPoints(), 250);
        }

        //Check that we can limit points by radius
        {
           AirwayGeneration generation(5);
           generation.SetDistributionRadius(0.51);

           generation.AddApex(0, loc0, current_direction, parent_direction);
           generation.AddApex(1, loc1, current_direction, parent_direction);
           generation.AddApex(2, loc2, current_direction, parent_direction);
           generation.AddApex(3, loc3, current_direction, parent_direction);

           std::set<unsigned> no_invalid;
           generation.DistributeGrowthPoints(CreatePointCube(), no_invalid);

           TS_ASSERT_EQUALS(generation.GetApices()[0].mPointCloud->GetNumberOfPoints(), 134);
           TS_ASSERT_EQUALS(generation.GetApices()[1].mPointCloud->GetNumberOfPoints(), 134);
           TS_ASSERT_EQUALS(generation.GetApices()[2].mPointCloud->GetNumberOfPoints(), 134);
           TS_ASSERT_EQUALS(generation.GetApices()[3].mPointCloud->GetNumberOfPoints(), 134);
        }

        //Check that growth points can be invalidated
        {
            AirwayGeneration generation(5);

            generation.AddApex(0, loc0, current_direction, parent_direction);
            generation.AddApex(1, loc1, current_direction, parent_direction);
            generation.AddApex(2, loc2, current_direction, parent_direction);
            generation.AddApex(3, loc3, current_direction, parent_direction);

            std::set<unsigned> invalid;
            invalid.insert(0);
            invalid.insert(1);
            invalid.insert(2);
            invalid.insert(3);
            invalid.insert(4);
            invalid.insert(997);
            invalid.insert(998);
            invalid.insert(999);
            generation.DistributeGrowthPoints(CreatePointCube(), invalid);

            TS_ASSERT_EQUALS(generation.GetApices().size(), 4u);
            TS_ASSERT_EQUALS(generation.GetApices()[0].mPointCloud->GetNumberOfPoints(), 245);
            TS_ASSERT_EQUALS(generation.GetApices()[1].mPointCloud->GetNumberOfPoints(), 250);
            TS_ASSERT_EQUALS(generation.GetApices()[2].mPointCloud->GetNumberOfPoints(), 250);
            TS_ASSERT_EQUALS(generation.GetApices()[3].mPointCloud->GetNumberOfPoints(), 247);
        }
#endif
    }

    void TestDummyClassCoverage()
    {
#if !(defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6))
       EXIT_IF_PARALLEL;

       AirwayGeneration generation;

#endif
    }

private:
#if defined(CHASTE_VTK) && ( (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 6) || VTK_MAJOR_VERSION >= 6)
    vtkSmartPointer<vtkPolyData> CreatePointCube()
    {
        vtkSmartPointer<vtkPolyData> cube = vtkSmartPointer<vtkPolyData>::New();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        cube->SetPoints(points);

        double side_length = 1.0;
        unsigned points_per_dimension = 10;
        unsigned points_per_dimension_minus_one = points_per_dimension - 1;

        for (unsigned i = 0; i <= points_per_dimension_minus_one; ++i)
        {
            for (unsigned j = 0; j <= points_per_dimension_minus_one; ++j)
            {
                for (unsigned k = 0; k <= points_per_dimension_minus_one; ++k)
                {
                    double point[3];
                    point[0] = (static_cast<double>(i)/points_per_dimension_minus_one)*side_length;
                    point[1] = (static_cast<double>(j)/points_per_dimension_minus_one)*side_length;
                    point[2] = (static_cast<double>(k)/points_per_dimension_minus_one)*side_length;

                    points->InsertNextPoint(point);
                }
            }
        }

        assert((unsigned)cube->GetNumberOfPoints() == points_per_dimension*points_per_dimension*points_per_dimension);

        return cube;
    }
#endif
};

#endif /* TESTAIRWAYGENERATOR_HPP_ */
