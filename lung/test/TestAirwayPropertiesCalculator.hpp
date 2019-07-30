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

#ifndef _TESTAIRWAYPROPERTIESCALCULATOR_HPP_
#define _TESTAIRWAYPROPERTIESCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <queue>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "AirwayPropertiesCalculator.hpp"

class TestAirwayPropertiesCalculator : public CxxTest::TestSuite
{
public:


    void TestBranchProperties()
    {
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("lung/test/data/TestSubject002");
        mesh.ConstructFromMeshReader(mesh_reader);

        AirwayPropertiesCalculator properties_calculator(mesh, 0u);

        TS_ASSERT_EQUALS(properties_calculator.GetBranches().size(), 123893u);

        properties_calculator.CalculateBranchProperties();
        TS_ASSERT_DELTA(properties_calculator.GetLengthOneOverLengthTwoMean(), 0.65, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetLengthOverDiameterMajorChildMean(), 10.097, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetLengthOverDiameterMean(), 10.364, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetLengthOverDiameterMinorChildMean(), 10.4807, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetLengthOverLengthParentMean(), 0.877, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetDiameterOverParentDiameterMean(), 0.8145, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetMajorDiameterOverParentDiameterMean(), 0.842716, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetMinorDiameterOverMajorDiameterMean(), 0.959747, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetMinorDiameterOverParentDiameterMean(), 0.802192, 1e-3);
        TS_ASSERT_DELTA(properties_calculator.GetPercentageLengthOverParentLengthLessThanOne(), 0.698899, 1e-3);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetPhiMean(), 90, 1e-1);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetThetaMajorBranches(), 39.3648, 1e-3);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetThetaMean(), 43.2006, 1e-3);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetThetaMinorBranches(), 44.8781, 1e-3);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetThetaParentDiameter2mmTo1mm(), 35.9925, 1e-3);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetThetaParentDiameter3mmTo2mm(), 30.4503, 1e-3);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetThetaParentDiameter4mmTo3mm(), 25.2783, 1e-3);
        TS_ASSERT_DELTA(180/M_PI*properties_calculator.GetThetaParentDiameterGreaterThan4mm(), 29.4599, 1e-3);

        // Test that an index has been assigned to each branch sequentially
        std::vector<AirwayBranch*> branches = properties_calculator.GetBranches();
        TS_ASSERT_EQUALS(branches.size(), 123893u);
        for (unsigned branch_idx = 0 ; branch_idx < branches.size() ; branch_idx++)
        {
            TS_ASSERT_EQUALS(branch_idx, branches[branch_idx]->GetIndex());
        }

        TS_ASSERT_EQUALS(properties_calculator.GetMaximumTerminalGeneration(), 29u);
        TS_ASSERT_EQUALS(properties_calculator.GetMinimumTerminalGeneration(), 10u);
        TS_ASSERT_EQUALS(properties_calculator.GetMeanTerminalGeneration(), 17u);
    }


    void TestOrders()
    {
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/three_generation_branch_mesh_refined");
        mesh.ConstructFromMeshReader(mesh_reader);

        AirwayPropertiesCalculator properties_calculator(mesh, 0u);
        AirwayTreeWalker walker(mesh, 0u);

        std::vector<AirwayBranch*> branches = properties_calculator.GetBranches();

        TS_ASSERT_EQUALS(branches.size(), 7u);

        TS_ASSERT_EQUALS(properties_calculator.GetBranchStrahlerOrder(branches[0]), walker.GetMaxElementStrahlerOrder());
        TS_ASSERT_EQUALS(properties_calculator.GetBranchStrahlerOrder(branches[1]), 2u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchStrahlerOrder(branches[2]), 2u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchStrahlerOrder(branches[3]), 1u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchStrahlerOrder(branches[4]), 1u);

        TS_ASSERT_EQUALS(properties_calculator.GetBranchHorsfieldOrder(branches[0]), walker.GetMaxElementHorsfieldOrder());
        TS_ASSERT_EQUALS(properties_calculator.GetBranchHorsfieldOrder(branches[1]), 2u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchHorsfieldOrder(branches[2]), 2u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchHorsfieldOrder(branches[3]), 1u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchHorsfieldOrder(branches[4]), 1u);

        TS_ASSERT_EQUALS(properties_calculator.GetBranchGeneration(branches[0]), 0u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchGeneration(branches[1]), 1u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchGeneration(branches[2]), 1u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchGeneration(branches[3]), 2u);
        TS_ASSERT_EQUALS(properties_calculator.GetBranchGeneration(branches[4]), 2u);
    }

    void TestSubtreeProperties()
    {
        // Load test mesh
        TetrahedralMesh<1, 3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("lung/test/data/TestSubtreeProperties");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check num nodes and elements are as expected
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 9u);

        AirwayPropertiesCalculator properties_calculator(mesh, 0u);
        properties_calculator.CalculateSubtreeProperties();

        std::vector<double> branch_lengths = properties_calculator.GetSubtreeBranchLengths();
        std::vector<double> branch_volumes = properties_calculator.GetSubtreeBranchVolumes();
        std::vector<double> branch_areas = properties_calculator.GetSubtreeBranchLateralSurfaceAreas();
        std::vector<double> branch_resistances = properties_calculator.GetSubtreePoiseuilleResistances();
        std::vector<c_vector<double, 3> > branch_centroids = properties_calculator.GetSubtreeCentroids();

        // Check correct number of branches found
        TS_ASSERT_EQUALS(properties_calculator.GetBranches().size(), 3u);

        // Check length properties
        TS_ASSERT_DELTA(branch_lengths[0], 6.0 + 2 * sqrt(72), 1e-6);
        TS_ASSERT_DELTA(branch_lengths[1], branch_lengths[2], 1e-6);
        TS_ASSERT_DELTA(branch_lengths[1], sqrt(72), 1e-6);

        // Check volume properties
        TS_ASSERT_DELTA(branch_volumes[0], (6.0 + 2 * sqrt(72)) * M_PI, 1e-6);
        TS_ASSERT_DELTA(branch_volumes[1], branch_volumes[2], 1e-6);
        TS_ASSERT_DELTA(branch_volumes[1], sqrt(72)*M_PI, 1e-6);

        // Check area properties
        TS_ASSERT_DELTA(branch_areas[0], (6.0 + 2 * sqrt(72)) * 2 * M_PI, 1e-6);
        TS_ASSERT_DELTA(branch_areas[1], branch_areas[2], 1e-6);
        TS_ASSERT_DELTA(branch_areas[1], sqrt(72) * 2 * M_PI, 1e-6);

        // Check Poiseuille resistance properties
        TS_ASSERT_DELTA(branch_resistances[0], 6.0 + 3 * sqrt(2), 1e-6);
        TS_ASSERT_DELTA(branch_resistances[1], branch_lengths[1], 1e-6);
        TS_ASSERT_DELTA(branch_resistances[2], branch_lengths[2], 1e-6);

        // Check centroid properties
        TS_ASSERT_DELTA(branch_centroids[0][0], (2 * sqrt(72) - 27)/7, 1e-6);
        TS_ASSERT_DELTA(branch_centroids[0][1], 0.0, 1e-6);
        TS_ASSERT_DELTA(branch_centroids[0][2], 0.0, 1e-6);
        TS_ASSERT_DELTA(branch_centroids[1][0], branch_centroids[2][0], 1e-6);
        TS_ASSERT_DELTA(branch_centroids[1][1], -branch_centroids[2][1], 1e-6);
        TS_ASSERT_DELTA(branch_centroids[1][2], branch_centroids[2][2], 1e-6);
        TS_ASSERT_DELTA(branch_centroids[1][0], -3.0, 1e-6);
        TS_ASSERT_DELTA(branch_centroids[1][1], 3.0, 1e-6);
        TS_ASSERT_DELTA(branch_centroids[1][2], 0.0, 1e-6);
    }

    void TestUpstreamProperties()
    {
        // Load test mesh
        TetrahedralMesh<1, 3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("lung/test/data/TestSubtreeProperties");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check num nodes and elements are as expected
        TS_ASSERT_EQUALS(mesh.GetNumAllNodes(), 10u);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), 9u);

        AirwayPropertiesCalculator properties_calculator(mesh, 0u);
        properties_calculator.CalculateUpstreamProperties();

        std::vector<double> branch_lengths = properties_calculator.GetUpstreamBranchLengths();
        std::vector<double> branch_volumes = properties_calculator.GetUpstreamBranchVolumes();
        std::vector<double> branch_areas = properties_calculator.GetUpstreamBranchLateralSurfaceAreas();
        std::vector<double> branch_resistances = properties_calculator.GetUpstreamPoiseuilleResistances();

        // Check correct number of branches found
        TS_ASSERT_EQUALS(properties_calculator.GetBranches().size(), 3u);

        // Check length properties
        TS_ASSERT_DELTA(branch_lengths[0], 6.0, 1e-6);
        TS_ASSERT_DELTA(branch_lengths[1], branch_lengths[2], 1e-6);
        TS_ASSERT_DELTA(branch_lengths[1], 6.0 + sqrt(72), 1e-6);

        // Check volume properties
        TS_ASSERT_DELTA(branch_volumes[0], 6.0 * M_PI, 1e-6);
        TS_ASSERT_DELTA(branch_volumes[1], branch_volumes[2], 1e-6);
        TS_ASSERT_DELTA(branch_volumes[1], (6.0 + sqrt(72)) * M_PI, 1e-6);

        // Check area properties
        TS_ASSERT_DELTA(branch_areas[0], 6.0 * 2 * M_PI, 1e-6);
        TS_ASSERT_DELTA(branch_areas[1], branch_areas[2], 1e-6);
        TS_ASSERT_DELTA(branch_areas[1], (6.0 + sqrt(72)) * 2 * M_PI, 1e-6);

        // Check Poiseuille resistance properties
        TS_ASSERT_DELTA(branch_resistances[0], 6.0, 1e-6);
        TS_ASSERT_DELTA(branch_resistances[1], branch_lengths[2], 1e-6);
        TS_ASSERT_DELTA(branch_resistances[1], 6.0 + sqrt(72), 1e-6);

    }
};

#endif /*_TESTAIRWAYPROPERTIESCALCULATOR_HPP_*/
