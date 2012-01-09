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

#ifndef TESTPOSTPROCESSINGWRITER_HPP_
#define TESTPOSTPROCESSINGWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "DistributedTetrahedralMesh.hpp"
#include "PostProcessingWriter.hpp"
#include "Hdf5DataReader.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistanceMapCalculator.hpp"
#include "HeartConfig.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "NumericFileComparison.hpp"
#include "TetrahedralMesh.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "BidomainProblem.hpp"
#include "LuoRudy1991.hpp"

class TestPostProcessingWriter : public CxxTest::TestSuite
{

public:
// These things are available in the XML output requests
//    <ActionPotentialDurationMap threshold="-30.0" threshold_unit="mV" repolarisation_percentage="90"/>
//    <UpstrokeTimeMap threshold="-30.0" threshold_unit="mV"/>
//    <MaxUpstrokeVelocityMap/>
//    <ConductionVelocityMap origin_node="10"/>
//    <ConductionVelocityMap origin_node="20"/>

    void TestWriterMethods() throw(Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        PostProcessingWriter<1,1> writer(mesh, "heart/test/data", "postprocessingapd", false);

        writer.WriteApdMapFile(60.0, -30.0);

        std::string file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd_60_-30_Map.dat";
        std::string file2 = "heart/test/data/PostProcessorWriter/good_apd_postprocessing.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));

        writer.WriteUpstrokeTimeMap(-30.0);

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/UpstrokeTimeMap_-30.dat";
        file2 = "heart/test/data/PostProcessorWriter/good_upstroke_time_postprocessing.dat";
        NumericFileComparison comp2(file1, file2);
        TS_ASSERT(comp2.CompareFiles(1e-12));

        writer.WriteMaxUpstrokeVelocityMap(-30.0);

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/MaxUpstrokeVelocityMap_-30.dat";
        file2 = "heart/test/data/PostProcessorWriter/good_upstroke_velocity_postprocessing.dat";
        NumericFileComparison comp3(file1, file2);
        TS_ASSERT(comp3.CompareFiles(1e-12));

        DistanceMapCalculator<1,1> dist_calculator(mesh);

        std::vector<unsigned> origin_node;
        origin_node.push_back(0);
        std::vector<double> distance_map_from_0;

        dist_calculator.ComputeDistanceMap(origin_node, distance_map_from_0);

        writer.WriteConductionVelocityMap(0u, distance_map_from_0);

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/ConductionVelocityFromNode0.dat";
        file2 = "heart/test/data/PostProcessorWriter/conduction_velocity_10_nodes_from_node_0.dat";
        NumericFileComparison comp4(file1, file2);
        TS_ASSERT(comp4.CompareFiles(1e-12));
    }

    void TestApdWritingWithNoApdsPresent() throw(Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig
        PostProcessingWriter<1,1> writer(mesh, "heart/test/data/Monodomain1d", "MonodomainLR91_1d", false);

        writer.WriteApdMapFile(90.0, -30.0);

        std::string file1 = OutputFileHandler::GetChasteTestOutputDirectory()  + output_dir + "/Apd_90_-30_Map.dat";
        std::string file2 = "heart/test/data/101_zeroes.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));
    }

    void TestPostProcessWriting() throw (Exception)
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_10_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<std::pair<double,double> > apd_maps;
        apd_maps.push_back(std::pair<double, double>(80,-30));//repolarisation percentage first, as per schema
        apd_maps.push_back(std::pair<double, double>(90,-20));//repolarisation percentage first, as per schema
        HeartConfig::Instance()->SetApdMaps(apd_maps);

        std::vector<double> upstroke_time_map;
        upstroke_time_map.push_back(-70.0);
        upstroke_time_map.push_back( 20.0);
        HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);

        std::vector<double> upstroke_velocity_map;
        upstroke_velocity_map.push_back(-50.0);
        upstroke_velocity_map.push_back(50.0);
        HeartConfig::Instance()->SetMaxUpstrokeVelocityMaps(upstroke_velocity_map);

        std::vector<unsigned> conduction_velocity_map;
        conduction_velocity_map.push_back(0u);
        HeartConfig::Instance()->SetConductionVelocityMaps(conduction_velocity_map);

        //test the mtehod that extrapolates nodal traces
        std::vector<unsigned> nodes_to_extrapolate;//1D test, does not cover node permutation case
        nodes_to_extrapolate.push_back(1u);
        nodes_to_extrapolate.push_back(99u);
        HeartConfig::Instance()->SetRequestedNodalTimeTraces(nodes_to_extrapolate);

        std::vector<ChastePoint<1> > pseudo_ecg_electrodes;
        pseudo_ecg_electrodes.push_back(ChastePoint<1>(11.0));
        pseudo_ecg_electrodes.push_back(ChastePoint<1>(-1.0));
        HeartConfig::Instance()->SetPseudoEcgElectrodePositions(pseudo_ecg_electrodes);

        PostProcessingWriter<1,1> writer(mesh, "heart/test/data/Monodomain1d", "MonodomainLR91_1d", false);

        writer.WritePostProcessingFiles();

        writer.WriteAboveThresholdDepolarisationFile(-40.0);

        std::string output_dir = "ChasteResults/output"; // default given by HeartConfig

        std::string file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd_80_-30_Map.dat";
        std::string file2 = "heart/test/data/101_zeroes.dat";
        NumericFileComparison comp1(file1, file2);
        TS_ASSERT(comp1.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/Apd_90_-20_Map.dat";
        file2 = "heart/test/data/101_zeroes.dat";
        NumericFileComparison comp2(file1, file2);
        TS_ASSERT(comp2.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/UpstrokeTimeMap_-70.dat";
        file2 = "heart/test/data/PostProcessorWriter/UpstrokeTimeMap_-70.dat";
        NumericFileComparison comp3(file1, file2);
        TS_ASSERT(comp3.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/UpstrokeTimeMap_20.dat";
        file2 = "heart/test/data/PostProcessorWriter/UpstrokeTimeMap_20.dat";
        NumericFileComparison comp4(file1, file2);
        TS_ASSERT(comp4.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/MaxUpstrokeVelocityMap_-50.dat";
        file2 = "heart/test/data/PostProcessorWriter/MaxUpstrokeVelocityMap_-50.dat";
        NumericFileComparison comp5(file1, file2);
        TS_ASSERT(comp5.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/MaxUpstrokeVelocityMap_50.dat";
        file2 = "heart/test/data/PostProcessorWriter/MaxUpstrokeVelocityMap_50.dat";
        NumericFileComparison comp6(file1, file2);
        TS_ASSERT(comp6.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/ConductionVelocityFromNode0.dat";
        file2 = "heart/test/data/PostProcessorWriter/conduction_velocity_100_nodes_from_node_0.dat";
        NumericFileComparison comp7(file1, file2);
        TS_ASSERT(comp7.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/PseudoEcgFromElectrodeAt_11_0_0.dat";
        file2 = "heart/test/data/PostProcessorWriter/PseudoEcgFromElectrodeAt_11_0_0.dat";
        NumericFileComparison comp8(file1, file2);
        TS_ASSERT(comp8.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/PseudoEcgFromElectrodeAt_-1_0_0.dat";
        file2 = "heart/test/data/PostProcessorWriter/PseudoEcgFromElectrodeAt_-1_0_0.dat";
        NumericFileComparison comp9(file1, file2);
        TS_ASSERT(comp9.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/AboveThresholdDepolarisations-40.dat";
        file2 = "heart/test/data/PostProcessorWriter/AboveThresholdDepolarisations-40.dat";
        NumericFileComparison comp10(file1, file2);
        TS_ASSERT(comp10.CompareFiles(1e-12));

        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/NodalTraces_V.dat";
        file2 = "heart/test/data/PostProcessorWriter/NodalTrace_V_Valid.dat";
        NumericFileComparison comp_nodes(file1, file2);
        TS_ASSERT(comp_nodes.CompareFiles(1e-12));

    }

    void TestExtractNodeTracesWithNodePermutation() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        HeartConfig::Instance()->Reset();
        //the point of this test is to check the method handles permutation properly...
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), false);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetSimulationDuration(2); //ms

        std::string output_dir = "TestingWriteNodesMethod";
        HeartConfig::Instance()->SetOutputDirectory(output_dir);
        HeartConfig::Instance()->SetOutputFilenamePrefix("NodalTracesTest");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;
        BidomainProblem<2> problem( &cell_factory );
        problem.SetMesh(&mesh);

        std::vector<unsigned> nodes;
        nodes.push_back(0u);
        nodes.push_back(13u);
        HeartConfig::Instance()->SetRequestedNodalTimeTraces(nodes);

        problem.Initialise();
        problem.Solve();

        // compare for two nodes (0 and 13) that I visually checked after running a no-permutation simulation (with Tetrahedral mesh)
        // and comparing the hdf file (with hdfview) with the postprocessing output
        // in order to consider the file "valid".

        //check the file with V (assuming default variable name was V)
        std::string file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/output/NodalTraces_V.dat";
        std::string file2 = "heart/test/data/PostProcessorWriter/NodalTrace_V_WithPermutationValid.dat";

        NumericFileComparison comp_V(file1, file2);
        TS_ASSERT(comp_V.CompareFiles(1e-3));

        //check file with Phi_e (assuming default variable name was Phi_e)
        file1 = OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/output/NodalTraces_Phi_e.dat";
        file2 = "heart/test/data/PostProcessorWriter/NodalTrace_Phi_e_WithPermutationValid.dat";

        NumericFileComparison comp_phie(file1, file2);
        TS_ASSERT(comp_phie.CompareFiles(1e-3));
    }

    void TestWritingEads() throw (Exception)
    {
        HeartConfig::Instance()->Reset();

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_10_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PostProcessingWriter<1,1> writer(mesh, "heart/test/data", "Ead", false);

        writer.WriteAboveThresholdDepolarisationFile(-30.0);

        //assuming output was written in the default directory ChasteResults
        std::string file1 = OutputFileHandler::GetChasteTestOutputDirectory() + "ChasteResults/output/AboveThresholdDepolarisations-30.dat";
        std::string file2 = "heart/test/data/PostProcessorWriter/AboveThresholdDepolarisations-30.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));
    }
};


#endif /*TESTPOSTPROCESSINGWRITER_HPP_*/
