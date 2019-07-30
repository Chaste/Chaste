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

#ifndef TESTPOSTPROCESSINGWRITER_HPP_
#define TESTPOSTPROCESSINGWRITER_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "PostProcessingWriter.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "DistanceMapCalculator.hpp"
#include "HeartConfig.hpp"
#include "NumericFileComparison.hpp"
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToCmguiConverter.hpp"

#include "PetscSetupAndFinalize.hpp"
//#include "VtkMeshReader.hpp" //Needed for commented out test, see #1660


class TestPostProcessingWriter : public CxxTest::TestSuite
{
    /**
     * Set the test output directory in HeartConfig, and return it.
     *
     * @param rOutputDirName  the leaf folder name
     * @return  FileFinder for the output folder within this
     */
    FileFinder GetPath(const std::string& rOutputDirName)
    {
        HeartConfig::Instance()->SetOutputDirectory(rOutputDirName);
        OutputFileHandler handler(rOutputDirName, false);
        return handler.FindFile("");
    }

    /**
     * Wipe a relevant test output directory and
     * move the results.h5 file from the repository
     * to the relevant cleaned test output directory.
     *
     * The h5 file will be used as the basis for postprocessing for the tests.
     */
    void CopyTestDataHdf5ToCleanTestOutputFolder(const FileFinder& rTestOutputDirectory, const std::string& rHdf5BaseName)
    {
        OutputFileHandler handler(rTestOutputDirectory);
        FileFinder hdf5_file("heart/test/data/" + rHdf5BaseName + ".h5", RelativeTo::ChasteSourceRoot);
        handler.CopyFileTo(hdf5_file);
        PetscTools::Barrier("CopyTestDataHdf5ToCleanTestOutputFolder");
    }

public:
    void TestWriterMethods()
    {
        FileFinder test_dir = GetPath("TestPostProcessingWriter_WriterMethods");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        CopyTestDataHdf5ToCleanTestOutputFolder(test_dir, "PostProcessingWriter/postprocessingapd");

        ///\todo #2359 - it isn't nice that the Postprocessing writer
        /// and HDF5 converter can't be in the same scope.
        {
            // Call post processing writer (output to HDF5)

            PostProcessingWriter<1,1> writer(mesh, test_dir, "postprocessingapd");
            writer.WriteApdMapFile(60.0, -30.0);
            writer.WriteUpstrokeTimeMap(-30.0);
            writer.WriteMaxUpstrokeVelocityMap(-30.0);

            DistanceMapCalculator<1,1> dist_calculator(mesh);
            std::vector<unsigned> origin_node;
            origin_node.push_back(0);
            std::vector<double> distance_map_from_0;
            dist_calculator.ComputeDistanceMap(origin_node, distance_map_from_0);
            writer.WriteConductionVelocityMap(0u, distance_map_from_0);
        }

        // Convert from HDF5 to meshalyzer format.
        Hdf5ToMeshalyzerConverter<1,1> converter(test_dir,
                                                 "postprocessingapd",
                                                 &mesh,
                                                 HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());

        // Check the meshalyzer files are good
        FileFinder output_file("output/Apd_60_minus_30_Map.dat", test_dir);
        TS_ASSERT(output_file.Exists());

        std::string file1 = output_file.GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessingWriter/good_apd_postprocessing.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));

        file1 = FileFinder("output/UpstrokeTimeMap_minus_30.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/good_upstroke_time_postprocessing.dat";
        NumericFileComparison comp2(file1, file2);
        TS_ASSERT(comp2.CompareFiles(1e-12));

        file1 = FileFinder("output/MaxUpstrokeVelocityMap_minus_30.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/good_upstroke_velocity_postprocessing.dat";
        NumericFileComparison comp3(file1, file2);
        TS_ASSERT(comp3.CompareFiles(1e-12));

        file1 = FileFinder("output/ConductionVelocityFromNode0.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/conduction_velocity_10_nodes_from_node_0.dat";
        NumericFileComparison comp4(file1, file2);
        TS_ASSERT(comp4.CompareFiles(1e-12));
    }

    void TestApdWritingWithNoApdsPresent()
    {
        FileFinder output_dir = GetPath("TestPostProcessingWriter_ApdWritingWithNoApdsPresent");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        CopyTestDataHdf5ToCleanTestOutputFolder(output_dir, "Monodomain1d/MonodomainLR91_1d");

        ///\todo #2359 allow PostProcessingWriter and Hdf5 converters to be in same scope (conflicting Hdf5DataReaders).
        {
            PostProcessingWriter<1,1> writer(mesh, output_dir, "MonodomainLR91_1d");
            writer.WriteApdMapFile(90.0, -30.0);
        }

        // Now (as part of #1660) call the converter in a separate step.
        Hdf5ToMeshalyzerConverter<1,1> converter(output_dir,
                                                 "MonodomainLR91_1d",
                                                 &mesh,
                                                 HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());

        std::string file1 = FileFinder("output/Apd_90_minus_30_Map.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessingWriter/101_zeroes.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));
    }

    void TestPostProcessWriting()
    {
        FileFinder test_dir = GetPath("TestPostProcessingWriter_PostProcessWriting");
        CopyTestDataHdf5ToCleanTestOutputFolder(test_dir, "Monodomain1d/MonodomainLR91_1d");

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
        upstroke_time_map.push_back(  0.4);
        upstroke_time_map.push_back( -0.4);
        HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);

        std::vector<double> upstroke_velocity_map;
        upstroke_velocity_map.push_back(-50.0);
        upstroke_velocity_map.push_back(50.0);
        HeartConfig::Instance()->SetMaxUpstrokeVelocityMaps(upstroke_velocity_map);

        std::vector<unsigned> conduction_velocity_map;
        conduction_velocity_map.push_back(0u);
        HeartConfig::Instance()->SetConductionVelocityMaps(conduction_velocity_map);

        //test the method that extrapolates nodal traces
        std::vector<unsigned> nodes_to_extrapolate;//1D test, does not cover node permutation case
        nodes_to_extrapolate.push_back(1u);
        nodes_to_extrapolate.push_back(99u);
        HeartConfig::Instance()->SetRequestedNodalTimeTraces(nodes_to_extrapolate);

        std::vector<ChastePoint<1> > pseudo_ecg_electrodes;
        pseudo_ecg_electrodes.push_back(ChastePoint<1>(11.0));
        pseudo_ecg_electrodes.push_back(ChastePoint<1>(-1.0));
        HeartConfig::Instance()->SetPseudoEcgElectrodePositions(pseudo_ecg_electrodes);

        ///\todo #2359 allow PostProcessingWriter and Hdf5 converters to be in same scope (conflicting Hdf5DataReaders).
        {
            PostProcessingWriter<1,1> writer(mesh, test_dir, "MonodomainLR91_1d");
            writer.WritePostProcessingFiles();
            writer.WriteAboveThresholdDepolarisationFile(-40.0);
        }

        // Now (as part of #1660) call the converter in a separate step.
        Hdf5ToMeshalyzerConverter<1,1> converter(test_dir,
                                                 "MonodomainLR91_1d",
                                                 &mesh,
                                                 HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());

        std::string file1 = FileFinder("output/Apd_80_minus_30_Map.dat", test_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessingWriter/101_zeroes.dat";
        NumericFileComparison comp1(file1, file2);
        TS_ASSERT(comp1.CompareFiles(1e-12));

        file1 = FileFinder("output/Apd_90_minus_20_Map.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/101_zeroes.dat";
        NumericFileComparison comp2(file1, file2);
        TS_ASSERT(comp2.CompareFiles(1e-12));

        file1 = FileFinder("output/UpstrokeTimeMap_minus_70.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/UpstrokeTimeMap_-70.dat";
        NumericFileComparison comp3(file1, file2);
        TS_ASSERT(comp3.CompareFiles(1e-12));

        file1 = FileFinder("output/UpstrokeTimeMap_20.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/UpstrokeTimeMap_20.dat";
        NumericFileComparison comp4(file1, file2);
        TS_ASSERT(comp4.CompareFiles(1e-12));

        file1 = FileFinder("output/MaxUpstrokeVelocityMap_minus_50.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/MaxUpstrokeVelocityMap_-50.dat";
        NumericFileComparison comp5(file1, file2);
        TS_ASSERT(comp5.CompareFiles(1e-12));

        file1 = FileFinder("output/MaxUpstrokeVelocityMap_50.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/MaxUpstrokeVelocityMap_50.dat";
        NumericFileComparison comp6(file1, file2);
        TS_ASSERT(comp6.CompareFiles(1e-12));

        file1 = FileFinder("output/ConductionVelocityFromNode0.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/conduction_velocity_100_nodes_from_node_0.dat";
        NumericFileComparison comp7(file1, file2);
        TS_ASSERT(comp7.CompareFiles(1e-12));

        file1 = FileFinder("output/PseudoEcgFromElectrodeAt_11_0_0.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/PseudoEcgFromElectrodeAt_11_0_0.dat";
        NumericFileComparison comp8(file1, file2);
        TS_ASSERT(comp8.CompareFiles(1e-12));

        file1 = FileFinder("output/PseudoEcgFromElectrodeAt_-1_0_0.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/PseudoEcgFromElectrodeAt_-1_0_0.dat";
        NumericFileComparison comp9(file1, file2);
        TS_ASSERT(comp9.CompareFiles(1e-12));

        file1 = FileFinder("AboveThresholdDepolarisations_minus_40.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/AboveThresholdDepolarisations-40.dat";
        NumericFileComparison comp10(file1, file2);
        TS_ASSERT(comp10.CompareFiles(1e-12));

        file1 = FileFinder("NodalTraces_V.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/NodalTrace_V_Valid.dat";
        NumericFileComparison comp_nodes(file1, file2);
        TS_ASSERT(comp_nodes.CompareFiles(1e-12));

    }

    void TestExtractNodeTracesWithNodePermutation()
    {
        HeartConfig::Instance()->Reset();
        FileFinder output_dir = GetPath("TestPostProcessingWriter_ExtractNodeTracesWithNodePermutation");
        HeartConfig::Instance()->SetOutputFilenamePrefix("NodalTracesTest");

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // The point of this test is to check the method handles permutation properly...
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), false);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetSimulationDuration(2); //ms

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
        std::string file1 = FileFinder("NodalTraces_V.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessingWriter/NodalTrace_V_WithPermutationValid.dat";

        NumericFileComparison comp_V(file1, file2);
        TS_ASSERT(comp_V.CompareFiles(1e-3));

        //check file with Phi_e (assuming default variable name was Phi_e)
        file1 = FileFinder("NodalTraces_Phi_e.dat", output_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/NodalTrace_Phi_e_WithPermutationValid.dat";

        NumericFileComparison comp_phie(file1, file2);
        TS_ASSERT(comp_phie.CompareFiles(1e-3));
    }

    void TestWritingEads()
    {
        HeartConfig::Instance()->Reset();
        FileFinder output_dir = GetPath("TestPostProcessingWriter_WritingEads");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_10_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        CopyTestDataHdf5ToCleanTestOutputFolder(output_dir, "PostProcessingWriter/Ead");

        PostProcessingWriter<1,1> writer(mesh, output_dir, "Ead");
        writer.WriteAboveThresholdDepolarisationFile(-30.0);

        std::string file1 = FileFinder("AboveThresholdDepolarisations_minus_30.dat", output_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessingWriter/AboveThresholdDepolarisations-30.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));
    }

    void TestSwitchingOutputFormat()
    {
        HeartConfig::Instance()->Reset();
        FileFinder test_dir = GetPath("TestPostProcessingWriter_SwitchingOutputFormat");
        HeartConfig::Instance()->SetVisualizeWithCmgui();

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        CopyTestDataHdf5ToCleanTestOutputFolder(test_dir, "PostProcessingWriter/postprocessingapd");

        ///\todo #2359 allow PostProcessingWriter and Hdf5 converters to be in same scope (conflicting Hdf5DataReaders).
        {
            PostProcessingWriter<1,1> writer(mesh, test_dir, "postprocessingapd");
            double upstroke_threshold = -30.0;
            writer.WriteApdMapFile(60.0, upstroke_threshold);
            writer.WriteApdMapFile(90.0, upstroke_threshold);
            writer.WriteAboveThresholdDepolarisationFile(upstroke_threshold); // This isn't written into a meshalyzer or cmgui format (two columns) - so keep exporting manually.
        }

        // Now (as part of #1660) call the converter in a separate step.
        Hdf5ToMeshalyzerConverter<1,1> converter(test_dir,
                                                 "postprocessingapd",
                                                 &mesh,
                                                 HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());

        // Now (as part of #1660) call the converter in a separate step.
        Hdf5ToCmguiConverter<1,1> converter2(test_dir,
                                            "postprocessingapd",
                                            &mesh);

//        // Now (as part of #1660) call the converter in a separate step.
//        Hdf5ToVtkConverter<1,1> converter3(test_dir,
//                                           "postprocessingapd",
//                                           &mesh,
//                                           false,
//                                           true);

        std::string file1 = FileFinder("output/Apd_60_minus_30_Map.dat", test_dir).GetAbsolutePath();
        std::string file2 = "heart/test/data/PostProcessingWriter/good_apd_postprocessing.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-12));

        file1 = FileFinder("AboveThresholdDepolarisations_minus_30.dat", test_dir).GetAbsolutePath();
        file2 = "heart/test/data/PostProcessingWriter/good_ead_postprocessing.dat";
        NumericFileComparison comp2(file1, file2);
        TS_ASSERT(comp2.CompareFiles(1e-12));

        for (unsigned i=0; i<451u; i++)
        {
            std::stringstream filename;
            filename << "postprocessingapd_" << i << ".exnode";
            FileFinder file("cmgui_output/" + filename.str(), test_dir);
            TS_ASSERT(file.IsFile());
        }
    }

    // This test checks that the APD map is put into the HDF5 file.
    // NB The dataset and variable name are Apd_60_minus_30_Map here...
    void TestHdfOutput()
    {
        HeartConfig::Instance()->Reset();
        FileFinder output_dir("TestPostProcessingWriter_AddingToHdf5", RelativeTo::ChasteTestOutput);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        CopyTestDataHdf5ToCleanTestOutputFolder(output_dir, "PostProcessingWriter/postprocessingapd");

        ///\todo #2359 See comment below for the reason for the block.
        {
            PostProcessingWriter<1,1> writer(mesh, output_dir, "postprocessingapd");
            writer.WriteApdMapFile(60.0, -30.0);
        }
        // This needs to be after the PostProcessingWriter has disappeared. This could perhaps
        // be because PostProcessingWriter::mpDataReader exists and needs to let go of the file
        // before someone can read the updated one? Not sure, but having the above lines in scope
        // causes this line to say it can't find the "Apd_60_minus_30_Map" dataset.

        Hdf5DataReader reader(output_dir, "postprocessingapd", "Apd_60_minus_30_Map"); // dataset name

        // Check that the dataset contains just the new APD map we requested
        std::vector<std::string> variable_names = reader.GetVariableNames();
        TS_ASSERT_EQUALS(variable_names.size(), 1u);
        TS_ASSERT_EQUALS(variable_names[0], "Apd_60_minus_30_Map"); // variable name

        // Get an error if you ask for a timestep (pace number) that isn't there...
        DistributedVectorFactory factory(reader.GetNumberOfRows());
        Vec data = factory.CreateVec();
        reader.GetVariableOverNodes(data, "Apd_60_minus_30_Map", 0/*timestep*/);

        ReplicatableVector rep_vec;
        rep_vec.ReplicatePetscVector(data);

        if (PetscTools::AmMaster())
        {
            TS_ASSERT_DELTA(rep_vec[0],  300.965, 1e-3);
            TS_ASSERT_DELTA(rep_vec[1],  301.023, 1e-3);
            TS_ASSERT_DELTA(rep_vec[2],  301.703, 1e-3);
            TS_ASSERT_DELTA(rep_vec[3],  303.247, 1e-3);
            TS_ASSERT_DELTA(rep_vec[4],  302.629, 1e-3);
            TS_ASSERT_DELTA(rep_vec[5],  301.575, 1e-3);
            TS_ASSERT_DELTA(rep_vec[6],  300.816, 1e-3);
            TS_ASSERT_DELTA(rep_vec[7],  299.559, 1e-3);
            TS_ASSERT_DELTA(rep_vec[8],  297.915, 1e-3);
            TS_ASSERT_DELTA(rep_vec[9],  296.589, 1e-3);
            TS_ASSERT_DELTA(rep_vec[10], 295.969, 1e-3);
        }

        PetscTools::Destroy(data);
    }

    void TestVtkOutput()
    {
#ifdef CHASTE_VTK
        HeartConfig::Instance()->Reset();

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetVisualizeWithVtk();
        // Switch off meshalyzer to make sure we don't need the "/output" directory.
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetOutputDirectory("TestPostProcessingWriter_VtkOutput");

        std::vector<std::pair<double,double> > apd_maps;
        apd_maps.push_back(std::pair<double, double>(90,-30)); // Repolarisation percentage and upstroke threshold.
        apd_maps.push_back(std::pair<double, double>(50,-30));
        HeartConfig::Instance()->SetApdMaps(apd_maps);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;
        BidomainProblem<2> problem( &cell_factory );
        problem.SetMesh(&mesh);

        problem.Initialise();
        problem.Solve();

        OutputFileHandler handler("TestPostProcessingWriter_VtkOutput", false);
        FileFinder output_dir = handler.FindFile("vtk_output");

        TS_ASSERT(output_dir.IsDir());

        FileFinder vtu_file = handler.FindFile("vtk_output/SimulationResults.vtu");
        TS_ASSERT(vtu_file.IsFile());

        // Read VTK file & check it doesn't cause any problems.
        VtkMeshReader<2,2> vtk_reader(output_dir.GetAbsolutePath() + "/SimulationResults.vtu");

        // Load up and check APD data.
        std::vector<double> apd_data;

        // APD90 map - all zeros since nothing recorded after such a short simulation.
        // See TestLongPostprocessing for one with some data in it.
        vtk_reader.GetPointData("Apd_90_minus_30_Map_000000", apd_data);
        TS_ASSERT_EQUALS(apd_data.size(), 221u);
        for (unsigned i=0; i<apd_data.size(); i++)
        {
            TS_ASSERT_DELTA(apd_data[i], 0.0, 1e-9);
        }

        // and APD50 map
        vtk_reader.GetPointData("Apd_50_minus_30_Map_000000", apd_data);
        TS_ASSERT_EQUALS(apd_data.size(), 221u);
        for (unsigned i=0; i<apd_data.size(); i++)
        {
            TS_ASSERT_DELTA(apd_data[i], 0.0, 1e-9);
        }
#else
        std::cout << "VTK is not installed / Chaste is not configured to use it, this test didn't do anything.\n";
#endif
    }

    void TestDifferentNumberOfPaces()
    {
        HeartConfig::Instance()->Reset();
        FileFinder output_dir("TestPostProcessingWriter_DifferentNumberOfPaces", RelativeTo::ChasteTestOutput);
        {
            DistributedTetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(1.0, 1.0);

            CopyTestDataHdf5ToCleanTestOutputFolder(output_dir, "PostProcessingWriter/DifferentNumberOfPaces");

            std::vector<double> upstroke_time_map;
            upstroke_time_map.push_back(0.0);
            HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_time_map);

            PostProcessingWriter<1,1> writer(mesh, output_dir, "DifferentNumberOfPaces");
            writer.WritePostProcessingFiles();
        }
        {
            Hdf5DataReader reader(output_dir, "DifferentNumberOfPaces", "UpstrokeTimeMap_0");

            std::vector<double> reference;
            reference.push_back(1.0);
            reference.push_back(1.0);
            reference.push_back(501.0);
            reference.push_back(-999.0);     // padded value from PostProcessingWriter.cpp

            DistributedVectorFactory factory(reader.GetNumberOfRows());
            Vec data = factory.CreateVec();
            reader.GetVariableOverNodes(data, "UpstrokeTimeMap_0");
            DistributedVector distributed_vector = factory.CreateDistributedVector(data);

            for (DistributedVector::Iterator index = distributed_vector.Begin();
                                 index!= distributed_vector.End();
                                 ++index)
            {
                TS_ASSERT_DELTA(distributed_vector[index], reference[index.Global], 1e-9);
            }

            PetscTools::Destroy(data);

        }
    }
};


#endif /*TESTPOSTPROCESSINGWRITER_HPP_*/
