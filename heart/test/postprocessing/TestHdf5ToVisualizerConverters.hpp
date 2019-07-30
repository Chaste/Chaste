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

#ifndef TESTHDF5TOVISUALIZERCONVERTERS_HPP_
#define TESTHDF5TOVISUALIZERCONVERTERS_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"

#include "Hdf5ToCmguiConverter.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "HeartConfig.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "FileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

typedef Hdf5ToCmguiConverter<3,3> CMGUI_3D;

/* HOW_TO_TAG Cardiac/Post-processing
 * Convert already generated simulation (HDF5) results to Cmgui format.
 */

class TestHdf5ToVisualizerConverters : public CxxTest::TestSuite
{
private:
    // Copies a file (relative to Chaste home to CHASTE_TEST_OUTPUT/dir
    void CopyToTestOutputDirectory(std::string file, std::string dir)
    {
        OutputFileHandler handler(dir);
        FileFinder file_finder(file, RelativeTo::ChasteSourceRoot);

        FileFinder copied_file = handler.CopyFileTo(file_finder);
        TS_ASSERT(copied_file.IsFile());
    }

public:

    // This test covers the case when the hdf5 file contains more than 3 variables
    void TestCmguiConversionLotsOfVariables()
    {
        std::string output_dir = "TestHdf5ToCmguiConversionManyVariables";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("heart/test/data/many_variables/many_variables.h5",
                                  output_dir);

        TrianglesMeshReader<1,1> mesh_reader("heart/test/data/many_variables/1D_65_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToCmguiConverter<1,1> converter(FileFinder(output_dir, RelativeTo::ChasteTestOutput),
                                            "many_variables", &mesh, false, HeartConfig::Instance()->GetVisualizerOutputPrecision());

        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        // Compare the results files
        FileComparison(test_output_directory + output_dir + "/cmgui_output/many_variables_0.exnode",
                       "heart/test/data/many_variables/many_variables_0.exnode").CompareFiles();

        // Check validity of cmgui script
        // Note that FileComparison ignored *all* comment lines (including the ones which are informative to the end-user)
        FileComparison(test_output_directory + output_dir + "/cmgui_output/LoadSolutions.com",
                       "heart/test/data/many_variables/CmguiValidScript.com").CompareFiles();
    }

    void TestMonodomainCmguiConversion3D()
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_monodomain";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements"); //Not used in the test for exceptions
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToCmguiConverter<3,3> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                            "cube_2mm_12_elements", &mesh);

        // Compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        FileComparison(test_output_directory + working_directory + "/cmgui_output/cube_2mm_12_elements_0.exnode",
                       "heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_0.exnode").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/cube_2mm_12_elements_1.exnode",
                       "heart/test/data/CmguiData/monodomain/cube_2mm_12_elements_1.exnode").CompareFiles();

        // Check validity of cmgui script
        FileComparison(test_output_directory + working_directory + "/cmgui_output/LoadSolutions.com",
                               "heart/test/data/CmguiData/monodomain/monodomain3dValidScript.com").CompareFiles();
    }

    void TestBidomainCmguiConversion3D()
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_bidomain";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_bidomain,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements"); //Not used in the test for exceptions
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToCmguiConverter<3,3> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                            "cube_2mm_12_elements", &mesh, false, HeartConfig::Instance()->GetVisualizerOutputPrecision());

        // Compare the voltage file with a correct version that is known to visualize correctly in Cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        FileComparison(test_output_directory + working_directory + "/cmgui_output/cube_2mm_12_elements_0.exnode",
                               "heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_0.exnode").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/cube_2mm_12_elements_1.exnode",
                               "heart/test/data/CmguiData/bidomain/cube_2mm_12_elements_1.exnode").CompareFiles();
    }

    void TestBidomainWithBathCmguiConversion1D()
    {
        std::string working_directory = "TestBidomainWithBathCmguiConversion1D";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestBidomainWithBathCmguiConversion1D,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d.h5",
                                  working_directory);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_two_attributes"); //Not used in the test for exceptions
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_with_bath_1d");
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<1,1> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                            "bidomain_with_bath_1d", &mesh, true, HeartConfig::Instance()->GetVisualizerOutputPrecision());

        // Compare the voltage file with a correct version that is known to visualize correctly in Cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        // Mesh file first, one exnode, one for bath and one for tissue
        FileComparison(test_output_directory + working_directory + "/cmgui_output/bidomain_with_bath_1d.exnode",
                                       "heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d.exnode").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/tissue.exelem",
                                               "heart/test/data/CmguiData/bidomain_with_bath/tissue.exelem").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/bath.exelem",
                                               "heart/test/data/CmguiData/bidomain_with_bath/bath.exelem").CompareFiles();

        // Then the data file
        FileComparison(test_output_directory + working_directory + "/cmgui_output/bidomain_with_bath_1d_0.exnode",
                                               "heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d_0.exnode").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/bidomain_with_bath_1d_1.exnode",
                                               "heart/test/data/CmguiData/bidomain_with_bath/bidomain_with_bath_1d_1.exnode").CompareFiles();
    }

    void TestMonodomainCmguiConversion2D()
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_monodomain2D";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_monodomain2D,
         * as that is where the reader reads from. This data file was generated on this mesh by
         * TestMonodomainProblem2DWithPointStimulusInTheVeryCentreOfTheMesh.
         */
        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements.h5",
                                  working_directory);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<2,2> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                            "2D_0_to_1mm_400_elements", &mesh);

        // Compare the voltage file with a correct version that visualizes Vm correctly in cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        FileComparison(test_output_directory + working_directory + "/cmgui_output/2D_0_to_1mm_400_elements_0.exnode",
                                               "heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements_0.exnode").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/2D_0_to_1mm_400_elements_1.exnode",
                                               "heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements_1.exnode").CompareFiles();
    }

    void TestBidomainCmguiConversion1D()
    {
        std::string working_directory = "TestHdf5ToCmguiConverter_bidomain1D";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter_bidomain1D,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements.h5",
                                  working_directory);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<1,1> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                            "1D_0_to_1_100_elements", &mesh);

        // Compare the voltage file with a correct version that visualizes both Vm and Phie correctly in cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        FileComparison(test_output_directory + working_directory + "/cmgui_output/1D_0_to_1_100_elements_0.exnode",
                                               "heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements_0.exnode").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/1D_0_to_1_100_elements_1.exnode",
                                               "heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements_1.exnode").CompareFiles();
    }

    void TestCmguiConversion1DWith3Variables()
    {
        std::string working_directory = "TestHdf5ToCmguiConverter3Vars";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToCmguiConverter3Vars,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("heart/test/data/three_variables/3_vars.h5",
                                  working_directory);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        HeartConfig::Instance()->SetOutputFilenamePrefix("3_vars");
        HeartConfig::Instance()->SetOutputDirectory(working_directory);
        Hdf5ToCmguiConverter<1,1> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                            "3_vars", &mesh);

        // Compare the voltage file with a correct version that visualizes both Vs and Phie correctly in cmgui
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/3_vars.exnode",
                                               "heart/test/data/CmguiData/extended_bidomain/3_vars.exnode").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/3_vars.exelem",
                                               "heart/test/data/CmguiData/extended_bidomain/3_vars.exelem").CompareFiles();

        FileComparison(test_output_directory + working_directory + "/cmgui_output/3_vars_25.exnode",
                                               "heart/test/data/CmguiData/extended_bidomain/3_vars_25.exnode").CompareFiles();
    }

    void TestExceptions()
    {
        std::string directory = "TestHdf5ConverterExceptions";

        // Not used until the number of nodes is checked
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        CopyToTestOutputDirectory("heart/test/data/CmguiData/monodomain/2D_0_to_1mm_400_elements.h5", directory);

        TS_ASSERT_THROWS_THIS(CMGUI_3D converter(FileFinder(directory, RelativeTo::ChasteTestOutput),
                                                         "2D_0_to_1mm_400_elements", &mesh),
                              "Mesh and HDF5 file have a different number of nodes");
    }
};

#endif /*TESTHDF5TOVISUALIZERCONVERTERS_HPP_*/
