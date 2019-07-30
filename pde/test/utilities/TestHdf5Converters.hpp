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

#ifndef TESTHDF5CONVERTERS_HPP_
#define TESTHDF5CONVERTERS_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5ToVtkConverter.hpp"
#include "Hdf5ToTxtConverter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToXdmfConverter.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "FileFinder.hpp"
#include "FileComparison.hpp"
#include "NumericFileComparison.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"


#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

typedef Hdf5ToVtkConverter<3,3> VTK_3D;
typedef Hdf5ToMeshalyzerConverter<3,3> MESHA_3D;

/* HOW_TO_TAG Cardiac/Post-processing
 * Convert already generated simulation (HDF5) results to text, VTK or Meshalyzer format.
 */

class TestHdf5Converters : public CxxTest::TestSuite
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

    /**
     * This tests the HDF5 to VTK converter using a 3D example
     * taken from a bidomain simulation.
     */
    void TestBidomainVtkConversion3D()
    {
#ifdef CHASTE_VTK // Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string working_directory = "TestHdf5ToVtkConverter_bidomain";
        std::string working_directory2 = "TestHdf5ToVtkConverter_bidomain2";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToVtkConverter_bidomain,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("pde/test/data/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToVtkConverter<3,3> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                          "cube_2mm_12_elements", &mesh, false, true);
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        // Now do something else before reading back...

        {
            //NOTE: Interleaved test

            // Show that trying to write .pvtu files from a TetrahedralMesh gives a warning (but writes anyway)
            CopyToTestOutputDirectory("pde/test/data/cube_2mm_12_elements.h5",
                                  working_directory2);
            Hdf5ToVtkConverter<3,3> converter2(FileFinder(working_directory2, RelativeTo::ChasteTestOutput),
                                               "cube_2mm_12_elements", &mesh, true, true);

            //The reading part of this test is below
        }

        /*
         * Note that VTK is not thread-safe. The master process has spawned
         * a child to write the mesh and may still be writing! This barrier
         * just slows things down a bit.
         */
        PetscTools::Barrier();


        VtkMeshReader<3,3> vtk_mesh_reader(test_output_directory + working_directory
                                           + "/vtk_output/cube_2mm_12_elements.vtu");
        TS_ASSERT_EQUALS(vtk_mesh_reader.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vtk_mesh_reader.GetNumElements(), 12u);

        std::vector<double> first_node = vtk_mesh_reader.GetNextNode();
        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[2], 0.0, 1e-6);

        std::vector<double> next_node = vtk_mesh_reader.GetNextNode();
        TS_ASSERT_DELTA(next_node[0], 0.2, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(next_node[2], 0.0, 1e-6);

        // V_m and phi_e samples
        std::vector<double> v_at_last, phi_at_last;
        vtk_mesh_reader.GetPointData("V_000001", v_at_last);
        TS_ASSERT_DELTA(v_at_last[0],  -46.3761, 1e-3);
        TS_ASSERT_DELTA(v_at_last[6],  -46.3761, 1e-3);
        TS_ASSERT_DELTA(v_at_last[11], -46.3760, 1e-3);
        vtk_mesh_reader.GetPointData("Phi_e_000001", phi_at_last);
        TS_ASSERT_DELTA(phi_at_last[0],  0.0, 1e-3);
        TS_ASSERT_DELTA(phi_at_last[6],  0.0, 1e-3);
        TS_ASSERT_DELTA(phi_at_last[11], 0.0, 1e-3);

        {
            //NOTE: Interleaved test
            //The writing part of this test is above
            VtkMeshReader<3,3> vtk_mesh_reader2(test_output_directory + working_directory2
                                                + "/vtk_output/cube_2mm_12_elements.vtu");
            TS_ASSERT_EQUALS(vtk_mesh_reader2.GetNumNodes(), 12u);
        }

#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * This tests the HDF5 to VTK converter in parallel using a 2D example
     * taken from a monodomain simulation.
     */
    void TestMonodomainParallelVtkConversion2D()
    {
#ifdef CHASTE_VTK // Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string working_directory = "TestHdf5ToVtkConverter_monodomain2D";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToVtkConverter_monodomain2D,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("pde/test/data/2D_0_to_1mm_400_elements.h5",
                                  working_directory);

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2D_0_to_1mm_400_elements");
        DistributedTetrahedralMesh<2,2> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToVtkConverter<2,2> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                          "2D_0_to_1mm_400_elements", &mesh, true, false);

        /*
         * Note that VTK is not thread-safe. The master process has spawned
         * a child to write the mesh and may still be writing! This barrier
         * just slows things down a bit.
         */
        PetscTools::Barrier();

        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::stringstream filepath;
        filepath << test_output_directory << working_directory << "/vtk_output/2D_0_to_1mm_400_elements";
        if (!PetscTools::IsSequential())
        {
            filepath << "_" << PetscTools::GetMyRank();
        }
        filepath <<  ".vtu";

        VtkMeshReader<2,2> vtk_mesh_reader(filepath.str());
        TS_ASSERT_EQUALS(vtk_mesh_reader.GetNumNodes(), mesh.GetNumLocalNodes() + mesh.GetNumHaloNodes()); // 221 in total
        TS_ASSERT_EQUALS(vtk_mesh_reader.GetNumElements(), mesh.GetNumLocalElements()); // 400 in total

        if (PetscTools::IsSequential())
        {
            std::vector<double> first_node = vtk_mesh_reader.GetNextNode();
            TS_ASSERT_DELTA(first_node[0], 0.0 , 1e-6);
            TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(first_node[2], 0.0 , 1e-6); // 2d VTK files still carry z-coordinate

            std::vector<double> next_node = vtk_mesh_reader.GetNextNode();
            TS_ASSERT_DELTA(next_node[0], 0.01, 1e-6);
            TS_ASSERT_DELTA(next_node[1], 0.0 , 1e-6);
            TS_ASSERT_DELTA(next_node[2], 0.0 , 1e-6); // 2d VTK files still carry z-coordinate
        }

        // V_m samples
        std::vector<double> v_at_last;
        vtk_mesh_reader.GetPointData("V_000020", v_at_last);
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_DELTA(v_at_last[0],   -83.8534, 1e-3);
            TS_ASSERT_DELTA(v_at_last[110], -83.8534, 1e-3);
            TS_ASSERT_DELTA(v_at_last[220], -83.8530, 1e-3);
        }

        // Show that trying to write .pvtu files with original node ordering gives a warning (but writes anyway)
        Hdf5ToVtkConverter<2,2> converter2(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                           "2D_0_to_1mm_400_elements", &mesh, true, true);

        /*
         * Note that VTK is not thread-safe. The master process has spawned
         * a child to write the mesh and may still be writing! This barrier
         * just slows things down a bit.
         */
        PetscTools::Barrier();

        VtkMeshReader<2,2> vtk_mesh_reader2(test_output_directory + working_directory
                                            + "/vtk_output/2D_0_to_1mm_400_elements.vtu");
        TS_ASSERT_EQUALS(vtk_mesh_reader2.GetNumNodes(), 221u);
#else
        std::cout << "This test was not run, as VTK is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }

    /**
     * This tests the HDF5 to .txt converter using a 3D example
     * taken from a bidomain simulation.
     */
    void TestBidomainTxtConversion3D()
    {
        std::string working_directory = "TestHdf5ToTxtConverter_bidomain";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToTxtConverter_bidomain,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("pde/test/data/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToTxtConverter<3,3> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                          "cube_2mm_12_elements", &mesh);

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("cube_2mm_12_elements_V_0.txt");
        files_to_compare.push_back("cube_2mm_12_elements_V_1.txt");
        files_to_compare.push_back("cube_2mm_12_elements_Phi_e_0.txt");
        files_to_compare.push_back("cube_2mm_12_elements_Phi_e_1.txt");

        for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            std::cout << "Comparing generated and reference " << files_to_compare[i] << std::endl;
            FileFinder generated_file(working_directory +"/txt_output/" + files_to_compare[i], RelativeTo::ChasteTestOutput);
            FileFinder reference_file("pde/test/data/" + files_to_compare[i], RelativeTo::ChasteSourceRoot);
            NumericFileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestMonodomainMeshalyzerConversion()
    {
        // Firstly, copy ./heart/test/data/MonoDg01d/*.h5 to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
        // as that is where the reader reads from.
        std::string output_folder("TestHdf5Converters_TestMonodomainMeshalyzerConversion");
        CopyToTestOutputDirectory("heart/test/data/Monodomain1d/MonodomainLR91_1d.h5", output_folder);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToMeshalyzerConverter<1,1> converter(FileFinder(output_folder,RelativeTo::ChasteTestOutput),
                                                 "MonodomainLR91_1d", &mesh, true, 10 /* precision specified for coverage */);

        // Compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        FileComparison(test_output_directory + output_folder + "/output/MonodomainLR91_1d_V.dat",
                       "heart/test/data/Monodomain1d/MonodomainLR91_1d_V.dat").CompareFiles();
        FileComparison(test_output_directory + output_folder + "/output/MonodomainLR91_1d_times.info",
                       "heart/test/data/Monodomain1d/MonodomainLR91_1d_times.info").CompareFiles();
    }

    void TestBidomainMeshalyzerConversion()
    {
        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
         * as that is where the reader reads from.
         */
        std::string output_folder("TestHdf5Converters_TestBidomainMeshalyzerConversion");
        CopyToTestOutputDirectory("heart/test/data/Bidomain1d/bidomain.h5", output_folder);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToMeshalyzerConverter<1,1> converter(FileFinder(output_folder, RelativeTo::ChasteTestOutput),
                                                 "bidomain", &mesh, true);

        // Compare the voltage file
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        FileComparison(test_output_directory + output_folder + "/output/bidomain_V.dat",
                       "heart/test/data/Bidomain1d/bidomain_V.dat").CompareFiles();

        // Compare the Phi_e file
        FileComparison(test_output_directory + output_folder + "/output/bidomain_Phi_e.dat",
                       "heart/test/data/Bidomain1d/bidomain_Phi_e.dat").CompareFiles();

        // Compare the time information file
        FileComparison(test_output_directory + output_folder + "/output/bidomain_times.info",
                       "heart/test/data/Bidomain1d/bidomain_times.info").CompareFiles();
    }

    // This test covers the case when the hdf5 file contains 3 variables (e.g., after solving a problem with PROBLEM_DIM=3)
    void TestMeshalyzerConversion3Variables()
    {
        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
         * as that is where the reader reads from.
         */
        std::string output_folder("TestHdf5Converters_TestMeshalyzerConversion3Variables");
        CopyToTestOutputDirectory("heart/test/data/three_variables/3_vars.h5", output_folder);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToMeshalyzerConverter<1,1> converter(FileFinder(output_folder, RelativeTo::ChasteTestOutput),
                                                 "3_vars", &mesh, true);

        // Compare the first voltage file
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        NumericFileComparison(test_output_directory + output_folder + "/output/3_vars_Vm_1.dat",
                              "heart/test/data/three_variables/extended_bidomain_Vm_1.dat").CompareFiles();

        // Compare the second voltage file
        NumericFileComparison(test_output_directory + output_folder +"/output/3_vars_Vm_2.dat",
                              "heart/test/data/three_variables/extended_bidomain_Vm_2.dat").CompareFiles();

        // Compare the Phi_e file
        NumericFileComparison(test_output_directory + output_folder + "/output/3_vars_Phi_e.dat",
                              "heart/test/data/three_variables/extended_bidomain_Phi_e.dat").CompareFiles();

        // Compare the time information file
        FileComparison(test_output_directory + output_folder + "/output/3_vars_times.info",
                       "heart/test/data/three_variables/extended_bidomain_times.info").CompareFiles();
    }

    // This test covers the case when the hdf5 file contains more than 3 variables
    void TestMeshalyzerConversionLotsOfVariables()
    {
        std::string output_dir = "TestHdf5Converters_TestMeshalyzerConversionLotsOfVariables";

        /*
         * Firstly, copy the .h5 file to CHASTE_TEST_OUTPUT/TestHdf5ToMeshalyzerConverter,
         * as that is where the reader reads from.
         */
        CopyToTestOutputDirectory("heart/test/data/many_variables/many_variables.h5", output_dir);

        TrianglesMeshReader<1,1> mesh_reader("heart/test/data/many_variables/1D_65_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToMeshalyzerConverter<1,1> converter(FileFinder(output_dir, RelativeTo::ChasteTestOutput),
                                                 "many_variables", &mesh, true);

        std::vector<std::string> variable_names;
        variable_names.push_back("V");
        variable_names.push_back("I_ks");
        variable_names.push_back("I_kr");
        variable_names.push_back("I_Ca_tot");
        variable_names.push_back("I_tot");
        variable_names.push_back("I_Na_tot");

        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        for (unsigned i=0; i<variable_names.size(); i++)
        {
            // Compare the results files
            FileComparison(test_output_directory + "/" + output_dir + "/output/many_variables_"
                           + variable_names[i] + ".dat",
                           "heart/test/data/many_variables/many_variables_"
                           + variable_names[i] + ".dat").CompareFiles();
        }

        // Compare the time information file
        FileComparison(test_output_directory + output_dir + "/output/many_variables_times.info",
                       "heart/test/data/many_variables/many_variables_times.info").CompareFiles();
    }

    /**
     * This tests the HDF5 to XDMF converter
     */
    void TestHdf5ToXdmfConverter()
    {
#ifndef _MSC_VER
        std::string working_directory = "TestHdf5Converters_TestHdf5ToXdmfConverter";

        CopyToTestOutputDirectory("pde/test/data/cube_2mm_12_elements.h5",
                                  working_directory);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Convert
        Hdf5ToXdmfConverter<3,3> converter(FileFinder(working_directory, RelativeTo::ChasteTestOutput),
                                           "cube_2mm_12_elements",
                                           &mesh);

        std::vector<std::string> files_to_compare;
        files_to_compare.push_back("cube_2mm_12_elements.xdmf");
        files_to_compare.push_back("cube_2mm_12_elements_geometry_0.xml");
        files_to_compare.push_back("cube_2mm_12_elements_topology_0.xml");

        for (unsigned i=0; i<files_to_compare.size(); i++)
        {
            std::cout << "Comparing generated and reference " << files_to_compare[i] << std::endl;
            FileFinder generated_file(working_directory +"/xdmf_output/" + files_to_compare[i], RelativeTo::ChasteTestOutput);
            FileFinder reference_file("pde/test/data/xdmf_output/" + files_to_compare[i], RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
#endif // _MSC_VER
    }
};

#endif /*TESTHDF5CONVERTERS_HPP_*/
