/*

Copyright (c) 2005-2012, University of Oxford.
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
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

typedef Hdf5ToVtkConverter<3,3> VTK_3D;

class TestHdf5Converters : public CxxTest::TestSuite
{
private:
    // Copies a file (relative to Chaste home to CHASTE_TEST_OUTPUT/dir
    void CopyToTestOutputDirectory(std::string file, std::string dir)
    {
        if (PetscTools::AmMaster())
        {
            std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
            std::string command = "mkdir -p " + test_output_directory + dir;
            int return_value;
            return_value = system(command.c_str());
            assert(return_value == 0);
            command = "cp " + file + " " + test_output_directory + dir+"/";
            return_value = system(command.c_str());
            assert(return_value == 0);
        }
        PetscTools::Barrier();
    }

public:

    /**
     * This tests the HDF5 to VTK converter using a 3D example
     * taken from a bidomain simulation.
     */
    void TestBidomainVtkConversion3D() throw(Exception)
    {
#ifdef CHASTE_VTK // Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string working_directory = "TestHdf5ToVtkConverter_bidomain";
        OutputFileHandler handler(working_directory);

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
        Hdf5ToVtkConverter<3,3> converter(working_directory, "cube_2mm_12_elements", &mesh, false, true);

        /*
         * Note that VTK is not thread-safe. The master process has spawned
         * a child to write the mesh and may still be writing! This barrier
         * just slows things down a bit.
         */
        PetscTools::Barrier();

        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

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

        // Show that trying to write .pvtu files from a TetrahedralMesh gives a warning (but writes anyway)
        Hdf5ToVtkConverter<3,3> converter2(working_directory, "cube_2mm_12_elements", &mesh, true, true);
        VtkMeshReader<3,3> vtk_mesh_reader2(test_output_directory + working_directory
                                            + "/vtk_output/cube_2mm_12_elements.vtu");
        TS_ASSERT_EQUALS(vtk_mesh_reader2.GetNumNodes(), 12u);
#endif //CHASTE_VTK
    }

    /**
     * This tests the HDF5 to VTK converter in parallel using a 2D example
     * taken from a monodomain simulation.
     */
    void TestMonodomainParallelVtkConversion2D() throw(Exception)
    {
#ifdef CHASTE_VTK // Requires  "sudo aptitude install libvtk5-dev" or similar
        std::string working_directory = "TestHdf5ToVtkConverter_monodomain2D";
        OutputFileHandler handler(working_directory);

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
        Hdf5ToVtkConverter<2,2> converter(working_directory, "2D_0_to_1mm_400_elements", &mesh, true, false);

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
        Hdf5ToVtkConverter<2,2> converter2(working_directory, "2D_0_to_1mm_400_elements", &mesh, true, true);

        /*
         * Note that VTK is not thread-safe. The master process has spawned
         * a child to write the mesh and may still be writing! This barrier
         * just slows things down a bit.
         */
        PetscTools::Barrier();

        VtkMeshReader<2,2> vtk_mesh_reader2(test_output_directory + working_directory
                                            + "/vtk_output/2D_0_to_1mm_400_elements.vtu");
        TS_ASSERT_EQUALS(vtk_mesh_reader2.GetNumNodes(), 221u);
#endif //CHASTE_VTK
    }

    /**
     * This tests the HDF5 to .txt converter using a 3D example
     * taken from a bidomain simulation.
     */
    void TestBidomainTxtConversion3D() throw(Exception)
    {
        std::string working_directory = "TestHdf5ToTxtConverter_bidomain";
        OutputFileHandler handler(working_directory);

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
        Hdf5ToTxtConverter<3,3> converter(working_directory, "cube_2mm_12_elements", &mesh);

        // Compare the voltage file with a correct version
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();
        std::string command = "diff -a -I \"Created by Chaste\" " + test_output_directory
                                              + working_directory +"/txt_output/cube_2mm_12_elements_V_0.txt"
                                              + " pde/test/data/cube_2mm_12_elements_V_0.txt";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        command = "diff -a -I \"Created by Chaste\" " + test_output_directory
                                              + working_directory +"/txt_output/cube_2mm_12_elements_V_1.txt"
                                              + " pde/test/data/cube_2mm_12_elements_V_1.txt";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        command = "diff -a -I \"Created by Chaste\" " + test_output_directory
                                              + working_directory +"/txt_output/cube_2mm_12_elements_Phi_e_0.txt"
                                              + " pde/test/data/cube_2mm_12_elements_Phi_e_0.txt";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        command = "diff -a -I \"Created by Chaste\" " + test_output_directory
                                              + working_directory +"/txt_output/cube_2mm_12_elements_Phi_e_1.txt"
                                              + " pde/test/data/cube_2mm_12_elements_Phi_e_1.txt";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);
}
};

#endif /*TESTHDF5CONVERTERS_HPP_*/
