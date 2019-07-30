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

#ifndef TESTVOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_
#define TESTVOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include "VoltageInterpolaterOntoMechanicsMesh.hpp"
#include "FileFinder.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Hdf5DataReader.hpp"
#include "ReplicatableVector.hpp"

class TestVoltageInterpolaterOntoMechanicsMesh : public CxxTest::TestSuite
{
private :
    // copies a file (relative to Chaste home to CHASTE_TEST_OUTPUT/dir
    void CopyToTestOutputDirectory(std::string file, std::string dir)
    {
        OutputFileHandler handler(dir);
        FileFinder file_finder(file, RelativeTo::ChasteSourceRoot);

        FileFinder copied_file = handler.CopyFileTo(file_finder);
        TS_ASSERT(copied_file.IsFile());
    }
public:

    void TestWith1dData(void)
    {
        // firstly, copy ./heart/test/data/MonoDg01d/*.h5 to CHASTE_TEST_OUTPUT/TestVoltageInterpolater1d,
        // as that is where the interpolater reads and writes to
        CopyToTestOutputDirectory("heart/test/data/Monodomain1d/MonodomainLR91_1d.h5",
                                  "TestVoltageInterpolater1d");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        QuadraticMesh<1> mech_mesh;
        TrianglesMeshReader<1,1> mesh_reader2("mesh/test/data/1D_0_to_1_10_elements_quadratic",2,1,false);
        mech_mesh.ConstructFromMeshReader(mesh_reader2);

        // move the node at 0.1 to 0.105, so it is 30% between the electrics nodes
        // 0.10 and 0.11
        mech_mesh.GetNode(1)->rGetModifiableLocation()[0] = 0.103;

        std::vector<std::string> variable_names;
        variable_names.push_back("V");

        VoltageInterpolaterOntoMechanicsMesh<1> interpolater(mesh,
                                                             mech_mesh,
                                                             variable_names,
                                                             "TestVoltageInterpolater1d",
                                                             "MonodomainLR91_1d");

        Hdf5DataReader fine_reader("TestVoltageInterpolater1d","MonodomainLR91_1d");
        DistributedVectorFactory factory1(mesh.GetNumNodes());
        Vec voltage_fine = factory1.CreateVec();

        Hdf5DataReader coarse_reader("TestVoltageInterpolater1d","voltage_mechanics_mesh");
        DistributedVectorFactory factory2(mech_mesh.GetNumNodes());
        Vec voltage_coarse = factory2.CreateVec();

        TS_ASSERT_EQUALS(fine_reader.GetUnlimitedDimensionValues().size(), coarse_reader.GetUnlimitedDimensionValues().size());

        bool invalid = true;//see below

        // The number of times is 1000, no point going through them all
        // only go through every 100th (note the "i+=100" in the for loop)
        assert(fine_reader.GetUnlimitedDimensionValues().size()==1001u);
        for (unsigned i=1; i<fine_reader.GetUnlimitedDimensionValues().size(); i+=100)
        {
            fine_reader.GetVariableOverNodes(voltage_fine, "V", i);
            coarse_reader.GetVariableOverNodes(voltage_coarse, "V", i);

            ReplicatableVector fine_repl(voltage_fine);
            ReplicatableVector coarse_repl(voltage_coarse);

            TS_ASSERT_EQUALS(coarse_repl.GetSize(), mech_mesh.GetNumNodes());

            TS_ASSERT_DELTA(fine_repl[0], coarse_repl[0], 1e-9);
            TS_ASSERT_DELTA(fine_repl[fine_repl.GetSize()-1], coarse_repl[coarse_repl.GetSize()-1], 1e-9);

            // Check the value at the displaced node is interpolated correctly
            TS_ASSERT_DELTA(coarse_repl[1], fine_repl[10]*0.7 + 0.3*fine_repl[11], 1e-9);

            // check fine_repl[10] != fine_repl[11] *for some time*, else the above
            // check is invalid
            if (fabs(fine_repl[10] - fine_repl[11]) > 1e-6)
            {
                invalid = false;
            }
        }
        UNUSED_OPT(invalid);
        assert(!invalid);

        PetscTools::Destroy(voltage_coarse);
        PetscTools::Destroy(voltage_fine);
    }

    // the data in this test came from TestCardiacElectroMechanicsProblem::TestImplicitNhs2dOneMechanicsElement()
    void TestWith2dData(void)
    {
        // firstly, copy .h5 file to CHASTE_TEST_OUTPUT/TestVoltageInterpolater2d,
        // as that is where the interpolater reads and writes to
        CopyToTestOutputDirectory("heart/test/data/Monodomain2d.h5",
                                  "TestVoltageInterpolater2d");


        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.01, 0.05, 0.05);

        QuadraticMesh<2> mech_mesh(0.05, 0.05, 0.05);

        std::vector<std::string> variable_names;
        variable_names.push_back("V");
        VoltageInterpolaterOntoMechanicsMesh<2> interpolater(mesh,
                                                             mech_mesh,
                                                             variable_names,
                                                             "TestVoltageInterpolater2d",
                                                             "Monodomain2d");


        Hdf5DataReader fine_reader("TestVoltageInterpolater2d","Monodomain2d");
        DistributedVectorFactory factory1(mesh.GetNumNodes());
        Vec voltage_fine = factory1.CreateVec();
        fine_reader.GetVariableOverNodes(voltage_fine, "V", 1); // time = first printed time is when the V varies most with space
        ReplicatableVector fine_repl(voltage_fine);

        Hdf5DataReader coarse_reader("TestVoltageInterpolater2d","voltage_mechanics_mesh");
        DistributedVectorFactory factory2(mech_mesh.GetNumNodes());
        Vec voltage_coarse = factory2.CreateVec();
        coarse_reader.GetVariableOverNodes(voltage_coarse, "V", 1);
        ReplicatableVector coarse_repl(voltage_coarse);

        // the first four nodes in the coarse mesh mesh are the corners and correspond
        // to nodes 0,5,30 and 35 in the fine mesh
        TS_ASSERT_DELTA(coarse_repl[0], fine_repl[0], 1e-9);
        TS_ASSERT_DELTA(coarse_repl[1], fine_repl[5], 2e-9);
        TS_ASSERT_DELTA(coarse_repl[2], fine_repl[30], 1e-9);
        TS_ASSERT_DELTA(coarse_repl[3], fine_repl[35], 1e-9);

        PetscTools::Destroy(voltage_coarse);
        PetscTools::Destroy(voltage_fine);
    }

    void TestWithMultipleVariables1D()
    {
        // firstly, copy .h5 file to CHASTE_TEST_OUTPUT/TestWithMultipleVariables1D,
        // as that is where the interpolater reads and writes to
        CopyToTestOutputDirectory("heart/test/data/CmguiData/bidomain/1D_0_to_1_100_elements.h5",
                                  "TestWithMultipleVariables1D");


        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        QuadraticMesh<1> mech_mesh;
        TrianglesMeshReader<1,1> mesh_reader2("mesh/test/data/1D_0_to_1_10_elements_quadratic",2,1,false);
        mech_mesh.ConstructFromMeshReader(mesh_reader2);

        std::vector<std::string> variable_names;
        variable_names.push_back("V");
        variable_names.push_back("Phi_e");
        VoltageInterpolaterOntoMechanicsMesh<1> interpolater(mesh,
                                                             mech_mesh,
                                                             variable_names,
                                                             "TestWithMultipleVariables1D",
                                                             "1D_0_to_1_100_elements");


        Hdf5DataReader fine_reader("TestWithMultipleVariables1D","1D_0_to_1_100_elements");
        DistributedVectorFactory factory1(mesh.GetNumNodes());
        Vec voltage_fine = factory1.CreateVec();
        Vec phi_e_fine = factory1.CreateVec();

        Hdf5DataReader coarse_reader("TestWithMultipleVariables1D","voltage_mechanics_mesh");
        DistributedVectorFactory factory2(mech_mesh.GetNumNodes());
        Vec voltage_coarse = factory2.CreateVec();
        Vec phi_e_coarse = factory2.CreateVec();

        bool invalid_1 = true;//see below
        bool invalid_2 = true;//see below

        TS_ASSERT_EQUALS(fine_reader.GetUnlimitedDimensionValues().size(), 101u);
        for (unsigned i=1; i<fine_reader.GetUnlimitedDimensionValues().size(); i++)
        {
            fine_reader.GetVariableOverNodes(voltage_fine, "V", i);
            fine_reader.GetVariableOverNodes(phi_e_fine, "Phi_e", i);
            ReplicatableVector fine_repl(voltage_fine);
            ReplicatableVector fine_repl_phi_e(phi_e_fine);

            coarse_reader.GetVariableOverNodes(voltage_coarse, "V", i);
            coarse_reader.GetVariableOverNodes(phi_e_coarse, "Phi_e", i);
            ReplicatableVector coarse_repl(voltage_coarse);
            ReplicatableVector coarse_repl_phi_e(phi_e_coarse);


            TS_ASSERT_EQUALS(coarse_repl.GetSize(), mech_mesh.GetNumNodes());
            TS_ASSERT_EQUALS(coarse_repl_phi_e.GetSize(), mech_mesh.GetNumNodes());
            TS_ASSERT_EQUALS(fine_repl.GetSize(), mesh.GetNumNodes());
            TS_ASSERT_EQUALS(fine_repl_phi_e.GetSize(), mesh.GetNumNodes());

            //first node
            TS_ASSERT_DELTA(fine_repl[0], coarse_repl[0], 1e-4);
            TS_ASSERT_DELTA(fine_repl_phi_e[0], coarse_repl_phi_e[0], 1e-4);

            //last node
            TS_ASSERT_DELTA(fine_repl[fine_repl.GetSize()-1], coarse_repl[coarse_repl.GetSize()-1], 1e-4);
            TS_ASSERT_DELTA(fine_repl_phi_e[fine_repl_phi_e.GetSize()-1], coarse_repl_phi_e[coarse_repl_phi_e.GetSize()-1], 1e-4);

            // Check that the value at a node that is not matching is interpolated correctly
            TS_ASSERT_DELTA(coarse_repl[1], fine_repl[10]*0.7 + 0.3*fine_repl[11], 1e-4);
            TS_ASSERT_DELTA(coarse_repl_phi_e[1], fine_repl_phi_e[10]*0.7 + 0.3*fine_repl_phi_e[11], 1e-4);

            // check fine_repl[10] != fine_repl[11] *for some time*, else the above
            // check is invalid
            if ((fabs(fine_repl[10] - fine_repl[11]) > 1e-6))
            {
                invalid_1 = false;
            }
            if (fabs(fine_repl_phi_e[10] - fine_repl_phi_e[11]) > 1e-6)
            {
                invalid_2 = false;
            }
        }
        TS_ASSERT_EQUALS(invalid_1, false);
        TS_ASSERT_EQUALS(invalid_2, false);

        PetscTools::Destroy(voltage_coarse);
        PetscTools::Destroy(voltage_fine);
        PetscTools::Destroy(phi_e_coarse);
        PetscTools::Destroy(phi_e_fine);
    }
};


#endif /*TESTVOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_*/
