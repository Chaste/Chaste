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
#ifndef TESTCONDUCTIVITYTENSORS_HPP_
#define TESTCONDUCTIVITYTENSORS_HPP_

#include "UblasCustomFunctions.hpp"

#include <cxxtest/TestSuite.h>
#include "OrthotropicConductivityTensors.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

typedef AxisymmetricConductivityTensors<2,2> AXI_2D;

class TestConductivityTensors : public CxxTest::TestSuite
{
public:
    void TestConstantTensor3D()
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,1);

        OrthotropicConductivityTensors<3,3> ortho_tensors;
        ortho_tensors.SetConstantConductivities(Create_c_vector(2.1, 0.8, 0.135));
        ortho_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 2.1);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 0.8);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), 0.135);
        }

        AxisymmetricConductivityTensors<3,3> axi_tensors;
        axi_tensors.SetConstantConductivities(Create_c_vector(2.1, 0.8, 0.8));
        axi_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), 2.1);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), 0.8);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), 0.8);
        }
    }

    void TestTensorException()
    {
        c_vector<double, 1> constant_conductivities_1d(Create_c_vector(2.1));
        c_vector<double, 2> constant_conductivities_2d(Create_c_vector(2.1, 0.8));
        c_vector<double, 3> constant_conductivities_3d(Create_c_vector(2.1, 0.8, 0.0));

        OrthotropicConductivityTensors<1,1> ortho_1d_tensors;
        TS_ASSERT_THROWS_THIS(ortho_1d_tensors.SetConstantConductivities(constant_conductivities_2d),"Wrong number of conductivities provided");
        TS_ASSERT_THROWS_THIS(ortho_1d_tensors.SetConstantConductivities(constant_conductivities_3d),"Wrong number of conductivities provided");

        OrthotropicConductivityTensors<3,3> ortho_3d_tensors;
        TS_ASSERT_THROWS_THIS(ortho_3d_tensors.SetConstantConductivities(constant_conductivities_1d),"Wrong number of conductivities provided");
        TS_ASSERT_THROWS_THIS(ortho_3d_tensors.SetConstantConductivities(constant_conductivities_2d),"Wrong number of conductivities provided");

        // AxisymmetricConductivityTensors only makes sense in 3D problems
        TS_ASSERT_THROWS_THIS(AXI_2D axi_tensor,"Axisymmetric anisotropic conductivity only makes sense in 3D");

        // Transversal and longitudinal conductivities should have the same value
        AxisymmetricConductivityTensors<3,3> axi_3d_tensor;
        TS_ASSERT_THROWS_THIS( axi_3d_tensor.SetConstantConductivities(Create_c_vector(0.5,0.25,0.15)),
                "Axisymmetric media defined: transversal and normal conductivities should have the same value" );

    }

    void TestFibreOrientationFileExceptions()
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,1);

        {
            OrthotropicConductivityTensors<3,3> ortho_tensors;
            ortho_tensors.SetConstantConductivities( Create_c_vector(2.1, 0.8, 0.135) );
            FileFinder file("non_existing_file.ortho", RelativeTo::CWD);
            ortho_tensors.SetFibreOrientationFile(file);
            TS_ASSERT_THROWS_CONTAINS(ortho_tensors.Init(&mesh),
                    "Failed to open fibre file"); // non existing file
        }

        {
            OrthotropicConductivityTensors<3,3> ortho_tensors;
            ortho_tensors.SetConstantConductivities( Create_c_vector(2.1, 0.8, 0.135) );
            FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic2D.ortho", RelativeTo::ChasteSourceRoot);
            ortho_tensors.SetFibreOrientationFile(file);
            TS_ASSERT_THROWS_CONTAINS(ortho_tensors.Init(&mesh),
                    "A line is incomplete in "); // mismatching SPACE_DIM and # vectors in file
        }

        {
            OrthotropicConductivityTensors<3,3> ortho_tensors;
            ortho_tensors.SetConstantConductivities( Create_c_vector(2.1, 0.8, 0.135) );
            FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic2DWrongFormat.ortho", RelativeTo::ChasteSourceRoot);
            ortho_tensors.SetFibreOrientationFile(file);
            TS_ASSERT_THROWS_THIS(ortho_tensors.Init(&mesh),
                    "First (non comment) line of the fibre orientation file should contain the number "
                    "of lines of data in the file (and possibly a BIN tag) at most"); // wrong file format
        }

        {
            OrthotropicConductivityTensors<3,3> ortho_tensors;
            ortho_tensors.SetConstantConductivities( Create_c_vector(2.1, 0.8, 0.135) );
            FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic3DShortFile.ortho", RelativeTo::ChasteSourceRoot);
            ortho_tensors.SetFibreOrientationFile(file);
            TS_ASSERT_THROWS_CONTAINS(ortho_tensors.Init(&mesh),"End of file"); // short file
        }

        {
            std::vector<c_vector<double, 3> > non_constant_conductivities;
            non_constant_conductivities.push_back(Create_c_vector(0,0,0));

            OrthotropicConductivityTensors<3,3> ortho_tensors;
            ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities); //1 element
            TS_ASSERT_THROWS_THIS(ortho_tensors.Init(&mesh),  "The size of the conductivities vector does not match the number of elements in the mesh"); //Mesh has 6 elements

            AxisymmetricConductivityTensors<3,3> axi_tensors;
            axi_tensors.SetNonConstantConductivities(&non_constant_conductivities); //1 element
            TS_ASSERT_THROWS_THIS(axi_tensors.Init(&mesh),  "The size of the conductivities vector does not match the number of elements in the mesh"); //Mesh has 6 elements
        }
        {

            OrthotropicConductivityTensors<3,3> ortho_tensors;
            FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic3D4Elements.ortho", RelativeTo::ChasteSourceRoot);
            ortho_tensors.SetFibreOrientationFile(file);
            TS_ASSERT_THROWS_THIS(ortho_tensors.Init(&mesh),  "The size of the fibre file does not match the number of elements in the mesh"); //Mesh has 6 elements

            AxisymmetricConductivityTensors<3,3> axi_tensors;
            axi_tensors.SetFibreOrientationFile(file);
            TS_ASSERT_THROWS_THIS(axi_tensors.Init(&mesh),  "The size of the fibre file does not match the number of elements in the mesh"); //Mesh has 6 elements
        }
    }

    void TestFibreOrientationTensor3D()
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,1);

        c_vector<double, 3> constant_conductivities(Create_c_vector(2.1,0.8,0.135));

        OrthotropicConductivityTensors<3,3> ortho_tensors;
        ortho_tensors.SetConstantConductivities(constant_conductivities);
        FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        ortho_tensors.SetFibreOrientationFile(file);
        ortho_tensors.Init(&mesh);

        /* In each element the SimpleOrthotropic3D tensor is just the coordinate axes.
         * The first direction x, which looks like the fibre direction, takes the first conductivity value etc.
         *
         */
        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), constant_conductivities[0]);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), constant_conductivities[1]);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), constant_conductivities[2]);
        }
    }

    void TestCompareOrthotropicAxisymmetricTensors()
    {
        c_vector<double, 3> constant_conductivities(Create_c_vector(7.0,3.5,3.5));
        // Note that the non-primary conductivities match, so that in-sheet and trans-sheet will
        // be equivalent in the orthotropic case - reducing it to an axisymmetric conductivity
        TetrahedralMesh<3,3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        mesh.ConstructFromMeshReader(mesh_reader);

        OrthotropicConductivityTensors<3,3> ortho_tensors;
        ortho_tensors.SetConstantConductivities(constant_conductivities);
        FileFinder ortho_file("heart/test/data/box_shaped_heart/box_heart.ortho", RelativeTo::ChasteSourceRoot);
        ortho_tensors.SetFibreOrientationFile(ortho_file);
        ortho_tensors.Init(&mesh);

        AxisymmetricConductivityTensors<3,3> axi_tensors;
        axi_tensors.SetConstantConductivities(constant_conductivities);
        FileFinder axi_file("heart/test/data/box_shaped_heart/box_heart.axi", RelativeTo::ChasteSourceRoot);
        axi_tensors.SetFibreOrientationFile(axi_file);
        axi_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](0,0), axi_tensors[tensor_index](0,0), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](0,1), axi_tensors[tensor_index](0,1), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](0,2), axi_tensors[tensor_index](0,2), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](1,0), axi_tensors[tensor_index](1,0), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](1,1), axi_tensors[tensor_index](1,1), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](1,2), axi_tensors[tensor_index](1,2), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](2,0), axi_tensors[tensor_index](2,0), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](2,1), axi_tensors[tensor_index](2,1), 1e-5);
            TS_ASSERT_DELTA(ortho_tensors[tensor_index](2,2), axi_tensors[tensor_index](2,2), 1e-5);
        }
    }

    void TestFibreOrientationAxisymmetric3D()
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,1);

        c_vector<double, 3> constant_conductivities(Create_c_vector(2.1,0.8,0.8));

        //OrthotropicConductivityTensors<3,3> ortho_tensors;
        //ortho_tensors.SetConstantConductivities(constant_conductivities);
        //FileFinder file1("heart/test/data/fibre_tests/SimpleAxisymmetric.axi", RelativeTo::ChasteSourceRoot);
        //ortho_tensors.SetFibreOrientationFile(file1);
        //TS_ASSERT_THROWS_CONTAINS(ortho_tensors.Init(&mesh), "Failed to open fibre file");
        AxisymmetricConductivityTensors<3,3> axi_tensors;
        axi_tensors.SetConstantConductivities(constant_conductivities);
        FileFinder file2("heart/test/data/fibre_tests/SimpleAxisymmetric.axi", RelativeTo::ChasteSourceRoot);
        axi_tensors.SetFibreOrientationFile(file2);
        axi_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), constant_conductivities[0]);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), constant_conductivities[1]);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), constant_conductivities[1]);
        }
    }

    void TestHeterogeneousConductivitiesTensor3D()
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,1);

        std::vector<c_vector<double, 3> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10,1));
        non_constant_conductivities.push_back(Create_c_vector(200,20,2));
        non_constant_conductivities.push_back(Create_c_vector(300,30,3));
        non_constant_conductivities.push_back(Create_c_vector(400,40,4));
        non_constant_conductivities.push_back(Create_c_vector(500,50,5));

        OrthotropicConductivityTensors<3,3> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        ortho_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<6; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), tensor_index);
        }

        AxisymmetricConductivityTensors<3,3> axi_tensors;
        axi_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        axi_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<6; tensor_index++)
        {
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), 10*tensor_index);
        }


    }

    void TestHeterogeneousConductivitiesTensorDistributed3D()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,2);

        std::vector<c_vector<double, 3> > non_constant_conductivities;
        //for (unsigned global_element_index=0; global_element_index<12; global_element_index++)
        for (DistributedTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            unsigned global_element_index=it->GetIndex();
            non_constant_conductivities.push_back(Create_c_vector(global_element_index*100,global_element_index*10,global_element_index));
        }

        OrthotropicConductivityTensors<3,3> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        ortho_tensors.Init(&mesh);

        for (DistributedTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            unsigned tensor_index=it->GetIndex();
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](2,2), tensor_index);
        }

        AxisymmetricConductivityTensors<3,3> axi_tensors;
        axi_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        axi_tensors.Init(&mesh);

        for (DistributedTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            unsigned tensor_index=it->GetIndex();
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](0,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,1), 10*tensor_index);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](1,2), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,0), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,1), 0.0);
            TS_ASSERT_EQUALS(axi_tensors[tensor_index](2,2), 10*tensor_index);
        }


    }
    void TestHeterogeneousCondPlusFibreOrientationTensor1D()
    {
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(4);

        std::vector<c_vector<double, 1> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0));
        non_constant_conductivities.push_back(Create_c_vector(100));
        non_constant_conductivities.push_back(Create_c_vector(200));
        non_constant_conductivities.push_back(Create_c_vector(300));

        OrthotropicConductivityTensors<1,1> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic1D.ortho", RelativeTo::ChasteSourceRoot);
        ortho_tensors.SetFibreOrientationFile(file);
        ortho_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<4; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
        }
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor2D()
    {
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(1,3);

        std::vector<c_vector<double, 2> > non_constant_conductivities;
        non_constant_conductivities.push_back(Create_c_vector(0,0));
        non_constant_conductivities.push_back(Create_c_vector(100,10));
        non_constant_conductivities.push_back(Create_c_vector(200,20));
        non_constant_conductivities.push_back(Create_c_vector(300,30));
        non_constant_conductivities.push_back(Create_c_vector(400,40));
        non_constant_conductivities.push_back(Create_c_vector(500,50));

        OrthotropicConductivityTensors<2,2> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        FileFinder file("heart/test/data/fibre_tests/SimpleOrthotropic2D.ortho", RelativeTo::ChasteSourceRoot);
        ortho_tensors.SetFibreOrientationFile(file);
        ortho_tensors.Init(&mesh);

        for (unsigned tensor_index=0; tensor_index<6; tensor_index++)
        {
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,0), 100*tensor_index);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](0,1), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,0), 0.0);
            TS_ASSERT_EQUALS(ortho_tensors[tensor_index](1,1), 10*tensor_index);
        }
    }

    void TestHeterogeneousCondPlusFibreOrientationTensor3DDistributedTetrahedralMesh()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(1,1,5); // This cube shape ensures that there exist elements not own by every processor (at least for p=2,3
        std::vector<c_vector<double, 3> > non_constant_conductivities;

        for (AbstractTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            unsigned element_index = it->GetIndex();
            non_constant_conductivities.push_back(element_index*Create_c_vector(100,10,1));
        }

        OrthotropicConductivityTensors<3,3> ortho_tensors;
        ortho_tensors.SetNonConstantConductivities(&non_constant_conductivities);
        FileFinder file("heart/test/data/fibre_tests/NonTrivialOrthotropic3D.ortho", RelativeTo::ChasteSourceRoot);
        ortho_tensors.SetFibreOrientationFile(file);
        ortho_tensors.Init(&mesh);

        for (AbstractTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            unsigned element_index = it->GetIndex();
            double g_f = element_index * 100;
            double g_l = element_index * 10;
            double g_n = element_index * 1;
            double v = element_index*2.0/30.0*M_PI;
            double tol = 1e-5;

            /*  The orthotropic tensors are calculated as T_i = F_i*G*F_i', where
             *
             *        [ cos(v_i) -sin(v_i) 0]         [ g_f   0   0 ]
             *  F_i = [ sin(v_i)  cos(v_i) 0] , G = [   0 g_l   0 ] and v_i = i*2.0/30.0*PI .
             *        [        0         0 1]       [   0   0 g_n ]
             *
             *  The lines below assert that this results in a tensor of the form
             *
             *        [ g_f*cos(v_i)*cos(v_i)+g_l*sin(v_i)*sin(v_i) g_f*sin(v_i)*cos(v_i)-g_l*cos(v_i)*sin(v_i)   0]
             *    T_i = [ g_f*cos(v_i)*sin(v_i)-g_l*sin(v_i)*cos(v_i) g_l*cos(v_i)*cos(v_i)+g_f*sin(v_i)*sin(v_i)   0]
             *          [                                           0                                              0 g_n]
             *
             *  for each element index i.
             *
             */

            TS_ASSERT_DELTA(ortho_tensors[element_index](0,0), g_f*cos(v)*cos(v) + g_l*sin(v)*sin(v), tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](0,1), g_f*sin(v)*cos(v) - g_l*cos(v)*sin(v), tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](0,2), 0.0, tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](1,0), g_f*cos(v)*sin(v) - g_l*sin(v)*cos(v), tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](1,1), g_l*cos(v)*cos(v) + g_f*sin(v)*sin(v), tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](1,2), 0.0, tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](2,0), 0.0, tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](2,1), 0.0, tol);
            TS_ASSERT_DELTA(ortho_tensors[element_index](2,2), g_n, tol);
        }

        AxisymmetricConductivityTensors<3,3> axi_tensors;
        axi_tensors.SetNonConstantConductivities(&non_constant_conductivities);

        axi_tensors.SetFibreOrientationFile(file);
        TS_ASSERT_THROWS_CONTAINS(axi_tensors.Init(&mesh),"Too many entries in a line in");

        FileFinder axi_file("heart/test/data/fibre_tests/NonTrivialAxisymmetric3D.axi", RelativeTo::ChasteSourceRoot);
        axi_tensors.SetFibreOrientationFile(axi_file);
        axi_tensors.Init(&mesh);

        for (AbstractTetrahedralMesh<3,3>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            unsigned element_index = it->GetIndex();
            double g_f = element_index * 100;
            double g_l = element_index * 10;
            double v = element_index*2.0/30.0*M_PI;
            double tol = 1e-5;

            /*  The axisymmetric tensors are calculated according to
             *
             *  T_i = g_l * I + (g_f - g_l) * f_i * f_i'
             *
             *  where I is the identity matrix, g_f is the conductivity in the fibre direction,
             *  g_l is the conductivity perpenticluar to the fibre direction and
             *
             *  f_i = [ cos(v_i) sin(v_i) 0 ]'
             *
             *  The lines below assert that this results in a tensor of the form
             *
             *        [ g_l + cos(v_i)*cos(v_i)*(g_f-g_l)      -cos(v_i)*sin(v_i)*(g_f-g_l)   0]
             *  T_i = [      -sin(v_i)*cos(v_i)*(g_f-g_l) g_l + sin(v_i)*sin(v_i)*(g_f-g_l)   0]
             *        [                                 0                                 0 g_l]
             *
             *  for each element index i.
             *
             */

            TS_ASSERT_DELTA(axi_tensors[element_index](0,0), g_l + cos(v)*cos(v)*(g_f - g_l), tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](0,1), cos(v)*sin(v)*(g_f - g_l), tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](0,2), 0.0, tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](1,0), sin(v)*cos(v)*(g_f - g_l), tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](1,1), g_l + sin(v)*sin(v)*(g_f - g_l), tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](1,2), 0.0, tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](2,0), 0.0, tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](2,1), 0.0, tol);
            TS_ASSERT_DELTA(axi_tensors[element_index](2,2), g_l, tol);
        }
    }
};

#endif /*TESTFIBREORIENTATIONTENSORS_HPP_*/
