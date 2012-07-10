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
#ifndef TESTSTREETERFIBREGENERATOR_HPP_
#define TESTSTREETERFIBREGENERATOR_HPP_

#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "NumericFileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

/*
 * HOW_TO_TAG Cardiac/Problem definition
 * Generate fibre field definitions for cardiac geometries using a mathematical rule approach
 *
 */
class TestStreeterFibreGenerator : public CxxTest::TestSuite
{
public:

    void TestSimpleOrthotropic() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        std::string epi_face_file = "heart/test/data/box_shaped_heart/epi.tri";
        std::string rv_face_file = "heart/test/data/box_shaped_heart/rv.tri";
        std::string lv_face_file = "heart/test/data/box_shaped_heart/lv.tri";

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("shorter_streeter", "box_heart.ortho", true);

        OutputFileHandler handler("shorter_streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "box_heart.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/box_shaped_heart/box_heart.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
    }

    void TestSimpleOrthotropicNotDistributed() throw (Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        std::string epi_face_file = "heart/test/data/box_shaped_heart/epi.tri";
        std::string rv_face_file = "heart/test/data/box_shaped_heart/rv.tri";
        std::string lv_face_file = "heart/test/data/box_shaped_heart/lv.tri";

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("shorter_streeter", "box_heart_not_dist.ortho", true);

        OutputFileHandler handler("shorter_streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "box_heart_not_dist.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/box_shaped_heart/box_heart.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
    }

    void TestExceptions()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);

        // No surfaces defined
        TS_ASSERT_THROWS_THIS(fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "file.fibres"),
                "Files defining the heart surfaces not set");

        // Wrong surface filename
        TS_ASSERT_THROWS_THIS(fibre_generator.SetSurfaceFiles("wrong_name", "wrong_name", "wrong_name", false),
                "Wrong surface definition file name wrong_name");

        std::string epi_face_file = "heart/test/data/box_shaped_heart/epi.tri";
        std::string rv_face_file = "heart/test/data/box_shaped_heart/rv.tri";
        std::string lv_face_file = "heart/test/data/box_shaped_heart/lv.tri";


        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        TS_ASSERT_THROWS_THIS(fibre_generator.GenerateOrthotropicFibreOrientation("shorter_streeter", "downsampled.ortho"),
            "Apex to base vector has not been set");
        TS_ASSERT_THROWS_THIS(fibre_generator.SetApexToBase(999),
            "Apex to base coordinate axis was out of range");

        c_vector<double, 3> axis;
        axis[0] = 0.0;
        axis[1] = 0.0;
        axis[2] = 0.0;
        TS_ASSERT_THROWS_THIS(fibre_generator.SetApexToBase(axis),
            "Apex to base vector should be non-zero");
        axis[1] = 42.0; //Will be normalised
        fibre_generator.SetApexToBase(axis);
    }

    void TestConstructStreeterOnLeftWedge() throw(Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/human_wedge_mesh/HumanWedgeMesh");
        std::string epi_face_file = "heart/test/data/human_wedge_mesh/epi.tri";
        std::string endo_face_file = "heart/test/data/human_wedge_mesh/endo.tri";

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);

        //Assume we are in the left ventricle
        fibre_generator.SetSurfaceFiles(epi_face_file, "", endo_face_file, true);

        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("human_wedge_mesh/", "HumanWedgeMeshLeft.ortho", true);

        OutputFileHandler handler("human_wedge_mesh", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "HumanWedgeMeshLeft.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/human_wedge_mesh/HumanWedgeMeshLeft.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
    }

    void TestConstructStreeterOnRightWedge() throw(Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/human_wedge_mesh/HumanWedgeMesh");
        std::string epi_face_file = "heart/test/data/human_wedge_mesh/epi.tri";
        std::string endo_face_file = "heart/test/data/human_wedge_mesh/endo.tri";

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);

        //Assume we are in the left ventricle
        fibre_generator.SetSurfaceFiles(epi_face_file, endo_face_file, "", true);

        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("human_wedge_mesh/", "HumanWedgeMeshRight.ortho", true);

        OutputFileHandler handler("human_wedge_mesh", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "HumanWedgeMeshRight.ortho";

        NumericFileComparison comp(fibre_file,"heart/test/data/human_wedge_mesh/HumanWedgeMeshRight.ortho");
        TS_ASSERT(comp.CompareFiles(1e-11));
    }

};

#endif /*TESTSTREETERFIBREGENERATOR_HPP_*/
