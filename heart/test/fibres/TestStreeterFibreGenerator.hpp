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
#ifndef TESTSTREETERFIBREGENERATOR_HPP_
#define TESTSTREETERFIBREGENERATOR_HPP_

#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "FibreReader.hpp"

#include "PetscSetupAndFinalize.hpp"

/*
 * HOW_TO_TAG Cardiac/Problem definition
 * Generate fibre field definitions for cardiac geometries using a mathematical rule approach
 *
 */
class TestStreeterFibreGenerator : public CxxTest::TestSuite
{
private:
    void CompareGeneratedWithReferenceFile(FileFinder& rGeneratedFile, FibreFileType generatedFileType, FileFinder& rReferenceFile, FibreFileType referenceFileType)
    {
        // Make sure that any parallel file writers have finished before attempting to open
        PetscTools::Barrier("Wait for WriteData()");
        FibreReader<3> fibre_reader1(rGeneratedFile, generatedFileType);
        FibreReader<3> fibre_reader2(rReferenceFile, referenceFileType);

        double tol = 1e-10;

        std::vector< c_vector<double, 3> > fibres1;
        std::vector< c_vector<double, 3> > second1;
        std::vector< c_vector<double, 3> > third1;

        std::vector< c_vector<double, 3> > fibres2;
        std::vector< c_vector<double, 3> > second2;
        std::vector< c_vector<double, 3> > third2;

        if (generatedFileType == ORTHO)
        {
            fibre_reader1.GetAllOrtho(fibres1, second1, third1);
        }
        else
        {
            fibre_reader1.GetAllAxi(fibres1);
        }

        if (referenceFileType == ORTHO)
        {
            fibre_reader2.GetAllOrtho(fibres2, second2, third2);
        }
        else
        {
            fibre_reader2.GetAllAxi(fibres2);
        }

        TS_ASSERT_EQUALS(fibres1.size(), fibres2.size());
        if (generatedFileType == ORTHO && referenceFileType == ORTHO)
        {
            TS_ASSERT_EQUALS(second1.size(), second2.size());
            TS_ASSERT_EQUALS(third1.size(), third2.size());
        }

        for (unsigned i = 0u; i < fibres1.size(); i++)
        {
            for (unsigned j=0u; j<3u; j++)
            {
                TS_ASSERT_DELTA( fibres1[i][j] , fibres2[i][j] , tol);
                if (generatedFileType == ORTHO && referenceFileType == ORTHO)
                {
                    TS_ASSERT_DELTA( second1[i][j] , second2[i][j] , tol);
                    TS_ASSERT_DELTA( third1[i][j] , third2[i][j] , tol);
                }
            }
        }
    }

public:

    void TestSimpleOrthotropic()
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

        OutputFileHandler handler("shorter_streeter", false);
        fibre_generator.WriteData(handler, "box_heart.ortho");

        FileFinder fibre_file_ascii = handler.FindFile("box_heart.ortho");
        FileFinder fibre_file_reference("heart/test/data/box_shaped_heart/box_heart.ortho", RelativeTo::ChasteSourceRoot);

        CompareGeneratedWithReferenceFile(fibre_file_ascii, ORTHO, fibre_file_reference, ORTHO);

        fibre_generator.SetWriteFileAsBinary();
        fibre_generator.WriteData(handler, "box_heart_binary.ortho");

        FileFinder fibre_file_binary = handler.FindFile("box_heart_binary.ortho");

        CompareGeneratedWithReferenceFile(fibre_file_binary, ORTHO, fibre_file_reference, ORTHO);
    }

    void TestSimpleOrthotropicNotDistributed()
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

        OutputFileHandler handler("shorter_streeter", false);
        fibre_generator.WriteData(handler, "box_heart_not_dist.ortho");

        FileFinder fibre_file1 = handler.FindFile("box_heart_not_dist.ortho");
        FileFinder fibre_file2("heart/test/data/box_shaped_heart/box_heart.ortho", RelativeTo::ChasteSourceRoot);

        CompareGeneratedWithReferenceFile(fibre_file1, ORTHO, fibre_file2, ORTHO);
    }

    void TestExceptions()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");

        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);

        // No surfaces defined
        OutputFileHandler handler("streeter", false);
        TS_ASSERT_THROWS_THIS(fibre_generator.WriteData(handler, "file.fibres"),
                "Files defining the heart surfaces not set");

        // Wrong surface filename
        TS_ASSERT_THROWS_THIS(fibre_generator.SetSurfaceFiles("wrong_name", "wrong_name", "wrong_name", false),
                "Wrong surface definition file name wrong_name");

        std::string epi_face_file = "heart/test/data/box_shaped_heart/epi.tri";
        std::string rv_face_file = "heart/test/data/box_shaped_heart/rv.tri";
        std::string lv_face_file = "heart/test/data/box_shaped_heart/lv.tri";


        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        OutputFileHandler shorter_handler("shorter_streeter", false);
        TS_ASSERT_THROWS_THIS(fibre_generator.WriteData(shorter_handler, "vector_not_set.ortho"),
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

    void TestConstructStreeterOnLeftWedge()
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
        fibre_generator.SetWriteFileAsBinary();

        OutputFileHandler handler("human_wedge_mesh", false);

        fibre_generator.WriteData(handler, "HumanWedgeMeshLeft.ortho");

        FileFinder fibre_file1 = handler.FindFile("HumanWedgeMeshLeft.ortho");
        FileFinder fibre_file2("heart/test/data/human_wedge_mesh/HumanWedgeMeshLeft.ortho", RelativeTo::ChasteSourceRoot);

        CompareGeneratedWithReferenceFile(fibre_file1, ORTHO, fibre_file2, ORTHO);
    }

    void TestConstructStreeterOnRightWedge()
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
        fibre_generator.SetWriteFileAsBinary();

        OutputFileHandler handler("human_wedge_mesh", false);

        fibre_generator.WriteData(handler, "HumanWedgeMeshRight.ortho");

        FileFinder fibre_file1 = handler.FindFile("HumanWedgeMeshRight.ortho");
        FileFinder fibre_file2("heart/test/data/human_wedge_mesh/HumanWedgeMeshRight.ortho", RelativeTo::ChasteSourceRoot);

        CompareGeneratedWithReferenceFile(fibre_file1, ORTHO, fibre_file2, ORTHO);
    }

    void TestSetLogInfo()
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

        OutputFileHandler handler("shorter_streeter_loginfo");

        fibre_generator.WriteData(handler, "box_heart.ortho");

        FileFinder node_regions_file = handler.FindFile("node_regions.data");
        FileFinder wall_thickness_file = handler.FindFile("wall_thickness.data");
        FileFinder averaged_thickness_file = handler.FindFile("averaged_thickness.data");

        TS_ASSERT_EQUALS(node_regions_file.IsFile(), false);
        TS_ASSERT_EQUALS(wall_thickness_file.IsFile(), false);
        TS_ASSERT_EQUALS(averaged_thickness_file.IsFile(), false);

        // Make the above parallel tests finish before the next part.
        // (When master process enters PreWriteCalculations it will create empty files for the node_regions etc.)
        PetscTools::Barrier("Wait for tests that files don't exist yet");
        fibre_generator.SetLogInfo(true);
        fibre_generator.WriteData(handler, "box_heart.ortho");

        TS_ASSERT_EQUALS(node_regions_file.IsFile(), true);
        TS_ASSERT_EQUALS(wall_thickness_file.IsFile(), true);
        TS_ASSERT_EQUALS(averaged_thickness_file.IsFile(), true);
    }
};

#endif /*TESTSTREETERFIBREGENERATOR_HPP_*/
