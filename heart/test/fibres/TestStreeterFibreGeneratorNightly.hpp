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
#ifndef TESTSTREETERFIBREGENERATORNIGHTLY_HPP_
#define TESTSTREETERFIBREGENERATORNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "NumericFileComparison.hpp"
#include "VtkMeshWriter.hpp"
#include "FibreReader.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestStreeterFibreGeneratorNightly : public CxxTest::TestSuite
{
private:
    void CompareOrthoFiles(std::string orthoFile1Absolute, std::string orthoFile2Relative)
    {
        //Read one
        FileFinder file_finder1(orthoFile1Absolute, RelativeTo::Absolute);
        FibreReader<3> fibre_reader1(file_finder1, ORTHO);
        std::vector< c_vector<double, 3> > fibres1;
        std::vector< c_vector<double, 3> > second1;
        std::vector< c_vector<double, 3> > third1;
        fibre_reader1.GetAllOrtho(fibres1, second1, third1);
        TS_ASSERT_EQUALS(fibres1.size(), second1.size());
        TS_ASSERT_EQUALS(fibres1.size(), third1.size());

        //Read two
        FileFinder file_finder2(orthoFile2Relative, RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader2(file_finder2, ORTHO);
        std::vector< c_vector<double, 3> > fibres2;
        std::vector< c_vector<double, 3> > second2;
        std::vector< c_vector<double, 3> > third2;
        fibre_reader2.GetAllOrtho(fibres2, second2, third2);
        TS_ASSERT_EQUALS(fibres2.size(), second2.size());
        TS_ASSERT_EQUALS(fibres2.size(), third2.size());

        //Compare data
        TS_ASSERT_EQUALS(fibres1.size(), fibres2.size());

        for (unsigned i = 0u; i< fibres1.size(); i++)
        {
            for (unsigned j = 0; j< 3u ; j++)
            {
                TS_ASSERT_DELTA(fibres1[i][j], fibres2[i][j], 1e-10);
                TS_ASSERT_DELTA(second1[i][j], second2[i][j], 1e-10);
                TS_ASSERT_DELTA(third1[i][j], third2[i][j], 1e-10);
            }
        }
    }

public:
    void TestSimpleOrthotropic()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50_bin");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";

        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);
        fibre_generator.SetLogInfo(true);

        OutputFileHandler handler("streeter_parallel", false);
        fibre_generator.WriteData(handler, "point50.ortho");

        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "point50.ortho";
        std::string wall_file = handler.GetOutputDirectoryFullPath() + "wall_thickness.data";

        CompareOrthoFiles(fibre_file, "heart/test/data/point50_heart_mesh/point50_bin.ortho");

        NumericFileComparison comp_wall(wall_file, "heart/test/data/point50_heart_mesh/wall_thickness.data");
        TS_ASSERT(comp_wall.CompareFiles(1e-11));
    }

    void TestSimpleOrthotropicNotDistributed()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50_bin");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);
        fibre_generator.SetLogInfo(true);

        OutputFileHandler handler("streeter", false);

        fibre_generator.WriteData(handler, "point50_not_dist.ortho");
        fibre_generator.SetWriteFileAsBinary();

        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "point50_not_dist.ortho";
        std::string wall_file = handler.GetOutputDirectoryFullPath() + "wall_thickness.data";

        CompareOrthoFiles(fibre_file, "heart/test/data/point50_heart_mesh/point50_bin.ortho");
        NumericFileComparison comp_wall(wall_file, "heart/test/data/point50_heart_mesh/wall_thickness.data");
        TS_ASSERT(comp_wall.CompareFiles(1e-11));
    }

    void TestDownSampledRabbit()
    {

        TrianglesMeshReader<3,3> mesh_reader("apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um");
        std::string epi_face_file = "apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um.epi";
        std::string rv_face_file = "apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um.rv";
        std::string lv_face_file = "apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um.lv";
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, true);
        fibre_generator.SetApexToBase(2);
        fibre_generator.SetLogInfo(true);
        fibre_generator.SetWriteFileAsBinary();

        OutputFileHandler handler("streeter", false);
        fibre_generator.WriteData(handler, "downsampled.ortho");


        //CompareOrthoFiles(fibre_file, "heart/test/data/fibre_tests/downsampled.ortho");
        // These elements came from random.randint(0, 431989)

        FileFinder file_finder( handler.GetOutputDirectoryFullPath() + "downsampled.ortho", RelativeTo::Absolute);
        FibreReader<3> fibre_reader1(file_finder, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader1.GetAllOrtho(fibres, second, third);
        TS_ASSERT_EQUALS(fibres.size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(second.size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(third.size(), mesh.GetNumElements());
        double tolerance = 1e-4;

        TS_ASSERT_DELTA(fibres[136708][0],  0.3942, tolerance);
        TS_ASSERT_DELTA(fibres[136708][1], -0.9189, tolerance);
        TS_ASSERT_DELTA(fibres[136708][2], -0.0011, tolerance);
        TS_ASSERT_DELTA(fibres[73380][0],  -0.2397, tolerance);
        TS_ASSERT_DELTA(fibres[73380][1],   0.7673, tolerance);
        TS_ASSERT_DELTA(fibres[73380][2],   0.5947, tolerance);
        TS_ASSERT_DELTA(fibres[127320][0], -0.8803, tolerance);
        TS_ASSERT_DELTA(fibres[127320][1], -0.4742, tolerance);
        TS_ASSERT_DELTA(fibres[127320][2],  0.0023, tolerance);


#ifdef CHASTE_VTK
        //Output to VTK.
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "downsampled_fibres", false);
        writer.AddCellData("OrthoFibres", fibres);
        writer.AddCellData("OrthoSecond", second);
        writer.AddCellData("OrthoThird", third);

        //Add debugging data
        std::string debug_files[3] = {"wall_thickness", "node_regions", "averaged_thickness"};
        for (unsigned file=0; file<3; file++)
        {
            std::ifstream ifs((handler.GetOutputDirectoryFullPath()+debug_files[file]+".data").c_str());
            TS_ASSERT(ifs.is_open());
            std::vector<double> payload;
            payload.reserve(mesh.GetNumNodes());
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double temp;
                ifs >> temp;
                payload.push_back(temp);
            }
            writer.AddPointData(debug_files[file], payload);
            ifs.close();
        }
        writer.WriteFilesUsingMesh(mesh);
#else
        std::cout << "This test ran, but did not test VTK-dependent methods." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on Chaste VTK support." << std::endl;
#endif //CHASTE_VTK
    }
};

#endif /*TESTSTREETERFIBREGENERATORNIGHTLY_HPP_*/
