/*

Copyright (C) University of Oxford, 2005-2011

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
#ifndef TESTSTREETERFIBREGENERATORNIGHTLY_HPP_
#define TESTSTREETERFIBREGENERATORNIGHTLY_HPP_

#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "StreeterFibreGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "MemfemMeshReader.hpp"
#include "NumericFileComparison.hpp"
#include "VtkMeshWriter.hpp"
#include "FibreReader.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestStreeterFibreGeneratorNightly : public CxxTest::TestSuite
{
private:
    void CompareOrthoFiles(std::string orthoFile1Absolute, std::string orthoFile2Relative) throw (Exception)
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

        for (unsigned i = 431980; i< fibres1.size(); i++)
        {
            for (unsigned j = 0; j< 3u ; j++)
            {
                TS_ASSERT_DELTA(fibres1[i][j], fibres2[i][j], 1e-16);
                TS_ASSERT_DELTA(second1[i][j], second2[i][j], 1e-16);
                TS_ASSERT_DELTA(third1[i][j], third2[i][j], 1e-16);
            }
        }
    }



public:
    void TestSimpleOrthotropic() throw (Exception)
    {
        MemfemMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";

        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("streeter_parallel", "point50.ortho", true);

        OutputFileHandler handler("streeter_parallel", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "point50.ortho";
        std::string wall_file = handler.GetOutputDirectoryFullPath() + "wall_thickness.data";

        CompareOrthoFiles(fibre_file, "heart/test/data/point50_heart_mesh/point50.ortho");

        NumericFileComparison comp_wall(wall_file,"heart/test/data/point50_heart_mesh/wall_thickness.data");
        TS_ASSERT(comp_wall.CompareFiles(1e-11));
    }

    void TestSimpleOrthotropicNotDistributed() throw (Exception)
    {
        MemfemMeshReader<3,3> mesh_reader("heart/test/data/point50_heart_mesh/point50");
        std::string epi_face_file = "heart/test/data/point50_heart_mesh/epi.tri";
        std::string rv_face_file = "heart/test/data/point50_heart_mesh/rv.tri";
        std::string lv_face_file = "heart/test/data/point50_heart_mesh/lv.tri";

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, false);
        fibre_generator.SetApexToBase(0);

        fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "point50_not_dist.ortho", true);

        OutputFileHandler handler("streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "point50_not_dist.ortho";
        std::string wall_file = handler.GetOutputDirectoryFullPath() + "wall_thickness.data";

        CompareOrthoFiles(fibre_file, "heart/test/data/point50_heart_mesh/point50.ortho");
        NumericFileComparison comp_wall(wall_file,"heart/test/data/point50_heart_mesh/wall_thickness.data");
        TS_ASSERT(comp_wall.CompareFiles(1e-11));
    }

    void TestDownSampledRabbit() throw (Exception)
    {

        TrianglesMeshReader<3,3> mesh_reader("apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles");
        std::string epi_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.epi";
        std::string rv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.rv";
        std::string lv_face_file = "apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles.lv";
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        StreeterFibreGenerator<3> fibre_generator(mesh);
        fibre_generator.SetSurfaceFiles(epi_face_file, rv_face_file, lv_face_file, true);
        fibre_generator.SetApexToBase(2);

        fibre_generator.GenerateOrthotropicFibreOrientation("streeter", "downsampled.ortho", true);

        OutputFileHandler handler("streeter", false);
        std::string fibre_file = handler.GetOutputDirectoryFullPath() + "downsampled.ortho";

        CompareOrthoFiles(fibre_file, "heart/test/data/fibre_tests/downsampled.ortho");
#ifdef CHASTE_VTK
        //Output to VTK.
        VtkMeshWriter<3,3> writer("TestVtkMeshWriter", "downsampled_fibres", false);
        FileFinder file("heart/test/data/fibre_tests/downsampled.ortho", RelativeTo::ChasteSourceRoot);
        FibreReader<3> fibre_reader(file, ORTHO);
        std::vector< c_vector<double, 3> > fibres;
        std::vector< c_vector<double, 3> > second;
        std::vector< c_vector<double, 3> > third;
        fibre_reader.GetAllOrtho(fibres, second, third);
        TS_ASSERT_EQUALS(fibres.size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(second.size(), mesh.GetNumElements());
        TS_ASSERT_EQUALS(third.size(), mesh.GetNumElements());
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
#endif
    }
};

#endif /*TESTSTREETERFIBREGENERATORNIGHTLY_HPP_*/
