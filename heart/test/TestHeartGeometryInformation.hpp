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
#ifndef TESTHEARTGEOMETRYINFORMATION_HPP_
#define TESTHEARTGEOMETRYINFORMATION_HPP_

#include "TrianglesMeshReader.hpp"
#include "HeartGeometryInformation.hpp"
#include "PetscTools.hpp"
#include "OutputFileHandler.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
class TestHeartGeometryInformation : public CxxTest::TestSuite
{
private:
    /* Helper method to write out 'fake' face files given the vectors containing the nodes at the surface*/
    void WriteFakeFaceFile(const std::string& rOutputDir, const std::string& rFilename,
                           const std::vector<unsigned>& rNodeLayers, unsigned spaceDim)
    {
        OutputFileHandler handler(rOutputDir, false);

        if (PetscTools::AmMaster())
        {
            out_stream p_file = handler.OpenOutputFile(rFilename);
            p_file->close();
        }
        //Each process may have a small number of the nodes
        for (unsigned proc_turn=0; proc_turn<PetscTools::GetNumProcs(); proc_turn++)
        {
            if (PetscTools::GetMyRank()==proc_turn)
            {
                out_stream p_file = handler.OpenOutputFile(rFilename, std::ios::app);

                for (unsigned i=0; i<rNodeLayers.size(); i++)
                {
                    for (unsigned j=0; j<spaceDim; j++)
                    {
                       * p_file << (j == 0 ? "" : "  ") << (1+rNodeLayers[i]);
                    }
                   * p_file << std::endl;
                }
                p_file->close();
            }
            PetscTools::Barrier();
        }
    }

public:
    void TestCalculateRelativeWallPositionSimple2dMesh()
    {
        DistributedTetrahedralMesh<2,2> mesh;
        //This mesh will have 6 nodes per face, spaced by 1
        mesh.ConstructRectangularMesh(5, 5);

        std::vector<unsigned> left_face;
        std::vector<unsigned> right_face;

        unsigned low_index=mesh.GetDistributedVectorFactory()->GetLow();
        unsigned high_index=mesh.GetDistributedVectorFactory()->GetHigh();

        for (unsigned index=low_index; index<high_index; index++)
        {
            // Get the nodes at the left face of the square
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
            {
                left_face.push_back(index);
            }
            // Get the nodes at the right face of the square
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                right_face.push_back(index);
            }

        }
        // Write our fake face files
        std::string output_dir = "HeartGeom2d";
        WriteFakeFaceFile(output_dir, "epi.tri", left_face, 2u);
        WriteFakeFaceFile(output_dir, "endo.tri", right_face, 2u);

        PetscTools::Barrier(); // Make sure files are written

        // Read in
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        //call the constructor that takes in the surface files...
        HeartGeometryInformation<2> info(mesh, dir_path + "/epi.tri", dir_path + "/endo.tri", false);


        //and then we test the method to evaluate the position in the wall (again)
        for (unsigned index=low_index; index<high_index; index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(5.0-x)/5.0);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEpicardium()[index],x);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEndocardium()[index],(5.0-x));
        }
    }

    void TestCalculateRelativeWallPositionSimple3dMesh()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        //This mesh will have 6 nodes per face, spaced by 1
        mesh.ConstructCuboid(5, 5, 5);

        std::vector<unsigned> left_face;
        std::vector<unsigned> right_face;

        unsigned low_index=mesh.GetDistributedVectorFactory()->GetLow();
        unsigned high_index=mesh.GetDistributedVectorFactory()->GetHigh();
        for (unsigned index=low_index; index<high_index; index++)
        {
            // Get the nodes at the left face of the cube
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
            {
                left_face.push_back(index);
            }
            // Get the nodes at the right face of the cube
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                right_face.push_back(index);
            }

        }
        std::string output_dir = "HeartGeom3d";
        WriteFakeFaceFile(output_dir, "epi.tri", left_face, 3u);
        WriteFakeFaceFile(output_dir, "endo.tri", right_face, 3u);
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        //call the constructor that takes in the surface files...
        HeartGeometryInformation<3> info(mesh, dir_path + "/epi.tri", dir_path + "/endo.tri", false);

        //Check get methods
        std::vector<unsigned> nodes_on_endo = info.rGetNodesOnEndoSurface();
        TS_ASSERT_EQUALS(nodes_on_endo.size(),  36u);
        std::vector<unsigned> nodes_on_epi = info.rGetNodesOnEpiSurface();
        TS_ASSERT_EQUALS(nodes_on_epi.size(),  36u);

        for (unsigned index=low_index; index<high_index; index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            TS_ASSERT_EQUALS(info.CalculateRelativeWallPosition(index),(5-x)/5);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEpicardium()[index],x);
            TS_ASSERT_EQUALS(info.rGetDistanceMapEndocardium()[index],(5.0-x));
        }
        ChasteCuboid<3> epi_bounding_box=info.CalculateBoundingBoxOfEpi();
        ChasteCuboid<3> endo_bounding_box=info.CalculateBoundingBoxOfEndo();
        for (unsigned i=0; i<3; i++)
        {
            if (i==0)
            {
                TS_ASSERT_DELTA(epi_bounding_box.rGetUpperCorner()[i], 0.0, 1e-10);
                TS_ASSERT_DELTA(endo_bounding_box.rGetLowerCorner()[i], 5.0, 1e-10);
            }
            else
            {
                TS_ASSERT_DELTA(epi_bounding_box.rGetUpperCorner()[i], 5.0, 1e-10);
                TS_ASSERT_DELTA(endo_bounding_box.rGetLowerCorner()[i], 0.0, 1e-10);
            }
            TS_ASSERT_DELTA(epi_bounding_box.rGetLowerCorner()[i], 0.0, 1e-10);
            TS_ASSERT_DELTA(endo_bounding_box.rGetUpperCorner()[i], 5.0, 1e-10);
        }
    }

    void TestCalculateRelativeWallPositionWithThreeSurfaces()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        //This mesh will have 9 nodes per side, spaced by 1, it is a cube
        mesh.ConstructCuboid(8, 8, 8);

        std::vector<unsigned> epi_face;
        std::vector<unsigned> lv_face;
        std::vector<unsigned> rv_face;
        //Define three surfaces, epi, lv and rv.
        unsigned low_index=mesh.GetDistributedVectorFactory()->GetLow();
        unsigned high_index=mesh.GetDistributedVectorFactory()->GetHigh();
        for (unsigned index=low_index; index<high_index; index++)
        {
            // Get the nodes at cube face considered to be epi (at both external faces)
            if ((fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
                ||(fabs(mesh.GetNode(index)->rGetLocation()[0]-8.0) < 1e-6))
            {
                epi_face.push_back(index);
            }
            // Get the nodes at cube face considered to be lv (at the plane defined by x=3)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-3.0) < 1e-6)
            {
                lv_face.push_back(index);
            }
            // Get the nodes at cube face considered to be rv (at the plane defined by x=5)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-5.0) < 1e-6)
            {
                rv_face.push_back(index);
            }

        }


        // Write our fake face files
        std::string output_dir = "HeartGeom3d";
        WriteFakeFaceFile(output_dir, "epi.tri", epi_face, 3u);
        WriteFakeFaceFile(output_dir, "lv.tri", lv_face, 3u);
        WriteFakeFaceFile(output_dir, "rv.tri", rv_face, 3u);

        PetscTools::Barrier(); // Make sure files are written

        // Read in
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        HeartGeometryInformation<3> info(mesh, dir_path + "/epi.tri", dir_path + "/lv.tri", dir_path + "/rv.tri", false);

        //first we test the get methods for the nodes on the surfaces
        std::vector<unsigned> nodes_on_lv = info.rGetNodesOnLVSurface();
        std::vector<unsigned> nodes_on_rv = info.rGetNodesOnRVSurface();
        std::vector<unsigned> nodes_on_epi = info.rGetNodesOnEpiSurface();
        TS_ASSERT_EQUALS(nodes_on_lv.size(),  81u);
        TS_ASSERT_EQUALS(nodes_on_rv.size(), 81u);
        TS_ASSERT_EQUALS(nodes_on_epi.size(), 162u);

        //and then we test the method to evaluate the position in the wall (again)
        for (unsigned index=low_index; index<high_index; index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];
            //in the lv...
            if (x<=3)
            {
                TS_ASSERT_DELTA(info.CalculateRelativeWallPosition(index),(3.0-x)/3.0,1e-12);
            }
            //..in the septum...
            else if ((x>3)&&(x<5))
            {
                TS_ASSERT_DELTA(info.CalculateRelativeWallPosition(index), 1.0/5.0, 1e-12);
            }
            //...and in the rv.
            else if (x>=5)
            {
                TS_ASSERT_DELTA(info.CalculateRelativeWallPosition(index),(x-5.0)/3.0,1e-12);
            }
        }

        ChasteCuboid<3> epi_bounding_box=info.CalculateBoundingBoxOfEpi();
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(epi_bounding_box.rGetUpperCorner()[i], 8.0, 1e-10);
            TS_ASSERT_DELTA(epi_bounding_box.rGetLowerCorner()[i], 0.0, 1e-10);
        }
        ChasteCuboid<3> lv_bounding_box=info.CalculateBoundingBoxOfLV();
        TS_ASSERT_DELTA(lv_bounding_box.rGetUpperCorner()[0], 3.0, 1e-10);
        TS_ASSERT_DELTA(lv_bounding_box.rGetLowerCorner()[0], 3.0, 1e-10);
        ChasteCuboid<3> rv_bounding_box=info.CalculateBoundingBoxOfRV();
        TS_ASSERT_DELTA(rv_bounding_box.rGetUpperCorner()[0], 5.0, 1e-10);
        TS_ASSERT_DELTA(rv_bounding_box.rGetLowerCorner()[0], 5.0, 1e-10);

    }

    void TestDetermineLayerForEachNodeWritingAndReading()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        //This mesh will have 31 nodes per side, spaced by 1, it is a cube
        mesh.ConstructCuboid(30, 30, 30);

        unsigned low_index=mesh.GetDistributedVectorFactory()->GetLow();
        unsigned high_index=mesh.GetDistributedVectorFactory()->GetHigh();
        std::vector<unsigned> epi_face;
        std::vector<unsigned> lv_face;
        std::vector<unsigned> rv_face;
        //Define three surfaces, epi, lv and rv.
        for (unsigned index=low_index; index<high_index; index++)
        {
            // Get the nodes at cube face considered to be epi (at both external faces)
            if ((fabs(mesh.GetNode(index)->rGetLocation()[0]) < 1e-6)
                ||(fabs(mesh.GetNode(index)->rGetLocation()[0]-30.0) < 1e-6))
            {
                epi_face.push_back(index);
            }
            // Get the nodes at cube face considered to be lv (at the plane defined by x=10)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-10.0) < 1e-6)
            {
                lv_face.push_back(index);
            }
            // Get the nodes at cube face considered to be rv (at the plane defined by x=20)
            if (fabs(mesh.GetNode(index)->rGetLocation()[0]-20.0) < 1e-6)
            {
                rv_face.push_back(index);
            }

        }

        // Write our fake face files
        std::string output_dir = "HeartGeom3d";
        WriteFakeFaceFile(output_dir, "epi.tri", epi_face, 3u);
        WriteFakeFaceFile(output_dir, "lv.tri", lv_face, 3u);
        WriteFakeFaceFile(output_dir, "rv.tri", rv_face, 3u);

        PetscTools::Barrier(); // Make sure files are written

        // Read in
        OutputFileHandler handler(output_dir, false);
        std::string dir_path = handler.GetOutputDirectoryFullPath();
        HeartGeometryInformation<3> info(mesh, dir_path + "/epi.tri", dir_path + "/lv.tri", dir_path + "/rv.tri", false);

        //covering exceptions
        TS_ASSERT_THROWS_THIS(info.DetermineLayerForEachNode(0.9, 0.9), "The sum of fractions of epicardial and endocardial layers must be lesser than 1");
        TS_ASSERT_THROWS_THIS(info.DetermineLayerForEachNode(0.9, -1.0), "A fraction of a layer must be positive");
        TS_ASSERT_THROWS_THIS(info.DetermineLayerForEachNode(-2.0, 1.0), "A fraction of a layer must be positive");


        info.DetermineLayerForEachNode(0.29, 0.51);
        info.WriteLayerForEachNode("TestHeartGeom","layers.het");

        for (unsigned index=low_index; index<high_index; index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];

            if (x < 3 || x > 27)
            {
                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],EPI);
            }
            else if (x < 5 || x > 25)
            {
                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],MID);
            }
            else
            {
                TS_ASSERT_EQUALS(info.rGetLayerForEachNode()[index],ENDO);
            }
        }

        //now we test the constructor that takes in the node heterogeneity file
        std::string file = OutputFileHandler::GetChasteTestOutputDirectory() + "TestHeartGeom/layers.het";
        HeartGeometryInformation<3> info_2(file);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), info_2.rGetLayerForEachNode().size());

        for (unsigned index=low_index; index<high_index; index++)
        {
            double x = mesh.GetNode(index)->rGetLocation()[0];

            if (x < 3 || x > 27)
            {
                TS_ASSERT_EQUALS(info_2.rGetLayerForEachNode()[index],EPI);
            }
            else if (x < 5 || x > 25)
            {
                TS_ASSERT_EQUALS(info_2.rGetLayerForEachNode()[index],MID);
            }
            else
            {
                TS_ASSERT_EQUALS(info_2.rGetLayerForEachNode()[index],ENDO);
            }
        }

        //covering exceptions
        TS_ASSERT_THROWS_THIS(HeartGeometryInformation<3> info_non_existent_file("rubbish"), "Could not open heterogeneities file (rubbish)");
        TS_ASSERT_THROWS_THIS(HeartGeometryInformation<3> info_bad_file("heart/test/data/ValidPseudoEcg1D.dat"), "A value in the heterogeneities file (heart/test/data/ValidPseudoEcg1D.dat) is out of range (or not an integer). It should be epi = 0, mid = 1, endo = 2");
    }

    void TestHeartGeometryTakingMeshFromFile()
    {
        //files containing list of nodes on each surface
        std::string epi_surface = "heart/test/data/box_shaped_heart/epi.tri";
        std::string lv_surface = "heart/test/data/box_shaped_heart/lv.tri";
        std::string rv_surface = "heart/test/data/box_shaped_heart/rv.tri";
        std::string bad_surface = "heart/test/data/box_shaped_heart/zero.tri";


        //read in the mesh
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
        //DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);
        DistributedTetrahedralMesh<3,3> mesh;
        //TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        //Check the indexing exception
        TS_ASSERT_THROWS_THIS(HeartGeometryInformation<3> info2(mesh, epi_surface, lv_surface, bad_surface, false),
                              "Error when reading surface file.  It was assumed not to be indexed from zero, but zero appeared in the list.");

        //calculate the geometry information
        HeartGeometryInformation<3> info(mesh, epi_surface, lv_surface, rv_surface, false);
        info.DetermineLayerForEachNode(0.25,0.375);
        //and write them out to file
        OutputFileHandler results_handler("BoxShaped", false);


        std::vector<double> relative_wall_position;
        for (unsigned index=0; index<mesh.GetNumNodes(); index++)
        {
            relative_wall_position.push_back(info.CalculateRelativeWallPosition(index));
        }

//        //This only needs to be done when we are regenerating the output
//        if (PetscTools::IsSequential())
//        {
//            out_stream p_file = results_handler.OpenOutputFile("heart_geometry_layers.dat");
//
//            for (unsigned index=0; index<mesh.GetNumNodes(); index++)
//            {
//                (*p_file)<<info.rGetLayerForEachNode()[index]<<std::endl;
//            }
//            p_file->close();
//        }
//
//        //This only needs to be done when we are regenerating the output
//        if (PetscTools::IsSequential())
//        {
//            out_stream p_file = results_handler.OpenOutputFile("heart_relative_wall_position.dat");
//            *p_file << std::setprecision(20);//Slightly more than machine precision...
//            for (unsigned index=0; index<mesh.GetNumNodes(); index++)
//            {
//                (*p_file)<<relative_wall_position[index]<<std::endl;
//            }
//            p_file->close();
//        }


        //This is the actual test - read the file into a vector
        std::vector<HeartLayerType> sequential_layers;
        {
            std::ifstream file_stream("heart/test/data/heart_geometry_layers.dat");
            unsigned layer;
            while (file_stream >> layer)
            {
                sequential_layers.push_back((HeartLayerType)layer);
            }
            file_stream.close();
        }
        std::vector<double> sequential_relative_wall_position;
        {
            std::ifstream file_stream("heart/test/data/heart_relative_wall_position.dat");
            double dist;
            while (file_stream >> dist)
            {
                sequential_relative_wall_position.push_back(dist);
            }
            file_stream.close();
        }

        TS_ASSERT_EQUALS(sequential_layers.size(), info.rGetLayerForEachNode().size());
        TS_ASSERT_EQUALS(sequential_layers.size(), sequential_relative_wall_position.size());

//#ifdef CHASTE_VTK
//// Requires  "sudo aptitude install libvtk5-dev" or similar
//        if  (PetscTools::IsParallel())
//        {
//            VtkMeshWriter<3,3> writer("", "epi_distance_par", false);
//            // Add distance from origin into the node "point" data
//            writer.AddPointData("Distance from epi", info.rGetDistanceMapEpicardium());
//            writer.AddPointData("Relative wall distance", relative_wall_position);
//            std::vector<double> errors;
//            for (unsigned i=0;i<relative_wall_position.size();i++)
//            {
//                errors.push_back(sequential_relative_wall_position[i] - relative_wall_position[mesh.rGetNodePermutation()[i]]);
//            }
//            writer.AddPointData("Error compared to sequential", errors);
//            writer.WriteFilesUsingMesh(mesh);
//        }
//#endif //CHASTE_VTK


        if (PetscTools::IsSequential())
        {
            TS_ASSERT(mesh.rGetNodePermutation().empty());
            for (unsigned i=0;i<sequential_layers.size();i++)
            {
                TS_ASSERT_DELTA(sequential_relative_wall_position[i], relative_wall_position[i], 1e-15);
                TS_ASSERT_EQUALS(sequential_layers[ i ], info.rGetLayerForEachNode()[ i ]);
            }
        }
        else
        {
            //In parallel we need to apply the permutation to original data
            //Node i in the original data has moved to index mesh.rGetNodePermutation()[i]
            //therefore we will find its data in rGetLayerForEachNode()[ permutation...]

            TS_ASSERT(mesh.rGetNodePermutation().empty() == false );
            for (unsigned i=0;i<sequential_layers.size();i++)
            {
                TS_ASSERT_DELTA(sequential_relative_wall_position[i], relative_wall_position[mesh.rGetNodePermutation()[i]], 1e-15);
                TS_ASSERT_EQUALS(sequential_layers[ i ], info.rGetLayerForEachNode()[ mesh.rGetNodePermutation()[i] ]);
            }
        }
    }


    void TestHeartGeometryOnlyLeftVentricle()
    {
         //files containing list of nodes on each surface
         std::string epi_surface = "heart/test/data/box_shaped_heart/epi.tri";
         std::string lv_surface = "heart/test/data/box_shaped_heart/lv.tri";

         //read in the mesh
         TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart");
         DistributedTetrahedralMesh<3,3> mesh;
         mesh.ConstructFromMeshReader(mesh_reader);


         TS_ASSERT_THROWS_THIS(HeartGeometryInformation<3> no_info(mesh, epi_surface, "", "", false),
                                  "At least one of left ventricle or right ventricle files is required");


         //calculate the geometry information
         HeartGeometryInformation<3> info(mesh, epi_surface, lv_surface, "", false);

         unsigned low_index=mesh.GetDistributedVectorFactory()->GetLow();
         unsigned high_index=mesh.GetDistributedVectorFactory()->GetHigh();
         for (unsigned node_index=low_index; node_index<high_index; node_index++)
         {
             TS_ASSERT_EQUALS(info.GetHeartRegion(node_index), HeartGeometryInformation<3>::LEFT_VENTRICLE_WALL);
         }
    }
};

#endif /*TESTHEARTGEOMETRYINFORMATION_HPP_*/
