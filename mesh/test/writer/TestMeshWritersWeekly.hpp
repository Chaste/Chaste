/*

Copyright (c) 2005-2013, University of Oxford.
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


#ifndef _TESTMESHWRITERSWEEKLY_HPP_
#define _TESTMESHWRITERSWEEKLY_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "MathsCustomFunctions.hpp"

#include "PetscSetupAndFinalize.hpp"


class TestMeshWritersWeekly : public CxxTest::TestSuite
{
private:
    void RunTest(unsigned numStepsInEachDimension)
    {
        MeshEventHandler::Reset();
        unsigned nodes = (unsigned) SmallPow(numStepsInEachDimension+1, 3);
        if (PetscTools::AmMaster())
        {
            std::cout<<"Number steps per dimension = " << std::setw(12) << numStepsInEachDimension<<std::endl;
            std::cout<<"Number nodes               = " << std::setw(12) << nodes<<std::endl;
        }
        std::string file_name = "cuboid";
        std::stringstream directory_stream;
        directory_stream << "TestMeshWritersWeekly/steps_" << std::setw(12) << std::setfill('0')<< numStepsInEachDimension;
        std::string directory_name = directory_stream.str();

        DistributedTetrahedralMesh<3,3> cuboid_mesh;
        cuboid_mesh.ConstructCuboid(numStepsInEachDimension, numStepsInEachDimension, numStepsInEachDimension);
        TS_ASSERT_EQUALS(cuboid_mesh.GetNumNodes(), nodes);


        MeshEventHandler::BeginEvent(MeshEventHandler::TRIANGLES);
        {
            TrianglesMeshWriter<3,3> mesh_writer(directory_name,  file_name, false);
            mesh_writer.WriteFilesUsingMesh(cuboid_mesh);
        }
        MeshEventHandler::EndEvent(MeshEventHandler::TRIANGLES);

//        MeshEventHandler::BeginEvent(MeshEventHandler::BINTRI);
//        {
//            TrianglesMeshWriter<3,3> bin_mesh_writer(directory_name,  file_name+"_bin", false);
//            bin_mesh_writer.SetWriteFilesAsBinary();
//            bin_mesh_writer.WriteFilesUsingMesh(cuboid_mesh);
//        }
//        MeshEventHandler::EndEvent(MeshEventHandler::BINTRI);


#ifdef CHASTE_VTK
//        MeshEventHandler::BeginEvent(MeshEventHandler::VTK);
//        {
//            VtkMeshWriter<3,3> vtk_writer(directory_name, file_name, false);
//            vtk_writer.WriteFilesUsingMesh(cuboid_mesh);
//        }
//        MeshEventHandler::EndEvent(MeshEventHandler::VTK);

//        MeshEventHandler::BeginEvent(MeshEventHandler::PVTK);
//        {
//            VtkMeshWriter<3,3> parallel_vtk_writer(directory_name, file_name+"par",false);
//            parallel_vtk_writer.SetParallelFiles(cuboid_mesh);
//            std::vector<double> dummy_data(cuboid_mesh.GetNumLocalNodes(), PetscTools::GetMyRank());
//            parallel_vtk_writer.AddPointData("Process", dummy_data);
//            parallel_vtk_writer.WriteFilesUsingMesh(cuboid_mesh);
//        }
//        MeshEventHandler::EndEvent(MeshEventHandler::PVTK);
#endif //CHASTE_VTK
        MeshEventHandler::Headings();
        MeshEventHandler::Report();
    }

public:


    void TestRunAllTests()
    {
        //for (unsigned pow = 0; pow<6; pow++)
        unsigned pow = 5;
        {
            if (PetscTools::AmMaster())
            {
                std::cout<<"Power                      = " << std::setw(12) << pow<<std::endl;
            }
            unsigned steps = (unsigned) SmallPow(2, pow);
            RunTest(steps);
        }
        /* See #2351
         * Rough times (sequential, 2-way, 3-way, 4-way GccOpt)
         *
            pow     nodes   time    time2   time3   time4
            4       4K      1s      1s      8s      16s
            5       36K     5s      22s     >5m     >5
            6       280K    45s     >5m
            7       2M      (uses more than 4Gb RAM)
         */
        /* Timing of pow=5 test for IntelProduction:
         *
        Power                      =            5
        Number steps per dimension =           32
        Number nodes               =        35937
              Construct        Tri write        Bin write        VTK write       PVTK write            Total
           1.794 ( 38%)     0.303 (  6%)     0.152 (  3%)     1.196 ( 25%)     1.330 ( 28%)     4.775 (100%)  (seconds)

        Proc       Construct        Tri write        Bin write        VTK write       PVTK write            Total
          0:    0.959 (  5%)     5.778 ( 29%)     5.604 ( 28%)     6.886 ( 35%)     0.727 (  4%)    19.955 (100%)  (seconds)
          1:    0.929 (  5%)     5.809 ( 29%)     5.604 ( 28%)     6.886 ( 35%)     0.559 (  3%)    19.956 (100%)  (seconds)
        avg:    0.944 (  5%)     5.794 ( 29%)     5.604 ( 28%)     6.886 ( 35%)     0.643 (  3%)    19.955 (100%)  (seconds)
        max:    0.959 (  5%)     5.809 ( 29%)     5.604 ( 28%)     6.886 ( 35%)     0.727 (  4%)    19.956 (100%)  (seconds)

        Proc       Construct        Tri write        Bin write        VTK write       PVTK write            Total
          0:    0.508 (  0%)   252.454 ( 31%)   231.030 ( 28%)   328.700 ( 40%)     0.314 (  0%)   813.054 (100%)  (seconds)
          1:    0.507 (  0%)   252.455 ( 31%)   231.030 ( 28%)   328.700 ( 40%)     0.362 (  0%)   813.054 (100%)  (seconds)
          2:    0.442 (  0%)   252.520 ( 31%)   231.030 ( 28%)   328.700 ( 40%)     0.315 (  0%)   813.054 (100%)  (seconds)
          3:    0.399 (  0%)   252.564 ( 31%)   231.030 ( 28%)   328.700 ( 40%)     0.281 (  0%)   813.054 (100%)  (seconds)
        avg:    0.464 (  0%)   252.498 ( 31%)   231.030 ( 28%)   328.700 ( 40%)     0.318 (  0%)   813.054 (100%)  (seconds)
        max:    0.508 (  0%)   252.564 ( 31%)   231.030 ( 28%)   328.700 ( 40%)     0.362 (  0%)   813.054 (100%)  (seconds)

        Proc       Construct        Tri write        Bin write        VTK write       PVTK write            Total
          0:    0.255 (  0%)   727.416 ( 25%)  1109.748 ( 38%)  1094.729 ( 37%)     0.190 (  0%)  2932.420 (100%)  (seconds)
          1:    0.254 (  0%)   727.416 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.273 (  0%)  2932.419 (100%)  (seconds)
          2:    0.254 (  0%)   727.417 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.198 (  0%)  2932.420 (100%)  (seconds)
          3:    0.254 (  0%)   727.417 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.236 (  0%)  2932.419 (100%)  (seconds)
          4:    0.255 (  0%)   727.416 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.189 (  0%)  2932.420 (100%)  (seconds)
          5:    0.258 (  0%)   727.413 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.185 (  0%)  2932.419 (100%)  (seconds)
          6:    0.253 (  0%)   727.418 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.185 (  0%)  2932.420 (100%)  (seconds)
          7:    0.209 (  0%)   727.463 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.149 (  0%)  2932.420 (100%)  (seconds)
        avg:    0.249 (  0%)   727.422 ( 25%)  1109.748 ( 38%)  1094.728 ( 37%)     0.201 (  0%)  2932.420 (100%)  (seconds)
        max:    0.258 (  0%)   727.463 ( 25%)  1109.748 ( 38%)  1094.729 ( 37%)     0.273 (  0%)  2932.420 (100%)  (seconds)
         *
         */
        // Changes at r18242
        /* Using finer timings and barrier synchronisation
         *
      Number nodes               =        35937
      Tri write       node write        ele write       face write            spare            Total
   0.299 (100%)     0.086 ( 29%)     0.202 ( 67%)     0.011 (  4%)     0.000 (  0%)     0.299 (100%)  (seconds)

Proc       Tri write       node write        ele write       face write            spare            Total
  0:    1.421 (100%)     0.095 (  7%)     1.254 ( 88%)     0.071 (  5%)     0.000 (  0%)     1.421 (100%)  (seconds)
  1:    1.480 (100%)     0.095 (  6%)     1.254 ( 85%)     0.071 (  5%)     0.000 (  0%)     1.480 (100%)  (seconds)
avg:    1.450 (100%)     0.095 (  7%)     1.254 ( 86%)     0.071 (  5%)     0.000 (  0%)     1.450 (100%)  (seconds)
max:    1.480 (100%)     0.095 (  6%)     1.254 ( 85%)     0.071 (  5%)     0.000 (  0%)     1.480 (100%)  (seconds)

Proc       Tri write       node write        ele write       face write            spare            Total
  0:  242.766 (100%)     1.768 (  1%)   240.746 ( 99%)     0.252 (  0%)     0.000 (  0%)   242.766 (100%)  (seconds)
  1:  242.772 (100%)     1.768 (  1%)   240.746 ( 99%)     0.252 (  0%)     0.000 (  0%)   242.772 (100%)  (seconds)
  2:  242.774 (100%)     1.768 (  1%)   240.746 ( 99%)     0.252 (  0%)     0.000 (  0%)   242.774 (100%)  (seconds)
  3:  242.815 (100%)     1.768 (  1%)   240.746 ( 99%)     0.252 (  0%)     0.000 (  0%)   242.815 (100%)  (seconds)
avg:  242.782 (100%)     1.768 (  1%)   240.746 ( 99%)     0.252 (  0%)     0.000 (  0%)   242.782 (100%)  (seconds)
max:  242.815 (100%)     1.768 (  1%)   240.746 ( 99%)     0.252 (  0%)     0.000 (  0%)   242.815 (100%)  (seconds)

Proc       Tri write       node write        ele write       face write            spare            Total
  0:  903.342 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.342 (100%)  (seconds)
  1:  903.348 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.348 (100%)  (seconds)
  2:  903.347 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.347 (100%)  (seconds)
  3:  903.349 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.349 (100%)  (seconds)
  4:  903.346 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.346 (100%)  (seconds)
  5:  903.350 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.350 (100%)  (seconds)
  6:  903.344 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.344 (100%)  (seconds)
  7:  903.394 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.394 (100%)  (seconds)
avg:  903.352 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.352 (100%)  (seconds)
max:  903.394 (100%)     2.740 (  0%)   900.387 (100%)     0.214 (  0%)     0.000 (  0%)   903.394 (100%)  (seconds)


         *
         */

    }

};

#endif //_TESTMESHWRITERSWEEKLY_HPP_
