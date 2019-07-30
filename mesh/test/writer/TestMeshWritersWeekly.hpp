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
        unsigned nodes = SmallPow(numStepsInEachDimension+1, 3);
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

        MeshEventHandler::BeginEvent(MeshEventHandler::BINTRI);
        {
            TrianglesMeshWriter<3,3> bin_mesh_writer(directory_name,  file_name+"_bin", false);
            bin_mesh_writer.SetWriteFilesAsBinary();
            bin_mesh_writer.WriteFilesUsingMesh(cuboid_mesh);
        }
        MeshEventHandler::EndEvent(MeshEventHandler::BINTRI);


#ifdef CHASTE_VTK
        MeshEventHandler::BeginEvent(MeshEventHandler::VTK);
        {
            VtkMeshWriter<3,3> vtk_writer(directory_name, file_name, false);
            vtk_writer.WriteFilesUsingMesh(cuboid_mesh);
        }
        MeshEventHandler::EndEvent(MeshEventHandler::VTK);

        MeshEventHandler::BeginEvent(MeshEventHandler::PVTK);
        {
            VtkMeshWriter<3,3> parallel_vtk_writer(directory_name, file_name+"par",false);
            parallel_vtk_writer.SetParallelFiles(cuboid_mesh);
            std::vector<double> dummy_data(cuboid_mesh.GetNumLocalNodes(), PetscTools::GetMyRank());
            parallel_vtk_writer.AddPointData("Process", dummy_data);
            parallel_vtk_writer.WriteFilesUsingMesh(cuboid_mesh);
        }
        MeshEventHandler::EndEvent(MeshEventHandler::PVTK);
#endif //CHASTE_VTK
        MeshEventHandler::Headings();
        MeshEventHandler::Report();
    }

public:


    void TestRunAllTests()
    {
        for (unsigned pow = 0; pow<6; pow++)
        //unsigned pow = 5;
        {
            if (PetscTools::AmMaster())
            {
                std::cout<<"Power                      = " << std::setw(12) << pow<<std::endl;
            }
            unsigned steps = SmallPow(2u, pow);
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


        // Changes at r18411 (or r18412 !)
        /* Using synchronised (blocking) sends
         *
Power                      =            5
Number steps per dimension =           32
Number nodes               =        35937
      Tri write     BinTri write        VTK write       PVTK write       node write        ele write       face write        ncl write            comm1            comm2            Total
   0.309 ( 11%)     0.163 (  6%)     1.203 ( 42%)     1.197 ( 42%)     0.090 (  3%)     0.311 ( 11%)     0.018 (  1%)     0.052 (  2%)     0.000 (  0%)     0.000 (  0%)     2.872 (100%)  (seconds)

Proc       Tri write     BinTri write        VTK write       PVTK write       node write        ele write       face write        ncl write            comm1            comm2            Total
  0:    1.525 ( 26%)     1.332 ( 22%)     2.404 ( 41%)     0.673 ( 11%)     0.114 (  2%)     2.534 ( 43%)     0.158 (  3%)     0.050 (  1%)     0.255 (  4%)     0.131 (  2%)     5.932 (100%)  (seconds)
  1:    1.560 ( 26%)     1.332 ( 22%)     2.403 ( 40%)     0.608 ( 10%)     0.129 (  2%)     4.736 ( 79%)     0.231 (  4%)     0.050 (  1%)     4.716 ( 79%)     0.160 (  3%)     5.968 (100%)  (seconds)
avg:    1.542 ( 26%)     1.332 ( 22%)     2.403 ( 40%)     0.640 ( 11%)     0.122 (  2%)     3.635 ( 61%)     0.194 (  3%)     0.050 (  1%)     2.486 ( 42%)     0.146 (  2%)     5.950 (100%)  (seconds)
max:    1.560 ( 26%)     1.332 ( 22%)     2.404 ( 40%)     0.673 ( 11%)     0.129 (  2%)     4.736 ( 79%)     0.231 (  4%)     0.050 (  1%)     4.716 ( 79%)     0.160 (  3%)     5.968 (100%)  (seconds)

Proc       Tri write     BinTri write        VTK write       PVTK write       node write        ele write       face write        ncl write            comm1            comm2            Total
  0:    2.215 ( 28%)     2.029 ( 26%)     3.290 ( 41%)     0.411 (  5%)     0.122 (  2%)     3.850 ( 48%)     0.212 (  3%)     0.060 (  1%)     0.457 (  6%)     0.308 (  4%)     7.945 (100%)  (seconds)
  1:    2.217 ( 28%)     2.029 ( 26%)     3.289 ( 41%)     0.385 (  5%)     0.062 (  1%)     2.630 ( 33%)     4.327 ( 54%)     0.060 (  1%)     6.823 ( 86%)     0.085 (  1%)     7.947 (100%)  (seconds)
  2:    2.219 ( 28%)     2.029 ( 26%)     3.289 ( 41%)     0.358 (  5%)     0.111 (  1%)     4.917 ( 62%)     2.074 ( 26%)     0.060 (  1%)     6.756 ( 85%)     0.187 (  2%)     7.949 (100%)  (seconds)
  3:    2.259 ( 28%)     2.029 ( 25%)     3.289 ( 41%)     0.301 (  4%)     0.146 (  2%)     6.787 ( 85%)     0.312 (  4%)     0.060 (  1%)     6.978 ( 87%)     0.076 (  1%)     7.989 (100%)  (seconds)
avg:    2.228 ( 28%)     2.029 ( 25%)     3.290 ( 41%)     0.364 (  5%)     0.111 (  1%)     4.546 ( 57%)     1.731 ( 22%)     0.060 (  1%)     5.253 ( 66%)     0.164 (  2%)     7.957 (100%)  (seconds)
max:    2.259 ( 28%)     2.029 ( 25%)     3.290 ( 41%)     0.411 (  5%)     0.146 (  2%)     6.787 ( 85%)     4.327 ( 54%)     0.060 (  1%)     6.978 ( 87%)     0.308 (  4%)     7.989 (100%)  (seconds)

Proc       Tri write     BinTri write        VTK write       PVTK write       node write        ele write       face write        ncl write            comm1            comm2            Total
  0:    2.648 ( 28%)     2.496 ( 27%)     3.898 ( 42%)     0.210 (  2%)     0.144 (  2%)     4.682 ( 50%)     0.251 (  3%)     0.065 (  1%)     0.586 (  6%)     0.475 (  5%)     9.312 (100%)  (seconds)
  1:    2.649 ( 28%)     2.496 ( 27%)     3.897 ( 42%)     0.271 (  3%)     0.035 (  0%)     1.490 ( 16%)     6.891 ( 74%)     0.065 (  1%)     8.308 ( 89%)     0.046 (  0%)     9.313 (100%)  (seconds)
  2:    2.652 ( 28%)     2.496 ( 27%)     3.897 ( 42%)     0.246 (  3%)     0.061 (  1%)     2.683 ( 29%)     5.716 ( 61%)     0.065 (  1%)     8.272 ( 89%)     0.100 (  1%)     9.316 (100%)  (seconds)
  3:    2.652 ( 28%)     2.496 ( 27%)     3.897 ( 42%)     0.203 (  2%)     0.087 (  1%)     3.885 ( 42%)     4.531 ( 49%)     0.065 (  1%)     8.291 ( 89%)     0.097 (  1%)     9.316 (100%)  (seconds)
  4:    2.647 ( 28%)     2.496 ( 27%)     3.897 ( 42%)     0.190 (  2%)     0.105 (  1%)     5.005 ( 54%)     3.431 ( 37%)     0.065 (  1%)     8.364 ( 90%)     0.046 (  0%)     9.311 (100%)  (seconds)
  5:    2.646 ( 28%)     2.496 ( 27%)     3.897 ( 42%)     0.192 (  2%)     0.124 (  1%)     6.124 ( 66%)     2.334 ( 25%)     0.065 (  1%)     8.384 ( 90%)     0.046 (  0%)     9.310 (100%)  (seconds)
  6:    2.651 ( 28%)     2.496 ( 27%)     3.897 ( 42%)     0.189 (  2%)     0.150 (  2%)     7.311 ( 78%)     1.164 ( 12%)     0.065 (  1%)     8.352 ( 90%)     0.095 (  1%)     9.315 (100%)  (seconds)
  7:    2.694 ( 29%)     2.496 ( 27%)     3.897 ( 42%)     0.155 (  2%)     0.176 (  2%)     8.193 ( 88%)     0.373 (  4%)     0.065 (  1%)     8.463 ( 90%)     0.080 (  1%)     9.358 (100%)  (seconds)
avg:    2.655 ( 28%)     2.496 ( 27%)     3.897 ( 42%)     0.207 (  2%)     0.110 (  1%)     4.921 ( 53%)     3.087 ( 33%)     0.065 (  1%)     7.378 ( 79%)     0.123 (  1%)     9.319 (100%)  (seconds)
max:    2.694 ( 29%)     2.496 ( 27%)     3.898 ( 42%)     0.271 (  3%)     0.176 (  2%)     8.193 ( 88%)     6.891 ( 74%)     0.065 (  1%)     8.463 ( 90%)     0.475 (  5%)     9.358 (100%)  (seconds)
         *
         */
    }
};

#endif //_TESTMESHWRITERSWEEKLY_HPP_
