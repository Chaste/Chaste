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
#include "GenericEventHandler.hpp"
#include "MathsCustomFunctions.hpp"

#include "PetscSetupAndFinalize.hpp"

class MeshEventHandler : public GenericEventHandler<6, MeshEventHandler>
{
public:
    static const char* EventName[6];

    typedef enum
    {
        CONSTR=0,
        TRIANGLES,
        BINTRI,
        VTK,
        PVTK
    } EventType;
};

const char* MeshEventHandler::EventName[] = { "Construct", "Tri write", "Bin write", "VTK write", "PVTK write", "Total"};


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

        MeshEventHandler::BeginEvent(MeshEventHandler::CONSTR);
        DistributedTetrahedralMesh<3,3> cuboid_mesh;
        cuboid_mesh.ConstructCuboid(numStepsInEachDimension, numStepsInEachDimension, numStepsInEachDimension);
        MeshEventHandler::EndEvent(MeshEventHandler::CONSTR);
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
        {
            if (PetscTools::AmMaster())
            {
                std::cout<<"Power                      = " << std::setw(12) << pow<<std::endl;
            }
            unsigned steps = SmallPow(2, pow);
            RunTest(steps);
        }
        /* Rough times (sequential, 2-way, 3-way, 4-way GccOpt)
         *
            pow     nodes   time    time2   time3   time4
            4       4K      1s      1s      8s      16s
            5       36K     5s      22s     >5m     >5
            6       280K    45s     >5m
            7       2M      (uses more than 4Gb RAM)
         */
    }

};

#endif //_TESTMESHWRITERSWEEKLY_HPP_
