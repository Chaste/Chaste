//
// Created by twin on 15/07/19.
//

#ifndef TESTTRIANGLEEDGEVERTEXMESHWRITER_HPP_
#define TESTTRIANGLEEDGEVERTEXMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "MutableMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "TriangleEdgeVertexMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "VertexMeshReader.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

class TestTriangleEdgeVertexMeshWriter : public CxxTest::TestSuite
{
public:

    void TestCellEdgeVertexMeshWriterWriteToFile()
    {
        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(2, 1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* We then write to file */
        TriangleEdgeVertexMeshWriter<2, 2> mesh_writer("", "results", false);
        mesh_writer.WriteVtkUsingMesh(*p_mesh, "0");

        FileFinder vtk_file("results_0.vtu", RelativeTo::ChasteTestOutput);
        TS_ASSERT(vtk_file.Exists());

    }

};

#endif /* TESTTRIANGLEEDGEVERTEXMESHWRITER_HPP_ */
