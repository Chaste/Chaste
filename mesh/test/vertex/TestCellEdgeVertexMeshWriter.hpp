//
// Created by twin on 15/07/19.
//

#ifndef TESTCELLEDGEVERTEXMESHWRITER_HPP_
#define TESTCELLEDGEVERTEXMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "MutableMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "CellEdgeVertexMeshWriter.hpp"
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

class TestCellEdgeVertexMeshWriter : public CxxTest::TestSuite
{
public:

    void TestCellEdgeVertexMeshWriterWriteToFile()
    {
        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(2, 1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* We then write to file */
        CellEdgeVertexMeshWriter<2, 2> mesh_writer("", "edges-results", false);
        mesh_writer.WriteVtkUsingMesh(*p_mesh, "test");

    }

};

#endif /* TESTCELLEDGEVERTEXMESHWRITER_HPP_ */
