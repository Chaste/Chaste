/*

Copyright (c) 2005-2026, University of Oxford.
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
#ifndef TESTSEMMESHWRITER_HPP_
#define TESTSEMMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "SemMeshWriter.hpp"
#include "FileComparison.hpp"
#include "OutputFileHandler.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif

class TestSemMeshWriter : public CxxTest::TestSuite
{
public:

    void TestSemMeshWriterIn2d()
    {
#ifdef CHASTE_VTK
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        std::vector<Node<2>*> element_nodes = nodes;
        std::vector<SemElement<2>*> elements;
        elements.push_back(new SemElement<2>(0, element_nodes));
        SemMesh<2> mesh(nodes, elements);
        TS_ASSERT(mesh.GetOutputElementSurfacesToVtk());

        SemMeshWriter<2> writer("TestSemMeshWriter", "surface_results", false);
        std::vector<double> point_data;
        point_data.push_back(1.0);
        point_data.push_back(2.0);
        point_data.push_back(3.0);
        point_data.push_back(4.0);
        writer.AddPointData("test point data", point_data);
        writer.WriteVtkUsingMesh(mesh);

        std::string results_file = OutputFileHandler::GetChasteTestOutputDirectory()
                                   + "TestSemMeshWriter/surface_results.vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridReader> p_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        p_reader->SetFileName(results_file.c_str());
        p_reader->Update();
        vtkUnstructuredGrid* p_grid = p_reader->GetOutput();

        TS_ASSERT_EQUALS(p_grid->GetNumberOfCells(), 5u);
        TS_ASSERT_EQUALS(p_grid->GetNumberOfPoints(), 8u);

        vtkDataArray* p_output_kind = p_grid->GetCellData()->GetArray("SemOutputKind");
        TS_ASSERT(p_output_kind != nullptr);
        TS_ASSERT_EQUALS(p_output_kind->GetNumberOfTuples(), 5u);
        TS_ASSERT_DELTA(p_output_kind->GetTuple1(0), 0.0, 1e-12);
        for (unsigned i = 1u; i < 5u; ++i)
        {
            TS_ASSERT_DELTA(p_output_kind->GetTuple1(i), 1.0, 1e-12);
        }

        vtkDataArray* p_point_data = p_grid->GetPointData()->GetArray("test point data");
        TS_ASSERT(p_point_data != nullptr);
        TS_ASSERT_EQUALS(p_point_data->GetNumberOfTuples(), 8u);
#else
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        std::vector<Node<2>*> element_nodes = nodes;
        std::vector<SemElement<2>*> elements;
        elements.push_back(new SemElement<2>(0, element_nodes));
        SemMesh<2> mesh(nodes, elements);

        SemMeshWriter<2> writer("TestSemMeshWriter", "surface_results", false);
        std::vector<double> data;
        data.push_back(1.0);

        TS_ASSERT_THROWS_CONTAINS(writer.AddPointData("test point data", data), "requires Chaste to be compiled with VTK");
        TS_ASSERT_THROWS_CONTAINS(writer.AddCellData("test cell data", data), "requires Chaste to be compiled with VTK");
        TS_ASSERT_THROWS_CONTAINS(writer.WriteVtkUsingMesh(mesh), "requires Chaste to be compiled with VTK");
#endif // CHASTE_VTK
    }

    void TestPerElementDataToVtk()
    {
#ifdef CHASTE_VTK
        // Two square elements (four nodes each), offset along x
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 3.0, 0.0));
        nodes.push_back(new Node<2>(5, false, 4.0, 0.0));
        nodes.push_back(new Node<2>(6, false, 4.0, 1.0));
        nodes.push_back(new Node<2>(7, false, 3.0, 1.0));

        std::vector<Node<2>*> element_0_nodes(nodes.begin(), nodes.begin() + 4);
        std::vector<Node<2>*> element_1_nodes(nodes.begin() + 4, nodes.end());
        std::vector<SemElement<2>*> elements;
        elements.push_back(new SemElement<2>(0, element_0_nodes));
        elements.push_back(new SemElement<2>(1, element_1_nodes));
        SemMesh<2> mesh(nodes, elements);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        SemMeshWriter<2> writer("TestSemMeshWriterElementData", "element_data_results", false);
        std::vector<double> element_values;
        element_values.push_back(10.0);
        element_values.push_back(20.0);
        writer.AddElementData("element colour", element_values);
        writer.WriteVtkUsingMesh(mesh);

        std::string results_file = OutputFileHandler::GetChasteTestOutputDirectory()
                                   + "TestSemMeshWriterElementData/element_data_results.vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridReader> p_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        p_reader->SetFileName(results_file.c_str());
        p_reader->Update();
        vtkUnstructuredGrid* p_grid = p_reader->GetOutput();

        vtkDataArray* p_element_index = p_grid->GetCellData()->GetArray("SemElementIndex");
        vtkDataArray* p_colour = p_grid->GetCellData()->GetArray("element colour");
        TS_ASSERT(p_element_index != nullptr);
        TS_ASSERT(p_colour != nullptr);
        TS_ASSERT_EQUALS(p_colour->GetNumberOfTuples(), p_grid->GetNumberOfCells());

        // Every VTK cell (point-cloud and surface) must carry its element's value
        for (unsigned cell = 0u; cell < static_cast<unsigned>(p_grid->GetNumberOfCells()); ++cell)
        {
            unsigned element_index = static_cast<unsigned>(p_element_index->GetTuple1(cell));
            double expected = (element_index == 0u) ? 10.0 : 20.0;
            TS_ASSERT_DELTA(p_colour->GetTuple1(cell), expected, 1e-12);
        }

        // A per-element payload of the wrong length is rejected when the VTK mesh is built
        SemMeshWriter<2> bad_writer("TestSemMeshWriterElementDataBad", "bad_results", false);
        std::vector<double> wrong_length_values;
        wrong_length_values.push_back(1.0);
        bad_writer.AddElementData("bad", wrong_length_values);
        TS_ASSERT_THROWS_CONTAINS(bad_writer.WriteVtkUsingMesh(mesh), "must have one entry per element");
#endif // CHASTE_VTK
    }
};

#endif /*TESTSEMMESHWRITER_HPP_*/
