/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "TrapezoidEdgeVertexMeshWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrapezoidEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::TrapezoidEdgeVertexMeshWriter(const std::string& rDirectory,
                                                                                     const std::string& rBaseName,
                                                                                     const bool clearOutputDir)
        : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir)
{
#ifdef CHASTE_VTK
    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
#endif // CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrapezoidEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::~TrapezoidEdgeVertexMeshWriter()
{
#ifdef CHASTE_VTK
    mpVtkUnstructedMesh->Delete();
#endif // CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrapezoidEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                                                              const std::string& stamp)
{
#ifdef CHASTE_VTK
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE

    // Create VTK mesh
    MakeVtkMesh(rMesh);

    // Now write VTK mesh to file
    assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
    vtkXMLUnstructuredGridWriter* p_writer = vtkXMLUnstructuredGridWriter::New();
#if VTK_MAJOR_VERSION >= 6
    p_writer->SetInputData(mpVtkUnstructedMesh);
#else
    p_writer->SetInput(mpVtkUnstructedMesh);
#endif
    // Uninitialised stuff arises (see #1079), but you can remove valgrind problems by removing compression:
    // **** REMOVE WITH CAUTION *****
    p_writer->SetCompressor(nullptr);
    // **** REMOVE WITH CAUTION *****

    std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName;
    if (!stamp.empty())
    {
        vtk_file_name += "_" + stamp;
    }
    vtk_file_name += ".vtu";

    p_writer->SetFileName(vtk_file_name.c_str());
    // p_writer->PrintSelf(std::cout, vtkIndent());
    p_writer->Write();
    p_writer->Delete(); // Reference counted
#endif // CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrapezoidEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    // Only 2D version is supported at the moment
    assert(SPACE_DIM == 2); // LCOV_EXCL_LINE
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    /*
     * Populating points. First, outer points of elements
     */
    const unsigned n_vertices = rMesh.GetNumNodes();
    for (unsigned node_num = 0; node_num < rMesh.GetNumNodes(); node_num++)
    {
        c_vector<double, 2> position;
        position = rMesh.GetNode(node_num)->rGetLocation();
        p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
    }
    /*
     * Populating inner points.
     * [_________________][_____][_________]....[_____]
     *      ^^^^^^^^^^^    ^^^^^  ^^^^^^^^        ^^^^
     *  Outer points       Cell_1  Cell_2        Cell_{num_elements}
     *                              Inner Points
     * cell_offset_dist stores the distance from the beginning of p_pts array for each element
     * Note that the number of inner points equals to the number of nodes of each element
     */
    const unsigned num_elements = rMesh.GetNumElements();
    std::vector<unsigned> cell_offset_dist(rMesh.GetNumElements());
    cell_offset_dist[0] = n_vertices;
    for (unsigned i = 1; i < num_elements; ++i)
        cell_offset_dist[i] = cell_offset_dist[i - 1] + rMesh.GetElement(i - 1)->GetNumNodes();
    // Coefficient 0<=alpha<=1 represents how thin the trapezoid is
    const double alpha = 0.8;
    typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator elem_end = rMesh.GetElementIteratorEnd();
    for (auto elem = rMesh.GetElementIteratorBegin(); elem != elem_end; ++elem)
    {
        const unsigned num_elem_nodes = elem->GetNumNodes();
        const c_vector<double, SPACE_DIM> elem_centroid = rMesh.GetCentroidOfElement(elem->GetIndex());
        for (unsigned elem_node_num = 0; elem_node_num < num_elem_nodes; elem_node_num++)
        {
            c_vector<double, SPACE_DIM> node_position(2);
            node_position = elem->GetNode(elem_node_num)->rGetLocation();
            const double new_x = (node_position[0] - elem_centroid[0]) * alpha + elem_centroid[0];
            const double new_y = (node_position[1] - elem_centroid[1]) * alpha + elem_centroid[1];
            p_pts->InsertPoint(cell_offset_dist[elem->GetIndex()] + elem_node_num,
                               new_x, new_y, 0.0);
        }
    }
    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted
    unsigned total_num_edges = 0;
    for (auto elem = rMesh.GetElementIteratorBegin(); elem != elem_end; ++elem)
    {
        // First do the trapezoids for each edge
        for (unsigned edge_index = 0; edge_index < elem->GetNumEdges(); ++edge_index)
        {
            vtkCell* p_cell;
            p_cell = vtkQuad::New();
            const unsigned num_trap_nodes = p_cell->GetNumberOfEdges(); // 4 in 2D, 8 in 3D
            assert(num_trap_nodes == 4);
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(num_trap_nodes);
            auto p_edge = elem->GetEdge(edge_index);
            assert(p_edge->GetNumNodes() == 2);

            // See the diagram above for storing pattern
            std::array<unsigned, 2> base_ids = { { p_edge->GetNode(0)->GetIndex(), p_edge->GetNode(1)->GetIndex() } };
            std::array<unsigned, 2> top_ids = { { elem->GetNodeLocalIndex(base_ids[0])
                                                      + cell_offset_dist[elem->GetIndex()],
                                                  elem->GetNodeLocalIndex(base_ids[1]) + cell_offset_dist[elem->GetIndex()] } };

            // Assuming counter-clockwise ordering
            p_cell_id_list->SetId(0, base_ids[0]);
            p_cell_id_list->SetId(1, base_ids[1]);
            p_cell_id_list->SetId(2, top_ids[1]);
            p_cell_id_list->SetId(3, top_ids[0]);
            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
            p_cell->Delete(); // Reference counted
            total_num_edges++;
        }

        // Now do the internal cell
        vtkCell* p_cell;
        p_cell = vtkPolygon::New();
        const unsigned num_elem_nodes = elem->GetNumNodes();
        vtkIdList* p_cell_id_list = p_cell->GetPointIds();
        p_cell_id_list->SetNumberOfIds(num_elem_nodes);
        for (unsigned j = 0; j < num_elem_nodes; ++j)
        {
            p_cell_id_list->SetId(j, cell_offset_dist[elem->GetIndex()] + j);
        }

        mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);

        p_cell->Delete(); // Reference counted
    }

    // For 2D case. For 3D, we should sum the total number of faces + num_elements
    assert(total_num_edges + num_elements == mpVtkUnstructedMesh->GetNumberOfCells());
#endif // CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrapezoidEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i = 0; i < dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif // CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrapezoidEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    /*
     * Blank as we're only using the class for VTK at the moment, i.e. we don't 
     * write mesh information for the case when we have trapezoid edge elements, 
     * since trapezoids associated with each edge are only there for 
     * visualisation purposes and do not represent the actual mesh. The actual 
     * mesh can be written by VertexMeshWriter class and reconstructed by 
     * VertexMeshReader class. This method needs to be overriden here, as it is 
     * declared as a pure virtual function in the parent class.
     */
}

// Explicit instantiation
template class TrapezoidEdgeVertexMeshWriter<1, 1>;
template class TrapezoidEdgeVertexMeshWriter<1, 2>;
template class TrapezoidEdgeVertexMeshWriter<1, 3>;
template class TrapezoidEdgeVertexMeshWriter<2, 2>;
template class TrapezoidEdgeVertexMeshWriter<2, 3>;
template class TrapezoidEdgeVertexMeshWriter<3, 3>;