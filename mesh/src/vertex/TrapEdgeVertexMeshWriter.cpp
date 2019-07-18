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

#include "TrapEdgeVertexMeshWriter.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::TrapEdgeVertexMeshWriter(const std::string &rDirectory,
                                                                           const std::string &rBaseName,
                                                                           const bool clearOutputDir)
        : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
        mpMesh(nullptr)
{

#ifdef CHASTE_VTK
    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
#endif //CHASTE_VTK

}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::~TrapEdgeVertexMeshWriter()
{

}


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
                                                                         std::string stamp)
{
#ifdef CHASTE_VTK
    assert(SPACE_DIM==3 || SPACE_DIM == 2);    // LCOV_EXCL_LINE

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
    if (stamp != "")
    {
        vtk_file_name += "_" + stamp;
    }
    vtk_file_name += ".vtu";

    p_writer->SetFileName(vtk_file_name.c_str());
    //p_writer->PrintSelf(std::cout, vtkIndent());
    p_writer->Write();
    p_writer->Delete(); // Reference counted
#endif //CHASTE_VTK

}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM> &rMesh)
{
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    //Vertex here means vertex of vtkCells, not the real vertices
    p_pts->GetData()->SetName("Vertex positions");
    typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem, elem_end;
    elem = rMesh.GetElementIteratorBegin();
    elem_end = rMesh.GetElementIteratorEnd();
    unsigned node_num = 0;
    const unsigned total_n_vertices = rMesh.GetNumNodes();
    //Coefficient 0<=alpha<=1 represents how close to the trapezoid vertex
    //should be close to the cell's vertex.
    const double alpha = 0.8;
    for (; elem != elem_end; ++elem)
    {
        const unsigned n_elem_nodes = elem->GetNumNodes();
        const c_vector<double, SPACE_DIM> elem_centroid = rMesh.GetCentroidOfElement(elem->GetIndex());
        for (unsigned node_num = 0; node_num < n_elem_nodes; node_num++)
        {
            const c_vector<double, SPACE_DIM> node_position = elem->GetNode(node_num)->rGetLocation();
            const double new_x = (node_position[0]-elem_centroid[0])*alpha + elem_centroid[0];
            const double new_y = (node_position[1]-elem_centroid[1])*alpha + elem_centroid[1];
            c_vector<double, SPACE_DIM> trap_node(new_x, new_y, 0.0);
            if (SPACE_DIM == 2)
            {
                p_pts->InsertPoint(node_num, node_position[0], node_position[1], 0.0);
                p_pts->InsertPoint(node_num+1, trap_node[0], trap_node[1], 0.0);
            }
            else
            {
                const double new_z = (node_position[2]-elem_centroid[2])*alpha + elem_centroid[2];
                trap_node[2] = new_z;
                p_pts->InsertPoint(node_num, node_position[0], node_position[1], node_position[2]);
                p_pts->InsertPoint(node_num+1, trap_node[0], trap_node[1], trap_node[2]);
            }
            node_num += 2;
        }
    }
    //The number of nodes in the new mesh must be double of the number of the vertices
    assert(node_num == 2*total_n_vertices);

    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted
    elem = rMesh.GetElementIteratorBegin();
    /*Only 2D is fully supported.
     *For 3D, if "edge" should be synonymous with face
     *then, in principle, below should work
     */
    unsigned n_trap_nodes; //4 in 2D, 8 in 3D
    for (; elem != elem_end; ++elem)
    {
        //First do the trapezoids for each edge
        if (SPACE_DIM==2)
        {
            for (unsigned edge_index = 0; edge_index <elem->GetNumEdges(); ++edge_index)
            {
                vtkCell* p_cell;
                p_cell = vtkQuad::New();
                n_trap_nodes = p_cell->GetNumberOfEdges;
                assert(n_trap_nodes == 4);
                // I *think* vtkIdList provides mapping between local cell node Ids
                // with global node indices, where the vector of global nodes was filled in p_pts
                vtkIdList* p_cell_id_list = p_cell->GetPointIds();
                p_cell_id_list->SetNumberOfIds(n_trap_nodes);
                auto p_edge = elem->GetEdge(edge_index);
                assert(p_edge->GetNumNodes()==2);
                //When filling p_pts array, following the real vertex the internal
                //trapezoid node is added immediately after.
                //Thus, the original mesh ordering is (should be) shifted by one
                c_vector<unsigned, 2> base_ids(p_edge->GetNode(0)->GetIndex(),
                                               p_edge->GetNode(1)->GetIndex()+1);
                c_vector<unsigned, 2> top_ids(base_ids[0]+1, base_ids[1]+1);
                p_cell_id_list->SetId(0, base_ids[0]);
                p_cell_id_list->SetId(1, base_ids[1]);
                p_cell_id_list->SetId(2, top_ids[0]);
                p_cell_id_list->SetId(3, top_ids[1]);
            }
        }
        else
        {
            /*
             * \todo For each face ...
             */
            vtkCell* p_cell;
            p_cell = vtkHexahedron::New();
            n_trap_nodes = p_cell->GetNumberOfFaces;
            assert(n_trap_nodes == 6);
            p_cell->Delete();
        }
        //Now do the internal cell
    }
    for (unsigned elem_index = 0; elem_index < rMesh.GetNumElements(); ++elem_index)
    {
        auto p_element = rMesh.GetElement(elem_index);

        for (unsigned edge_index = 0; edge_index < p_element->GetNumEdges(); ++edge_index)
        {
            auto p_edge = p_element->GetEdge(edge_index);

            vtkCell* p_cell;
            p_cell = vtkTriangle::New();

            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(3);

            assert(p_edge->GetNumNodes() == 2);
            for (unsigned j=0; j<p_edge->GetNumNodes(); ++j)
            {
                p_cell_id_list->SetId(j, p_edge->GetNode(j)->GetIndex());
            }
            p_cell_id_list->SetId(2, num_nodes + elem_index); //The 3rd point is the centre location


            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
            p_cell->Delete(); // Reference counted
        }
    }
#endif //CHASTE_VTK
}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void
TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif //CHASTE_VTK

}

template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void TrapEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    //Blank as we're only using the class for VTK at the moment
}

///////// Explicit instantiation///////

template class TrapEdgeVertexMeshWriter<1,1>;
template class TrapEdgeVertexMeshWriter<1,2>;
template class TrapEdgeVertexMeshWriter<1,3>;
template class TrapEdgeVertexMeshWriter<2,2>;
template class TrapEdgeVertexMeshWriter<2,3>;
template class TrapEdgeVertexMeshWriter<3,3>;
