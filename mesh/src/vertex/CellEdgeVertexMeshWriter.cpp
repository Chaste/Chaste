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

#include "CellEdgeVertexMeshWriter.hpp"


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
CellEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::CellEdgeVertexMeshWriter(const std::string &rDirectory,
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
CellEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::~CellEdgeVertexMeshWriter()
{

}


template<unsigned int ELEMENT_DIM, unsigned int SPACE_DIM>
void CellEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM> &rMesh,
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
void CellEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM> &rMesh)
{
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");

    unsigned num_nodes = rMesh.GetNumNodes();

    for (unsigned node_num=0; node_num<rMesh.GetNumNodes(); node_num++)
    {
        c_vector<double, SPACE_DIM> position = rMesh.GetNode(node_num)->rGetLocation();
        if (SPACE_DIM==2)
        {
            p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
        }
        else
        {
            p_pts->InsertPoint(node_num, position[0], position[1], position[2]);
        }
    }

    //Inserts the centre points into the vertex index
    unsigned centre_point_counter = num_nodes;
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator iter = rMesh.GetElementIteratorBegin();
         iter != rMesh.GetElementIteratorEnd();
         ++iter)
    {
        auto position = rMesh.GetCentroidOfElement(iter->GetIndex());
        if (SPACE_DIM==2)
        {
            p_pts->InsertPoint(centre_point_counter, position[0], position[1], 0.0);
        }
        else
        {
            p_pts->InsertPoint(centre_point_counter, position[0], position[1], position[2]);
        }
        ++centre_point_counter;
    }


    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted
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
CellEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
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
void CellEdgeVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    //Blank as we're only using the class for VTK at the moment
}

///////// Explicit instantiation///////

template class CellEdgeVertexMeshWriter<1,1>;
template class CellEdgeVertexMeshWriter<1,2>;
template class CellEdgeVertexMeshWriter<1,3>;
template class CellEdgeVertexMeshWriter<2,2>;
template class CellEdgeVertexMeshWriter<2,3>;
template class CellEdgeVertexMeshWriter<3,3>;
