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

#include "VertexMeshWriter.hpp"
#include "Version.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "Toroidal2dVertexMesh.hpp"

/**
 * Convenience collection of iterators, primarily to get compilation to happen.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MeshWriterIterators
{
    /** Iterator over nodes */
    typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator* pNodeIter;
    /** Iterator over vertex elements */
    typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator* pElemIter;
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::VertexMeshWriter(const std::string& rDirectory,
                                                           const std::string& rBaseName,
                                                           const bool clearOutputDir)
    : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
      mpMesh(nullptr),
      mpIters(new MeshWriterIterators<ELEMENT_DIM, SPACE_DIM>),
      mpNodeMap(nullptr),
      mNodeMapCurrentIndex(0)
{
    mpIters->pNodeIter = nullptr;
    mpIters->pElemIter = nullptr;

#ifdef CHASTE_VTK
     // Dubious, since we shouldn't yet know what any details of the mesh are.
     mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::~VertexMeshWriter()
{
    if (mpIters->pNodeIter)
    {
        delete mpIters->pNodeIter;
        delete mpIters->pElemIter;
    }

    delete mpIters;

    if (mpNodeMap)
    {
        delete mpNodeMap;
    }

#ifdef CHASTE_VTK
     // Dubious, since we shouldn't yet know what any details of the mesh are.
     mpVtkUnstructedMesh->Delete(); // Reference counted
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    if (mpMesh)
    {
        // Sanity check
        assert(this->mNumNodes == mpMesh->GetNumNodes());

        std::vector<double> coordinates(SPACE_DIM+1);

        // Get the node coordinates using the node iterator (thus skipping deleted nodes)
        for (unsigned j=0; j<SPACE_DIM; j++)
        {
            coordinates[j] = (*(mpIters->pNodeIter))->GetPoint()[j];
        }
        coordinates[SPACE_DIM] = (*(mpIters->pNodeIter))->IsBoundaryNode();

        ++(*(mpIters->pNodeIter));

        return coordinates;
    }
    else
    {
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextNode();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElementData VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElementWithFaces()
{
    // This method should only be called in 3D
    assert(SPACE_DIM == 3);    // LCOV_EXCL_LINE
    assert(mpMesh);
    assert(this->mNumElements == mpMesh->GetNumElements());

    // Create data structure for this element
    VertexElementData elem_data;

    // Store node indices owned by this element
    elem_data.NodeIndices.resize((*(mpIters->pElemIter))->GetNumNodes());
    for (unsigned j=0; j<elem_data.NodeIndices.size(); j++)
    {
        unsigned old_index = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
        elem_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
    }

    // Store faces owned by this element
    elem_data.Faces.resize((*(mpIters->pElemIter))->GetNumFaces());
    for (unsigned i=0; i<elem_data.Faces.size(); i++)
    {
        // Get pointer to this face
        VertexElement<ELEMENT_DIM-1, SPACE_DIM>* p_face = (*(mpIters->pElemIter))->GetFace(i);

        // Create data structure for this face
        ElementData face_data;

        // Store this face's index
        face_data.AttributeValue = p_face->GetIndex();

        // Store node indices owned by this face
        face_data.NodeIndices.resize(p_face->GetNumNodes());
        for (unsigned j=0; j<face_data.NodeIndices.size(); j++)
        {
            unsigned old_index = p_face->GetNodeGlobalIndex(j);
            face_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
        }

        // Store this face
        elem_data.Faces[i] = face_data;

        ///\todo Store face orientations? (#1076/#1377)
    }

    // Set attribute
    elem_data.AttributeValue = (unsigned)(*(mpIters->pElemIter))->GetAttribute();
    ++(*(mpIters->pElemIter));

    return elem_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement()
{
    ///\todo Assert this method should only be called in 2D? (#1076/#1377)

    if (mpMesh)
    {
        assert(this->mNumElements == mpMesh->GetNumElements());

        ElementData elem_data;
        elem_data.NodeIndices.resize((*(mpIters->pElemIter))->GetNumNodes());
        for (unsigned j=0; j<elem_data.NodeIndices.size(); j++)
        {
            unsigned old_index = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
            elem_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
        }

        // Set attribute
        elem_data.AttributeValue = (*(mpIters->pElemIter))->GetAttribute();
        ++(*(mpIters->pElemIter));

        return elem_data;
    }
    else
    {
        return AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string stamp)
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

/**
 * Write VTK file using a mesh.
 *
 * @param rMesh reference to the vertex-based mesh
 * @param stamp is an optional stamp (like a time-stamp) to put into the name of the file
 */
template<>
void VertexMeshWriter<2, 2>::WriteVtkUsingMesh(VertexMesh<2, 2>& rMesh, std::string stamp)
{
#ifdef CHASTE_VTK
    // Create VTK mesh
    VertexMesh<2, 2>* p_mesh_for_vtk = rMesh.GetMeshForVtk();
    MakeVtkMesh(*p_mesh_for_vtk);

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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
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

    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator iter = rMesh.GetElementIteratorBegin();
         iter != rMesh.GetElementIteratorEnd();
         ++iter)
    {
        vtkCell* p_cell;
        if (SPACE_DIM == 2)
        {
            p_cell = vtkPolygon::New();
        }
        else
        {
            p_cell = vtkConvexPointSet::New();
        }
        vtkIdList* p_cell_id_list = p_cell->GetPointIds();
        p_cell_id_list->SetNumberOfIds(iter->GetNumNodes());
        for (unsigned j=0; j<iter->GetNumNodes(); ++j)
        {
            p_cell_id_list->SetId(j, iter->GetNodeGlobalIndex(j));
        }
        mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
        p_cell->Delete(); // Reference counted
    }
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i=0; i<dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif //CHASTE_VTK
}

///\todo Mesh should be const (#1076)
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(VertexMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    this->mpMeshReader = nullptr;
    mpMesh = &rMesh;

    this->mNumNodes = mpMesh->GetNumNodes();
    this->mNumElements = mpMesh->GetNumElements();

    typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator ElemIterType;
    mpIters->pElemIter = new ElemIterType(mpMesh->GetElementIteratorBegin());

    // Set up node map if we might have deleted nodes
    mNodeMapCurrentIndex = 0;
    if (mpMesh->IsMeshChanging())
    {
        mpNodeMap = new NodeMap(mpMesh->GetNumAllNodes());
        for (NodeIterType it = mpMesh->GetNodeIteratorBegin(); it != mpMesh->GetNodeIteratorEnd(); ++it)
        {
            mpNodeMap->SetNewIndex(it->GetIndex(), mNodeMapCurrentIndex++);
        }
    }
    WriteFiles();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();

    // Write node file
    std::string node_file_name = this->mBaseName + ".node";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    // Write the node header
    unsigned num_attr = 0;
    unsigned max_bdy_marker = 1; // as we include boundary node information in the node file
    unsigned num_nodes = this->GetNumNodes();

    *p_node_file << num_nodes << "\t";
    *p_node_file << SPACE_DIM << "\t";
    *p_node_file << num_attr << "\t";
    *p_node_file << max_bdy_marker << "\n";
    *p_node_file << std::setprecision(6);

    // Write each node's data
    for (unsigned item_num=0; item_num<num_nodes; item_num++)
    {
        std::vector<double> current_item = this->GetNextNode();
        *p_node_file << item_num;
        for (unsigned i=0; i<SPACE_DIM+1; i++)
        {
            *p_node_file << "\t" << current_item[i];
        }
        *p_node_file << "\n";
    }
    *p_node_file << comment << "\n";
    p_node_file->close();

    // Write element file
    std::string element_file_name = this->mBaseName + ".cell";
    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name);

    // Write the element header
    num_attr = 1; //Always write element attributes
    unsigned num_elements = this->GetNumElements();
    *p_element_file << num_elements << "\t" << num_attr << "\n";

    // Write each element's data
    for (unsigned item_num=0; item_num<num_elements; item_num++)
    {
        if (SPACE_DIM == 2) // In 2D, write the node indices owned by this element
        {
            // Get data for this element
            ElementData elem_data = this->GetNextElement();

            // Get the node indices owned by this element
            std::vector<unsigned> node_indices = elem_data.NodeIndices;

            // Write this element's index and the number of nodes owned by it to file
            *p_element_file << item_num <<  "\t" << node_indices.size();

            // Write the node indices owned by this element to file
            for (unsigned i=0; i<node_indices.size(); i++)
            {
                *p_element_file << "\t" << node_indices[i];
            }

            *p_element_file << "\t" << elem_data.AttributeValue;

            // New line
            *p_element_file << "\n";
        }
        else // 3D
        {
            assert(SPACE_DIM == 3);    // LCOV_EXCL_LINE

            // Get data for this element
            VertexElementData elem_data_with_faces = this->GetNextElementWithFaces();

            // Get the node indices owned by this element
            std::vector<unsigned> node_indices = elem_data_with_faces.NodeIndices;

            // Write this element's index and the number of nodes owned by it to file
            *p_element_file << item_num <<  "\t" << node_indices.size();

            // Write the node indices owned by this element to file
            for (unsigned i=0; i<node_indices.size(); i++)
            {
                *p_element_file << "\t" << node_indices[i];
            }

            // Get the faces owned by this element (if any)
            std::vector<ElementData> faces = elem_data_with_faces.Faces;

            // Write the number of faces owned by this element to file
            *p_element_file << "\t" << faces.size();

            for (unsigned j=0; j<faces.size(); j++)
            {
                // Get the node indices owned by this face
                std::vector<unsigned> face_node_indices = faces[j].NodeIndices;

                // Write this face's index and the number of nodes owned by it to file
                *p_element_file << "\t" << faces[j].AttributeValue <<  "\t" << face_node_indices.size();

                // Write the node indices owned by this element to file
                for (unsigned i=0; i<face_node_indices.size(); i++)
                {
                    *p_element_file << "\t" << face_node_indices[i];
                }

                ///\todo Store face orientations? (#1076/#1377)
            }

            *p_element_file << "\t" << elem_data_with_faces.AttributeValue;

            // New line
            *p_element_file << "\n";
        }
    }

    *p_element_file << comment << "\n";
    p_element_file->close();
}

///////// Explicit instantiation///////

template class VertexMeshWriter<1,1>;
template class VertexMeshWriter<1,2>;
template class VertexMeshWriter<1,3>;
template class VertexMeshWriter<2,2>;
template class VertexMeshWriter<2,3>;
template class VertexMeshWriter<3,3>;
