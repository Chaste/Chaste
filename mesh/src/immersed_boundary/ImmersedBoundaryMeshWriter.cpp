/*

Copyright (c) 2005-2018, University of Oxford.
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

#include "ImmersedBoundaryMeshWriter.hpp"
#include "Version.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include <boost/multi_array.hpp>

/**
 * Convenience collection of iterators, primarily to get compilation to happen.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MeshWriterIterators
{
    /** Iterator over nodes */
    typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator* pNodeIter;
    /** Iterator over immersed boundary elements */
    typename ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>::ImmersedBoundaryElementIterator* pElemIter;
    /** Iterator over immersed boundary laminas */
    typename ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>::ImmersedBoundaryLaminaIterator* pLamIter;
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryMeshWriter(const std::string& rDirectory,
                                                                               const std::string& rBaseName,
                                                                               const bool clearOutputDir)
        : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
          mpMesh(NULL),
          mpIters(new MeshWriterIterators<ELEMENT_DIM, SPACE_DIM>)
{
    mpIters->pNodeIter = NULL;
    mpIters->pElemIter = NULL;
    mpIters->pLamIter = NULL;

#ifdef CHASTE_VTK
    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::~ImmersedBoundaryMeshWriter()
{
    if (mpIters->pNodeIter)
    {
        delete mpIters->pNodeIter;
        delete mpIters->pElemIter;
        delete mpIters->pLamIter;
    }

    delete mpIters;

#ifdef CHASTE_VTK
// Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh->Delete(); // Reference counted
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
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
ImmersedBoundaryElementData ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextImmersedBoundaryElement()
{
    ///\todo Assert this method should only be called in 2D? (#1076/#1377)

    assert(this->mNumElements == mpMesh->GetNumElements());

    ImmersedBoundaryElementData elem_data;
    elem_data.NodeIndices.resize((*(mpIters->pElemIter))->GetNumNodes());
    for (unsigned j=0; j<elem_data.NodeIndices.size(); j++)
    {
        elem_data.NodeIndices[j] = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
    }

    // Set attribute
    elem_data.AttributeValue = (*(mpIters->pElemIter))->GetAttribute();

    ++(*(mpIters->pElemIter));

    return elem_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElementData ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextImmersedBoundaryLamina()
{
    ///\todo Assert this method should only be called in 2D? (#1076/#1377)

    assert(mNumLaminas == mpMesh->GetNumLaminas());

    ImmersedBoundaryElementData lamina_data;
    lamina_data.NodeIndices.resize((*(mpIters->pLamIter))->GetNumNodes());
    for (unsigned j=0; j<lamina_data.NodeIndices.size(); j++)
    {
        lamina_data.NodeIndices[j] = (*(mpIters->pLamIter))->GetNodeGlobalIndex(j);
    }

    // Set attribute
    lamina_data.AttributeValue = (*(mpIters->pLamIter))->GetAttribute();

    ++(*(mpIters->pLamIter));

    return lamina_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string stamp)
{
#ifdef CHASTE_VTK
    assert(SPACE_DIM == 2);
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
    p_writer->SetCompressor(NULL);
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
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    /**
     * To allow viewing in Paraview, we have to treat differently cells which overlap the boundaries, as there is no
     * support for periodicity in Paraview.
     *
     * We overcome this by first identifying which cells overlap, and breaking them in to pieces as necessary, so that
     * each cell when displayed stays contiguous.
     *
     * Cell overlaps should have already been calculated - this is done by a call to CalculateCellOverlaps() from
     * ImmersedBoundaryCellPopulation::WriteVtkResultsToFile().
     *
     * Because no node can be present in more than one cell, it is safe to add points to mpVtkUnstructuredMesh as we go
     * through each cell, rather than having to add nodes before cells.
     */
//#ifdef CHASTE_VTK
    // Assert mesh is 2D to simplify vtk creation
    assert(SPACE_DIM == 2);

    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");

    // Next, we decide how to output the VTK data for elements depending on the type of overlap
    for (typename ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>::ImmersedBoundaryElementIterator iter = rMesh.GetElementIteratorBegin();
            iter != rMesh.GetElementIteratorEnd();
            ++iter)
    {
        unsigned elem_idx = iter->GetIndex();
        unsigned num_nodes = iter->GetNumNodes();

        // Case 1:  no overlap at all
        if ( !mHOverlaps[elem_idx] && !mVOverlaps[elem_idx] )
        {
            // Double check no points of overlap were found
            assert( (mHOverlapPoints[elem_idx].size() == 0) && (mVOverlapPoints[elem_idx].size() == 0) );

            vtkCell* p_cell = vtkPolygon::New();
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(iter->GetNumNodes());

            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                // Get node, index and location
                Node<SPACE_DIM>* p_node = iter->GetNode(node_idx);
                unsigned global_idx = p_node->GetIndex();
                c_vector<double, SPACE_DIM> position = p_node->rGetLocation();

                p_pts->InsertPoint(global_idx, position[0], position[1], 0.0);

                p_cell_id_list->SetId(node_idx, global_idx);
            }

            mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
            p_cell->Delete(); // Reference counted
        }

        // Case 2:  only horizontal OR vertical overlap (exclusive)
        else if ( ( mHOverlaps[elem_idx] && !mVOverlaps[elem_idx] ) ||
                  ( mVOverlaps[elem_idx] && !mHOverlaps[elem_idx] ) )
        {
            // There should be exactly two points of overlap found - if not, there is likely to be some weird geometry
            // which has not been considered
            assert( ( (mHOverlapPoints[elem_idx].size() == 2) && (mVOverlapPoints[elem_idx].size() == 0) ) ||
                    ( (mVOverlapPoints[elem_idx].size() == 2) && (mHOverlapPoints[elem_idx].size() == 0) ) );

            // Decide whether we need the horizontal or vertical overlap points
            std::vector<unsigned> overlap = (mHOverlapPoints[elem_idx].size() == 2) ? mHOverlapPoints[elem_idx] : mVOverlapPoints[elem_idx];

            // We will need to create two 'cells' for output.  First find number of nodes in each of these half-cells
            //std::sort(overlap.begin(), overlap.end());
            unsigned num_nodes_a = overlap[1] - overlap[0];
            unsigned num_nodes_b = num_nodes - num_nodes_a;

            vtkCell* p_cell_a = vtkPolygon::New();
            vtkIdList* p_cell_id_list_a = p_cell_a->GetPointIds();
            p_cell_id_list_a->SetNumberOfIds(num_nodes_a);

            for (unsigned node_idx = 0; node_idx < num_nodes_a; node_idx++)
            {
                // Get node, index and location
                Node<SPACE_DIM>* p_node = iter->GetNode(node_idx + overlap[0]);
                unsigned global_idx = p_node->GetIndex();
                c_vector<double, SPACE_DIM> position = p_node->rGetLocation();

                p_pts->InsertPoint(global_idx, position[0], position[1], 0.0);

                p_cell_id_list_a->SetId(node_idx, global_idx);
            }

            mpVtkUnstructedMesh->InsertNextCell(p_cell_a->GetCellType(), p_cell_id_list_a);
            p_cell_a->Delete(); // Reference counted


            vtkCell* p_cell_b = vtkPolygon::New();
            vtkIdList* p_cell_id_list_b = p_cell_b->GetPointIds();
            p_cell_id_list_b->SetNumberOfIds(num_nodes_b);

            for (unsigned node_idx = 0; node_idx < num_nodes_b; node_idx++)
            {
                // Get node, index and location
                Node<SPACE_DIM>* p_node = iter->GetNode((node_idx + overlap[1]) % num_nodes);
                unsigned global_idx = p_node->GetIndex();
                c_vector<double, SPACE_DIM> position = p_node->rGetLocation();

                p_pts->InsertPoint(global_idx, position[0], position[1], 0.0);

                p_cell_id_list_b->SetId(node_idx, global_idx);
            }

            mpVtkUnstructedMesh->InsertNextCell(p_cell_b->GetCellType(), p_cell_id_list_b);
            p_cell_b->Delete(); // Reference counted
        }


        // Case 3:  multiple overlaps
        else //( (h_overlaps[elem_idx] == true) && (v_overlaps[elem_idx] == true)
        {
            //\todo: implement this

            // There should be exactly two points of overlap found - if not, there is likely to be some weird geometry
            // which has not been considered
            assert( (mHOverlapPoints[elem_idx].size() == 2) && (mVOverlapPoints[elem_idx].size() == 2) );

            NEVER_REACHED;
        }
    }

    // Finally, output the VTK data for laminas
    for (typename ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>::ImmersedBoundaryLaminaIterator iter = rMesh.GetLaminaIteratorBegin();
         iter != rMesh.GetLaminaIteratorEnd();
         ++iter)
    {
        unsigned num_nodes = iter->GetNumNodes();

        vtkCell* p_cell = vtkPolygon::New();
        vtkIdList* p_cell_id_list = p_cell->GetPointIds();
        p_cell_id_list->SetNumberOfIds(num_nodes);

        for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
        {
            // Get node, index and location
            Node<SPACE_DIM>* p_node = iter->GetNode(node_idx);
            unsigned global_idx = p_node->GetIndex();
            c_vector<double, SPACE_DIM> position = p_node->rGetLocation();

            p_pts->InsertPoint(global_idx, position[0], position[1], 0.0);

            p_cell_id_list->SetId(node_idx, global_idx);
        }

        mpVtkUnstructedMesh->InsertNextCell(3, p_cell_id_list);
        p_cell->Delete(); // Reference counted
    }

    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted

//#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
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
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
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
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    this->mpMeshReader = NULL;
    mpMesh = &rMesh;

    this->mNumNodes = mpMesh->GetNumNodes();
    this->mNumElements = mpMesh->GetNumElements();
    mNumLaminas = mpMesh->GetNumLaminas();

    typedef typename AbstractMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>::ImmersedBoundaryElementIterator ElemIterType;
    mpIters->pElemIter = new ElemIterType(mpMesh->GetElementIteratorBegin());

    typedef typename ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>::ImmersedBoundaryLaminaIterator LamIterType;
    mpIters->pLamIter = new LamIterType(mpMesh->GetLaminaIteratorBegin());

    WriteFiles();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
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
    *p_node_file << mpMesh->GetCharacteristicNodeSpacing() << "\t";
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
    std::string element_file_name = this->mBaseName + ".elem";
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
            ImmersedBoundaryElementData elem_data = this->GetNextImmersedBoundaryElement();

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

            *p_element_file << "\t" << elem_data.SpringConstant;

            *p_element_file << "\t" << elem_data.RestLength;

            // New line
            *p_element_file << "\n";
        }
        else // 3D
        {
        }
    }
    *p_element_file << comment << "\n";
    p_element_file->close();

    // Write lamina file
    std::string lamina_file_name = this->mBaseName + ".lam";
    out_stream p_lamina_file = this->mpOutputFileHandler->OpenOutputFile(lamina_file_name);

    // Write the lamina header
    num_attr = 1; //Always write element attributes
    unsigned num_laminas = mpMesh->GetNumLaminas();
    *p_lamina_file << num_laminas << "\t" << num_attr << "\n";

    // Write each lamina's data
    for (unsigned item_num=0; item_num<num_laminas; item_num++)
    {
        if (SPACE_DIM == 2) // In 2D, write the node indices owned by this element
        {
            // Get data for this element
            ImmersedBoundaryElementData lamina_data = this->GetNextImmersedBoundaryLamina();

            // Get the node indices owned by this element
            std::vector<unsigned> node_indices = lamina_data.NodeIndices;

            // Write this element's index and the number of nodes owned by it to file
            *p_lamina_file << item_num <<  "\t" << node_indices.size();

            // Write the node indices owned by this element to file
            for (unsigned i=0; i<node_indices.size(); i++)
            {
                *p_lamina_file << "\t" << node_indices[i];
            }

            *p_lamina_file << "\t" << lamina_data.AttributeValue;

            *p_lamina_file << "\t" << lamina_data.SpringConstant;

            *p_lamina_file << "\t" << lamina_data.RestLength;

            // New line
            *p_lamina_file << "\n";
        }
        else // 3D
        {
        }
    }
    *p_lamina_file << comment << "\n";
    p_lamina_file->close();

    // Write grid file
    std::string grid_file_name = this->mBaseName + ".grid";
    out_stream p_grid_file = this->mpOutputFileHandler->OpenOutputFile(grid_file_name);

    // Write the element header
    unsigned num_gridpts_x = mpMesh->GetNumGridPtsX();
    unsigned num_gridpts_y = mpMesh->GetNumGridPtsY();

    *p_grid_file << num_gridpts_x << "\t" << num_gridpts_y << "\n";

    // Write grid data
    const multi_array<double, 3>& vel_grids = mpMesh->rGet2dVelocityGrids();

    for (unsigned y_idx = 0; y_idx < num_gridpts_y; y_idx ++)
    {
        for (unsigned x_idx = 0; x_idx < num_gridpts_x; x_idx ++)
        {
            *p_grid_file << vel_grids[0][x_idx][y_idx] << "\t";
        }
        *p_grid_file << "\n";
    }

    for (unsigned y_idx = 0; y_idx < num_gridpts_y; y_idx ++)
    {
        for (unsigned x_idx = 0; x_idx < num_gridpts_x; x_idx ++)
        {
            *p_grid_file << vel_grids[1][x_idx][y_idx] << "\t";
        }
        *p_grid_file << "\n";
    }

    *p_grid_file << comment << "\n";
    p_grid_file->close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::CalculateCellOverlaps(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    assert(SPACE_DIM == 2);

    // Initialise all vectors to the correct length
    unsigned num_elem = rMesh.GetNumAllElements();
    mHOverlaps.resize(num_elem);
    mVOverlaps.resize(num_elem);
    mHOverlapPoints.resize(num_elem);
    mVOverlapPoints.resize(num_elem);
    mNumCellParts.resize(num_elem);

    // Helper variables
    c_vector<double, SPACE_DIM> prev_location;
    c_vector<double, SPACE_DIM> curr_location;

    // We loop first over each element to work out which overlap due to periodic boundaries
    for (typename ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>::ImmersedBoundaryElementIterator iter = rMesh.GetElementIteratorBegin();
         iter != rMesh.GetElementIteratorEnd();
         ++iter)
    {
        unsigned elem_idx = iter->GetIndex();

        unsigned num_nodes = iter->GetNumNodes();
        assert(num_nodes > 1);
        prev_location = iter->GetNode(num_nodes - 1)->rGetLocation();

        for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
        {
            // Get node, index and location
            Node<SPACE_DIM>* p_node = iter->GetNode(node_idx);
            c_vector<double, SPACE_DIM> curr_location = p_node->rGetLocation();

            if ( fabs (curr_location[0] - prev_location[0]) > 0.5)
            {
                mHOverlaps[elem_idx] = true;
                mHOverlapPoints[elem_idx].push_back(node_idx);
            }

            if ( fabs (curr_location[1] - prev_location[1]) > 0.5)
            {
                mVOverlaps[elem_idx] = true;
                mVOverlapPoints[elem_idx].push_back(node_idx);
            }
            prev_location = curr_location;
        }
    }

    // We then loop over each element again to determine how many 'cells' we need to split each output element in to
    for (unsigned elem_idx = 0; elem_idx < num_elem; elem_idx++)
    {
        // If no overlap, we only need one cell
        if (!mHOverlaps[elem_idx] && !mVOverlaps[elem_idx])
        {
            mNumCellParts[elem_idx] = 1;
        }
        // If only one overlap is false, then we need two cells
        else if ( mHOverlaps[elem_idx] * mVOverlaps[elem_idx] == 0)
        {
            mNumCellParts[elem_idx] = 2;
        }
        // The remaining case is there is horizontal and vertical overlap.  We may need either 3 or 4 cells
        else
        {
            mNumCellParts[elem_idx] = 3;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned>& ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::rGetNumCellParts() const
{
    return mNumCellParts;
}

// Explicit instantiation
template class ImmersedBoundaryMeshWriter<1,1>;
template class ImmersedBoundaryMeshWriter<1,2>;
template class ImmersedBoundaryMeshWriter<1,3>;
template class ImmersedBoundaryMeshWriter<2,2>;
template class ImmersedBoundaryMeshWriter<2,3>;
template class ImmersedBoundaryMeshWriter<3,3>;
