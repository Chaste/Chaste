/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "ImmersedBoundaryArray.hpp"
#include "MathsCustomFunctions.hpp"
#include "UblasCustomFunctions.hpp"
#include "Version.hpp"

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
          mpMesh(nullptr),
          mpIters(new MeshWriterIterators<ELEMENT_DIM, SPACE_DIM>)
{
    mpIters->pNodeIter = nullptr;
    mpIters->pElemIter = nullptr;
    mpIters->pLamIter = nullptr;

    switch (SPACE_DIM)
    {
        case 2:
        {
            geom_point corner_0(0.0, 0.0);
            geom_point corner_1(1.0, 0.0);
            geom_point corner_2(0.0, 1.0);
            geom_point corner_3(1.0, 1.0);

            mBoundaryEdges[0] = geom_segment(corner_0, corner_2);  // left edge
            mBoundaryEdges[1] = geom_segment(corner_1, corner_3);  // right edge
            mBoundaryEdges[2] = geom_segment(corner_0, corner_1);  // bottom edge
            mBoundaryEdges[3] = geom_segment(corner_2, corner_3);  // top edge

            break;
        }

        default:
            NEVER_REACHED;
    }

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
        return AbstractMeshWriter<ELEMENT_DIM,SPACE_DIM>::GetNextNode(); //LCOV_EXCL_LINE fairly sure this is unreachable
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElementData ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextImmersedBoundaryElement()
{
    // Can only be constructed in 2D currently
    //LCOV_EXCL_START
    if constexpr (SPACE_DIM != 2) {
        EXCEPTION("ImmersedBoundaryMeshWriter::GetNextImmersedBoundaryElement is not yet implemetned in 3D!");
    }
    //LCOV_EXCL_STOP

    assert(this->mNumElements == mpMesh->GetNumElements());

    ImmersedBoundaryElementData elem_data;
    elem_data.NodeIndices.resize((*(mpIters->pElemIter))->GetNumNodes());
    for (unsigned j=0; j<elem_data.NodeIndices.size(); j++)
    {
        elem_data.NodeIndices[j] = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
    }

    // Set attribute
    elem_data.AttributeValue = (*(mpIters->pElemIter))->GetAttribute();

    // Immersed boundary specific data
    auto& element = *(*mpIters->pElemIter);

    if (element.GetFluidSource() != nullptr) {
        elem_data.hasFluidSource = true;
        elem_data.fluidSourceIndex = element.GetFluidSource()->GetIndex();
    } else {
        elem_data.hasFluidSource = false;
    }

    for (auto cornerNode : element.rGetCornerNodes()) {
        elem_data.cornerNodeIndices.push_back(cornerNode->GetIndex());
    }

    elem_data.averageNodeSpacing = element.GetAverageNodeSpacing();

    elem_data.isBoundaryElement = element.IsElementOnBoundary();

    ++(*(mpIters->pElemIter));

    return elem_data;
} // LCOV_EXCL_LINE

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElementData ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextImmersedBoundaryLamina()
{
    //LCOV_EXCL_START
    if constexpr (SPACE_DIM != 2) {
        EXCEPTION("ImmersedBoundaryMeshWriter::GetNextImmersedBoundaryLamina is not yet implemetned in 3D!");
    }
    //LCOV_EXCL_STOP

    assert(mNumLaminas == mpMesh->GetNumLaminas());

    ImmersedBoundaryElementData lamina_data;
    lamina_data.NodeIndices.resize((*(mpIters->pLamIter))->GetNumNodes());
    for (unsigned j=0; j<lamina_data.NodeIndices.size(); j++)
    {
        lamina_data.NodeIndices[j] = (*(mpIters->pLamIter))->GetNodeGlobalIndex(j);
    }

    // Set attribute
    lamina_data.AttributeValue = (*(mpIters->pLamIter))->GetAttribute();

    // Immersed boundary specific data
    auto& lamina = *(*mpIters->pLamIter);

    try {
        lamina.GetFluidSource();
        // Can't really support 2D laminas at the moment
        //LCOV_EXCL_START
        if (lamina.GetFluidSource() != nullptr) {
            lamina_data.hasFluidSource = true;
            lamina_data.fluidSourceIndex = lamina.GetFluidSource()->GetIndex();
        } else {
            lamina_data.hasFluidSource = false;
        }
        //LCOV_EXCL_STOP
    } catch(Exception& e) {
        lamina_data.hasFluidSource = false;
    }

    for (auto cornerNode : lamina.rGetCornerNodes()) {
        lamina_data.cornerNodeIndices.push_back(cornerNode->GetIndex()); // LCOV_EXCL_LINE
    }

    lamina_data.averageNodeSpacing = lamina.GetAverageNodeSpacing();

    lamina_data.isBoundaryElement = lamina.IsElementOnBoundary();

    ++(*(mpIters->pLamIter));

    return lamina_data;
} // LCOV_EXCL_LINE

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string stamp)
{
#ifdef CHASTE_VTK
    if constexpr (SPACE_DIM == 2)
    {
        // Create VTK mesh
        MakeVtkMesh(rMesh);
        // Now write VTK mesh to file
        //assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
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
        p_writer->Write();
        p_writer->Delete(); // Reference counted
    }
    else
    {
        NEVER_REACHED;
    }
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
#ifdef CHASTE_VTK
    /**
     * To allow viewing in Paraview, we have to treat differently cells which overlap the boundaries, as there is no
     * support for periodicity in Paraview.
     *
     * We overcome this by first identifying which cells overlap, and breaking them in to pieces as necessary, so that
     * each cell when displayed stays contiguous.
     *
     * Cell overlaps should have already been calculated - this is done by a call to FindElementOverlaps() from
     * ImmersedBoundaryCellPopulation::WriteVtkResultsToFile().
     *
     * Because no node can be present in more than one cell, it is safe to add points to mpVtkUnstructuredMesh as we go
     * through each cell, rather than having to add nodes before cells.
     */
    // Assert mesh is 2D to simplify vtk creation
    if constexpr (SPACE_DIM == 2)
    {
        FindElementOverlaps(rMesh);

        // Make the Vtk mesh
        vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
        p_pts->GetData()->SetName("Vertex positions");

        // Keep a track of the number of extra points that have been added to the purpose of visualisation
        unsigned num_pts_added = 0;

        // Next, we decide how to output the VTK data for elements depending on the type of overlap
        for (auto iter = rMesh.GetElementIteratorBegin(); iter != rMesh.GetElementIteratorEnd(); ++iter)
        {
            unsigned elem_idx = iter->GetIndex();
            unsigned num_nodes = iter->GetNumNodes();

            // Case 1: no overlap
            if (mElementParts.size() == 0 || mElementParts[elem_idx].empty())
            {
                vtkCell* p_cell = vtkPolygon::New();
                vtkIdList* p_cell_id_list = p_cell->GetPointIds();
                p_cell_id_list->SetNumberOfIds(num_nodes);

                for (unsigned node_local_idx = 0; node_local_idx < num_nodes; ++node_local_idx)
                {
                    // Get node, index and location
                    unsigned global_idx = iter->GetNode(node_local_idx)->GetIndex();
                    const auto& r_location = iter->GetNode(node_local_idx)->rGetLocation();

                    p_pts->InsertPoint(global_idx, r_location[0], r_location[1], 0.0);
                    p_cell_id_list->SetId(node_local_idx, global_idx);
                }

                mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
                p_cell->Delete(); // Reference counted
            }
            else  // Case 2: at least one overlap
            {
                // Get the node index at the start of each part
                const std::vector<unsigned>& start_idx_each_part = mElementParts[iter->GetIndex()];

                // Get the number of parts, and put the first index onto the back of the vector for convenience
                const auto num_parts = start_idx_each_part.size();

                for (unsigned part = 0; part < num_parts; ++part)
                {
                    const long this_start = start_idx_each_part[part];
                    const long next_start = start_idx_each_part[AdvanceMod(part, 1, num_parts)];

                    const long num_nodes_this_part = next_start > this_start ? next_start - this_start : num_nodes + next_start - this_start;

                    // Identify the extra points that need to be added
                    std::vector<c_vector<double, SPACE_DIM>> extra_locations;

                    // Start of the part
                    {
                        const auto& r_start_pos = iter->GetNode(this_start)->rGetLocation();
                        const auto& r_end_pos = iter->GetNode(AdvanceMod(this_start, -1, num_nodes))->rGetLocation();

                        const c_vector<double, SPACE_DIM> vec_a2b = rMesh.GetVectorFromAtoB(r_start_pos, r_end_pos);
                        extra_locations.emplace_back(GetIntersectionOfEdgeWithBoundary(r_start_pos, r_start_pos + vec_a2b));
                    }

                    // End of the part
                    {
                        const auto& r_start_pos = iter->GetNode(AdvanceMod(next_start, -1, num_nodes))->rGetLocation();
                        const auto& r_end_pos = iter->GetNode(next_start)->rGetLocation();

                        const c_vector<double, SPACE_DIM> vec_a2b = rMesh.GetVectorFromAtoB(r_start_pos, r_end_pos);
                        extra_locations.emplace_back(GetIntersectionOfEdgeWithBoundary(r_start_pos, r_start_pos + vec_a2b));
                    }

                    // If the two additional points are not on the same edge of the boundary, we also need to add a corner
                    if (std::fabs(extra_locations.front()[0] - extra_locations.back()[0]) > DBL_EPSILON &&
                        std::fabs(extra_locations.front()[1] - extra_locations.back()[1]) > DBL_EPSILON)
                    {
                        extra_locations.emplace_back(GetNearestCorner(extra_locations.front(), extra_locations.back()));
                    }

                    vtkCell* p_cell = vtkPolygon::New();
                    vtkIdList* p_cell_id_list = p_cell->GetPointIds();
                    p_cell_id_list->SetNumberOfIds(num_nodes_this_part + extra_locations.size());

                    for (long idx = 0; idx < num_nodes_this_part; ++idx)
                    {
                        unsigned node_idx = AdvanceMod(idx, this_start, num_nodes);

                        // Get index and location
                        unsigned global_idx = iter->GetNode(node_idx)->GetIndex();
                        const auto& r_location = iter->GetNode(node_idx)->rGetLocation();

                        p_pts->InsertPoint(global_idx, r_location[0], r_location[1], 0.0);
                        p_cell_id_list->SetId(idx, global_idx);
                    }

                    // Now add the extra locations
                    if (extra_locations.size() == 2)
                    {
                        const unsigned global_index = rMesh.GetNumNodes() + num_pts_added;

                        p_pts->InsertPoint(global_index, extra_locations[1][0], extra_locations[1][1], 0.0);      // end
                        p_pts->InsertPoint(global_index + 1, extra_locations[0][0], extra_locations[0][1], 0.0);  // start

                        p_cell_id_list->SetId(num_nodes_this_part, global_index);
                        p_cell_id_list->SetId(num_nodes_this_part + 1, global_index + 1);

                        num_pts_added += 2;
                    }
                    else if (extra_locations.size() == 3)
                    {
                        const unsigned global_index = rMesh.GetNumNodes() + num_pts_added;

                        p_pts->InsertPoint(global_index, extra_locations[1][0], extra_locations[1][1], 0.0);      // end
                        p_pts->InsertPoint(global_index + 1, extra_locations[2][0], extra_locations[2][1], 0.0);  // corner
                        p_pts->InsertPoint(global_index + 2, extra_locations[0][0], extra_locations[0][1], 0.0);  // start

                        p_cell_id_list->SetId(num_nodes_this_part, global_index);
                        p_cell_id_list->SetId(num_nodes_this_part + 1, global_index + 1);
                        p_cell_id_list->SetId(num_nodes_this_part + 2, global_index + 2);

                        num_pts_added += 3;
                    }
                    else
                    {
                        NEVER_REACHED;
                    }

                    mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
                    p_cell->Delete(); // Reference counted
                }
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
    }
    else
    {
        NEVER_REACHED;
    }
#endif //CHASTE_VTK
}

//LCOV_EXCL_START
/**
 * Template instnantiation for unused code path
 * @param rMesh the mesh
 */
template <>
void ImmersedBoundaryMeshWriter<1, 1>::MakeVtkMesh(ImmersedBoundaryMesh<1, 1>& rMesh)
{
}
//LCOV_EXCL_STOP

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
    for (double scalar : dataPayload)
    {
        p_scalars->InsertNextValue(scalar);
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(ImmersedBoundaryMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
{
    this->mpMeshReader = nullptr;
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
            for (unsigned node_index : node_indices)
            {
                *p_element_file << "\t" << node_index;
            }

            *p_element_file << "\t" << elem_data.AttributeValue;

            // Has fluid source
            *p_element_file << "\t" << elem_data.hasFluidSource;

            // Fluid source index
            if (elem_data.hasFluidSource) {
                *p_element_file << "\t" << elem_data.fluidSourceIndex;
            }

            // Corner node indices
            *p_element_file << "\t" << elem_data.cornerNodeIndices.size();

            for (auto nodeIndex : elem_data.cornerNodeIndices) {
                *p_element_file << "\t" << nodeIndex;
            }

            // Average node spacing
            *p_element_file << "\t" << elem_data.averageNodeSpacing;

            // Is boundary element
            *p_element_file << "\t" << elem_data.isBoundaryElement;

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

            // Has fluid source
            *p_lamina_file << "\t" << lamina_data.hasFluidSource;

            // Fluid source index
            if (lamina_data.hasFluidSource) {
                *p_lamina_file << "\t" << lamina_data.fluidSourceIndex; //LCOV_EXCL_LINE
            }

            // Corner node indices
            *p_lamina_file << "\t" << lamina_data.cornerNodeIndices.size();

            for (auto nodeIndex : lamina_data.cornerNodeIndices) {
                *p_lamina_file << "\t" << nodeIndex; //LCOV_EXCL_LINE
            }

            // Average node spacing
            *p_lamina_file << "\t" << lamina_data.averageNodeSpacing;

            // Is boundary lamina
            *p_lamina_file << "\t" << lamina_data.isBoundaryElement;


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
void ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::FindElementOverlaps(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    if constexpr (SPACE_DIM == 2)
    {
        // Resize and initialise the vector of overlaps (bools; whether each element has an overlap)
        mElementParts.resize(rMesh.GetNumAllElements());
        for (auto& parts : mElementParts)
        {
            parts.clear();
        }

        // We loop over each element and the node index at each discontinuity due to periodic boundaries
        for (auto iter = rMesh.GetElementIteratorBegin(); iter != rMesh.GetElementIteratorEnd(); ++iter)
        {
            for (unsigned node_idx = 0; node_idx < iter->GetNumNodes(); ++node_idx)
            {
                const unsigned prev_idx = AdvanceMod(node_idx, -1, iter->GetNumNodes());

                const auto& this_location = iter->GetNode(node_idx)->rGetLocation();
                const auto& prev_location = iter->GetNode(prev_idx)->rGetLocation();

                if (norm_inf(this_location - prev_location) > 0.5)
                {
                    mElementParts[iter->GetIndex()].emplace_back(node_idx);
                }
            }
        }
    }
    else
    {
        NEVER_REACHED;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetIntersectionOfEdgeWithBoundary(
        const c_vector<double, SPACE_DIM>& rStart,
        const c_vector<double, SPACE_DIM>& rEnd)
{
    // Note that this function relies on SPACE_DIM == 2

    // Turn the c_vector start and end into a boost::geometry segment object
    geom_point start(rStart[0], rStart[1]);
    geom_point end(rEnd[0], rEnd[1]);
    geom_segment edge(start, end);

    // Identify which boundary edge the intersection is with
    std::vector<geom_point> intersections;
    for (const auto& boundary_edge : mBoundaryEdges)
    {
        if (boost::geometry::intersects(edge, boundary_edge))
        {
            boost::geometry::intersection(edge, boundary_edge, intersections);
            break;
        }
    }

    // There should be exactly one intersection
    if (intersections.size() != 1)
    {
        NEVER_REACHED;
    }

    return Create_c_vector(intersections.front().get<0>(), intersections.front().get<1>());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNearestCorner(
        const c_vector<double, SPACE_DIM>& rA, const c_vector<double, SPACE_DIM>& rB) const
{
    // Identify the nearest corner to the average location of the two vectors
    const double x = 0.5 * (rA[0] + rB[0]) < 0.5 ? 0.0 : 1.0;
    const double y = 0.5 * (rA[1] + rB[1]) < 0.5 ? 0.0 : 1.0;

    return Create_c_vector(x, y);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<std::vector<unsigned>>& ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM>::rGetElementParts() const
{
    return mElementParts;
}

// Explicit instantiation
template class ImmersedBoundaryMeshWriter<1,1>;
template class ImmersedBoundaryMeshWriter<1,2>;
template class ImmersedBoundaryMeshWriter<1,3>;
template class ImmersedBoundaryMeshWriter<2,2>;
template class ImmersedBoundaryMeshWriter<2,3>;
template class ImmersedBoundaryMeshWriter<3,3>;
