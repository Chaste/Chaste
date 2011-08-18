/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "PottsMeshWriter.hpp"
#include "Version.hpp"

/**
 * Convenience collection of iterators, primarily to get compilation to happen.
 */
template<unsigned SPACE_DIM>
struct MeshPottsWriterIterators
{
    /** Iterator over nodes */
    typename AbstractMesh<SPACE_DIM,SPACE_DIM>::NodeIterator* pNodeIter;
    /** Iterator over potts elements */
    typename PottsMesh<SPACE_DIM>::PottsElementIterator* pElemIter;
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned SPACE_DIM>
PottsMeshWriter<SPACE_DIM>::PottsMeshWriter(const std::string& rDirectory,
                                                           const std::string& rBaseName,
                                                           const bool clearOutputDir)
    : AbstractMeshWriter<SPACE_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
      mpMesh(NULL),
      mpIters(new MeshPottsWriterIterators<SPACE_DIM>),
      mpNodeMap(NULL),
      mNodeMapCurrentIndex(0)
{
    mpIters->pNodeIter = NULL;
    mpIters->pElemIter = NULL;
}

template<unsigned SPACE_DIM>
PottsMeshWriter<SPACE_DIM>::~PottsMeshWriter()
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
}

template<unsigned SPACE_DIM>
std::vector<double> PottsMeshWriter<SPACE_DIM>::GetNextNode()
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
        return AbstractMeshWriter<SPACE_DIM,SPACE_DIM>::GetNextNode();
    }
}

template<unsigned SPACE_DIM>
ElementData PottsMeshWriter<SPACE_DIM>::GetNextElement()
{
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
        elem_data.AttributeValue = (*(mpIters->pElemIter))->GetRegion();
        ++(*(mpIters->pElemIter));

        return elem_data;
    }
    else
    {
        return AbstractMeshWriter<SPACE_DIM, SPACE_DIM>::GetNextElement();
    }
}

///\todo Mesh should be const (#1663)
template<unsigned SPACE_DIM>
void PottsMeshWriter<SPACE_DIM>::WriteFilesUsingMesh(PottsMesh<SPACE_DIM>& rMesh)
{
    this->mpMeshReader = NULL;
    mpMesh = &rMesh;

    this->mNumNodes = mpMesh->GetNumNodes();
    this->mNumElements = mpMesh->GetNumElements();

    typedef typename AbstractMesh<SPACE_DIM,SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename PottsMesh<SPACE_DIM>::PottsElementIterator ElemIterType;
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

template<unsigned SPACE_DIM>
void PottsMeshWriter<SPACE_DIM>::WriteFiles()
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
    unsigned num_elements = this->GetNumElements();

    unsigned first_elem_attribute_value = (*(mpIters->pElemIter))->GetRegion();
    if (first_elem_attribute_value != 0)
    {
        num_attr = 1;
    }

    *p_element_file << num_elements << "\t";
    *p_element_file << num_attr << "\n";

    // Write each element's data
    for (unsigned item_num=0; item_num<num_elements; item_num++)
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

        // Write the element attribute if necessary
        if (elem_data.AttributeValue != 0)
        {
            *p_element_file << "\t" << elem_data.AttributeValue;
        }

        // New line
        *p_element_file << "\n";
    }

    *p_element_file << comment << "\n";
    p_element_file->close();
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class PottsMeshWriter<1>;
template class PottsMeshWriter<2>;
template class PottsMeshWriter<3>;
