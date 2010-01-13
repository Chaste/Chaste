/*

Copyright (C) University of Oxford, 2005-2010

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
#include "VertexMeshReader.hpp"
#include "Exception.hpp"

#include <sstream>


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::VertexMeshReader(std::string pathBaseName)
    : mFilesBaseName(pathBaseName),
      mIndexFromZero(false), // initially assume that nodes are not numbered from zero
      mNumNodes(0),
      mNumElements(0),
      mNodesRead(0),
      mElementsRead(0),
      mNumElementAttributes(0)
{
    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes() const
{
    return mNumElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    /// \todo Implement this method
    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{
    /// \todo Implement this method
    ElementData ret;
    ret.NodeIndices = std::vector<unsigned>();
    ret.AttributeValue = 0;
    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumEdges() const
{
    /// \todo Implement this method
    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
    CloseFiles();
    OpenFiles();
    ReadHeaders();

    mNodesRead = 0;
    mElementsRead = 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    std::vector<double> node_data;

    std::string buffer;
    GetNextLineFromStream(mNodesFile, buffer);

    std::stringstream buffer_stream(buffer);

    unsigned index;
    buffer_stream >> index;

    unsigned offset = mIndexFromZero ? 0 : 1;
    if (index != mNodesRead + offset)
    {
        std::stringstream error;
        error << "Data for node " << mNodesRead << " missing";
        EXCEPTION(error.str());
    }

    double node_value;
    for (unsigned i=0; i<3; i++)
    {
        buffer_stream >> node_value;
        node_data.push_back(node_value);
    }

    mNodesRead++;
    return node_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
    ElementData element_data;

    std::string buffer;
    GetNextLineFromStream(mElementsFile, buffer);

    std::stringstream buffer_stream(buffer);

    unsigned element_index;
    buffer_stream >> element_index;

    unsigned offset = mIndexFromZero ? 0 : 1;
    if (element_index != mElementsRead + offset)
    {
        std::stringstream error;
        error << "Data for element " << mElementsRead << " missing";
        EXCEPTION(error.str());
    }

    unsigned num_nodes_in_element;
    buffer_stream >> num_nodes_in_element;

    unsigned node_index;
    for (unsigned i=0; i<num_nodes_in_element; i++)
    {
        buffer_stream >> node_index;
        element_data.NodeIndices.push_back(node_index - offset);
    }

    if (mNumElementAttributes > 0)
    {
        assert(mNumElementAttributes==1);

        unsigned attribute_value;
        buffer_stream >> attribute_value;
        element_data.AttributeValue = attribute_value;
    }
    else
    {
        element_data.AttributeValue = 0;
    }

    mElementsRead++;
    return element_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenFiles()
{
    OpenNodeFile();
    OpenElementsFile();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenNodeFile()
{
    // Nodes definition
    std::string file_name = mFilesBaseName + ".node";
    mNodesFile.open(file_name.c_str());
    if (!mNodesFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenElementsFile()
{
    // Elements definition
    std::string file_name;
    file_name = mFilesBaseName + ".cell";

    mElementsFile.open(file_name.c_str());
    if (!mElementsFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadHeaders()
{
    std::string buffer;

    GetNextLineFromStream(mNodesFile, buffer);
    std::stringstream buffer_stream(buffer);
    buffer_stream >> mNumNodes >> mNumNodeAttributes;

    // Get the next line to see if nodes are indexed from zero or not
    GetNextLineFromStream(mNodesFile, buffer);
    std::stringstream node_buffer_stream(buffer);

    unsigned first_index;
    node_buffer_stream >> first_index;
    assert(first_index == 0 || first_index == 1);
    mIndexFromZero = (first_index == 0);

    // Close, reopen, skip header
    mNodesFile.close();
    OpenNodeFile();
    GetNextLineFromStream(mNodesFile, buffer);

    GetNextLineFromStream(mElementsFile, buffer);
    std::stringstream element_buffer_stream(buffer);

    element_buffer_stream >> mNumElements >> mNumElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::CloseFiles()
{
    mNodesFile.close();
    mElementsFile.close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextLineFromStream(std::ifstream& fileStream, std::string& rawLine)
{
    bool line_is_blank;

    do
    {
        getline(fileStream, rawLine);

        if (fileStream.eof())
        {
            EXCEPTION("Cannot get the next line from node or element file due to incomplete data");
        }

        // Get rid of any comment
        rawLine = rawLine.substr(0,rawLine.find('#', 0));

        line_is_blank = (rawLine.find_first_not_of(" \t", 0) == std::string::npos);
    }
    while (line_is_blank);
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VertexMeshReader<1,1>;
template class VertexMeshReader<1,2>;
template class VertexMeshReader<1,3>;
template class VertexMeshReader<2,2>;
template class VertexMeshReader<2,3>;
template class VertexMeshReader<3,3>;
