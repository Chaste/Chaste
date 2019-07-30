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
    /// \todo Implement this method (#1076)
    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{
    /// \todo Implement this method (#1076, #1377)
    ElementData ret;
    ret.NodeIndices = std::vector<unsigned>();
    ret.AttributeValue = 0;
    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumEdges() const
{
    /// \todo Implement this method (#1076)
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
        EXCEPTION("Data for node " << mNodesRead << " missing");
    }

    double node_value;
    for (unsigned i=0; i<SPACE_DIM+1; i++)
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
    // Create data structure for this element
    ElementData element_data;

    std::string buffer;
    GetNextLineFromStream(mElementsFile, buffer);

    std::stringstream buffer_stream(buffer);

    unsigned element_index;
    buffer_stream >> element_index;

    unsigned offset = mIndexFromZero ? 0 : 1;
    if (element_index != mElementsRead + offset)
    {
        EXCEPTION("Data for element " << mElementsRead << " missing");
    }

    unsigned num_nodes_in_element;
    buffer_stream >> num_nodes_in_element;

    // Store node indices owned by this element
    unsigned node_index;
    for (unsigned i=0; i<num_nodes_in_element; i++)
    {
        buffer_stream >> node_index;
        element_data.NodeIndices.push_back(node_index - offset);
    }

    if (mNumElementAttributes > 0)
    {
        assert(mNumElementAttributes == 1);

        buffer_stream >>  element_data.AttributeValue;
    }
    else
    {
        element_data.AttributeValue = 0;
    }

    mElementsRead++;
    return element_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElementData VertexMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementDataWithFaces()
{
    // Create data structure for this element
    VertexElementData element_data;

    std::string buffer;
    GetNextLineFromStream(mElementsFile, buffer);

    std::stringstream buffer_stream(buffer);

    unsigned element_index;
    buffer_stream >> element_index;

    unsigned offset = mIndexFromZero ? 0 : 1;
    if (element_index != mElementsRead + offset)
    {
        EXCEPTION("Data for element " << mElementsRead << " missing");
    }

    // Get number of nodes owned by this element
    unsigned num_nodes_in_element;
    buffer_stream >> num_nodes_in_element;

    // Store node indices owned by this element
    for (unsigned i=0; i<num_nodes_in_element; i++)
    {
        unsigned node_index;
        buffer_stream >> node_index;
        element_data.NodeIndices.push_back(node_index - offset);
    }

    // Get number of faces owned by this element
    unsigned num_faces_in_element;
    buffer_stream >> num_faces_in_element;

    element_data.Faces.resize(num_faces_in_element);
    for (unsigned j=0; j<num_faces_in_element; j++)
    {
        // Create data structure for this face
        ElementData face_data;

        // Get face index
        unsigned face_index;
        buffer_stream >> face_index;
        face_data.AttributeValue = face_index;

        // Get number of nodes owned by this face
        unsigned num_nodes_in_face;
        buffer_stream >> num_nodes_in_face;

        // Store node indices owned by this face
        for (unsigned i=0; i<num_nodes_in_face; i++)
        {
            unsigned node_index_face;
            buffer_stream >> node_index_face;
            face_data.NodeIndices.push_back(node_index_face - offset);
        }

        ///\todo Store face orientations? (#1076/#1377)

        element_data.Faces[j] = face_data;
    }

    //For back compatibility (we always store attributes on elements now)
    element_data.AttributeValue = 0;
    if (mNumElementAttributes > 0)
    {
        assert(mNumElementAttributes==1);

        buffer_stream >> element_data.AttributeValue;
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

///////// Explicit instantiation///////

template class VertexMeshReader<1,1>;
template class VertexMeshReader<1,2>;
template class VertexMeshReader<1,3>;
template class VertexMeshReader<2,2>;
template class VertexMeshReader<2,3>;
template class VertexMeshReader<3,3>;
