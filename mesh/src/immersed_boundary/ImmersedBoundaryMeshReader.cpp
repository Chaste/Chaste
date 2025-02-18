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

#include "ImmersedBoundaryMeshReader.hpp"
#include "Exception.hpp"

#include <sstream>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryMeshReader(std::string pathBaseName)
    : mFilesBaseName(pathBaseName),
      mIndexFromZero(false), // initially assume that nodes are not numbered from zero
      mNumNodes(0),
      mNumElements(0),
      mNumLaminas(0),
      mNodesRead(0),
      mElementsRead(0),
      mLaminasRead(0),
      mNumElementAttributes(0)
{
    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumLaminas() const
{
    return mNumLaminas;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumGridPtsX() const
{
    return mNumGridPtsX;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumGridPtsY() const
{
    return mNumGridPtsY;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes() const
{
    return mNumElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumLaminaAttributes() const
{
    return mNumLaminaAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetCharacteristicNodeSpacing()
{
    return mCharacteristicNodeSpacing;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
    CloseFiles();
    OpenFiles();
    ReadHeaders();

    mNodesRead = 0;
    mElementsRead = 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    std::vector<double> node_data;

    std::string buffer = "";
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
std::vector<double> ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextGridRow()
{
    std::vector<double> grid_row;
    grid_row.resize(mNumGridPtsX);

    std::string buffer = "";
    GetNextLineFromStream(mGridFile, buffer);

    std::stringstream buffer_stream(buffer);

    double value;
    for (unsigned i=0; i<mNumGridPtsX; i++)
    {
        buffer_stream >> value;
        grid_row[i] = value;
    }

    return grid_row;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElementData ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextImmersedBoundaryElementData()
{
    // Create data structure for this element
    ImmersedBoundaryElementData element_data;

    std::string buffer = "";
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

        buffer_stream >> element_data.AttributeValue;
    }
    else
    {
        element_data.AttributeValue = 0;
    }

    // Immersed boundary specific data
    buffer_stream >> element_data.hasFluidSource;

    if (element_data.hasFluidSource) {
        buffer_stream >> element_data.fluidSourceIndex;
    }

    unsigned num_corner_nodes;
    buffer_stream >> num_corner_nodes;
    for (unsigned i = 0; i < num_corner_nodes; i++) {
        buffer_stream >> node_index;
        element_data.cornerNodeIndices.push_back(node_index);
    }

    buffer_stream >> element_data.averageNodeSpacing;

    buffer_stream >> element_data.isBoundaryElement;

    mElementsRead++;
    return element_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElementData ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextImmersedBoundaryLaminaData()
{
    // Create data structure for this lamina
    ImmersedBoundaryElementData lamina_data;

    std::string buffer = "";
    GetNextLineFromStream(mLaminasFile, buffer);

    std::stringstream buffer_stream(buffer);

    unsigned lamina_idx;
    buffer_stream >> lamina_idx;

    unsigned offset = mIndexFromZero ? 0 : 1;
    if (lamina_idx != mLaminasRead + offset)
    {
        EXCEPTION("Data for lamina " << mLaminasRead << " missing");
    }

    unsigned num_nodes_in_lamina;
    buffer_stream >> num_nodes_in_lamina;

    // Store node indices owned by this element
    unsigned node_index;
    for (unsigned i=0; i<num_nodes_in_lamina; i++)
    {
        buffer_stream >> node_index;
        lamina_data.NodeIndices.push_back(node_index - offset);
    }

    if (mNumLaminaAttributes > 0)
    {
        assert(mNumLaminaAttributes == 1);

        buffer_stream >> lamina_data.AttributeValue;
    }
    else
    {
        lamina_data.AttributeValue = 0;
    }

    // Immersed boundary specific data
    buffer_stream >> lamina_data.hasFluidSource;

    if (lamina_data.hasFluidSource) {
        buffer_stream >> lamina_data.fluidSourceIndex; //LCOV_EXCL_LINE
    }

    unsigned num_corner_nodes;
    buffer_stream >> num_corner_nodes;
    for (unsigned i = 0; i < num_corner_nodes; i++) {
        buffer_stream >> node_index;
        lamina_data.cornerNodeIndices.push_back(node_index);
    }

    buffer_stream >> lamina_data.averageNodeSpacing;

    buffer_stream >> lamina_data.isBoundaryElement;

    mLaminasRead++;
    return lamina_data;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenFiles()
{
    OpenNodeFile();
    OpenElementsFile();
    OpenGridFile();
    OpenLaminasFile();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenNodeFile()
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
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenElementsFile()
{
    // Elements definition
    std::string file_name;
    file_name = mFilesBaseName + ".elem";

    mElementsFile.open(file_name.c_str());
    if (!mElementsFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenLaminasFile()
{
    // Elements definition
    std::string file_name;
    file_name = mFilesBaseName + ".lam";

    mLaminasFile.open(file_name.c_str());
    if (!mLaminasFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenGridFile()
{
    // Grid definition
    std::string file_name;
    file_name = mFilesBaseName + ".grid";

    mGridFile.open(file_name.c_str());
    if (!mGridFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadHeaders()
{
    std::string buffer = "";
    GetNextLineFromStream(mNodesFile, buffer);

    unsigned local_space_dim;

    std::stringstream buffer_stream(buffer);
    buffer_stream >> mNumNodes >> local_space_dim >> mNumNodeAttributes >> mCharacteristicNodeSpacing;

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

    // Element header
    GetNextLineFromStream(mElementsFile, buffer);
    std::stringstream element_buffer_stream(buffer);
    element_buffer_stream >> mNumElements >> mNumElementAttributes;

    // Lamina header
    GetNextLineFromStream(mLaminasFile, buffer);
    std::stringstream lamina_buffer_stream(buffer);
    lamina_buffer_stream >> mNumLaminas >> mNumLaminaAttributes;

    // Grid header
    GetNextLineFromStream(mGridFile, buffer);
    std::stringstream grid_buffer_stream(buffer);
    grid_buffer_stream >> mNumGridPtsX >> mNumGridPtsY;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::CloseFiles()
{
    mNodesFile.close();
    mElementsFile.close();
    mLaminasFile.close();
    mGridFile.close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextLineFromStream(std::ifstream& fileStream, std::string& rawLine)
{
    bool line_is_blank;

    do
    {
        getline(fileStream, rawLine);

        if (fileStream.eof())
        {
            EXCEPTION("Cannot get the next line from node, element, lamina or grid file due to incomplete data");
        }

        // Get rid of any comment
        rawLine = rawLine.substr(0,rawLine.find('#', 0));

        line_is_blank = (rawLine.find_first_not_of(" \t", 0) == std::string::npos);
    }
    while (line_is_blank);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return 0;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
    ElementData ret;
    ret.NodeIndices = std::vector<unsigned>();
    ret.AttributeValue = 0;
    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData ImmersedBoundaryMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{
    ElementData ret;
    ret.NodeIndices = std::vector<unsigned>();
    ret.AttributeValue = 0;
    return ret;
}

// Explicit instantiation
template class ImmersedBoundaryMeshReader<1,1>;
template class ImmersedBoundaryMeshReader<1,2>;
template class ImmersedBoundaryMeshReader<1,3>;
template class ImmersedBoundaryMeshReader<2,2>;
template class ImmersedBoundaryMeshReader<2,3>;
template class ImmersedBoundaryMeshReader<3,3>;
