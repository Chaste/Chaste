/*

Copyright (C) University of Oxford, 2005-2012

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

#include "AbstractCachedMeshReader.hpp"
#include "Exception.hpp"

#include <fstream>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::AbstractCachedMeshReader()
    : mNumNodeAttributes(0),
      mMaxNodeBdyMarker(0),
      mNumElementNodes(0),
      mNumElementAttributes(0),
      mMaxFaceBdyMarker(0),
      mIndexFromZero(false) // Initially assume that nodes are not numbered from zero
{
    // We have initialized all numeric variables to zero
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::string> AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetRawDataFromFile(
        const std::string& rFileName)
{
    // Open raw data file

    std::vector<std::string> raw_data;
    std::ifstream data_file(rFileName.c_str());

    // Checks that input file has been opened correctly. If not throws an
    // exception that should be caught by the user.
    if (!data_file.is_open())
    {
        EXCEPTION("Could not open data file " + rFileName);
    }

    // Read each line in turn
    std::string raw_line;
    getline(data_file, raw_line);

    while (data_file)
    {
        // Remove comments (everything from a hash to the end of the line)
        // If there is no hash, then hashLocation = string::npos = -1 = 4294967295 = UINT_MAX
        // (so it works with unsigneds but is a little nasty)
        long hash_location = raw_line.find('#', 0);
        if (hash_location >= 0)
        {
            raw_line = raw_line.substr(0, hash_location);
        }
        // Remove blank lines.  This is unnecessary, since the tokenizer will
        // ignore blank lines anyway.
        long not_blank_location = raw_line.find_first_not_of(" \t", 0);
        if (not_blank_location >= 0)
        {
            raw_data.push_back(raw_line);
        }

        // Move onto next line
        getline(data_file, raw_line);
    }

    data_file.close(); // Closes the data file
    return raw_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetMaxNodeIndex()
{
    // Initialize an interator for the vector of nodes
    std::vector<std::vector<unsigned> >::iterator the_iterator;

    unsigned max_node_index = 0; // Nice if it were negative

    for (the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
    {
        std::vector<unsigned> indices = *the_iterator; // the_iterator points at each line in turn

        for (unsigned i = 0; i < ELEMENT_DIM+1; i++)
        {
            if ( indices[i] >  max_node_index)
            {
                max_node_index = indices[i];
            }
        }
    }

    return max_node_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetMinNodeIndex()
{
    // Initialize an interator for the vector of nodes
    std::vector<std::vector<unsigned> >::iterator the_iterator;

    unsigned min_node_index = UINT_MAX; // A large integer

    for (the_iterator = mElementData.begin(); the_iterator < mElementData.end(); the_iterator++)
    {
        std::vector<unsigned> indices = *the_iterator; // the_iterator points at each line in turn

        for (unsigned i = 0; i < ELEMENT_DIM+1; i++)
        {
            if (indices[i] < min_node_index)
            {
                min_node_index = indices[i];
            }
        }
    }

    return min_node_index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    // Checks that there are still some nodes left to read. If not throws an
    // exception that must be caught by the user.
    if (mpNodeIterator == mNodeData.end())
    {
        EXCEPTION("All nodes already got");
    }

    std::vector<double> next_node = *mpNodeIterator;

    mpNodeIterator++;

    return next_node;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
    // Checks that there are still some elements left to read. If not throws an
    // exception that must be caught by the user.
    if (mpElementIterator == mElementData.end())
    {
        EXCEPTION("All elements already got");
    }

    ElementData ret;
    ret.NodeIndices = *mpElementIterator;
    ret.AttributeValue = 0;

    mpElementIterator++;

    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
    mpElementIterator = mElementData.begin();
    mpFaceIterator = mFaceData.begin();
    mpNodeIterator = mNodeData.begin();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{
    // Checks that there are still some faces left to read. If not throws an
    // exception that must be caught by the user.
    if (mpFaceIterator == mFaceData.end())
    {
        EXCEPTION("All faces (or edges) already got");
    }

    ElementData ret;
    ret.NodeIndices = *mpFaceIterator;
    ret.AttributeValue = 0;

    mpFaceIterator++;

    return ret;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElementData.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNodeData.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCachedMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaceData.size();
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractCachedMeshReader<1,1>;
template class AbstractCachedMeshReader<1,2>;
template class AbstractCachedMeshReader<1,3>;
template class AbstractCachedMeshReader<2,2>;
template class AbstractCachedMeshReader<2,3>;
template class AbstractCachedMeshReader<3,3>;
