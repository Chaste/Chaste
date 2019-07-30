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

#include "MemfemMeshReader.hpp"
#include "Exception.hpp"

#include <sstream>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::MemfemMeshReader(const std::string& rPathBaseName)
{
    // Open node file and store the lines as a vector of strings (minus the comments)
    std::string node_file_name = rPathBaseName + ".pts";
    this->mNodeRawData = this->GetRawDataFromFile(node_file_name);

    // Read single line header which is the number of nodes */
    std::stringstream node_header_stream(this->mNodeRawData[0]);
    unsigned num_nodes;
    node_header_stream >> num_nodes;

    // All Memfem data is in 3D
    if (SPACE_DIM != 3  || ELEMENT_DIM != 3)
    {
        EXCEPTION("You have asked to read non-3D data. All Memfem data is in 3D.");
    }

    // Read the rest of the node data using TokenizeStringsToDoubles method
    this->mNodeData = TokenizeStringsToDoubles(this->mNodeRawData);

    // Initialise iterator for public GetNextNode method
    this->mpNodeIterator = this->mNodeData.begin();

    // Check that the size of the data matches the information in the header
    if (num_nodes != this->mNodeData.size())
    {
        // ignored from coverage because otherwise would have to create files
        // for a bad mesh just to test this line
// LCOV_EXCL_START
        EXCEPTION("Number of nodes does not match expected number declared in header");
// LCOV_EXCL_STOP
    }

    // Open element file and store the lines as a vector of strings (minus the comments)
    std::string element_file_name = rPathBaseName + ".tetras";
    this->mElementRawData = this->GetRawDataFromFile(element_file_name);

    // Read single line header which is the number of elements
    std::stringstream element_header_stream(this->mElementRawData[0]);
    unsigned num_elements;
    element_header_stream >> num_elements;

    // Read the rest of the element data using TokenizeStringsToInts method
    this->mElementData = TokenizeStringsToInts(this->mElementRawData, SPACE_DIM+1, true);
    this->mpElementIterator = this->mElementData.begin();

    // Check that the size of the data matches the information in the header
    if (num_elements != this->mElementData.size())
    {
        // Ignored from coverage because otherwise we would have to create files
        // for a bad mesh just to test this line
// LCOV_EXCL_START
        EXCEPTION("Number of elements does not match expected number declared in header");
// LCOV_EXCL_STOP
    }

    // Open boundary face file and store the lines as a vector of strings (minus the comments)
    std::string face_file_name = rPathBaseName + ".tri";
    this->mFaceRawData = this->GetRawDataFromFile(face_file_name);

    // There is no header

    // Read the face/edge data using TokenizeStringsToInts method
    this->mFaceData = TokenizeStringsToInts(this->mFaceRawData, SPACE_DIM, false);
    this->mpFaceIterator = this->mFaceData.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::~MemfemMeshReader()
{}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<double> > MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToDoubles(
    const std::vector<std::string>& rRawData)
{
    std::vector<std::vector<double> > tokenized_data; // Output

    // Iterate over the lines of input
    std::vector<std::string>::const_iterator the_iterator;
    for (the_iterator = rRawData.begin(); the_iterator != rRawData.end(); the_iterator++ )
    {
        const std::string& r_line_of_data = *the_iterator;
        std::stringstream line_stream(r_line_of_data);

        if (the_iterator != rRawData.begin()) // Ignore the header string
        {
            std::vector<double> current_coords;

            // Form the vector which represents the position of this item
            for (unsigned i=0; i<SPACE_DIM; i++)
            {
                double item_coord;
                line_stream >> item_coord;
                current_coords.push_back(item_coord);
            }

            // Put item onto main output vector
            tokenized_data.push_back(current_coords);
        }
    }

    return tokenized_data;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<std::vector<unsigned> > MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts(
    const std::vector<std::string>& rRawData,
    unsigned dimensionOfObject,
    bool readHeader)
{
    std::vector<std::vector<unsigned> > tokenized_data;

    std::vector<std::string>::const_iterator the_iterator;
    for (the_iterator = rRawData.begin(); the_iterator != rRawData.end(); the_iterator++ )
    {
        const std::string& r_line_of_data = *the_iterator;
        std::stringstream line_stream(r_line_of_data);

        if (readHeader == false || the_iterator != rRawData.begin())
        {
            std::vector<unsigned> current_indices;

            for (unsigned i=0; i<dimensionOfObject; i++)
            {
                unsigned item_index;
                line_stream >> item_index;
                // The nodes have been indexed from one so we need to shift the indices
                item_index -= 1;
                current_indices.push_back(item_index);
            }

            tokenized_data.push_back(current_indices);
        }
    }

    return tokenized_data;
}

// Explicit instantiation
template class MemfemMeshReader<1,1>;
template class MemfemMeshReader<1,2>;
template class MemfemMeshReader<1,3>;
template class MemfemMeshReader<2,2>;
template class MemfemMeshReader<2,3>;
template class MemfemMeshReader<3,3>;
