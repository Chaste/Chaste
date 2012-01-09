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
#define COVERAGE_IGNORE
        EXCEPTION("Number of nodes does not match expected number declared in header");
#undef COVERAGE_IGNORE
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
#define COVERAGE_IGNORE
        EXCEPTION("Number of elements does not match expected number declared in header");
#undef COVERAGE_IGNORE
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
        const std::string& line_of_data = *the_iterator;
        std::stringstream line_stream(line_of_data);

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
        const std::string& line_of_data = *the_iterator;
        std::stringstream line_stream(line_of_data);

        if ( readHeader == false || the_iterator != rRawData.begin() )
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

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class MemfemMeshReader<1,1>;
template class MemfemMeshReader<1,2>;
template class MemfemMeshReader<1,3>;
template class MemfemMeshReader<2,2>;
template class MemfemMeshReader<2,3>;
template class MemfemMeshReader<3,3>;
