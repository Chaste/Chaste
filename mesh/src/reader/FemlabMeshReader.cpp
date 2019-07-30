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

#include "FemlabMeshReader.hpp"
#include "Exception.hpp"

#include <sstream>

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::FemlabMeshReader(const std::string& rPathBaseName,
                                                           const std::string& rNodeFileName,
                                                           const std::string& rElementFileName,
                                                           const std::string& rEdgeFileName)
{

    // Open node file and store the lines as a vector of strings (minus the comments)
    std::string node_file_name = rPathBaseName + rNodeFileName;
    this->mNodeRawData = this->GetRawDataFromFile(node_file_name);

    // Read the node data using TokenizeStringsToDoubles method
    this->mNodeData = TokenizeStringsToDoubles(this->mNodeRawData);

    // Initialise iterator for public GetNextNode method
    this->mpNodeIterator = this->mNodeData.begin();


    // Open element file and store the lines as a vector of strings (minus the comments)
    std::string element_file_name = rPathBaseName + rElementFileName;
    this->mElementRawData = this->GetRawDataFromFile(element_file_name);

    // Read the rest of the element data using TokenizeStringsToInts method
    this->mElementData = TokenizeStringsToInts(this->mElementRawData, SPACE_DIM + 1);
    this->mpElementIterator = this->mElementData.begin();

    /*
     * Open edge file and store the lines as a vector of strings (minus the comments)
     * We store edges as "faces" but the superclass
     * provides a GetNextEdgeData method which queries this data.
     */

    std::string edge_file_name = rPathBaseName + rEdgeFileName;
    this->mFaceRawData = this->GetRawDataFromFile(edge_file_name);

    // Read the rest of the face/edge data using TokenizeStringsToInts method
    this->mFaceData = TokenizeStringsToInts(this->mFaceRawData, SPACE_DIM);
    this->mpFaceIterator = this->mFaceData.begin();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::~FemlabMeshReader()
{}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector < std::vector<double> >
    FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToDoubles(const std::vector<std::string>& rRawData)
{
    std::vector < std::vector < double > >tokenized_data; // Output

    // Iterate over the lines of input
    unsigned dimension_count = 0;
    std::vector < std::string >::const_iterator the_iterator;
    for (the_iterator = rRawData.begin(); the_iterator != rRawData.end();
         the_iterator++)
    {
        const std::string& r_line_of_data = *the_iterator;
        std::stringstream line_stream (r_line_of_data);

        if (dimension_count == 0)
        {
            // First iteration, build the tokenized_data vector and push in x coordinates
            while (!line_stream.eof())
            {
                double item_coord;

                std::vector < double >x_coord;
                line_stream >> item_coord;
                x_coord.push_back (item_coord);
                tokenized_data.push_back (x_coord);
            }
        }
        else
        {
            unsigned current_node = 0;

            // Other iterations, push in coordinates other than x
            while (!line_stream.eof())
            {
                double item_coord;
                line_stream >> item_coord;
                tokenized_data[current_node].push_back (item_coord);
                current_node++;
            }
        }
        // Dimension of mesh is the same as the line of rawData
        dimension_count++;
    }

    if (SPACE_DIM != dimension_count)
    {
        EXCEPTION("SPACE_DIM  != dimension read from file");
    }
    return (tokenized_data);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector < std::vector < unsigned > >
    FemlabMeshReader<ELEMENT_DIM, SPACE_DIM>::TokenizeStringsToInts(const std::vector<std::string>& rRawData, unsigned dimensionOfObject)
{
    std::vector < std::vector < unsigned > >tokenized_data;

    // There are dimensionOfObject lines to be read
    for (unsigned i = 0; i < dimensionOfObject; i++)
    {
        const std::string& r_line_of_data = rRawData[i];
        std::stringstream line_stream (r_line_of_data);

        if (i == 0)
        {
            // First iteration, build the tokenized_data vector and push in x coordinates
            while (!line_stream.eof())
            {
                double item_index;

                std::vector < unsigned >first_index;
                line_stream >> item_index;
                first_index.push_back ((unsigned) (item_index - 0.5)); // item indices should be minus 1
                tokenized_data.push_back (first_index);
            }
        }
        else
        {
            unsigned current_node = 0;

            // Other iterations, push in coordinates other than x.
            while (!line_stream.eof())
            {
                double item_index;
                line_stream >> item_index;
                tokenized_data[current_node].
                push_back ((unsigned) (item_index - 0.5));
                current_node++;
            }
        }
    }
    return (tokenized_data);
}


// Explicit instantiation
template class FemlabMeshReader<1,1>;
template class FemlabMeshReader<1,2>;
template class FemlabMeshReader<1,3>;
template class FemlabMeshReader<2,2>;
template class FemlabMeshReader<2,3>;
template class FemlabMeshReader<3,3>;
