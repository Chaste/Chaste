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
#include <cassert>
#include <sstream>
#include <iostream>

#include "GmshMeshReader.hpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GmshMeshReader(std::string pathBaseName,
                                                       unsigned orderOfElements,
                                                       unsigned orderOfBoundaryElements) :
       mFileName(pathBaseName),
       mOrderOfElements(orderOfElements),
       mOrderOfBoundaryElements(orderOfBoundaryElements)
{
    // Only linear and quadratic elements
    assert(mOrderOfElements==1 || mOrderOfElements==2);

    if (mOrderOfElements==1)
    {
        mNodesPerElement = ELEMENT_DIM+1;
    }
    else
    {
        mNodesPerElement = (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2;
    }

    if (mOrderOfBoundaryElements==1)
    {
        mNodesPerBoundaryElement = ELEMENT_DIM;
    }
    else
    {
        mNodesPerBoundaryElement = ELEMENT_DIM*(ELEMENT_DIM+1)/2;
    }

    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::~GmshMeshReader()
{
    CloseFiles();
}

template<unsigned  ELEMENT_DIM, unsigned  SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenFiles()
{
    // Open mesh file
    mNodeFile.open(mFileName.c_str());
    mElementFile.open(mFileName.c_str());
    mFaceFile.open(mFileName.c_str());
    if (!mNodeFile.is_open() || !mElementFile.is_open() || !mFaceFile.is_open() )
    {
        EXCEPTION("Could not open data file: " + mFileName);
    }
}

template<unsigned  ELEMENT_DIM, unsigned  SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::CloseFiles()
{
    mNodeFile.close();
    mElementFile.close();
    mFaceFile.close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadHeaders()
{
    /*
     * Read mesh format information from the file header
     */
    std::string this_line;
    getline(mNodeFile, this_line);

    assert(this_line == "$MeshFormat");

    //Read the version no.
    getline(mNodeFile, this_line);
    std::stringstream line(this_line);

    line >> mVersionNumber >> mFileType >> mDataSize;

    if (mVersionNumber != 2.2)
    {
        EXCEPTION("Only .msh version 2.2 files are supported.");
    }
    assert(mFileType == 0);

    //Check mesh format close string
    getline(mNodeFile, this_line);
    assert(this_line == "$EndMeshFormat");

    ReadNodeHeader();
    ReadElementHeader(); // This reads the total number of elements in the file into mTotalNumElementsAndFaces
    ReadFaceHeader();

    if (mTotalNumElementsAndFaces != mNumElements + mNumFaces)
    {
        EXCEPTION("Unrecognised element types present in the .msh file: check mesh generation settings in gmsh.");
        // If you hit this, then the .msh file lists elements that are not simply
        // in 2D:
        //   linear/quadratic lines for faces, triangles for elements
        // in 3D:
        //   linear/quadratic traingles for faces, tets for elements.
        // which probably means something funny happened in the mesh generation.
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadNodeHeader()
{
    // Search for the start of the node section
    std::string this_line;
    std::stringstream line(this_line);

    while (this_line != "$Nodes")
    {
        getline(mNodeFile, this_line);
    }
    getline(mNodeFile, this_line);

    line.clear();
    line.str(this_line);
    line >> mNumNodes; // mNodesFile should now be pointing at the start of the node lines in the file.
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadElementHeader()
{
    // Search for the start of the elements section
    std::string this_line;
    std::stringstream line(this_line);

    getline(mElementFile, this_line);
    while (this_line != "$Elements")
    {
        getline(mElementFile, this_line);
    }
    getline(mElementFile, this_line); //Throw away the line that says the number of elements specified in the file
    line.clear();
    line.str(this_line);
    line >> mTotalNumElementsAndFaces;
    int ele_start = mElementFile.tellg(); //Pointer to the start of the element block.

    mNumElements = 0u;
    getline(mElementFile, this_line);

    while (this_line != "$EndElements")
    {
        line.clear();
        line.str(this_line);

        unsigned ele_index;
        unsigned ele_type;
        line >> ele_index >> ele_type;

        if (ELEMENT_DIM == 2 && (ele_type == GmshTypes::TRIANGLE || ele_type == GmshTypes::QUADRATIC_TRIANGLE))
        {
            mNumElements++;
        }
        else if (ELEMENT_DIM == 3 && (ele_type == GmshTypes::TETRAHEDRON || ele_type == GmshTypes::QUADRATIC_TETRAHEDRON))
        {
            mNumElements++;
        }

        getline(mElementFile, this_line);
    }

    mElementFile.seekg(ele_start); //mElementFile should now be pointing at the start of the node lines in the file.
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadFaceHeader()
{
    // Search for the start of the elements section
    std::string this_line;
    std::stringstream line(this_line);

    getline(mFaceFile, this_line);
    while (this_line != "$Elements")
    {
        getline(mFaceFile, this_line);
    }
    getline(mFaceFile, this_line); // Throw away the line that says the number of elements specified in the file
    int face_start = mFaceFile.tellg(); // Pointer to the start of the element block.

    mNumFaces = 0u;
    getline(mFaceFile, this_line);

    while (this_line != "$EndElements")
    {
        line.clear();
        line.str(this_line);

        unsigned ele_index;
        unsigned ele_type;
        line >> ele_index >> ele_type;

        if (ELEMENT_DIM == 2 && (ele_type == GmshTypes::LINE || ele_type == GmshTypes::QUADRATIC_LINE))
        {
            mNumFaces++;
        }
        else if (ELEMENT_DIM == 3 && (ele_type == GmshTypes::TRIANGLE || ele_type == GmshTypes::QUADRATIC_TRIANGLE))
        {
            mNumFaces++;
        }

        getline(mFaceFile, this_line);
    }

    mFaceFile.seekg(face_start); //mFacesFile should now be pointing at the start of the node lines in the file.
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mNumFaces;
}

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements() const
{
    NEVER_REACHED;
    //return mNumCableElements;
}
// LCOV_EXCL_STOP

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetOrderOfElements()
{
    return mOrderOfElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetOrderOfBoundaryElements()
{
    return mOrderOfBoundaryElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes() const
{
    return mNumElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaceAttributes() const
{
    return mNumFaceAttributes;
}

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumCableElementAttributes() const
{
    NEVER_REACHED;
    //return mNumCableElementAttributes;
}
// LCOV_EXCL_STOP

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
    CloseFiles();

    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    std::vector<double> ret_coords(SPACE_DIM);

    std::string this_line;
    getline(mNodeFile, this_line);

    std::stringstream line(this_line);

    unsigned node_index;
    line >> node_index >> ret_coords[0] >> ret_coords[1];

    if (SPACE_DIM == 3)
    {
        line >> ret_coords[2];
    }

    return ret_coords;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNodeAttributes()
{
    return std::vector<double>(0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
    ElementData element_data;
    element_data.NodeIndices.resize(mNodesPerElement);
    element_data.AttributeValue = 0.0; // If an attribute is not read this stays as zero, otherwise overwritten.

    std::string this_line;
    std::stringstream line(this_line);

    unsigned ele_index;
    unsigned ele_type;
    unsigned ele_attributes = 0;
    bool volume_element_found = false;

    while (!volume_element_found)
    {
        getline(mElementFile, this_line);
        line.clear();
        line.str(this_line);

        line >> ele_index >> ele_type >> ele_attributes;

        if ((ELEMENT_DIM == 2 && (ele_type == GmshTypes::TRIANGLE || ele_type == GmshTypes::QUADRATIC_TRIANGLE) ) ||
           (ELEMENT_DIM == 3 && (ele_type == GmshTypes::TETRAHEDRON || ele_type == GmshTypes::QUADRATIC_TETRAHEDRON)))
        {
            volume_element_found = true;
        }
    }

    //Gmsh can have arbitrary numbers of element attributes, but Chaste only handles one. We pick the first attribute and throw
    //away the remainder.
    if (ele_attributes > 0)
    {
        mNumElementAttributes = 1u;
        line >> element_data.AttributeValue;
    }
    unsigned unused_attr;
    for (unsigned attr_index = 0; attr_index < (ele_attributes-1); ++attr_index)
    {
        line >> unused_attr;
    }

    //Read the node indices
    for (unsigned node_index = 0; node_index < mNodesPerElement; ++node_index)
    {
        line >> element_data.NodeIndices[node_index];

        element_data.NodeIndices[node_index]--; //Gmsh *always* indexes from 1, we index from 0
    }

    return element_data;
}

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextCableElementData()
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{
    ElementData face_data;
    face_data.NodeIndices.resize(mNodesPerBoundaryElement);
    face_data.AttributeValue = 0.0; // If an attribute is not read this stays as zero, otherwise overwritten.

    std::string this_line;
    std::stringstream line(this_line);

    unsigned face_index;
    unsigned face_type;
    unsigned face_attributes=0;
    bool surface_element_found = false;

    while (!surface_element_found)
    {
        getline(mFaceFile, this_line);
        line.clear();
        line.str(this_line);

        line >> face_index >> face_type >> face_attributes;

        if ((ELEMENT_DIM == 2 && (face_type == GmshTypes::LINE || face_type == GmshTypes::QUADRATIC_LINE) ) ||
           (ELEMENT_DIM == 3 && (face_type == GmshTypes::TRIANGLE || face_type == GmshTypes::QUADRATIC_TRIANGLE)))
        {
            surface_element_found = true;
        }
    }

    //Gmsh can have arbitrary numbers of element attributes, but Chaste only handles one. We pick the first attribute and throw
    //away the remainder.
    if (face_attributes > 0)
    {
        mNumFaceAttributes = 1u;
        line >> face_data.AttributeValue;
    }
    unsigned unused_attr;
    for (unsigned attr_index = 0; attr_index < (face_attributes-1); ++attr_index)
    {
        line >> unused_attr;
    }

    //Read the node indices
    for (unsigned node_index = 0; node_index < mNodesPerBoundaryElement; ++node_index)
    {
        line >> face_data.NodeIndices[node_index];

        face_data.NodeIndices[node_index]--; //Gmsh *always* indexes from 1, we index from 0
    }

    return face_data;
}

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetElementData(unsigned index)
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetFaceData(unsigned index)
{
    NEVER_REACHED;
}
// LCOV_EXCL_STOP

// Explicit instantiation
template class GmshMeshReader<0,1>;
template class GmshMeshReader<1,1>;
template class GmshMeshReader<1,2>;
template class GmshMeshReader<1,3>;
template class GmshMeshReader<2,2>;
template class GmshMeshReader<2,3>;
template class GmshMeshReader<3,3>;
