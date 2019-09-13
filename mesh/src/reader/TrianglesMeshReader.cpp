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

#include "TrianglesMeshReader.hpp"
#include "Exception.hpp"

static const char* NODES_FILE_EXTENSION = ".node";
static const char* ELEMENTS_FILE_EXTENSION = ".ele";
static const char* FACES_FILE_EXTENSION = ".face";
static const char* EDGES_FILE_EXTENSION = ".edge";
static const char* NCL_FILE_EXTENSION = ".ncl";
static const char* CABLE_FILE_EXTENSION = ".cable";

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::TrianglesMeshReader(std::string pathBaseName,
                                                                 unsigned orderOfElements,
                                                                 unsigned orderOfBoundaryElements,
                                                                 bool readContainingElementForBoundaryElements)
    : mFilesBaseName(pathBaseName),
      mNodeItemWidth(0),
      mElementItemWidth(0),
      mFaceItemWidth(0),
      mNumNodes(0),
      mNumElements(0),
      mNumFaces(0),
      mNumCableElements(0),
      mNodesRead(0),
      mElementsRead(0),
      mCableElementsRead(0),
      mFacesRead(0),
      mBoundaryFacesRead(0),
      mNclItemsRead(0),
      mNumNodeAttributes(0),
      mNumElementAttributes(0),
      mNumFaceAttributes(0),
      mNumCableElementAttributes(0),
      mOrderOfElements(orderOfElements),
      mOrderOfBoundaryElements(orderOfBoundaryElements),
      mEofException(false),
      mReadContainingElementOfBoundaryElement(readContainingElementForBoundaryElements),
      mFilesAreBinary(false),
      mMeshIsHexahedral(false),
      mNodeFileReadBuffer(nullptr),
      mElementFileReadBuffer(nullptr),
      mFaceFileReadBuffer(nullptr),
      mNodePermutationDefined(false)
{
    // Only linear and quadratic elements
    assert(orderOfElements==1 || orderOfElements==2);
    if (mOrderOfBoundaryElements == 2 &&  mReadContainingElementOfBoundaryElement)
    {
        EXCEPTION("Boundary element file should not have containing element info if it is quadratic");
    }
    if (mOrderOfElements==1)
    {
        mNodesPerElement = ELEMENT_DIM+1;
    }
    else
    {
        assert(SPACE_DIM==ELEMENT_DIM); // LCOV_EXCL_LINE
        mNodesPerElement = (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2;
    }

    if (mOrderOfBoundaryElements==1)
    {
        mNodesPerBoundaryElement = ELEMENT_DIM;
    }
    else
    {
        assert(SPACE_DIM==ELEMENT_DIM); // LCOV_EXCL_LINE
        mNodesPerBoundaryElement = ELEMENT_DIM*(ELEMENT_DIM+1)/2;
    }

    mIndexFromZero = false; // Initially assume that nodes are not numbered from zero

    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::~TrianglesMeshReader()
{
    delete[] mNodeFileReadBuffer;
    delete[] mElementFileReadBuffer;
    delete[] mFaceFileReadBuffer;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return mNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mNumFaces;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements() const
{
    return mNumCableElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes() const
{
    return mNumElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaceAttributes() const
{
    return mNumFaceAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumCableElementAttributes() const
{
    return mNumCableElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
    CloseFiles();

    mNodesRead = 0;
    mElementsRead = 0;
    mFacesRead = 0;
    mBoundaryFacesRead = 0;
    mCableElementsRead = 0;
    mNclItemsRead = 0;
    mEofException = false;

    OpenFiles();
    ReadHeaders();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    std::vector<double> ret_coords(SPACE_DIM);

    mNodeAttributes.clear(); // clear attributes for this node
    GetNextItemFromStream(mNodesFile, mNodesRead, ret_coords, mNumNodeAttributes, mNodeAttributes);

    mNodesRead++;
    return ret_coords;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
    ElementData element_data;
    element_data.NodeIndices.resize(mNodesPerElement);
    element_data.AttributeValue = 0.0; // If an attribute is not read this stays as zero, otherwise overwritten.

    std::vector<double> element_attributes;
    GetNextItemFromStream(mElementsFile, mElementsRead, element_data.NodeIndices, mNumElementAttributes, element_attributes);

    if (mNumElementAttributes > 0)
    {
        element_data.AttributeValue = element_attributes[0];///only one element attribute registered for the moment
    }

    EnsureIndexingFromZero(element_data.NodeIndices);

    mElementsRead++;

    if (mNodePermutationDefined)
    {
        for (std::vector<unsigned>::iterator node_it = element_data.NodeIndices.begin();
             node_it != element_data.NodeIndices.end();
             ++ node_it)
        {
            assert(*node_it < mPermutationVector.size());
            *node_it =  mPermutationVector[*node_it];
        }
    }

    return element_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextCableElementData()
{
    ElementData element_data;
    element_data.NodeIndices.resize(2u);
    element_data.AttributeValue = 0; // If an attribute is not read this stays as zero, otherwise overwritten.

    std::vector<double> cable_element_attributes;
    GetNextItemFromStream(mCableElementsFile, mCableElementsRead, element_data.NodeIndices, mNumCableElementAttributes, cable_element_attributes);
    if (mNumCableElementAttributes > 0)
    {
        element_data.AttributeValue = cable_element_attributes[0];///only one element attribute registered for the moment
    }

    EnsureIndexingFromZero(element_data.NodeIndices);

    mCableElementsRead++;

    // Node permutation can only be done with binary data...
//    if (mNodePermutationDefined)
//    {
//        for (std::vector<unsigned>::iterator node_it = element_data.NodeIndices.begin();
//             node_it != element_data.NodeIndices.end();
//             ++ node_it)
//        {
//            assert(*node_it < mPermutationVector.size());
//            *node_it =  mPermutationVector[*node_it];
//        }
//    }

    return element_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{
    ElementData face_data;
    std::vector<unsigned> ret_indices;

    // In the first case there's no file, all the nodes are set as faces
    if (ELEMENT_DIM == 1)
    {
        ret_indices.push_back( mOneDimBoundary[mBoundaryFacesRead] );
    }
    else
    {
        ret_indices.resize(mNodesPerBoundaryElement);

        assert(ELEMENT_DIM != 0); // LCOV_EXCL_LINE //Covered in earlier exception, but needed in loop guard here.
        do
        {
            face_data.AttributeValue = 1.0; // If an attribute is not read this stays as one, otherwise overwritten.

            std::vector<double> face_attributes; //will store face attributes, if any
            if (mReadContainingElementOfBoundaryElement)
            {
                assert(mNumFaceAttributes == 0);
                GetNextItemFromStream(mFacesFile, mFacesRead, ret_indices, 1, face_attributes);

                if (face_attributes.size() > 0)
                {
                    face_data.ContainingElement = (unsigned) face_attributes[0];// only one face attribute registered for the moment
                }

            }
            else
            {
                GetNextItemFromStream(mFacesFile, mFacesRead, ret_indices, mNumFaceAttributes,
                                      face_attributes);

                if (mNumFaceAttributes > 0)
                {
                    face_data.AttributeValue = face_attributes[0]; //only one face attribute registered for the moment
                }
            }

            EnsureIndexingFromZero(ret_indices);

            mFacesRead++;
        }
        while (ELEMENT_DIM==2 && face_data.AttributeValue==0.0); //In triangles format we ignore internal edges (which are marked with attribute 0)
    }

    mBoundaryFacesRead++;

    if (mNodePermutationDefined)
    {
        for (std::vector<unsigned>::iterator node_it = ret_indices.begin();
             node_it != ret_indices.end();
             ++ node_it)
        {
            assert(*node_it < mPermutationVector.size());
            *node_it =  mPermutationVector[*node_it];
        }
    }

    face_data.NodeIndices = ret_indices;

    return face_data;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
    if (!mFilesAreBinary)
    {
        EXCEPTION("Random access is only implemented in mesh readers for binary mesh files.");
    }
    if (index >= mNumNodes)
    {
        EXCEPTION("Node does not exist - not enough nodes.");
    }

    if (mNodePermutationDefined)
    {
        assert(index<mInversePermutationVector.size());
        index = mInversePermutationVector[index];
    }

    // Put the file stream pointer to the right location
    if (index > mNodesRead)
    {
        // This is a monotonic (but non-contiguous) read.  Let's assume that it's more efficient
        // to seek from the current position rather than from the start of the file
        mNodesFile.seekg( mNodeItemWidth*(index-mNodesRead), std::ios_base::cur);
    }
    else if (mNodesRead != index)
    {
        mNodesFile.seekg(mNodeFileDataStart + mNodeItemWidth*index, std::ios_base::beg);
    }

    mNodesRead = index; // Allow GetNextNode() to note the position of the item after this one
    // Read the next item.
    return GetNextNode();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetElementData(unsigned index)
{
    if (!mFilesAreBinary)
    {
        EXCEPTION("Random access is only implemented in mesh readers for binary mesh files.");
    }
    if (index >= mNumElements)
    {
        EXCEPTION("Element " << index << " does not exist - not enough elements (only " << mNumElements << ").");
    }

    // Put the file stream pointer to the right location
    if (index > mElementsRead)
    {
        // This is a monotonic (but non-contiguous) read.  Let's assume that it's more efficient
        // to seek from the current position rather than from the start of the file
        mElementsFile.seekg( mElementItemWidth*(index-mElementsRead), std::ios_base::cur);
    }
    else if (mElementsRead != index)
    {
        mElementsFile.seekg(mElementFileDataStart + mElementItemWidth*index, std::ios_base::beg);
    }

    mElementsRead = index; // Allow GetNextElementData() to note the position of the item after this one
    // Read the next item.
    return GetNextElementData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetFaceData(unsigned index)
{
    if (!mFilesAreBinary)
    {
        EXCEPTION("Random access is only implemented in mesh readers for binary mesh files.");
    }
    if (index >=mNumFaces)
    {
        EXCEPTION("Face does not exist - not enough faces.");
    }

    /*
     *
    if (index > mFacesRead)
    {
        // This would be a monotonic (but non-contiguous) read. But we don't actually read faces with this access pattern.
        ///\todo Revisit #1930?
        mFacesFile.seekg( mFaceItemWidth*(index-mFacesRead), std::ios_base::cur);
    }
    else
    */
    // Put the file stream pointer to the right location
    if (mFacesRead != index)
    {
        mFacesFile.seekg(mFaceFileDataStart + mFaceItemWidth*index, std::ios_base::beg);
    }
    mFacesRead = index; // Allow next call to mark the position in the file stream

    // Read the next item
    return GetNextFaceData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetContainingElementIndices(unsigned index)
{
    if (!mFilesAreBinary)
    {
        EXCEPTION("NCL file functionality is only implemented in mesh readers for binary mesh files.");
    }

    if (!mNclFileAvailable)
    {
        EXCEPTION("No NCL file available for this mesh.");
    }
    if (index >= mNumNodes)
    {
        EXCEPTION("Connectivity list does not exist - not enough nodes.");
    }

    if (mNodePermutationDefined)
    {
        assert(index < mInversePermutationVector.size());
        index = mInversePermutationVector[index];
    }

    // Put the file stream pointer to the right location
    if (index > mNclItemsRead)
    {
        // This is a monotonic (but non-contiguous) read.  Let's assume that it's more efficient
        // to seek from the current position rather than from the start of the file
        mNclFile.seekg( mNclItemWidth*(index-mNclItemsRead), std::ios_base::cur);
    }
    else if  ( mNclItemsRead != index )
    {
        mNclFile.seekg(mNclFileDataStart + mNclItemWidth*index, std::ios_base::beg);
    }

    // Read the next item
    std::vector<unsigned> containing_element_indices;
    containing_element_indices.resize(mMaxContainingElements);

    std::vector<double> dummy; // unused here
    GetNextItemFromStream(mNclFile, index, containing_element_indices, 0, dummy);
    mNclItemsRead = index + 1; //Ready for the next call

    EnsureIndexingFromZero(containing_element_indices);

    unsigned num_containing_elements = mMaxContainingElements;
    while ( containing_element_indices[num_containing_elements-1] == UINT_MAX )
    {
        num_containing_elements--;
    }

    containing_element_indices.resize(num_containing_elements);

    return containing_element_indices;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenFiles()
{
    OpenNodeFile();
    OpenElementsFile();
    OpenFacesFile();
    OpenNclFile();
    OpenCableElementsFile();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenNodeFile()
{
    // Nodes definition
    std::string file_name = mFilesBaseName + NODES_FILE_EXTENSION;
    mNodesFile.open(file_name.c_str(), std::ios::binary);
    if (!mNodesFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenElementsFile()
{
    // Elements definition
    std::string file_name;
    if (ELEMENT_DIM == SPACE_DIM)
    {
        file_name = mFilesBaseName + ELEMENTS_FILE_EXTENSION;
    }
    else
    {
        if (ELEMENT_DIM == 1)
        {
            file_name = mFilesBaseName + EDGES_FILE_EXTENSION;
        }
        else if (ELEMENT_DIM == 2)
        {
            file_name = mFilesBaseName + FACES_FILE_EXTENSION;
        }
        else
        {
            EXCEPTION("Can't have a zero-dimensional mesh in a one-dimensional space");
        }
    }

    mElementsFile.open(file_name.c_str(), std::ios::binary);
    if (!mElementsFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenFacesFile()
{
    // Faces/edges definition
    std::string file_name;
    if (ELEMENT_DIM == 3)
    {
        file_name = mFilesBaseName + FACES_FILE_EXTENSION;
    }
    else if (ELEMENT_DIM == 2)
    {
        file_name = mFilesBaseName + EDGES_FILE_EXTENSION;
    }
    else //if (ELEMENT_DIM == 1)
    {
        // There is no file, data will be read from the node file (with boundaries marked)
        return;
    }

    mFacesFile.open(file_name.c_str(), std::ios::binary);
    if (!mFacesFile.is_open())
    {
        EXCEPTION("Could not open data file: " + file_name);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenNclFile()
{
    std::string file_name = mFilesBaseName + NCL_FILE_EXTENSION;
    mNclFile.open(file_name.c_str(), std::ios::binary);

    mNclFileAvailable = mNclFile.is_open();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::OpenCableElementsFile()
{
    std::string file_name = mFilesBaseName + CABLE_FILE_EXTENSION;
    mCableElementsFile.open(file_name.c_str(), std::ios::binary);
    if (!mCableElementsFile.is_open())
    {
        mNumCableElements = 0u;
        mNumCableElementAttributes = 0u;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNodeAttributes()
{
    return mNodeAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadHeaders()
{
    /*
     * Reading node file header
     */
    std::string buffer;
    GetNextLineFromStream(mNodesFile, buffer);
    std::stringstream node_header_line(buffer);
    unsigned dimension;
    node_header_line >> mNumNodes >> dimension >> mNumNodeAttributes >> mMaxNodeBdyMarker;
    if (SPACE_DIM != dimension)
    {
        EXCEPTION("SPACE_DIM  != dimension read from file ");
    }

    // Is there anything else on the header line?
    std::string extras;
    node_header_line >> extras;
    if (extras == "BIN")
    {
        mFilesAreBinary = true;
        mNodeFileDataStart = mNodesFile.tellg(); // Record the position of the first byte after the header.
        mNodeItemWidth = SPACE_DIM * sizeof(double);

        // We enforce that all binary files (written by Chaste) are indexed from zero
        mIndexFromZero = true;
    }
    else
    {
        // #1621 - ncl files are only supported in binary read mode.
        assert(!mNclFileAvailable);

        // Get the next line to see if it is indexed from zero or not
        GetNextLineFromStream(mNodesFile, buffer);
        std::stringstream node_first_line(buffer);
        unsigned first_index;
        node_first_line >> first_index;
        assert(first_index == 0 || first_index == 1);
        mIndexFromZero = (first_index == 0);

        // Close, reopen, skip header
        mNodesFile.close();
        OpenNodeFile();
        GetNextLineFromStream(mNodesFile, buffer);
    }

    /*
     * Reading element file header
     */
    GetNextLineFromStream(mElementsFile, buffer);
    std::stringstream element_header_line(buffer);

    unsigned extra_attributes = 0;

    if (ELEMENT_DIM == SPACE_DIM)
    {
        element_header_line >> mNumElements >> mNumElementNodes >> mNumElementAttributes;

        extra_attributes = mNumElementAttributes;

        // Is there anything else on the header line?
        std::string element_extras;
        element_header_line >> element_extras;
        if (element_extras == "BIN")
        {
            // Double check for binaryness
            assert (mFilesAreBinary);
        }
        else if (element_extras == "HEX")
        {
            mMeshIsHexahedral = true;
            if (ELEMENT_DIM == 2)
            {
                mNodesPerElement = 4;
                mNodesPerBoundaryElement = 2;
            }
            if (ELEMENT_DIM == 3)
            {
                mNodesPerElement = 8;
                mNodesPerBoundaryElement = 4;
            }
        }
        else
        {
            assert (element_extras == "");
        }

        // The condition below allows the writer to cope with a NodesOnlyMesh
        if (mNumElements != 0)
        {
            if (mNumElementNodes != mNodesPerElement)
            {
                EXCEPTION("Number of nodes per elem, " << mNumElementNodes << ", does not match "
                      << "expected number, " << mNodesPerElement << " (which is calculated given "
                      << "the order of elements chosen, " << mOrderOfElements << " (1=linear, 2=quadratics)");
            }
        }
    }
    else
    {
        // Note that .face files don't have the number of nodes in a face element in the header (its element dim +1)
        element_header_line >> mNumElements >> mNumFaceAttributes;

        extra_attributes = mNumFaceAttributes;

        if (ELEMENT_DIM == 1 || ELEMENT_DIM == 2)
        {
            mNumElementAttributes = mNumFaceAttributes;
        }

        // Is there anything else on the header line?
        std::string element_extras;
        element_header_line >> element_extras;
        if (element_extras == "BIN")
        {
            // Double check for binaryness
            assert (mFilesAreBinary);
        }

        mNodesPerElement = ELEMENT_DIM+1;
    }

    if (mFilesAreBinary)
    {
        mElementFileDataStart = mElementsFile.tellg(); // Record the position of the first byte after the header.
        mElementItemWidth = mNodesPerElement*sizeof(unsigned) +  extra_attributes*sizeof(double);
    }

    /*
     * Reading face/edge file header.
     * The condition below allows the writer to cope with a NodesOnlyMesh.
     */
    if (mNumElements != 0)
    {
        if (ELEMENT_DIM == 1)
        {
           GetOneDimBoundary();
           mNumFaces = mOneDimBoundary.size();
        }
        else
        {
            GetNextLineFromStream(mFacesFile, buffer);
            std::stringstream face_header_line(buffer);

            face_header_line >> mNumFaces >> mNumFaceAttributes;
            assert(mNumFaceAttributes==0 || mNumFaceAttributes==1);

            /*
             * If mNumFaceAttributes=1 then loop over and set mNumFaces to be
             * the number of faces which are marked as boundary faces.
             * Double check for binaryness.
             */
            std::string face_extras;
            face_header_line >> face_extras;
            assert (mFilesAreBinary == (face_extras == "BIN"));
            if (mNumFaceAttributes==1)
            {
                unsigned num_boundary_faces = 0;
                bool end_of_file=false;
                while (!end_of_file)
                {
                    try
                    {
                        GetNextFaceData();
                        num_boundary_faces++;
                    }
                    catch(Exception& e)
                    {
                        if (mEofException)
                        {
                            end_of_file = true;
                        }
                        else
                        {
                            throw e;
                        }
                    }
                }
                mNumFaces = num_boundary_faces;

    //// This exception would be helpful to have until #1116 is done, unfortunately some meshes do
    //// actually have no boundary elements (eg closed 2d meshes in 3d space).
    //            if (mNumFaces==0)
    //            {
    //                EXCEPTION("No boundary elements found. NOTE: elements in face/edge file with an attribute value of 0 are considered to be internal (non-boundary) elements");
    //            }

                // close the file, reopen, and skip the header again
                mFacesFile.close();
                mFacesFile.clear(); // Older versions of gcc don't explicitly reset "fail" and "eof" flags in std::ifstream after calling close()
                OpenFacesFile();
                GetNextLineFromStream(mFacesFile, buffer);
                mFacesRead = 0;
                mBoundaryFacesRead = 0;
            }
        }
    }

    if (mFilesAreBinary)
    {
        mFaceFileDataStart = mFacesFile.tellg(); // Record the position of the first byte after the header.
        mFaceItemWidth = ELEMENT_DIM*sizeof(unsigned) + mNumFaceAttributes*sizeof(double);
    }

    /*
     * Read NCL file (if one is available)
     */
    if (mNclFileAvailable)
    {
        GetNextLineFromStream(mNclFile, buffer);
        std::stringstream ncl_header_line(buffer);
        unsigned num_nodes_in_file;
        ncl_header_line >> num_nodes_in_file >> mMaxContainingElements;

        if (mNumNodes != num_nodes_in_file)
        {
            EXCEPTION("NCL file does not contain the correct number of nodes for mesh");
        }

        mNclFileDataStart = mNclFile.tellg(); // Record the position of the first byte after the header
        mNclItemWidth = mMaxContainingElements * sizeof(unsigned);
    }

    /*
     * Read cable file (if one is available)
     */
    if (mCableElementsFile.is_open())
    {
        GetNextLineFromStream(mCableElementsFile, buffer);
        std::stringstream cable_header_line(buffer);
        unsigned num_nodes_per_cable_element;
        cable_header_line >> mNumCableElements >> num_nodes_per_cable_element >> mNumCableElementAttributes;
        assert(num_nodes_per_cable_element == 2u);
        mCableElementsRead = 0u;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::CloseFiles()
{
    mNodesFile.close();
    mElementsFile.close();
    mFacesFile.close();
    mNclFile.close();
    mCableElementsFile.close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextLineFromStream(std::ifstream& rFileStream, std::string& rRawLine)
{
    bool line_is_blank;
    mEofException = false;
    do
    {
        getline(rFileStream, rRawLine, '\n');
        if (rFileStream.eof())
        {
            mEofException = true;
            EXCEPTION("File contains incomplete data: unexpected end of file.");
        }

        // Get rid of any comment
        rRawLine = rRawLine.substr(0, rRawLine.find('#',0));

        line_is_blank = (rRawLine.find_first_not_of(" \t",0) == std::string::npos);
    }
    while (line_is_blank);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
template<class T_DATA>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextItemFromStream(std::ifstream& rFileStream, unsigned expectedItemNumber,
                               std::vector<T_DATA>& rDataPacket, const unsigned& rNumAttributes, std::vector<double>& rAttributes)
{
    if (mFilesAreBinary)
    {
        if (!rDataPacket.empty()) // Avoid MSVC 10 assertion
        {
            rFileStream.read((char*)&rDataPacket[0], rDataPacket.size()*sizeof(T_DATA));
        }
        if (rNumAttributes > 0)
        {
            for (unsigned i = 0; i < rNumAttributes; i++)
            {
                double attribute;
                rFileStream.read((char*) &attribute, sizeof(double));
                rAttributes.push_back(attribute);
            }
        }
    }
    else
    {
        std::string buffer;
        GetNextLineFromStream(rFileStream, buffer);
        std::stringstream buffer_stream(buffer);

        unsigned item_index;
        buffer_stream >> item_index;

        // If we are indexing from zero our expected item number is one larger
        expectedItemNumber += mIndexFromZero ? 0 : 1;

        if (item_index != expectedItemNumber)
        {
            if (!mIndexFromZero)
            {
                // To fix the exception message to agree with file format
                expectedItemNumber--;
            }
            EXCEPTION("Data for item " << expectedItemNumber << " missing");
        }

        for (unsigned i=0; i<rDataPacket.size(); i++)
        {
            buffer_stream >> rDataPacket[i];
        }

        if (rNumAttributes > 0)
        {
            for (unsigned i = 0; i < rNumAttributes; i++)
            {
                double attribute;
                buffer_stream >> attribute;
                if (buffer_stream.fail())
                {
                    EXCEPTION("Error in reading attribute index " << i << " (out of " << rNumAttributes << ") in one of the files in " << mFilesBaseName);
                }
                rAttributes.push_back(attribute);
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetMeshFileBaseName()
{
    return mFilesBaseName;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::GetOneDimBoundary()
{
    assert(ELEMENT_DIM == 1);    // LCOV_EXCL_LINE
    mNumFaceAttributes = 0;
    if (!mOneDimBoundary.empty())
    {
        // We have already read this and have reset the reader (probably from the mesh class)
        return;
    }
    std::vector<unsigned> node_indices(2);
    std::vector<double> dummy_attribute; // unused

    // Count how many times we see each node
    std::vector<unsigned> node_count(mNumNodes); // Covers the case if it's indexed from 1
    for (unsigned element_index=0; element_index<mNumElements;element_index++)
    {
        GetNextItemFromStream(mElementsFile, element_index, node_indices, mNumElementAttributes, dummy_attribute);
        if (!mIndexFromZero)
        {
            // Adjust so we are indexing from zero
            node_indices[0]--;
            node_indices[1]--;
        }
        node_count[node_indices[0]]++;
        node_count[node_indices[1]]++;
    }

    // Find the ones which are terminals (only one mention)
    for (unsigned node_index=0; node_index<mNumNodes;node_index++)
    {
        if (node_count[node_index] == 1u)
        {
            mOneDimBoundary.push_back(node_index);
        }
    }

    // Close the file, reopen, and skip the header again
    mElementsFile.close();
    mElementsFile.clear(); // Older versions of gcc don't explicitly reset "fail" and "eof" flags in std::ifstream after calling close()
    OpenElementsFile();
    std::string buffer;
    GetNextLineFromStream(mElementsFile, buffer);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::EnsureIndexingFromZero(std::vector<unsigned>& rNodeIndices)
{
    if (!mIndexFromZero) // If node indices do not start at zero move them all down one so they do
    {
        for (unsigned i=0; i<rNodeIndices.size(); i++)
        {
            rNodeIndices[i]--;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::IsFileFormatBinary()
{
    return mFilesAreBinary;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::HasNclFile()
{
    return mNclFileAvailable;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::SetReadBufferSize(unsigned bufferSize)
{
    mNodeFileReadBuffer = new char[bufferSize];
    mElementFileReadBuffer = new char[bufferSize];
    mFaceFileReadBuffer = new char[bufferSize];

    mNodesFile.rdbuf()->pubsetbuf(mNodeFileReadBuffer, bufferSize);
    mElementsFile.rdbuf()->pubsetbuf(mElementFileReadBuffer, bufferSize);
    mFacesFile.rdbuf()->pubsetbuf(mFaceFileReadBuffer, bufferSize);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::SetNodePermutation(std::vector<unsigned>& rPermutationVector)
{
    if (!mFilesAreBinary)
    {
        // It would be too inefficient otherwise...
        EXCEPTION("Permuted read can only be used with binary files since it requires random access to the node file.");
    }

    mNodePermutationDefined = true;
    mPermutationVector = rPermutationVector;
    mInversePermutationVector.resize(mPermutationVector.size());
    for (unsigned index=0; index<mPermutationVector.size(); index++)
    {
        mInversePermutationVector[mPermutationVector[index]]=index;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::HasNodePermutation()
{
    return(mNodePermutationDefined);
}

/**
* @return the node permutation if a node permutation has been applied to this reader (or an empty permutation)
*/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned>& TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>::rGetNodePermutation()
{
    return mPermutationVector;
}

// Explicit instantiation
template class TrianglesMeshReader<0,1>;
template class TrianglesMeshReader<1,1>;
template class TrianglesMeshReader<1,2>;
template class TrianglesMeshReader<1,3>;
template class TrianglesMeshReader<2,2>;
template class TrianglesMeshReader<2,3>;
template class TrianglesMeshReader<3,3>;


/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template void TrianglesMeshReader<0,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<0,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<1,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<1,1>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<1,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<1,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<1,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<1,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<2,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<2,2>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<2,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<2,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<3,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<unsigned>&, const unsigned&, std::vector<double>&);
template void TrianglesMeshReader<3,3>::GetNextItemFromStream(std::ifstream&, unsigned, std::vector<double>  &, const unsigned&, std::vector<double>&);
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
