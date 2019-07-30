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


#ifndef _TRIANGLESMESHREADER_HPP_
#define _TRIANGLESMESHREADER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include "AbstractMeshReader.hpp"

/**
 * Concrete version of the AbstractCachedMeshReader class.
 * Once constructed the public methods of the AbstractCachedMeshReader
 * (std::vector<double> GetNextNode(); etc) can be called to interrogate the
 * data.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshReader : public AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>
{

    friend class TestTrianglesMeshReader;//for testing

private:

    bool mIndexFromZero;            /**< True if input data is numbered from zero, false otherwise */

    std::string mFilesBaseName;     /**< The base name for mesh files. */

    std::ifstream mNodesFile;       /**< The nodes file for the mesh. */
    std::ifstream mElementsFile;    /**< The elements file for the mesh. */
    std::ifstream mFacesFile;       /**< The faces (edges) file for the mesh. */
    std::ifstream mNclFile;         /**< The node connectivity list file for the mesh. */
    std::ifstream mCableElementsFile; /**< The elements file for the mesh. */

    std::streampos mNodeFileDataStart; /**< The start of the binary data*/
    std::streamoff mNodeItemWidth;  /**< The number of bytes in a line of the node file*/
    std::streampos mElementFileDataStart; /**< The start of the binary element data*/
    std::streamoff mElementItemWidth;  /**< The number of bytes in a line of the element file*/
    std::streampos mFaceFileDataStart; /**< The start of the binary face data*/
    std::streamoff mFaceItemWidth;  /**< The number of bytes in a line of the face file*/
    std::streampos mNclFileDataStart; /**< The start of the binary data*/
    std::streamoff mNclItemWidth;  /**< The number of bytes in a line of the node file*/

    unsigned mNumNodes;             /**< Number of nodes in the mesh. */
    unsigned mNumElements;          /**< Number of elements in the mesh. */
    unsigned mNumFaces;             /**< Number of faces in the mesh. */
    unsigned mNumCableElements;     /**< Number of cable elements in the mesh. */

    unsigned mNodesRead;            /**< Number of nodes read in (or the index of the last one read in + 1).  */
    unsigned mElementsRead;         /**< Number of elements read in (or the index of the last one read in + 1). */
    unsigned mCableElementsRead;    /**< Number of cable elements read in. */
    unsigned mFacesRead;            /**< Number of faces read in (or the index of the last one read in + 1).*/
    unsigned mBoundaryFacesRead;    /**< Number of boundary faces read in. */
    unsigned mNclItemsRead;         /**< Number of NCL element items read in */
     std::vector<unsigned> mOneDimBoundary; /**<Indices of nodes which are at the boundary of a 1D mesh*/

    unsigned mNumNodeAttributes;    /**< Is the number of attributes stored at each node. */
    std::vector<double> mNodeAttributes; /**<Will contain the nodal attributes at each node. Cleared and re-filled at each node*/
    unsigned mMaxNodeBdyMarker;     /**< Is the maximum node boundary marker. */
    unsigned mNumElementNodes;      /**< Is the number of nodes per element. */
    unsigned mNumElementAttributes; /**< Is the number of attributes stored for each element. */
    unsigned mNumFaceAttributes;    /**< Is the number of attributes stored for each face. */
    unsigned mNumCableElementAttributes; /**< Is the number of attributes stored for each cable element. */

    unsigned mOrderOfElements;      /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mOrderOfBoundaryElements; /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mNodesPerElement;      /**< The number of nodes contained in each element. */
    unsigned mNodesPerBoundaryElement; /**< The number of nodes in each boundary element. */

    unsigned mMaxContainingElements; /**< The maximum number of elements that any node is contained in. */

    bool mEofException; /**< Set to true when end-of-file exception is thrown (for use in a try-catch) */

    bool mReadContainingElementOfBoundaryElement; /**< Whether to read containing element info for each boundary element (obtaining by doing tetgen with the -nn flag) */
    bool mFilesAreBinary; /**< Whether to read all data as binary (determined by a magic number in the node file header)*/
    bool mMeshIsHexahedral; /**< Whether the mesh is hexahedral (determined by a magic number in the element file header) */
    bool mNclFileAvailable; /**< Whether a ncl file exists */

    char* mNodeFileReadBuffer; /**< Buffer for node file read with std::ifstream */
    char* mElementFileReadBuffer; /**< Buffer for element file read with std::ifstream */
    char* mFaceFileReadBuffer; /**< Buffer for face file read with std::ifstream */

    bool mNodePermutationDefined; /**< Whether to consider a user-defined node permutation when reading a mesh from file.*/
    std::vector<unsigned> mPermutationVector; /**< Permutation to be considered, i-th entry of the vector contains new index for original node i.*/
    std::vector<unsigned> mInversePermutationVector; /**< Permutation inverse, stored for performance reasons.*/

//    /** The containing element for each boundary element (obtaining by doing tetgen with the -nn flag).
//     *  In a std::vector rather than the struct to save space if not read.
//     */
//    std::vector<unsigned> mContainingElementsOfBoundaryElement;
//
//    unsigned mIndexIntoContainingElementsVector; /**< Which index to use when GetNextContainingElementOfBoundaryElement() is called */

public:

    /**
     * Constructor.
     *
     * @param pathBaseName  the base name of the files from which to read the mesh data
     *    (either absolute, or relative to the current directory)
     * @param orderOfElements  the order of each element: 1 for linear, 2 for quadratic (defaults to 1)
     * @param orderOfBoundaryElements the order of each boundary element: 1 for linear, 2 for quadratic (defaults to 1. May
     *  or may not be different to orderOfElements (Note tetgen with the -o2 flag creates quadratic elements but doesn't
     *  create quadratic faces, hence the need for this third parameter)
     * @param readContainingElementsForBoundaryElements Whether to read in the containing element information
     *  for each boundary element (in the .face file if tetgen was run with '-nn').
     */
    TrianglesMeshReader(std::string pathBaseName,
                        unsigned orderOfElements=1,
                        unsigned orderOfBoundaryElements=1,
                        bool readContainingElementsForBoundaryElements=false);

    /**
     * Destructor
     */
    ~TrianglesMeshReader();

    /** @return the number of elements in the mesh */
    unsigned GetNumElements() const;

    /** @return the number of nodes in the mesh */
    unsigned GetNumNodes() const;

    /** @return the number of faces in the mesh (synonym of GetNumEdges()) */
    unsigned GetNumFaces() const;

    /** @return the number of cable elements in the mesh */
    unsigned GetNumCableElements() const;

    /** @return the number of attributes in the mesh */
    unsigned GetNumElementAttributes() const;

    /** @return the number of attributes in the mesh */
    unsigned GetNumFaceAttributes() const;

    /** @return the number of cable element attributes in the mesh */
    unsigned GetNumCableElementAttributes() const;

    /** Resets pointers to beginning*/
    void Reset();

    /** @return a vector of the coordinates of each node in turn */
    std::vector<double> GetNextNode();

    /** @return a vector of the nodes of each element (and any attribute information, if there is any) in turn */
    ElementData GetNextElementData();

    /** @return a vector of the nodes of each face in turn (synonym of GetNextEdgeData()) */
    ElementData GetNextFaceData();

    /** @return a vector of the node indices of each cable element (and any attribute information, if there is any) in turn */
    ElementData GetNextCableElementData();


    /**
     * @return the expected order of the element file (1=linear, 2=quadratic)
     */
    unsigned GetOrderOfElements()
    {
        return mOrderOfElements;
    }
    /**
     * @return the expected order of the boundary element file (1=linear, 2=quadratic)
     */
    unsigned GetOrderOfBoundaryElements()
    {
        return mOrderOfBoundaryElements;
    }

    /**
     * @return true if the boundary element file is linear, but contains information about neighbouring elements
     */
    bool GetReadContainingElementOfBoundaryElement()
    {
        return mReadContainingElementOfBoundaryElement;
    }

    /**
     * @return the vector of node attributes
     */
    std::vector<double> GetNodeAttributes();

    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global node index
     * @return a vector of the coordinates of the node
     */
    std::vector<double> GetNode(unsigned index);

    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global element index
     * @return a vector of the node indices of the element (and any attribute information, if there is any)
     */
    ElementData GetElementData(unsigned index);

    /**
     *  Normally throws an exception.  Only implemented for tetrahedral mesh reader of binary files.
     *
     * @param index  The global face index
     * @return a vector of the node indices of the face (and any attribute/containment information, if there is any)
     */
    ElementData GetFaceData(unsigned index);

    /**
     *  Normally throws an exception.  When a NCL file is available, returns a list of the elements
     *  that contain the node (only available for binary files).
     *
     * @param index  The global node index
     * @return a vector of the node indices of the face (and any attribute/containment information, if there is any)
     */
    std::vector<unsigned> GetContainingElementIndices(unsigned index);


    /*** @return true if reading binary files, false if reading ascii files */
    bool IsFileFormatBinary();

    /**
     * @return true if there is a node connectivity list (NCL) file available.
     *
     * @return whether there is a node connectivity list (NCL) file available
     */
    bool HasNclFile();

    /**
     * Sets size of std:ifstream internal read buffer. Use it for tuning I/O.
     *
     * @param bufferSize The size of the read buffer in bytes.
     */
    void SetReadBufferSize(unsigned bufferSize);

    /**
     * Sets a node permutation to use when reading in node file.
     * \todo #2452 We need a way of propagating this back to the mesh
     * (and insuring that the mesh is NOT being re-partitioned)
     *
     * @param rPermutationVector Permutation vector
     */
    void SetNodePermutation(std::vector<unsigned>& rPermutationVector);

    /**
     * @return true if a node permutation has been applied by SetNodePermutation.
     */
    bool HasNodePermutation();

    /**
     * @return the node permutation if a node permutation has been applied to this reader (or an empty permutation)
     */
    const std::vector<unsigned>& rGetNodePermutation();


private:

    /** Open mesh files. */
    void OpenFiles();

    /** Open node file. \todo Change name to OpenNodesFile for consistency with OpenElementsFile and OpenFacesFile? (#991) */
    void OpenNodeFile();

    /** Open elements file. */
    void OpenElementsFile();

    /** Open faces file. */
    void OpenFacesFile();

    /** Open node connectivity list file. */
    void OpenNclFile();

    /** Open the cable elements definition file, if it exists. */
    void OpenCableElementsFile();

    /** Read the header from each mesh file. */
    void ReadHeaders();

    /** Close mesh files. */
    void CloseFiles();

    /**
     * Read in the next line.
     *
     * @param rFileStream  The file to read from
     * @param rRawLine  Will be filled in with the next line
     */
    void GetNextLineFromStream(std::ifstream& rFileStream, std::string& rRawLine);

    /**
     * @return the Item details from the next line via a call to GetNextLineFromStream()
     *
     * @param rFileStream  The file to read from
     * @param expectedItemNumber  To check file syntax, what item is expected to be on the next line.
     * @param rDataPacket  Assumed to be of the right size but is allowed to contain dirty data on entry.
     * @param rNumAttributes  The number of attributes per item that we expect to read. Either #mNumFaceAttributes or #mNumElementAttributes.
     * @param rAttributes  Will be filled with the attribute values if rNumAttributes > 0, otherwise empty.  Note that floating point attributes are now standard
     */
    template<class T_DATA>
    void GetNextItemFromStream(std::ifstream& rFileStream, unsigned expectedItemNumber,
                               std::vector<T_DATA>& rDataPacket, const unsigned& rNumAttributes,
                               std::vector<double>& rAttributes);

    /** @return #mFilesBaseName. */
    std::string GetMeshFileBaseName();

    /** Get method specialized to 1D meshes */
    void GetOneDimBoundary();

    /**
     * Helper method to ensure we are indexing the nodes from 0
     * (some files have them indexed from 1)
     * decides according to the property mIndexFromZero
     *
     * @param rNodeIndices  The nodes we have read in.
     */
    void EnsureIndexingFromZero(std::vector<unsigned>& rNodeIndices);
};

#endif //_TRIANGLESMESHREADER_HPP_
