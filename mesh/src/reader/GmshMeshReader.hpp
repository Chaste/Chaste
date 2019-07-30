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


#ifndef _GMSHMESHREADER_HPP_
#define _GMSHMESHREADER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include "AbstractMeshReader.hpp"

/**
 * Structure encapsulating the enumeration of
 */
struct GmshTypes
{
    /**
     * What things a path can be relative to.
     */
    enum Value
    {
        LINE = 1u,
        QUADRATIC_LINE = 8u,
        TRIANGLE = 2u,
        QUADRATIC_TRIANGLE = 9u,
        TETRAHEDRON = 4u,
        QUADRATIC_TETRAHEDRON = 11u
    };
};

/**
 * Class to enable reading of Gmsh format mesh files (see #2312).
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class GmshMeshReader : public AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>
{
public:

    /**
     * Constructor.
     *
     * @param pathBaseName  the base name of the files from which to read the mesh data
     *    (either absolute, or relative to the current directory)
     * @param orderOfElements  the order of each element: 1 for linear, 2 for quadratic (defaults to 1)
     * @param orderOfBoundaryElements the order of each boundary element: 1 for linear, 2 for quadratic (defaults to 1. May
     *  or may not be different to orderOfElements
     */
    GmshMeshReader(std::string pathBaseName,
                   unsigned orderOfElements = 1,
                   unsigned orderOfBoundaryElements = 1);

    /**
     * Destructor
     */
    ~GmshMeshReader();

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

    /** * @return the expected order of the element file (1=linear, 2=quadratic) */
    unsigned GetOrderOfElements();

    /** @return the expected order of the element file (1=linear, 2=quadratic) */
    unsigned GetOrderOfBoundaryElements();

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

private:
    /** Read all the header information from the mesh file. */
    void ReadHeaders();

    /** Read the node header from the mesh file. */
    void ReadNodeHeader();

    /** Read the element header from the mesh file. */
    void ReadElementHeader();

    /** Read the face header from the mesh file. */
    void ReadFaceHeader();

    /** Opens the .msh file descriptors */
    void OpenFiles();

    /** Closes the .msh file descriptors */
    void CloseFiles();

    std::string mFileName; /**< The name of the mesh file. */
    std::ifstream mNodeFile; /**< A file stream used to read the node (and header) part of the file. */
    std::ifstream mElementFile; /**< A file stream used to read the volume elements of the file. */
    std::ifstream mFaceFile; /**< A file stream used to read the boundary elements of the file. */
    double mVersionNumber; /**< The version number of the file. */
    unsigned mFileType; /**< The type of the mesh file being read (should always be 0) */
    unsigned mDataSize; /**< The number of floating point numbers in the file */
    unsigned mNumNodes; /**< Number of nodes in the mesh. */
    unsigned mNumElements; /**< Number of elements in the mesh. */
    unsigned mNumFaces; /**< Number of faces in the mesh. */
    unsigned mTotalNumElementsAndFaces; /**<Total number of elements and faces in the mesh. */
    unsigned mNumElementAttributes; /**< Is the number of attributes stored for each element. */
    unsigned mNumFaceAttributes; /**< Is the number of attributes stored for each face. */
    unsigned mOrderOfElements; /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mOrderOfBoundaryElements; /**< The order of each element (1 for linear, 2 for quadratic). */
    unsigned mNodesPerElement; /**< The number of nodes contained in each element. */
    unsigned mNodesPerBoundaryElement; /**< The number of nodes contained in each boundary element. */
};

#endif //_GMSHMESHREADER_HPP_
