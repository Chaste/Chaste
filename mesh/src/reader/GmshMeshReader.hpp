/*

Copyright (c) 2005-2013, University of Oxford.
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
     */
    GmshMeshReader(std::string pathBaseName);

    /**
     * Destructor
     */
     ~GmshMeshReader();

    /** Returns the number of elements in the mesh */
    unsigned GetNumElements() const;

    /** Returns the number of nodes in the mesh */
    unsigned GetNumNodes() const;

    /** Returns the number of faces in the mesh (synonym of GetNumEdges()) */
    unsigned GetNumFaces() const;

    /** Returns the number of cable elements in the mesh */
    unsigned GetNumCableElements() const;

    /** Returns the number of attributes in the mesh */
    unsigned GetNumElementAttributes() const;

    /** Returns the number of attributes in the mesh */
    unsigned GetNumFaceAttributes() const;

    /** Returns the number of cable element attributes in the mesh */
    unsigned GetNumCableElementAttributes() const;

    /** Resets pointers to beginning*/
    void Reset();

    /** Returns a vector of the coordinates of each node in turn */
    std::vector<double> GetNextNode();

    /** Returns a vector of the nodes of each element (and any attribute information, if there is any) in turn */
    ElementData GetNextElementData();

    /** Returns a vector of the nodes of each face in turn (synonym of GetNextEdgeData()) */
    ElementData GetNextFaceData();

    /** Returns a vector of the node indices of each cable element (and any attribute information, if there is any) in turn */
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

    /** Read the header from the mesh file. */
    void ReadHeader();

private:
    std::string mFileName;     /**< The name of the mesh file. */
    std::ifstream mFile;       /**< The file for the mesh. */
    double mVersionNumber;     /**< The version number of the file. */
    unsigned mFileType;  	   /**< The type of the mesh file being read (should always be 0) */
    unsigned mDataSize;		   /**< The number of floating point numbers in the file */

};

#endif //_GMSHMESHREADER_HPP_
