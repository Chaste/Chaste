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
#include <cassert>
#include <sstream>
#include <iostream>

#include "GmshMeshReader.hpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GmshMeshReader(std::string pathBaseName) : mFileName(pathBaseName)
{
	// Open mesh file
	mFile.open(mFileName.c_str());
	if (!mFile.is_open())
	{
		EXCEPTION("Could not open data file: " + mFileName);
	}

	ReadHeader();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::~GmshMeshReader()
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::ReadHeader()
{
    /*
     * Read mesh format information from the file header
     */
    std::string actual_line;

    mFile >> actual_line;
    assert(actual_line == "$MeshFormat");

    //Read the version no.
    mFile >> mVersionNumber >> mFileType >> mDataSize;

    if(mVersionNumber != 2.2)
    {
    	EXCEPTION("We can only read Gmsh version 2.2 files.");
    }
    assert(mFileType == 0);

    //Check mesh format close string
	mFile >> actual_line;
	assert(actual_line == "$EndMeshFormat");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
	NEVER_REACHED;
    //return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
	NEVER_REACHED;
    //return mNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
	NEVER_REACHED;
    //return mNumFaces;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements() const
{
	NEVER_REACHED;
    //return mNumCableElements;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumElementAttributes() const
{
	NEVER_REACHED;
    //return mNumElementAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumFaceAttributes() const
{
	NEVER_REACHED;
    //return mNumFaceAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNumCableElementAttributes() const
{
	NEVER_REACHED;
    //return mNumCableElementAttributes;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::Reset()
{
	NEVER_REACHED;
    //
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
	NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNodeAttributes()
{
	NEVER_REACHED;
    //return mNodeAttributes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextElementData()
{
	NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextCableElementData()
{
	NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNextFaceData()
{
	NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetNode(unsigned index)
{
	NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetElementData(unsigned index)
{
	NEVER_REACHED;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData GmshMeshReader<ELEMENT_DIM, SPACE_DIM>::GetFaceData(unsigned index)
{
	NEVER_REACHED;
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class GmshMeshReader<0,1>;
template class GmshMeshReader<1,1>;
template class GmshMeshReader<1,2>;
template class GmshMeshReader<1,3>;
template class GmshMeshReader<2,2>;
template class GmshMeshReader<2,3>;
template class GmshMeshReader<3,3>;
