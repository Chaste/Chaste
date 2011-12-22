/*

Copyright (C) University of Oxford, 2005-2011

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

#include <cassert>
#include "AbstractMeshWriter.hpp"
#include "PetscTools.hpp"

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::AbstractMeshWriter(const std::string &rDirectory,
                                                               const std::string &rBaseName,
                                                               const bool clearOutputDir)
    : mBaseName(rBaseName),
      mpMeshReader(NULL)
{
    mpOutputFileHandler = new OutputFileHandler(rDirectory, clearOutputDir);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::~AbstractMeshWriter()
{
    delete mpOutputFileHandler;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputDirectory()
{
    return mpOutputFileHandler->GetOutputDirectoryFullPath();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNumNodes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumBoundaryFaces()
{
    return mNumBoundaryElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNumCableElements()
{
    return mNumCableElements;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    assert(mpMeshReader!=NULL);
    return mpMeshReader->GetNextNode();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement()
{
    assert(mpMeshReader!=NULL);
    return mpMeshReader->GetNextElementData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextBoundaryElement()
{
    assert(mpMeshReader!=NULL);
    return mpMeshReader->GetNextFaceData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextCableElement()
{
    assert(mpMeshReader!=NULL);
    return mpMeshReader->GetNextCableElementData();
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(
        AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    mpMeshReader = &rMeshReader;
    mNumNodes = mpMeshReader->GetNumNodes();
    mNumElements = mpMeshReader->GetNumElements();
    mNumBoundaryElements = mpMeshReader->GetNumFaces();

    ///Only triangles mesh readers know about cable elements
    mNumCableElements = mpMeshReader->GetNumCableElements();

    if (!PetscTools::AmMaster())
    {
        return;
    }

    WriteFiles();
}

/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class AbstractMeshWriter<1,1>;
template class AbstractMeshWriter<1,2>;
template class AbstractMeshWriter<1,3>;
template class AbstractMeshWriter<2,2>;
template class AbstractMeshWriter<2,3>;
template class AbstractMeshWriter<3,3>;
