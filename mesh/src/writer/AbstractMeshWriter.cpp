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
#include "AbstractMeshWriter.hpp"
#include "PetscTools.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::AbstractMeshWriter(const std::string& rDirectory,
                                                               const std::string& rBaseName,
                                                               const bool clearOutputDir)
    : mBaseName(rBaseName),
      mpMeshReader(nullptr)
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
    assert(mpMeshReader != nullptr);
    return mpMeshReader->GetNextNode();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement()
{
    assert(mpMeshReader != nullptr);
    return mpMeshReader->GetNextElementData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextBoundaryElement()
{
    assert(mpMeshReader != nullptr);
    return mpMeshReader->GetNextFaceData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextCableElement()
{
    assert(mpMeshReader != nullptr);
    return mpMeshReader->GetNextCableElementData();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)
{
    mpMeshReader = &rMeshReader;
    mNumNodes = mpMeshReader->GetNumNodes();
    mNumElements = mpMeshReader->GetNumElements();
    mNumBoundaryElements = mpMeshReader->GetNumFaces();

    // Only triangles mesh readers know about cable elements
    mNumCableElements = mpMeshReader->GetNumCableElements();

    if (PetscTools::AmMaster())
    {
        WriteFiles();
    }
}

// Explicit instantiation
template class AbstractMeshWriter<1,1>;
template class AbstractMeshWriter<1,2>;
template class AbstractMeshWriter<1,3>;
template class AbstractMeshWriter<2,2>;
template class AbstractMeshWriter<2,3>;
template class AbstractMeshWriter<3,3>;
