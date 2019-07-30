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

#ifndef XDMFMESHWRITER_HPP_
#define XDMFMESHWRITER_HPP_

#include "AbstractTetrahedralMeshWriter.hpp"
// Xerces is currently not supported in the Windows port
#ifndef _MSC_VER

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xsd/cxx/xml/string.hxx>

#ifndef X //Also used in XmlTools in the heart component
/**
 * Convenience macro for transcoding C++ strings to Xerces' format.
 * @param str  the string to transcode
 */
#define X(str) xsd::cxx::xml::string(str).c_str()
#endif //X

#endif // _MSC_VER
/**
 * A class for writing from a Chaste mesh to the geometry/topology components of
 * an XDMF file.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class XdmfMeshWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
protected:
    unsigned mNumberOfTimePoints; /**< Defaults to 1, when we are writing geometry only.  Used in HDF5 converter which has "protected" access as a derived class.*/
    double mTimeStep; /**< Defaults to 1.0.*/

private:
    /**
     * Write the master file.  This just contains references to the geometry/topology files.
     * @param numberOfChunks  is the number of geometric pieces which is 1 for sequential code and for non-distributed meshes.
     */
    void WriteXdmfMasterFile(unsigned numberOfChunks=1u);

#ifndef _MSC_VER
    /**
     * Generate Attribute tags and append to the element.  Here this is a dummy class, but can be
     * overloaded with real variables elsewhere (see pde/src/postprocesssing/Hdf5toXdmfConverter).
     * @param pGridElement  Pointer to DOMElement to append Attribute tags to.
     * @param pDomDocument  Pointer to DOMDocument to generate new elements.
     * @param timeStep  Index of time point to write.
     */
    virtual void AddDataOnNodes(XERCES_CPP_NAMESPACE_QUALIFIER DOMElement* pGridElement,
                                XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument* pDomDocument,
                                unsigned timeStep)
    {
        //Empty body - implemented in derived classes
    }
#endif // _MSC_VER

public:
    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    XdmfMeshWriter(const std::string& rDirectory,
                   const std::string& rBaseName,
                   const bool clearOutputDir=true);
    /**
     * Write the files using a mesh reader.  Called from WriteFilesUsingMeshReader in the base class.
     */
    void WriteFiles();
    /**
     * Write the files using mesh -- this allows the parallel distributed case to write only local data.
     * @param rMesh the mesh
     * @param keepOriginalElementIndexing  Whether to write the mesh with the same element ordering.
     *                                     Ignored in this derived class
     */
    void WriteFilesUsingMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                             bool keepOriginalElementIndexing=true);
};

#endif /* XDMFMESHWRITER_HPP_ */
