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


#ifndef _TRIANGLESMESHWRITER_HPP_
#define _TRIANGLESMESHWRITER_HPP_

#include "AbstractTetrahedralMeshWriter.hpp"
#include "OutputFileHandler.hpp"


/**
 * A concrete mesh writer class that writes Triangle output files.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TrianglesMeshWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    TrianglesMeshWriter(const std::string& rDirectory,
                        const std::string& rBaseName,
                        const bool clearOutputDir=true);
    /**
     * Switch this mesh write to write binary files
     *
     * (set to write ascii files in the constructor)
     */
     void SetWriteFilesAsBinary();

    /**
     * Write mesh data to files.
     */
    void WriteFiles();

    /**
     * Write elements as faces (used in the case ELEMENT_DIM== SPACE_DIM-1)
     */
    void WriteElementsAsFaces();

    /**
     * Write faces as edges (used in the case ELEMENT_DIM==2, SPACE_DIM==3)
     */
    void WriteFacesAsEdges();

    /**
     * Write a line (ascii format) to a specific file stream
     * Templated over std::vector dataPacket contents of unsigned or doubles.
     * Templated over type of attribute.
     *
     * @param pFile Pointer to file stream
     * @param itemNumber Index of the element, node or face
     * @param dataPacket List of unsigneds (for node indices) or doubles (for node locations)
     * @param rAttributes A vector of attributes (double precision).
     */

     template<class  T_DATA>
     void WriteItem(out_stream &pFile, unsigned itemNumber, const std::vector<T_DATA> &dataPacket, const std::vector<double> &rAttributes);

     /**
      * Write a line (ascii format) to a specific file stream
      * Templated over std::vector dataPacket contents of unsigned or doubles.
      *
      * @param pFile Pointer to file stream
      * @param itemNumber Index of the element, node or face
      * @param dataPacket List of unsigneds (for node indices) or doubles (for node locations)
      */
      template<class  T_DATA>
      void WriteItem(out_stream &pFile, unsigned itemNumber, const std::vector<T_DATA> &dataPacket);


    /**
     * Destructor.
     */
    virtual ~TrianglesMeshWriter();
};

#endif //_TRIANGLESMESHWRITER_HPP_
