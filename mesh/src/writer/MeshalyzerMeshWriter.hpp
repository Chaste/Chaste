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


#ifndef _MESHALYZERMESHWRITER_HPP_
#define _MESHALYZERMESHWRITER_HPP_

#include "AbstractTetrahedralMeshWriter.hpp"

/**
 * A concrete Meshalyzer mesh writer class.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MeshalyzerMeshWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * Open the file node information is written to.
     *
     * @return file handler
     * @param append  whether to append to the file, or overwrite it
     */
    out_stream OpenNodeFile(bool append=false);

    /**
     * Open the file element information is written to.
     *
     * @return file handler
     * @param append  whether to append to the file, or overwrite it
     */
    out_stream OpenElementFile(bool append=false);

    /**
     * Open the file face information is written to.
     * @return file handler
     * @param append  whether to append to the file, or overwrite it
     */
    out_stream OpenFaceFile(bool append=false);

    /**
     * Write the meta information file.
     */
    void WriteMetaFile();

    /**
     * Append footers to output files.
     */
    void WriteFilesFooter();

    /**
     * @return the mode to use when opening files.
     *
     * @param append  whether to append to the file, or overwrite it
     */
    std::ios_base::openmode GetOpenMode(bool append);

protected:
    /**
     * Create output files and add headers.
     */
    void CreateFilesWithHeaders();

    /**
     * Append local mesh data to output files.
     */
    void AppendLocalDataToFiles();

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param rCleanDirectory  whether to clean the directory (defaults to true)
     * @param rSetCoolGraphics (defaults to false)
     */
    MeshalyzerMeshWriter(const std::string& rDirectory,
                         const std::string& rBaseName,
                         const bool& rCleanDirectory=true,
                         const bool& rSetCoolGraphics=false);

    /**
     * Write mesh data to files.
     */
    void WriteFiles();

    /**
     * Destructor.
     */
    virtual ~MeshalyzerMeshWriter();
};

#endif //_MESHALYZERMESHWRITER_HPP_
