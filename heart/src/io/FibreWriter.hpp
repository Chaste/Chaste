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

#ifndef FIBREWRITER_HPP_
#define FIBREWRITER_HPP_

#include <string>
#include <fstream>
#include <vector>

#include "OutputFileHandler.hpp"
#include "UblasIncludes.hpp"


/**
 * A class for writing .axi files (files which define the fibre direction
 * for each element) and .ortho files (files which define the fibre, sheet
 * and normal directions for each element.
 */
template<unsigned DIM>
class FibreWriter
{
private:
    OutputFileHandler* mpOutputFileHandler; /**< Output file handler */

    std::string mBaseName; /**< Base name for the input files */
    bool mFileIsBinary;  /**< Whether all data is to be written as binary*/

    /**
     * Open a fibre file for writing and write the header line.
     *
     * @param rFileName The path and name of the file to create
     * @param numItems The number of items (~ number of elements)
     * @return an outstream pointer to the open file.
     */
    out_stream OpenFileAndWriteHeader(const std::string& rFileName, unsigned numItems);

public:
    /**
     * Create a new FibreWriter.
     *
     * @param rDirectory  the directory in which to write the fibre to file
     * @param rBaseName  the base name of the files in which to write the fibre data
     * @param clearOutputDir  whether to clean the directory (defaults to true)
     */
    FibreWriter(const std::string& rDirectory,
                const std::string& rBaseName,
                const bool clearOutputDir=true);

    /**
     *  Destructor
     */
    ~FibreWriter();

     /**
      * Writes all axisymmetric vectors to the file.
      * @param fibres a vector of fibre direction vectors
      */
     void WriteAllAxi(const std::vector< c_vector<double, DIM> >& fibres);

     /**
      * Writes all orthonormal tensors to the file.
      * @param fibres a vector of longitudial fibre direction vectors
      * @param second a vector of transverse vectors (orthogonal to fibres)
      * @param third a vector of normal vectors (orthogonal to others)
      */
     void WriteAllOrtho(const std::vector< c_vector<double, DIM> >& fibres,
                        const std::vector< c_vector<double, DIM> >& second,
                        const std::vector< c_vector<double, DIM> >& third);
    /**
     * Switch to write binary fibre file
     *
     * (set to write ascii files in the constructor)
     */
     void SetWriteFileAsBinary();
};

#endif /*FIBREWRITER_HPP_*/
