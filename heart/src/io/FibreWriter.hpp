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
