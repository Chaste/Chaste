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

#ifndef FIBREREADER_HPP_
#define FIBREREADER_HPP_

#include <string>
#include <fstream>
#include <vector>

#include "UblasIncludes.hpp"
#include "FileFinder.hpp"

/**
 * Simple enumeration for use in FibreReader constructor
 */
typedef enum FibreFileType_
{
    AXISYM=0,
    ORTHO
} FibreFileType;

/**
 * A class for reading .axi files (files which define the fibre direction
 * for each element) and .ortho files (files which define the fibre, sheet
 * and normal directions for each element.
 */
template<unsigned DIM>
class FibreReader
{
private:
    /** File stream to use for GetTokensAtNextLine */
    std::ifstream mDataFile;

    /** Absolute path of the file being read */
    std::string mFilePath;

    /** Number of lines of data in the file, read from the first line of the file */
    unsigned mNumLinesOfData;

    bool mFileIsBinary;  /**< Whether the data file has binary entries*/

    /** How many items we expect to find per line: DIM for axisymmetric, DIM*DIM for orthotropic */
    unsigned mNumItemsPerLine;

    /** The next index we expect to read from the file. */
    unsigned mNextIndex;

    /** Vector which entries read from a line in a file is put into. */
    std::vector<double> mTokens;

    /**
     *  Read a line of numbers from #mDataFile.
     *  Sets up the member variable #mTokens with the data in the next line.
     *  @return the number of data entries put into #mTokens
     */
    unsigned GetTokensAtNextLine();

    /**
     *  Read number of elements from #mDataFile.
     *  Note: Must be called before GetTokensAtNextLine (it assumes that
     *  it's reading the first line).
     */
    void ReadNumLinesOfDataFromFile();

public:
    /**
     * Create a new FibreReader.
     *
     * @param rFileFinder  the path to the fibre direction file
     * @param fibreFileType AXISYM or ORTHO depending on type of file to be read
     */
    FibreReader(const FileFinder& rFileFinder, FibreFileType fibreFileType);

    /**
     *  Destructor closes file.
     */
    ~FibreReader();


    /**
     * Read a fibre direction matrix from the file.  Must only be used when
     * reading an orthotropic file.  These have lines of the form
     * \code
     *  fibre0 fibre1 fibre2 sheet0 sheet1 sheet2 normal0 normal1 normal2
     * \endcode
     * which are converted to the matrix
     * \code
     *     [ fibre0   sheet0   normal0  ]
     *     [ fibre1   sheet1   normal1  ]
     *     [ fibre2   sheet2   normal2  ]
     * \endcode
     *
     * @param fibreIndex  which fibre vector to read.  Note that vectors must be read
     *     in monotonically increasing order, so subsequent calls to this method must
     *     always pass a strictly greater index.  They may skip vectors, however.
     * @param rFibreMatrix  matrix to be filled in
     * @param checkOrthogonality  if true, checks if the matrix is orthogonal
     *    and throws an exception if not
     */
    void GetFibreSheetAndNormalMatrix(unsigned fibreIndex, c_matrix<double,DIM,DIM>& rFibreMatrix, bool checkOrthogonality=true);

    /**
     * Read a fibre direction vector from the file.  Must only be used when
     * reading an axisymmetric file.  These have lines of the form
     * \code
     *  fibre0 fibre1 fibre2
     * \endcode
     *
     * @param fibreIndex  which fibre vector to read.  Note that vectors must be read
     *     in monotonically increasing order, so subsequent calls to this method must
     *     always pass a strictly greater index.  They may skip vectors, however.
     * @param rFibreVector  vector to be filled in
     * @param checkNormalised  if true, checks if the read vector is normalised
     *   and throws an exception if not
     */
    void GetFibreVector(unsigned fibreIndex, c_vector<double,DIM>& rFibreVector, bool checkNormalised=true);

    /**
     *  @return the number of lines of data in the file - this is the value read from
     *  the first line.
     */
    unsigned GetNumLinesOfData()
    {
        return mNumLinesOfData;
    }

    /**
     * Get every line of a fibre file (axi-symmetric) in vector of vectors format.
     * This is useful for adding to a visualizer.
     * (Do not use with GetNext...)
     *
     * @param direction  an empty vector which will be filled with data from file
     */
    void GetAllAxi(std::vector< c_vector<double, DIM> >& direction);

    /**
     * Get every line of a fibre file (orthotropic) in vector of vectors format.
     * This is useful for adding to a visualizer.
     * (Do not use with GetNext...)
     *
     * @param first_direction  an empty vector which will be filled with data from file
     * @param second_direction  an empty vector which will be filled with data from file
     * @param third_direction  an empty vector which will be filled with data from file (or will remain empty in 2D)
     */
    void GetAllOrtho(std::vector< c_vector<double, DIM> >& first_direction,
                         std::vector< c_vector<double, DIM> >& second_direction,
                         std::vector< c_vector<double, DIM> >& third_direction);

    /**
     * @return Whether the fibre file contains binary data.
     */
    bool IsBinary()
    {
        return mFileIsBinary;
    }
};

#endif /*FIBREREADER_HPP_*/
