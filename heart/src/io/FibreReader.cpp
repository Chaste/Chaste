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

#include "FibreReader.hpp"

#include <sstream>
#include "Exception.hpp"

template<unsigned DIM>
FibreReader<DIM>::FibreReader(const FileFinder& rFileFinder, FibreFileType fibreFileType)
   : mFileIsBinary(false), // overwritten by ReadNumLinesOfDataFromFile() if applicable.
     mNextIndex(0u)
{
    if (fibreFileType == AXISYM)
    {
        mNumItemsPerLine = DIM;
    }
    else //(fibreFileType == ORTHO)
    {
        mNumItemsPerLine = DIM*DIM;
    }
    mTokens.resize(mNumItemsPerLine);

    mFilePath = rFileFinder.GetAbsolutePath();
    mDataFile.open(mFilePath.c_str());
    if (!mDataFile.is_open())
    {
        EXCEPTION("Failed to open fibre file " + rFileFinder.GetAbsolutePath());
    }

    // Note: this method will close the file on error
    ReadNumLinesOfDataFromFile();
}

template<unsigned DIM>
FibreReader<DIM>::~FibreReader()
{
    mDataFile.close();
}

template<unsigned DIM>
void FibreReader<DIM>::GetAllAxi(std::vector< c_vector<double, DIM> >& direction)
{
    assert(direction.empty());
    if (mNumItemsPerLine != DIM)
    {
        EXCEPTION("Use GetAllOrtho when reading orthotropic fibres");
    }
    direction.reserve(mNumLinesOfData);
    for (unsigned i=0; i<mNumLinesOfData; i++)
    {
        c_vector<double, DIM> temp_vector;
        GetFibreVector(i, temp_vector, false);
        direction.push_back(temp_vector);
    }
}

template<unsigned DIM>
void FibreReader<DIM>::GetAllOrtho(std::vector< c_vector<double, DIM> >& first_direction,
                                   std::vector< c_vector<double, DIM> >& second_direction,
                                   std::vector< c_vector<double, DIM> >& third_direction)
{
    assert(first_direction.empty());
    assert(second_direction.empty());
    assert(third_direction.empty());
    if (mNumItemsPerLine != DIM*DIM)
    {
        EXCEPTION("Use GetAllAxi when reading axisymmetric fibres");
    }
    for (unsigned i=0; i<mNumLinesOfData; i++)
    {
        c_matrix<double, DIM, DIM> temp_matrix;
        GetFibreSheetAndNormalMatrix(i, temp_matrix, true);

        //Note that although the matrix appears row-wise in the ascii .ortho file,
        //for convenience it is stored column-wise.
        matrix_column<c_matrix<double, DIM, DIM> > col0(temp_matrix, 0);
        first_direction.push_back(col0);
        if (DIM>=2)
        {
            matrix_column<c_matrix<double, DIM, DIM> > col1(temp_matrix, 1);
            second_direction.push_back(col1);
        }
        if (DIM==3)
        {
            matrix_column<c_matrix<double, DIM, DIM> > col2(temp_matrix, 2);
            third_direction.push_back(col2);
        }
    }
}

template<unsigned DIM>
void FibreReader<DIM>::GetFibreSheetAndNormalMatrix(unsigned fibreIndex,
                                                    c_matrix<double,DIM,DIM>& rFibreMatrix,
                                                    bool checkOrthogonality)
{
    if (mNumItemsPerLine != DIM*DIM)
    {
        EXCEPTION("Use GetFibreVector when reading axisymmetric fibres");
    }
    if (fibreIndex < mNextIndex)
    {
        EXCEPTION("Fibre reads must be monotonically increasing; " << fibreIndex
                << " is before expected next index " << mNextIndex);
    }
    if (mFileIsBinary)
    {

        // Skip to the desired index
        mDataFile.seekg((fibreIndex-mNextIndex)*mNumItemsPerLine*sizeof(double), std::ios::cur);
        // Take mNumItemsPerLine from the ifstream
        mDataFile.read((char*)&(rFibreMatrix(0,0)), mNumItemsPerLine*sizeof(double));
        mNextIndex = fibreIndex+1;
    }
    else
    {
        unsigned num_entries = 0u;
        while (fibreIndex >= mNextIndex)
        {
             num_entries = GetTokensAtNextLine();
             mNextIndex++;
        }
        if (num_entries < mNumItemsPerLine)
        {
            EXCEPTION("A line is incomplete in " << mFilePath
                          << " - each line should contain " << DIM*DIM << " entries");
        }
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                rFibreMatrix(i,j) = mTokens[DIM*i + j];
            }
        }
    }

    // The binary file and ascii file are row-major. However, we store column major matrices.
    rFibreMatrix = trans(rFibreMatrix);

    if (checkOrthogonality)
    {
        // Note that we define this matrix before setting it as otherwise the profiling build will break (see #2367)
        c_matrix<double,DIM,DIM> temp;
        temp = prod(trans(rFibreMatrix), rFibreMatrix);

        // Check temp is equal to the identity
        for (unsigned i=0; i<DIM; i++)
        {
            for (unsigned j=0; j<DIM; j++)
            {
                double val = (i==j ? 1.0 : 0.0);

                if (fabs(temp(i,j)-val) > 1e-4)
                {
                    EXCEPTION("Read fibre-sheet matrix, " << rFibreMatrix << " from file "
                                  << " which is not orthogonal (tolerance 1e-4)");
                }
            }
        }
    }
}

template<unsigned DIM>
void FibreReader<DIM>::GetFibreVector(unsigned fibreIndex,
                                      c_vector<double,DIM>& rFibreVector,
                                      bool checkNormalised)
{
    if (mNumItemsPerLine != DIM)
    {
        EXCEPTION("Use GetFibreSheetAndNormalMatrix when reading orthotropic fibres");
    }
    if (fibreIndex < mNextIndex)
    {
        EXCEPTION("Fibre reads must be monotonically increasing; " << fibreIndex
                  << " is before expected next index " << mNextIndex);
    }

    if (mFileIsBinary)
    {
        // Skip to the desired index
        mDataFile.seekg((fibreIndex-mNextIndex)*mNumItemsPerLine*sizeof(double), std::ios::cur);
        // Take mNumItemsPerLine from the ifstream
        mDataFile.read((char*)&rFibreVector[0], mNumItemsPerLine*sizeof(double));
        mNextIndex = fibreIndex+1;
    }
    else
    {
        unsigned num_entries = 0u;
        while (fibreIndex >= mNextIndex)
        {
            num_entries = GetTokensAtNextLine();
            mNextIndex++;
        }
        if (num_entries < mNumItemsPerLine)
        {
            EXCEPTION("A line is incomplete in " << mFilePath
                          << " - each line should contain " << DIM << " entries");
        }
        for (unsigned i=0; i<DIM; i++)
        {
            rFibreVector(i) = mTokens[i];
        }
    }


    if (checkNormalised && fabs(norm_2(rFibreVector)-1)>1e-4)
    {
        EXCEPTION("Read vector " << rFibreVector << " from file "
                      << mFilePath << " which is not normalised (tolerance 1e-4)");
    }
}

template<unsigned DIM>
unsigned FibreReader<DIM>::GetTokensAtNextLine()
{
    assert(mTokens.size() == mNumItemsPerLine);

    std::string line;
    bool blank_line;

    do
    {
        getline(mDataFile, line);

        if (line.empty() && mDataFile.eof())
        {
            mDataFile.close();
            std::string error =   "End of file " + mFilePath + " reached. Either file contains fewer "
                                + "definitions than defined in header, or one of the GetNext[..] methods "
                                + "has been called too often";
            EXCEPTION(error);
        }

        // Get rid of any comments
        line = line.substr(0, line.find('#'));

        blank_line = (line.find_first_not_of(" \t",0) == std::string::npos);
    }
    while (blank_line);

    // Get rid of any trailing whitespace
    std::string::iterator iter = line.end();
    iter--;
    unsigned nchars2delete = 0;
    while (*iter == ' ' || *iter == '\t')
    {
        nchars2delete++;
        iter--;
    }
    line.erase(line.length()-nchars2delete);

    std::stringstream line_stream(line);

    unsigned index = 0;
    while (!line_stream.eof())
    {
        double item;
        line_stream >> item;
        if (index >= mNumItemsPerLine)
        {
            EXCEPTION("Too many entries in a line in " + mFilePath);
        }
        mTokens[index++] = item;
    }

    return index; // the number of entries put into mTokens
}

template<unsigned DIM>
void FibreReader<DIM>::ReadNumLinesOfDataFromFile()
{
    std::string raw_line;
    bool blank_line = false;
    do
    {
        getline(mDataFile, raw_line);
        //Strip comments following a hash
        raw_line = raw_line.substr(0, raw_line.find('#'));
        //Check for blank line
        blank_line = (raw_line.find_first_not_of(" \t",0) == std::string::npos);
    }
    while (blank_line);

    std::stringstream header_line(raw_line);

    header_line >> mNumLinesOfData;

    std::string extras;
    header_line >> extras;

    if (extras == "BIN")
    {
        mFileIsBinary = true;
    }
    else if (extras!="")
    {
        mDataFile.close();
        EXCEPTION("First (non comment) line of the fibre orientation file should contain the number of lines of data in the file (and possibly a BIN tag) at most");
    }
}

template class FibreReader<1>;
template class FibreReader<2>;
template class FibreReader<3>;
