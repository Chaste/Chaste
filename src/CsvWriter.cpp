/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "CsvWriter.hpp"
#include "FileFinder.hpp"

CsvWriter::CsvWriter()
    : mDirectoryName(""),
      mFileName(""),
      mDataLength(0),
      mDataLengthSet(false),
      mHeader(false)
{
}

CsvWriter::~CsvWriter()
{
    mVecUnsigned.clear();
    mVecDoubles.clear();
    mVecStrings.clear();
}

void CsvWriter::AddHeaders(std::vector<std::string>& rHeaders)
{
    if (rHeaders.size() < 1)
    {
        EXCEPTION("No strings provided in vector of headers");
    }

    mHeaderStrings = rHeaders;
    mHeader = true;
}

void CsvWriter::AddData(std::vector<unsigned>& rData)
{
    this->ValidateNewData(rData.size());
    mVecUnsigned.push_back(rData);
}

void CsvWriter::AddData(std::vector<double>& rData)
{
    this->ValidateNewData(rData.size());
    mVecDoubles.push_back(rData);
}

void CsvWriter::AddData(std::vector<std::string>& rData)
{
    this->ValidateNewData(rData.size());
    mVecStrings.push_back(rData);
}

void CsvWriter::WriteDataToFile()
{
    if (mDataLength < 1)
    {
        EXCEPTION("There is no data to write");
    }

    if (mDirectoryName == "")
    {
        EXCEPTION("Output directory has not been specified");
    }

    if (mFileName == "")
    {
        EXCEPTION("File name has not been specified");
    }

    if (mHeader && mHeaderStrings.size() != mVecUnsigned.size() + mVecDoubles.size() + mVecStrings.size())
    {
        EXCEPTION("Expecting to write a header row, but header lenght does not match number of data rows");
    }

    // Open file for writing
    std::ofstream output;
    output.open(mDirectoryName.c_str());
    assert(output.is_open());

    // Output header, if present
    if (mHeader)
    {
        unsigned headers_printed = 0;
        std::string sep = ",";

        for (unsigned header_idx = 0; header_idx < mHeaderStrings.size(); header_idx++)
        {
            headers_printed++;
            if (headers_printed == mHeaderStrings.size())
            {
                sep.clear();
            }
            output << mHeaderStrings[header_idx] << sep.c_str();
        }

        output << "\n";
    }

    unsigned num_rows = (unsigned int) (mVecUnsigned.size() + mVecDoubles.size() + mVecStrings.size());

    // Output data: unsigned then doubles then strings
    for (unsigned data_idx = 0; data_idx < mDataLength; data_idx++)
    {
        unsigned rows_printed = 0;
        std::string sep = ",";

        for (unsigned vec_idx = 0; vec_idx < mVecUnsigned.size(); vec_idx++)
        {
            rows_printed++;
            if (rows_printed == num_rows)
            {
                sep.clear();
            }
            output << mVecUnsigned[vec_idx][data_idx] << sep.c_str();
        }

        for (unsigned vec_idx = 0; vec_idx < mVecDoubles.size(); vec_idx++)
        {
            rows_printed++;
            if (rows_printed == num_rows)
            {
                sep.clear();
            }
            output << std::scientific << mVecDoubles[vec_idx][data_idx] << sep.c_str();
        }

        for (unsigned vec_idx = 0; vec_idx < mVecStrings.size(); vec_idx++)
        {
            rows_printed++;
            if (rows_printed == num_rows)
            {
                sep.clear();
            }
            output << mVecStrings[vec_idx][data_idx] << sep.c_str();
        }

        output << "\n";
    }

    // Close file
    output.close();
}

void CsvWriter::ValidateNewData(unsigned dataLength)
{
    // The first vector to be added determines the data length
    if (!mDataLengthSet)
    {
        if (dataLength == 0)
        {
            EXCEPTION("Vectors passed to this writer must contain at least one element");
        }

        mDataLength = dataLength;
        mDataLengthSet = true;
    }

    // All subsequent vectors added must be the same length as the first
    if (dataLength != mDataLength)
    {
        EXCEPTION("All data vectors added must be the same size");
    }
}

const std::string& CsvWriter::GetDirectoryName() const
{
    return mDirectoryName;
}

const std::string& CsvWriter::GetFileName() const
{
    return mFileName;
}

void CsvWriter::SetDirectoryName(std::string directoryName)
{
    FileFinder file_finder(directoryName, RelativeTo::ChasteTestOutput);

    mDirectoryName = file_finder.GetAbsolutePath();

    ///\todo Why is this code commented out?
//    // Append "/" if necessary as later this will be concatenated with file name
//    if (*(directoryName.rbegin()) != char('/'))
//    {
//        directoryName += "/";
//    }
//
//    // Validate directory name by checking it exists
//    struct stat sb;
//    if (stat(directoryName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
//    {
//        mDirectoryName = directoryName;
//    }
//    else
//    {
//        EXCEPTION("Invalid directory provided");
//    }
}

void CsvWriter::SetFileName(std::string fileName)
{
    if (fileName.find_last_of(".") == std::string::npos)
    {
        fileName += ".csv";
    }

    if (fileName.substr(fileName.find_last_of(".")) != ".csv")
    {
        EXCEPTION("This class is designed to write a file with extension .csv");
    }

    mFileName = fileName;
}
