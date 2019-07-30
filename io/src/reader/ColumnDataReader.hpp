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

#ifndef _COLUMNDATAREADER_HPP_
#define _COLUMNDATAREADER_HPP_

#include "AbstractDataReader.hpp"

#include <string>
#include <vector>
#include <map>

#include "FileFinder.hpp"

/**
 * A concrete column data reader class.
 */
class ColumnDataReader : public AbstractDataReader
{
private:

    std::map<std::string, int> mVariablesToColumns;       /**< Map between variable names and data column numbers. \todo Change int to unsigned? (#991) */
    std::map<std::string, std::string> mVariablesToUnits; /**< Map between variable names and variable units. */
    int mNumFixedDimensions;                              /**< The number of fixed dimensions in data file. */
    bool mHasUnlimitedDimension;                          /**< Whether the data file has an unlimited dimension. */
    int mNumVariables;                                    /**< The number of variables in the data file. */
    std::string mInfoFilename;                            /**< The name of the info file.*/
    std::string mDataFilename;                            /**< The name of the data file.*/
    std::string mAncillaryFilename;                       /**< The name of the ancillary file.*/
    std::vector<double> mValues;                          /**< Vector to hold values for a variable.*/
    unsigned mFieldWidth;                                 /**< Width of each column in the text file (excludes column headers). Determined from the first data entry*/
    /**
     * Push back an entry from the data file into #mValues.
     *
     * @param rLine the line of the data file
     * @param col  the column number
     */
    void PushColumnEntryFromLine(const std::string& rLine, int col);

    /**
     * Read in a given column from a data file into #mValues.
     *
     * @param rFilename the file name
     * @param col  the column number
     */
    void ReadColumnFromFile(const std::string& rFilename, int col);

    /**
     * Push back an entry from a file into #mValues.
     *
     * @param rFilename the file name
     * @param col  the column number
     * @param row  the row number
     */
    void ReadValueFromFile(const std::string& rFilename, int col, int row);

    /**
     * Set up internal data structures based on file structure, checking that they
     * contain data in roughly the expected format.
     *
     * @param rDirectory  Absolute path of the directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     */
    void CheckFiles(const std::string& rDirectory, const std::string& rBaseName);

public:

    /**
     * Read data from the given files into memory.  The files should be formatted as if
     * written by ColumnDataWriter, with fixed width columns (except for the header line)
     * and fields in scientific notation.
     *
     * This will attempt to determine the field width from the input file.  However, it
     * needs some data to work with in order to do so.  Provided at least one correctly
     * formatted entry exists, it will be able to determine the field width, assuming
     * that space is allowed for 3 digits in the exponent.  Release 1.1 and earlier of
     * Chaste only allowed 2 digits for the exponent; we can cope with this provided that
     * the first data entry in the file has another entry immediately to the right of it.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     * @param makeAbsolute  Whether to convert directory to an absolute path using the
     *                      OutputFileHandler (defaults to true)
     */
    ColumnDataReader(const std::string& rDirectory,
                     const std::string& rBaseName,
                     bool makeAbsolute=true);

    /**
     * Alternative constructor using FileFinder to specify the directory files are stored
     * in.
     *
     * @param rDirectory  The directory the files are stored in
     * @param rBaseName  The base name of the files to read (i.e. without the extensions)
     */
    ColumnDataReader(const FileFinder& rDirectory,
                     const std::string& rBaseName);

    /**
     * @return the entries for a given variable.
     *
     * @param rVariableName
     */
    std::vector<double> GetValues(const std::string& rVariableName);

    /**
     * @return the entries for a given variable with fixed dimension.
     *
     * @param rVariableName
     * @param fixedDimension
     */
    std::vector<double> GetValues(const std::string& rVariableName, int fixedDimension);

    /**
     * @return the entries for a given variable with unlimited dimension.
     */
    std::vector<double> GetUnlimitedDimensionValues();

    /**
     * @return true if the data file has entries for a given variable.
     *
     * @param rVariableName
     */
    bool HasValues(const std::string& rVariableName);

    /**
     *  @return the field width (the number of characters (excl. preceding '+' or '-') printed for each data entry in the file).
     */
    unsigned GetFieldWidth();
};

#endif //_COLUMNDATAREADER_HPP_
