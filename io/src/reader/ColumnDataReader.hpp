/*

Copyright (C) University of Oxford, 2005-2012

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
    static const int SPACING = 2;                         /**< Space between columns (includes minus sign) */

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
     * Get the entries for a given variable.
     *
     * @param rVariableName
     */
    std::vector<double> GetValues(const std::string& rVariableName);

    /**
     * Get the entries for a given variable with fixed dimension.
     *
     * @param rVariableName
     * @param fixedDimension
     */
    std::vector<double> GetValues(const std::string& rVariableName, int fixedDimension);

    /**
     * Get the entries for a given variable with unlimited dimension.
     */
    std::vector<double> GetUnlimitedDimensionValues();

    /**
     * Determine whether the data file has entries for a given variable.
     *
     * @param rVariableName
     */
    bool HasValues(const std::string& rVariableName);

    /**
     *  Get the field width (the number of characters (excl. preceding '+' or '-') printed for each data entry in the file).
     */
    unsigned GetFieldWidth();
};

#endif //_COLUMNDATAREADER_HPP_
