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

/**
 * @file
 *
 * Implementation file for ColumnDataReader class.
 */

#include "ColumnDataReader.hpp"
#include "ColumnDataConstants.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <climits>
#include <cctype> //for isdigit
#include "OutputFileHandler.hpp"
#include "Exception.hpp"

/**
 * Variables read in from the data file are initialised to the
 * following constant so one can check if they were read correctly.
 */
const int NOT_READ = INT_UNSET;

ColumnDataReader::ColumnDataReader(const std::string& rDirectory,
                                   const std::string& rBaseName,
                                   bool makeAbsolute)
{
    // Find out where files are really stored
    std::string directory;
    if (makeAbsolute)
    {
        OutputFileHandler output_file_handler(rDirectory, false);
        directory = output_file_handler.GetOutputDirectoryFullPath();
    }
    else
    {
        // Add a trailing slash if needed
        if (!(*(rDirectory.end()-1) == '/'))
        {
            directory = rDirectory + "/";
        }
        else
        {
            directory = rDirectory;
        }
    }
    CheckFiles(directory, rBaseName);
}

ColumnDataReader::ColumnDataReader(const FileFinder& rDirectory,
                                   const std::string& rBaseName)
{
    if (!rDirectory.IsDir() || !rDirectory.Exists())
    {
        EXCEPTION("Directory does not exist: " + rDirectory.GetAbsolutePath());
    }
    CheckFiles(rDirectory.GetAbsolutePath(), rBaseName);
}

void ColumnDataReader::CheckFiles(const std::string& rDirectory, const std::string& rBaseName)
{
    // Read in info file
    mInfoFilename = rDirectory + rBaseName + ".info";
    std::ifstream infofile(mInfoFilename.c_str(), std::ios::in);

    // If it doesn't exist - throw exception
    if (!infofile.is_open())
    {
        EXCEPTION("Couldn't open info file: " + mInfoFilename);
    }
    std::string junk;
    mNumFixedDimensions = NOT_READ;
    mHasUnlimitedDimension = false;
    mNumVariables = NOT_READ;

    infofile >> junk;
    infofile >> mNumFixedDimensions >> junk;
    infofile >> mHasUnlimitedDimension >> junk;
    infofile >> mNumVariables;

    if (mNumFixedDimensions == NOT_READ || mNumVariables == NOT_READ)
    {
        infofile.close();
        EXCEPTION("Couldn't read info file correctly");
    }

    // Read in variables and associated them with a column number
    if (mHasUnlimitedDimension)
    {
        if (mNumFixedDimensions < 1)
        {
            mDataFilename = rDirectory + rBaseName + ".dat";
        }
        else
        {
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << 0;

            mDataFilename = rDirectory + rBaseName + "_" + suffix.str() + ".dat";

            /*
             * The ancillary path needs to come from a single place that is
             * used by both the reader & writer, otherwise all will be bad.
             */
            mAncillaryFilename = rDirectory + rBaseName + "_unlimited.dat";

            // Extract the units and place into map
            std::ifstream ancillaryfile(mAncillaryFilename.c_str(), std::ios::in);

            // If it doesn't exist - throw exception
            if (!ancillaryfile.is_open())
            {
                EXCEPTION("Couldn't open ancillary data file");
            }
            std::string dimension;
            std::getline(ancillaryfile, dimension);
            std::stringstream dimension_stream(dimension);
            std::string dimension_unit, dimension_name, header;
            dimension_stream >> header;

            // Separate into variable name and units
            int unitpos = header.find("(") + 1;

            dimension_name = header.substr(0, unitpos - 1);
            dimension_unit = header.substr(unitpos, header.length() - unitpos - 1);

            mVariablesToUnits[dimension_name] = dimension_unit;
            ancillaryfile.close();
        }
    }
    else
    {
        mDataFilename = rDirectory + rBaseName + ".dat";
    }

    std::ifstream datafile(mDataFilename.c_str(), std::ios::in);
    // If it doesn't exist - throw exception
    if (!datafile.is_open())
    {
        EXCEPTION("Couldn't open data file");
    }

    std::string variables;
    std::getline(datafile, variables);
    std::stringstream variable_stream(variables);
    std::string header, variable, unit;
    int column = 0;

    // Insert variables into map
    while (variable_stream >> header)
    {
        // Separate into variable name and units
        int unitpos = header.find("(") + 1;

        variable = header.substr(0, unitpos - 1);
        unit = header.substr(unitpos, header.length() - unitpos - 1);

        mVariablesToColumns[variable] = column;
        mVariablesToUnits[variable] = unit;

        column++;
    }

    /*
     * Now read the first line of proper data to determine the field width used when this
     * file was created. Do this by
     * 1. reading the first entry and measuring the distance from
     * the decimal point to the 'e'.  This gives the precision; the field width is then
     * precision + 7 (With MSVC on Windows, it's precision + 8).
     * e.g. if the first entry is
     *   6.3124e+01         => field width = 11 // chaste release 1 and 1.1
     *  -3.5124e+01         => field width = 11 // chaste release 1 and 1.1
     *  +1.00000000e+00     => field width = 15
     *  -1.20000000e+01     => field width = 15
     *  -1.12345678e-321    => field width = 15
     * 2. Because the first column has a varying number of spaces read a few columns and
     *    do some modular arithmetic to work out the correct width
     */
    std::string first_line;
    std::string first_entry;
    unsigned last_pos=0u;
    // Read the first entry of the line. If there is no first entry, move to the next line..
    while (first_entry.length()==0 && !datafile.eof())
    {
        std::getline(datafile, first_line);
        std::stringstream stream(first_line);
        stream >> first_entry;
        last_pos = stream.tellg(); // Where the first number ends (but it might be in the column 2 or 3)
        while (stream.good() && last_pos <170) //Avoid reading more than about 10 columns, because we want to avoid last_pos being divisible by too many factors
        {
            std::string last_entry;
            stream >> last_entry;
            if (stream.tellg() > 0)
            {
                last_pos = stream.tellg();
            }
        }
    }

    if (datafile.eof() && first_entry.length()==0)
    {
        EXCEPTION("Unable to determine field width from file as cannot find any data entries");
    }
    assert (last_pos > 0u);

    size_t dot_pos = first_entry.find(".");
    size_t e_pos = first_entry.find("e");
    if (dot_pos == std::string::npos || e_pos == std::string::npos)
    {
        EXCEPTION("Badly formatted scientific data field");
    }

    unsigned est_field_width = e_pos - dot_pos - 1 + 8; // = Precision + 8

    if (last_pos % est_field_width == 0)
    {
        mFieldWidth = est_field_width;
    }
    else
    {
        assert ( last_pos % (est_field_width+1) == 0  || (last_pos+1) % (est_field_width+1) == 0 );
        mFieldWidth = est_field_width+1;
    }
    infofile.close();
    datafile.close();
}

std::vector<double> ColumnDataReader::GetValues(const std::string& rVariableName)
{
    if (mNumFixedDimensions > 0)
    {
        EXCEPTION("Data file has fixed dimension which must be specified");
    }

    std::map<std::string, int>::iterator col = mVariablesToColumns.find(rVariableName);
    if (col == mVariablesToColumns.end())
    {
        std::stringstream variable_name;
        variable_name << rVariableName;
        EXCEPTION("'" + variable_name.str() + "' is an unknown variable.");
    }

    int column = (*col).second;
    ReadColumnFromFile(mDataFilename, column);

    return mValues;
}

std::vector<double> ColumnDataReader::GetValues(const std::string& rVariableName,
                                                int fixedDimension)
{
    if (mNumFixedDimensions < 1)
    {
        EXCEPTION("Data file has no fixed dimension");
    }

    mValues.clear();
    if (mHasUnlimitedDimension)
    {
        std::string datafile = mDataFilename;
        std::map<std::string, int>::iterator col = mVariablesToColumns.find(rVariableName);
        if (col == mVariablesToColumns.end())
        {
            EXCEPTION("Unknown variable");
        }
        int column = (*col).second;

        int counter = 1;
        while (true)
        {
            try
            {
                ReadValueFromFile(datafile, column, fixedDimension);
            }
            catch (Exception)
            {
                break;
            }

            // Advance counter
            std::string::size_type underscore_pos = datafile.rfind("_", datafile.length());
            std::stringstream suffix;

            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << counter;

            if (underscore_pos != std::string::npos)
            {
                datafile = datafile.substr(0, underscore_pos+1) + suffix.str() + ".dat";
            }
            counter++;
        }
    }
    else
    {
        int column = mVariablesToColumns[rVariableName];
        if (0 == column)
        {
            EXCEPTION("Unknown variable");
        }
        ReadValueFromFile(mDataFilename, column, fixedDimension);
    }

    return mValues;
}

std::vector<double> ColumnDataReader::GetUnlimitedDimensionValues()
{
    mValues.clear();
    if (!mHasUnlimitedDimension)
    {
        EXCEPTION("Data file has no unlimited dimension");
    }
    if (mNumFixedDimensions > 0)
    {
        // Read in from the ancillary file
        ReadColumnFromFile(mAncillaryFilename, 0);
    }
    else
    {
        // Read the first column
        ReadColumnFromFile(mDataFilename, 0);
    }
    return mValues;
}

void ColumnDataReader::ReadValueFromFile(const std::string& rFilename, int col, int row)
{
    std::ifstream datafile(rFilename.c_str(), std::ios::in);
    // If it doesn't exist - throw exception
    if (!datafile.is_open())
    {
        EXCEPTION("Couldn't open data file");
    }
    std::string variable_values;
    for (int i=0; i<row+1; i++)
    {
        std::getline(datafile, variable_values);
    }

    std::getline(datafile, variable_values);
    this->PushColumnEntryFromLine(variable_values, col);

    datafile.close();
}

void ColumnDataReader::ReadColumnFromFile(const std::string& rFilename, int col)
{
    // Empty the values vector
    mValues.clear();

    // Read in from the ancillary file
    std::ifstream datafile(rFilename.c_str(), std::ios::in);
    std::string value;

    // We should have already checked that this file can be opened.
    assert(datafile.is_open());

    // The current variable becomes true just after reading the last line
    bool end_of_file_reached = false;

    // Skip header line
    end_of_file_reached = std::getline(datafile, value).eof();

    while (!end_of_file_reached)
    {
        end_of_file_reached = std::getline(datafile, value).eof();
        this->PushColumnEntryFromLine(value, col);
    }
    datafile.close();
}

void ColumnDataReader::PushColumnEntryFromLine(const std::string& rLine, int col)
{
    std::string value;
    unsigned startpos = col * mFieldWidth;
    value = rLine.substr(startpos, mFieldWidth);

    std::stringstream variable_stream(value);
    double d_value;
    variable_stream >> d_value;
    if (variable_stream.fail())
    {
        if (variable_stream.eof()) //Missing data from column
        {
            d_value = DBL_MAX;
        }
        else
        {
// LCOV_EXCL_START
           //  Clang Objective C++ (on Mac OSX) treats reading very small numbers (<2e-308) as an error but other compilers just round to zero
           d_value = 0.0;
// LCOV_EXCL_STOP
        }
    }
    mValues.push_back(d_value);
}

bool ColumnDataReader::HasValues(const std::string& rVariableName)
{
    std::map<std::string, int>::iterator col = mVariablesToColumns.find(rVariableName);
    return !(col == mVariablesToColumns.end());
}

unsigned ColumnDataReader::GetFieldWidth()
{
    return mFieldWidth;
}
