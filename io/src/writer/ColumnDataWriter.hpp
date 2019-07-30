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

#ifndef COLUMNDATAWRITER_HPP_
#define COLUMNDATAWRITER_HPP_

#include <string>
#include <vector>

#include "AbstractDataWriter.hpp"
#include "DataWriterVariable.hpp"
#include "OutputFileHandler.hpp"

/**
 * A concrete column data writer class.  Writes grid-formatted data in
 * space separated column form.  Each file has a header row with names
 * and optional units for each column.
 */
class ColumnDataWriter : public AbstractDataWriter
{
protected:

    OutputFileHandler mOutputFileHandler; /**< For opening data files. */

    std::string mDirectory; /**< Directory output files will be stored in. */
    std::string mBaseName; /**< The base name for the output data files. */
    bool mIsInDefineMode; /**< Is the DataWriter in define mode or not */
    bool mIsFixedDimensionSet; /**< Is the fixed dimension set */
    bool mIsUnlimitedDimensionSet; /**< Is the unlimited dimension set */
    long mUnlimitedDimensionPosition; /**< The position along the unlimited dimension that writing of variables will take place */
    long mFixedDimensionSize; /**< The size of the fixed dimension */
    out_stream mpCurrentOutputFile; /**< Filestream currently being addressed */
    out_stream mpCurrentAncillaryFile; /**< Ancillary filestream currently being addressed (required for two dimensional output) eg. time file*/
    DataWriterVariable* mpUnlimitedDimensionVariable; /**< The variable corresponding to the unlimited dimension */
    DataWriterVariable* mpFixedDimensionVariable; /**< The variable corresponding to the fixed dimension */

    std::string mUnlimitedDimensionName; /**< The name of the unlimited dimension. */
    std::string mUnlimitedDimensionUnits; /**< The physical units of the unlimited dimension. */

    std::string mFixedDimensionName; /**< The name of the fixed dimension */
    std::string mFixedDimensionUnits; /**< The units of the fixed dimension */

    std::vector<DataWriterVariable> mVariables; /**< The data variables */

    const unsigned mFieldWidth; /**< Width of each column in the text file (excludes column headers) */
    const unsigned mPrecision; /**< Precision used in writing the data */
    static const int SPACING = 1; /**< Space between columns (includes minus sign) */
    static const int FIXED_DIMENSION_VAR_ID = -1; /**< id of fixed dimension variable */
    static const int UNLIMITED_DIMENSION_VAR_ID = -2; /**< id of unlimited dimension variable */

    std::string mFileExtension; /**< Extension of output files */

    int mRowStartPosition; /**< The position of the file pointer when it's at the beginning of the current row */
    int mRowWidth; /**< The width in characters of a row in the file */

    int mAncillaryRowStartPosition; /**< The position of the ancillary file pointer when it's at the beginning of the current row */
    int mAncillaryRowWidth; /**< The width in characters of a row in the ancillary file */

    bool mHasPutVariable; /**< Whether a variable value has been output to a file. */
    bool mNeedAdvanceAlongUnlimitedDimension; /**< Whether we need to advance along the unlimited dimension. */

    /**
     *  (Optional) comment that can be set and will be written to the info file
     *  when CreateInfoFile() is called in EndDefineMode()
     */
    std::string mCommentForInfoFile;

    /**
     * Create the output file and write out the header for it.
     *
     * @param rFileName  the name of the file to write to, relative to the output directory
     */
    void CreateFixedDimensionFile(const std::string& rFileName);

    /**
     * Create the info file.
     *
     * @param rFileName  the name of the file to create, relative to the output directory
     */
    void CreateInfoFile(const std::string& rFileName);

    /**
     * Check name of variable is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     *
     * @param rName variable name
     */
    void CheckVariableName(const std::string& rName);

    /**
     * Check name of unit is allowed, i.e. contains only alphanumeric & _, and isn't blank.
     *
     * @param rName unit name
     */
    void CheckUnitsName(const std::string& rName);

    /**
     * Advance along the unlimited dimension. Normally this will be called
     * when all variables in a row have been input.
     */
    void DoAdvanceAlongUnlimitedDimension();

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the data to file
     * @param rBaseName  the name of the file in which to write the data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     * @param precision the precision with which to write the data (i.e. exactly
     *    how many digits to display after the decimal point).  Defaults to 8.
     *    Must be between 2 and 20 (inclusive).
     */
    ColumnDataWriter(const std::string& rDirectory,
                     const std::string& rBaseName,
                     bool cleanDirectory=true,
                     unsigned precision=8);

    /**
     * Destructor. Closes any open files.
     */
    virtual ~ColumnDataWriter();

    /**
     * Define the unlimited dimension, i.e. the dimension that increases as the simulation progresses.
     *
     * @param rDimensionName The name of the unlimited dimension
     * @param rDimensionUnits The physical units of the unlimited dimension
     *
     * @return The identifier of the variable
     */
    int DefineUnlimitedDimension(const std::string& rDimensionName,
                                 const std::string& rDimensionUnits);

    /**
     * Define the fixed dimension.
     *
     * @param rDimensionName The name of the dimension
     * @param rDimensionUnits The physical units of the dimension
     * @param dimensionSize The size of the dimension
     *
     * @return The identifier of the variable
     */
    int DefineFixedDimension(const std::string& rDimensionName,
                             const std::string& rDimensionUnits,
                             long dimensionSize);

    /**
     * Define a variable.
     *
     * @param rVariableName The name of the variable
     * @param rVariableUnits The physical units of the variable
     *
     * @return The identifier of the variable
     */
    int DefineVariable(const std::string& rVariableName,
                       const std::string& rVariableUnits);

    /**
     * Set a comment to be written in the info file (optional).
     * This needs to be called before EndDefineMode().
     *
     * @param comment  the comment
     */
    void SetCommentForInfoFile(std::string comment)
    {
        mCommentForInfoFile = comment;
    }

    /**
     * End the define mode of the DataWriter.
     */
    virtual void EndDefineMode();

    /**
     *  Dummy function for DoAdvanceAlongUnlimitedDimension.
     */
    virtual void AdvanceAlongUnlimitedDimension();

    /**
     * Input the variable value to the output file or ancillary file.
     *
     * @param variableID
     * @param variableValue
     * @param dimensionPosition  The position in column (defaults to -1). This is required if
     *      there is a fixed dimension, and will be the position along that dimension
     */
    virtual void PutVariable(int variableID, double variableValue, long dimensionPosition = -1);

    /**
     * Close any open files.
     */
    virtual void Close();

    /**
     * @return the full pathname of the directory where we're writing files.
     */
    std::string GetOutputDirectory();
};

#endif //COLUMNDATAWRITER_HPP_
