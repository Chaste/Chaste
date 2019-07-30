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
* Implementation file for ColumnDataWriter class.
*
*/

#include <ctype.h>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "ColumnDataWriter.hpp"
#include "ColumnDataConstants.hpp"
#include "Exception.hpp"
#include "Version.hpp"


ColumnDataWriter::ColumnDataWriter(const std::string& rDirectory,
                                   const std::string& rBaseName,
                                   bool cleanDirectory,
                                   unsigned precision)
    : mOutputFileHandler(rDirectory, cleanDirectory),
      mDirectory(rDirectory),
      mBaseName(rBaseName),
      mIsInDefineMode(true),
      mIsFixedDimensionSet(false),
      mIsUnlimitedDimensionSet(false),
      mUnlimitedDimensionPosition(0),
      mFixedDimensionSize(-1),
      mpCurrentOutputFile(nullptr),
      mpCurrentAncillaryFile(nullptr),
      mpUnlimitedDimensionVariable(nullptr),
      mpFixedDimensionVariable(nullptr),
      mFieldWidth(precision+8),
      mPrecision(precision),
      mHasPutVariable(false),
      mNeedAdvanceAlongUnlimitedDimension(false),
      mCommentForInfoFile("")
{
    if (mPrecision<2 || mPrecision>20)
    {
        EXCEPTION("Precision must be between 2 and 20 (inclusive)");
    }
}

ColumnDataWriter::~ColumnDataWriter()
{
    // Close any open output files
    Close();

    // Delete memory allocated for variables
    if (mpUnlimitedDimensionVariable != nullptr)
    {
        delete mpUnlimitedDimensionVariable;
    }
    if (mpFixedDimensionVariable != nullptr)
    {
        delete mpFixedDimensionVariable;
    }
}

std::string ColumnDataWriter::GetOutputDirectory()
{
    return mOutputFileHandler.GetOutputDirectoryFullPath();
}

void ColumnDataWriter::Close()
{
    if (mpCurrentOutputFile.get() != nullptr)
    {
        mpCurrentOutputFile->close();
        mpCurrentOutputFile = out_stream(nullptr);
    }

    if (mpCurrentAncillaryFile.get() != nullptr)
    {
        mpCurrentAncillaryFile->close();
        mpCurrentAncillaryFile = out_stream(nullptr);
    }
}

void ColumnDataWriter::CheckVariableName(const std::string& rName)
{
    if (rName.length() == 0)
    {
        EXCEPTION("Variable name not allowed: may not be blank.");
    }
    CheckUnitsName(rName);
}

void ColumnDataWriter::CheckUnitsName(const std::string& rName)
{
    for (unsigned i=0; i<rName.length(); i++)
    {
        if (!isalnum(rName[i]) && !(rName[i]=='_'))
        {
            std::string error = "Variable name/units '" + rName + "' not allowed: may only contain alphanumeric characters or '_'.";
            EXCEPTION(error);
        }
    }
}

int ColumnDataWriter::DefineUnlimitedDimension(const std::string& rDimensionName,
                                               const std::string& rDimensionUnits)
{
    if (mIsUnlimitedDimensionSet)
    {
        EXCEPTION("Unlimited dimension already set. Cannot be defined twice");
    }

    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }

    CheckVariableName(rDimensionName);
    CheckUnitsName(rDimensionUnits);

    mUnlimitedDimensionName = rDimensionName;
    mUnlimitedDimensionUnits = rDimensionUnits;

    mpUnlimitedDimensionVariable = new DataWriterVariable;
    mpUnlimitedDimensionVariable->mVariableName = rDimensionName;
    mpUnlimitedDimensionVariable->mVariableUnits = rDimensionUnits;

    mIsUnlimitedDimensionSet = true;

    return UNLIMITED_DIMENSION_VAR_ID;
}

int ColumnDataWriter::DefineFixedDimension(const std::string& rDimensionName,
                                           const std::string& rDimensionUnits,
                                           long dimensionSize)
{
    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }
    if (dimensionSize < 1)
    {
        EXCEPTION("Fixed dimension must be at least 1 long");
    }

    CheckVariableName(rDimensionName);
    CheckUnitsName(rDimensionUnits);

    mFixedDimensionName = rDimensionName;
    mFixedDimensionUnits = rDimensionUnits;
    mFixedDimensionSize = dimensionSize;

    mIsFixedDimensionSet = true;

    mpFixedDimensionVariable = new DataWriterVariable;
    mpFixedDimensionVariable->mVariableName = rDimensionName;
    mpFixedDimensionVariable->mVariableUnits = rDimensionUnits;
    return FIXED_DIMENSION_VAR_ID;
}

int ColumnDataWriter::DefineVariable(const std::string& rVariableName,
                                     const std::string& rVariableUnits)
{
    if (!mIsInDefineMode)
    {
        EXCEPTION("Cannot define variables when not in Define mode");
    }

    CheckVariableName(rVariableName);
    CheckUnitsName(rVariableUnits);

    int variable_id;

    if (rVariableName == mUnlimitedDimensionName)
    {
        EXCEPTION("Variable name: " + rVariableName + " already in use as unlimited dimension");
    }
    else if (rVariableName == mFixedDimensionName)
    {
        EXCEPTION("Variable name: " + rVariableName + " already in use as fixed dimension");
    }
    else // ordinary variable
    {
        // Add the variable to the variable vector
        DataWriterVariable new_variable;
        new_variable.mVariableName = rVariableName;
        new_variable.mVariableUnits = rVariableUnits;
        mVariables.push_back(new_variable);

        // Use the index of the variable vector as the variable ID.
        // This is ok since there is no way to remove variables.
        variable_id = mVariables.size()-1;
    }

    return variable_id;
}

void ColumnDataWriter::EndDefineMode()
{
    // Check that a dimension has been defined
    if (mIsFixedDimensionSet == false && mIsUnlimitedDimensionSet == false)
    {
        EXCEPTION("Cannot end define mode. No dimensions have been defined.");
    }
    // Check that at least one variable has been defined
    if (mVariables.size() < 1)
    {
        EXCEPTION("Cannot end define mode. No variables have been defined.");
    }
    // Calculate the width of each row
    int unlimited_dimension_variable = (mpUnlimitedDimensionVariable != nullptr);
    int fixed_dimension_variable = (mpFixedDimensionVariable != nullptr);
    if (mIsUnlimitedDimensionSet)
    {
        if (mIsFixedDimensionSet)
        {
            mRowWidth = (mVariables.size() + fixed_dimension_variable) * (mFieldWidth + SPACING);
            mAncillaryRowWidth = mFieldWidth + SPACING;

            // Write out the headers for the first position along the unlimited dimension
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << mUnlimitedDimensionPosition;

            if (mpUnlimitedDimensionVariable != nullptr)
            {
                std::string ancillary_filename = mBaseName + "_unlimited.dat";
                mpCurrentAncillaryFile = mOutputFileHandler.OpenOutputFile(ancillary_filename, std::ios::out | std::ios::binary);
                (*mpCurrentAncillaryFile) << std::setiosflags(std::ios::scientific);
                (*mpCurrentAncillaryFile) << std::setprecision(mPrecision);
                if (mpUnlimitedDimensionVariable != nullptr)
                {
                    (*mpCurrentAncillaryFile) << mpUnlimitedDimensionVariable->mVariableName
                                              << "(" << mpUnlimitedDimensionVariable->mVariableUnits << ") ";
                }
            }
            mAncillaryRowStartPosition = mpCurrentAncillaryFile->tellp();
            std::string filename = mBaseName + "_" + suffix.str() + ".dat";
            this->CreateFixedDimensionFile(filename);
        }
        else
        {
            mRowWidth = (mVariables.size() + unlimited_dimension_variable) * (mFieldWidth + SPACING);

            // Write out the column headers
            std::string filename = mBaseName + ".dat";
            mpCurrentOutputFile = mOutputFileHandler.OpenOutputFile(filename, std::ios::out);
            (*mpCurrentOutputFile) << std::setiosflags(std::ios::scientific);
            (*mpCurrentOutputFile) << std::setprecision(mPrecision);
            if (mpUnlimitedDimensionVariable != nullptr)
            {
                (*mpCurrentOutputFile) << mpUnlimitedDimensionVariable->mVariableName
                                       << "(" << mpUnlimitedDimensionVariable->mVariableUnits << ") ";
            }
            /*
             * Write out header(which may contain several variables) for output file.
             * In this scope the method "CreateFixedDimensionFile" has not been invoked,
             * because there is no mFixedDimensionSize available.
             */
            for (unsigned i=0; i<mVariables.size(); i++)
            {
                (*mpCurrentOutputFile) << mVariables[i].mVariableName << "(" << mVariables[i].mVariableUnits << ")";
                if (i < mVariables.size()-1)
                {
                    (*mpCurrentOutputFile) << " ";
                }
            }
            (*mpCurrentOutputFile) << std::endl;
            mRowStartPosition = mpCurrentOutputFile->tellp();

            // Write out a line of blank space which is #variables * (mFieldWidth + 1) -1
            std::string blank_line(mRowWidth, ' ');
            (*mpCurrentOutputFile) << blank_line;
        }
    }
    else
    {
        // The fixed dimension must be set at this point or we wouldn't be here
        mRowWidth = (mVariables.size() + fixed_dimension_variable) * (mFieldWidth + SPACING);
        std::string filename = mBaseName + ".dat";
        this->CreateFixedDimensionFile(filename);
    }

    // Write info file
    std::string infoname = mBaseName + ".info";
    this->CreateInfoFile(infoname);

    mIsInDefineMode = false;
}

void ColumnDataWriter::CreateFixedDimensionFile(const std::string& rFileName)
{
    // Create new data file
    mpCurrentOutputFile = mOutputFileHandler.OpenOutputFile(rFileName, std::ios::out | std::ios::binary);
    (*mpCurrentOutputFile) << std::setiosflags(std::ios::scientific);
    (*mpCurrentOutputFile) << std::setprecision(mPrecision);
    if (mpFixedDimensionVariable != nullptr)
    {
        (*mpCurrentOutputFile) << mpFixedDimensionVariable->mVariableName
                               << "(" << mpFixedDimensionVariable->mVariableUnits << ") ";
    }
    // Write out the column headers and spaces for the rest of the file
    for (unsigned i = 0; i < mVariables.size(); i++)
    {
        (*mpCurrentOutputFile) << mVariables[i].mVariableName << "(" << mVariables[i].mVariableUnits << ")";
        if (i < mVariables.size()-1)
        {
            (*mpCurrentOutputFile) << " ";
        }
    }
    (*mpCurrentOutputFile) << std::endl;
    mRowStartPosition = mpCurrentOutputFile->tellp();
    std::string blank_line(mRowWidth, ' ');
    for (int i = 0; i < mFixedDimensionSize; i++)
    {
        (*mpCurrentOutputFile) << blank_line << std::endl;
    }
}

void ColumnDataWriter::CreateInfoFile(const std::string& rFileName)
{
    // Create new info file
    out_stream p_info_file = mOutputFileHandler.OpenOutputFile(rFileName, std::ios::out | std::ios::binary);
    (*p_info_file) << "FIXED " << mFixedDimensionSize << std::endl;
    (*p_info_file) << "UNLIMITED " << mIsUnlimitedDimensionSet << std::endl;
    (*p_info_file) << "VARIABLES " << mVariables.size() << std::endl;
    if (mCommentForInfoFile != "")
    {
        *p_info_file << mCommentForInfoFile << std::endl;
    }
    *p_info_file << ChasteBuildInfo::GetProvenanceString();
    p_info_file->close();
}

void ColumnDataWriter::DoAdvanceAlongUnlimitedDimension()
{
    mHasPutVariable = false;
    mNeedAdvanceAlongUnlimitedDimension = false;

    if (mIsUnlimitedDimensionSet)
    {
        if (mIsFixedDimensionSet)
        {
            //first close the current file before creating another one
            mpCurrentOutputFile->close();
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << mUnlimitedDimensionPosition + 1;

            std::string filename = mBaseName + "_" + suffix.str() + ".dat";
            this->CreateFixedDimensionFile(filename);
        }
        else
        {
            //go to the end of the current line
            mpCurrentOutputFile->seekp(mRowStartPosition+mRowWidth);
            (*mpCurrentOutputFile) << std::endl;
            mRowStartPosition = mpCurrentOutputFile->tellp();
            std::string blank_line(mRowWidth,' ');
            (*mpCurrentOutputFile) << blank_line;
        }
    }
    else
    {
        EXCEPTION("Cannot advance along unlimited dimension if it is not defined");
    }
    mUnlimitedDimensionPosition++;
}

void ColumnDataWriter::AdvanceAlongUnlimitedDimension()
{
    if (mHasPutVariable)
    {
        mNeedAdvanceAlongUnlimitedDimension = true;
    }
}

void ColumnDataWriter::PutVariable(int variableID, double variableValue, long dimensionPosition)
{

    if (mNeedAdvanceAlongUnlimitedDimension)
    {
        DoAdvanceAlongUnlimitedDimension();
    }

    // Check that we are not in define mode
    if (mIsInDefineMode)
    {
        EXCEPTION("Cannot put variables when in Define mode");
    }
    // Check that variableID is in range (exception)
    if (variableID > (int)mVariables.size() ||
        (variableID != UNLIMITED_DIMENSION_VAR_ID &&
         variableID != FIXED_DIMENSION_VAR_ID &&
         variableID < 0))
    {
        EXCEPTION("variableID unknown");
    }

    if (mIsFixedDimensionSet)
    {
        if (dimensionPosition == -1 && variableID != UNLIMITED_DIMENSION_VAR_ID)
        {
            EXCEPTION("Dimension position not supplied");
        }
        if (dimensionPosition < -1 || dimensionPosition >= mFixedDimensionSize)
        {
            EXCEPTION("Dimension position out of range");
        }
        if (dimensionPosition != -1 && variableID == UNLIMITED_DIMENSION_VAR_ID)
        {
            EXCEPTION("Dimension position supplied, but not required");
        }
    }

    if (mIsUnlimitedDimensionSet)
    {
        if (mIsFixedDimensionSet)
        {
            // Go to the correct position in the file
            if (variableID == UNLIMITED_DIMENSION_VAR_ID)
            {
                (*mpCurrentAncillaryFile) << std::endl << " ";
                mpCurrentAncillaryFile->width(mFieldWidth);
                (*mpCurrentAncillaryFile) << variableValue;
            }
            else
            {
                int position;
                if (variableID == FIXED_DIMENSION_VAR_ID)
                {
                    position = mRowStartPosition + (mRowWidth+SPACING) * dimensionPosition;
                }
                else
                {
                    // ordinary variables
                    position = mRowStartPosition + (mRowWidth+SPACING) * dimensionPosition +
                               ((variableID + (mpFixedDimensionVariable != nullptr)) * (mFieldWidth + SPACING));
                }

                mpCurrentOutputFile->seekp(position);
                mpCurrentOutputFile->width(mFieldWidth);
                (*mpCurrentOutputFile) << variableValue;
            }
        }
        else
        {
            // Go to the correct position in the file
            int position;
            if (variableID == UNLIMITED_DIMENSION_VAR_ID)
            {
                position = mRowStartPosition;
            }
            else
            {
                position = (variableID + (mpUnlimitedDimensionVariable != nullptr)) * (mFieldWidth + SPACING) +
                           mRowStartPosition;
            }

            mpCurrentOutputFile->seekp(position);
            mpCurrentOutputFile->width(mFieldWidth);
            (*mpCurrentOutputFile) << variableValue;

        }
    }
    else
    {
        // Go to the correct position in the file
        int position;
        if (variableID == FIXED_DIMENSION_VAR_ID)
        {
            position = mRowStartPosition + (mRowWidth+SPACING) * dimensionPosition;
        }
        else
        {
            position = mRowStartPosition + (mRowWidth+SPACING) * dimensionPosition +
                       ((variableID + (mpFixedDimensionVariable != nullptr)) * (mFieldWidth + SPACING));
        }
        mpCurrentOutputFile->seekp(position);
        mpCurrentOutputFile->width(mFieldWidth);
        (*mpCurrentOutputFile) << variableValue;
    }

    mHasPutVariable = true;
}
