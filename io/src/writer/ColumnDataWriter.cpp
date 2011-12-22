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

/**
* Implementation file for ColumnDataWriter class.
*
*/

#include "ColumnDataWriter.hpp"
#include "ColumnDataConstants.hpp"
#include "Exception.hpp"
#include "Version.hpp"

#include <ctype.h>
#include <sstream>
#include <iomanip>
#include <fstream>

//#include <sys/stat.h> // For chmod()

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
      mpCurrentOutputFile(NULL),
      mpCurrentAncillaryFile(NULL),
      mpUnlimitedDimensionVariable(NULL),
      mpFixedDimensionVariable(NULL),
      mFieldWidth(precision+7), // Allow for numbers like -1.111e-321 (where precision=3)
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
    if (mpUnlimitedDimensionVariable != NULL)
    {
        delete mpUnlimitedDimensionVariable;
    }
    if (mpFixedDimensionVariable != NULL)
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
    if (mpCurrentOutputFile.get() != NULL)
    {
        mpCurrentOutputFile->close();
        mpCurrentOutputFile = out_stream(NULL);
    }

    if (mpCurrentAncillaryFile.get() != NULL)
    {
        mpCurrentAncillaryFile->close();
        mpCurrentAncillaryFile = out_stream(NULL);
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
    int unlimited_dimension_variable = (mpUnlimitedDimensionVariable != NULL);
    int fixed_dimension_variable = (mpFixedDimensionVariable != NULL);
    if (mIsUnlimitedDimensionSet)
    {
        if (mIsFixedDimensionSet)
        {
            mRowWidth = (mVariables.size() + fixed_dimension_variable) * (mFieldWidth + SPACING);
            mAncillaryRowWidth = mFieldWidth + SPACING;

            // Write out the headers for the first position along the unlimited dimension
            std::stringstream suffix;
            suffix << std::setfill('0') << std::setw(FILE_SUFFIX_WIDTH) << mUnlimitedDimensionPosition;

            if (mpUnlimitedDimensionVariable != NULL)
            {
                std::string ancillary_filename = mBaseName + "_unlimited.dat";
                mpCurrentAncillaryFile = mOutputFileHandler.OpenOutputFile(ancillary_filename, std::ios::out);
                (*mpCurrentAncillaryFile) << std::setiosflags(std::ios::scientific);
                (*mpCurrentAncillaryFile) << std::setprecision(mPrecision);
                if (mpUnlimitedDimensionVariable != NULL)
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
            if (mpUnlimitedDimensionVariable != NULL)
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
    mpCurrentOutputFile = mOutputFileHandler.OpenOutputFile(rFileName, std::ios::out);
    (*mpCurrentOutputFile) << std::setiosflags(std::ios::scientific);
    (*mpCurrentOutputFile) << std::setprecision(mPrecision);
    if (mpFixedDimensionVariable != NULL)
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
    out_stream p_info_file = mOutputFileHandler.OpenOutputFile(rFileName, std::ios::out);
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
                    position = mRowStartPosition + (mRowWidth+1) * dimensionPosition + SPACING - 1;
                }
                else
                {
                    // ordinary variables
                    position = mRowStartPosition + (mRowWidth+1) * dimensionPosition +
                               ((variableID + (mpFixedDimensionVariable != NULL)) * (mFieldWidth + SPACING)) + SPACING - 1;
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
                position = mRowStartPosition + SPACING - 1;
            }
            else
            {
                position = (variableID + (mpUnlimitedDimensionVariable != NULL)) * (mFieldWidth + SPACING) +
                           mRowStartPosition + SPACING - 1;
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
            position = mRowStartPosition + (mRowWidth+1) * dimensionPosition + SPACING - 1;
        }
        else
        {
            position = mRowStartPosition + (mRowWidth+1) * dimensionPosition +
                       ((variableID + (mpFixedDimensionVariable != NULL)) * (mFieldWidth + SPACING)) + SPACING - 1;
        }
        mpCurrentOutputFile->seekp(position);
        mpCurrentOutputFile->width(mFieldWidth);
        (*mpCurrentOutputFile) << variableValue;
    }

    mHasPutVariable = true;
}
