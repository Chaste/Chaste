/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "Exception.hpp"
#include "AbstractHdf5Access.hpp"

bool AbstractHdf5Access::DoesDatasetExist(const std::string& rDatasetName)
{
#if H5_VERS_MAJOR>=1 && H5_VERS_MINOR>=8
    // This is a nice method for testing existence, introduced in HDF5 1.8.0
    htri_t dataset_status = H5Lexists(mFileId, rDatasetName.c_str(), H5P_DEFAULT);
    return (dataset_status>0);
#else
    bool result=false;
    // This is not a nice way of doing it because the error stack produces a load of 'HDF failed' output.
    // The "TRY" macros are a convenient way to temporarily turn the error stack off.
    H5E_BEGIN_TRY
    {
        hid_t dataset_id = H5Dopen(mFileId, rDatasetName.c_str());
        if (dataset_id>0)
        {
            H5Dclose(dataset_id);
            result = true;
        }
    }
    H5E_END_TRY;

    return result;
#endif
}

void AbstractHdf5Access::SetTimeDatasetId()
{
    // Now deal with time
    assert(mUnlimitedDimensionName=="Time");
    // Files being extended/read are always assumed to have an unlimited dimension
    // called "Time" (even though this writer can vary this!).
    mUnlimitedDimensionUnit = "ms"; // Assumed by Chaste...

    // Files pre - r16738 use "Time" for "Data"'s unlimited variable.
    // Files post - r16738 use "<DatasetName>_Time" for "<DatasetName>"'s unlimited variable, if this is missing we look
    // for simply "Time".
    if (DoesDatasetExist(mDatasetName + "_" + mUnlimitedDimensionName))
    {
        mTimeDatasetId = H5Dopen(mFileId, (mDatasetName + "_" + mUnlimitedDimensionName).c_str());
    }
    else if (DoesDatasetExist(mUnlimitedDimensionName))
    {
        mTimeDatasetId = H5Dopen(mFileId, mUnlimitedDimensionName.c_str());
    }
    else
    {
        NEVER_REACHED;
    }
}

AbstractHdf5Access::AbstractHdf5Access(const std::string& rDirectory,
                   const std::string& rBaseName,
                   const std::string& rDatasetName,
                   bool makeAbsolute)
 : mBaseName(rBaseName),
   mDatasetName(rDatasetName),
   mIsDataComplete(true),
   mUnlimitedDimensionName("Time"),
   mUnlimitedDimensionUnit("ms"),
   mIsUnlimitedDimensionSet(false)
{
    RelativeTo::Value relative_to;
    if (makeAbsolute)
    {
        relative_to = RelativeTo::ChasteTestOutput;
    }
    else
    {
        relative_to = RelativeTo::Absolute;
    }
    mDirectory.SetPath(rDirectory, relative_to);
}

AbstractHdf5Access::AbstractHdf5Access(const FileFinder& rDirectory,
                   const std::string& rBaseName,
                   const std::string& rDatasetName)
 : mBaseName(rBaseName),
   mDatasetName(rDatasetName),
   mDirectory(rDirectory),
   mIsDataComplete(true),
   mUnlimitedDimensionName("Time"),
   mUnlimitedDimensionUnit("ms"),
   mIsUnlimitedDimensionSet(false)
{
}

AbstractHdf5Access::~AbstractHdf5Access()
{
}


bool AbstractHdf5Access::IsDataComplete()
{
    return mIsDataComplete;
}


std::vector<unsigned> AbstractHdf5Access::GetIncompleteNodeMap()
{
    return mIncompleteNodeIndices;
}



