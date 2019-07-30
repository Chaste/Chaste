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

#include <cassert>

#include "ProgressReporter.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

ProgressReporter::ProgressReporter(std::string outputDirectory, double startTime, double endTime)
    : mStartTime(startTime),
      mEndTime(endTime),
      mLastPercentage(UINT_MAX)
{
    assert(startTime < endTime);

    // Note we make sure we don't delete anything in the output directory
    OutputFileHandler handler(outputDirectory, false);

    // Open the file on the master process only
    if (PetscTools::AmMaster())
    {
        mpFile = handler.OpenOutputFile("progress_status.txt");
    }
}

ProgressReporter::~ProgressReporter()
{
    if (PetscTools::AmMaster())
    {
        if (mLastPercentage!=100)
        {
            *mpFile << "100% completed" << std::endl;
        }
        *mpFile << "..done!" << std::endl;
        mpFile->close();
    }
}

void ProgressReporter::Update(double currentTime)
{
    unsigned percentage = (unsigned)( (currentTime - mStartTime)/(mEndTime - mStartTime)*100 );
    if (mLastPercentage==UINT_MAX || percentage > mLastPercentage)
    {
        if (PetscTools::AmMaster())
        {
            *mpFile << percentage << "% completed" << std::endl;
        }
        mLastPercentage = percentage;
    }
}

void ProgressReporter::PrintFinalising()
{
    if (PetscTools::AmMaster())
    {
        *mpFile << "Finalising.." << std::endl;
    }
}

void ProgressReporter::PrintInitialising()
{
    if (PetscTools::AmMaster())
    {
        *mpFile << "Initialising.." << std::endl;
    }
}
