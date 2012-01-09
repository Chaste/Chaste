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
    if ( PetscTools::AmMaster() )
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
        if ( PetscTools::AmMaster() )
        {
            *mpFile << percentage << "% completed" << std::endl;
        }
        mLastPercentage = percentage;
    }
}

void ProgressReporter::PrintFinalising()
{
    if ( PetscTools::AmMaster() )
    {
        *mpFile << "Finalising.." << std::endl;
    }
}

void ProgressReporter::PrintInitialising()
{
    if ( PetscTools::AmMaster() )
    {
        *mpFile << "Initialising.." << std::endl;
    }
}
