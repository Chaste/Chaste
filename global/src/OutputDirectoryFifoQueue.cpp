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

#include "OutputDirectoryFifoQueue.hpp"
#include "OutputFileHandler.hpp"
#include <cassert>
#include <iostream>
#include "Exception.hpp"
#include "PetscTools.hpp"

OutputDirectoryFifoQueue::OutputDirectoryFifoQueue(const std::string& rBaseDirectory, unsigned queueMaxSize) :
    mBaseDirectory(rBaseDirectory),
    mQueueMaxSize(queueMaxSize)
{
    // Create the base directory
    OutputFileHandler handler(mBaseDirectory);

    // Paranoia
    assert(mQueue.empty());
}

std::string OutputDirectoryFifoQueue::CreateNextDir(const std::string& rSubdirectoryName)
{
    std::string subdirectory_full_name = mBaseDirectory + "/" + rSubdirectoryName;

    if (mQueue.size() == mQueueMaxSize)
    {
        std::string directory_to_remove =  OutputFileHandler::GetChasteTestOutputDirectory() + mBaseDirectory + "/" + mQueue.front();
        if (PetscTools::AmMaster())
        {
            ABORT_IF_NON0(system, "rm -rf "+directory_to_remove);
        }
        PetscTools::Barrier("OutputDirectoryFifoQueue::CreateNextDir");

        mQueue.pop();
    }

    mQueue.push(rSubdirectoryName);
    OutputFileHandler handler(subdirectory_full_name);

    assert(mQueue.size() <= mQueueMaxSize);

    return subdirectory_full_name;
}
