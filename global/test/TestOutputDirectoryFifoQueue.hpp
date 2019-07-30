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

#ifndef TESTOUTPUTDIRECTORYFIFOQUEUE_HPP_
#define TESTOUTPUTDIRECTORYFIFOQUEUE_HPP_

#include <cxxtest/TestSuite.h>
#include "OutputFileHandler.hpp"
#include "OutputDirectoryFifoQueue.hpp"
#include "FileFinder.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestOutputDirectoryFifoQueue : public CxxTest::TestSuite
{
public:

    void TestQueueCreatesDirectories()
    {
        FileFinder checkpoints("checkpoints", RelativeTo::ChasteTestOutput);
        // Remove directory in case it was there from previous executions.
        if (PetscTools::AmMaster())
        {
            ABORT_IF_THROWS(checkpoints.Remove());
        }
        PetscTools::Barrier("TestQueueCreatesDirectories-1");
        TS_ASSERT(!checkpoints.Exists());
        PetscTools::Barrier("TestQueueCreatesDirectories-2");

        OutputDirectoryFifoQueue fifo_queue("checkpoints", 2);
        TS_ASSERT(checkpoints.IsDir());

        fifo_queue.CreateNextDir("0.1");
        FileFinder dir1("0.1", checkpoints);
        TS_ASSERT(dir1.IsDir());

        fifo_queue.CreateNextDir("0.2");
        FileFinder dir2("0.2", checkpoints);
        TS_ASSERT(dir2.IsDir());
    }

    void TestQueueRemovesAndCreatesDirectories()
    {
        FileFinder checkpoints("checkpoints2", RelativeTo::ChasteTestOutput);
        // Remove directory in case it was there from previous executions.
        PetscTools::Barrier("TestQueueRemovesAndCreatesDirectories-0");
        if (PetscTools::AmMaster())
        {
            ABORT_IF_THROWS(checkpoints.Remove());
        }
        PetscTools::Barrier("TestQueueRemovesAndCreatesDirectories-1");
        TS_ASSERT(!checkpoints.Exists());
        PetscTools::Barrier("TestQueueRemovesAndCreatesDirectories-2");

        OutputDirectoryFifoQueue fifo_queue("checkpoints2", 2);
        TS_ASSERT(checkpoints.IsDir());

        fifo_queue.CreateNextDir("0.1");
        FileFinder dir1("0.1", checkpoints);
        TS_ASSERT(dir1.IsDir());

        fifo_queue.CreateNextDir("0.2");
        FileFinder dir2("0.2", checkpoints);
        TS_ASSERT(dir2.IsDir());
        TS_ASSERT(dir1.IsDir());

        PetscTools::Barrier("TestQueueRemovesAndCreatesDirectories-3");
        fifo_queue.CreateNextDir("0.3");
        FileFinder dir3("0.3", checkpoints);
        TS_ASSERT(dir3.IsDir());
        TS_ASSERT(dir2.IsDir());
        TS_ASSERT(!dir1.Exists());

        PetscTools::Barrier("TestQueueRemovesAndCreatesDirectories-4");
        fifo_queue.CreateNextDir("0.4");
        FileFinder dir4("0.4", checkpoints);
        TS_ASSERT(dir4.IsDir());
        TS_ASSERT(dir3.IsDir());
        TS_ASSERT(!dir2.Exists());
        TS_ASSERT(!dir1.Exists());
    }
};

#endif /*TESTOUTPUTDIRECTORYFIFOQUEUE_HPP_*/
