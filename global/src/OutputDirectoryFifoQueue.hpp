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

#ifndef OUTPUTDIRECTORYFIFOQUEUE_HPP_
#define OUTPUTDIRECTORYFIFOQUEUE_HPP_

#include <string>
#include <queue>

/**
 * This is a helper class to handle a FIFO collection of subdirectories.
 *
 * All the subdirectories will be created inside a base directory provided
 * in the constructor. The maximum number of concurrent subdirectories is
 * specified in the the constructor. Once this number is reached, the next
 * call to CreateNextDir() will delete the oldest directory as a side effect.
 */
class OutputDirectoryFifoQueue
{
private:

    std::string mBaseDirectory; /**< Base directory for all the subdirectories to be created. */
    unsigned mQueueMaxSize; /**< Maximum number of subdirectories*/
    std::queue<std::string> mQueue;  /**<The queue of names of subdirectories currently on the disk*/

public:

    /**
     * Constructor.
     *
     * @param rBaseDirectory base directory for all the subdirectories to be created
     * @param queueMaxSize maximum number of subdirectories
     */
    OutputDirectoryFifoQueue(const std::string& rBaseDirectory, unsigned queueMaxSize);

    /**
     * Creates a subdirectory called rSubdirectoryName deleting the oldest subdirectory
     * if the maximum number has been reached.
     *
     * @note Must be called collectively.
     *
     * @param rSubdirectoryName subdirectory name
     * @return new directory name relative to the base directory
     */
    std::string CreateNextDir(const std::string& rSubdirectoryName);
};

#endif /*OUTPUTDIRECTORYFIFOQUEUE_HPP_*/
