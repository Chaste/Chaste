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
