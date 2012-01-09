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
#ifndef PROGRESSREPORTER_HPP_
#define PROGRESSREPORTER_HPP_

#include "OutputFileHandler.hpp"

/**
 * This class creates the file 'progress_status.txt' in the specified directory
 * and writes "n% completed" etc in the file when n% of a simulation has been
 * done, for integer n.
 *
 * You can watch the progress of your simulation by doing one of
 * a) watch tail \<outputDirectory>/progress_status.txt
 * b) tail -f \<outputDirectory>/progress_status.txt
 */
class ProgressReporter
{
private:

    double mStartTime;        /**< Start time of the simulation */
    double mEndTime;          /**< End time of the simulation */
    out_stream mpFile;        /**< Progress status file */
    unsigned mLastPercentage; /**< Last percentage that was written */

public:

    /**
     * Constuctor saves times and opens output file ('progress_status.txt').
     *
     * @param outputDirectory where to open the output file
     * @param startTime the start time
     * @param endTime the end time
     */
    ProgressReporter(std::string outputDirectory, double startTime, double endTime);

    /**
     * Destructor.
     */
    ~ProgressReporter();

    /**
     * Calculates the percentage completed using the time given and the start and end
     * time and prints to file if another percent has been done.
     *
     * @param currentTime the given time
     */
    void Update(double currentTime);

    /**
     * Print finalising message to file.
     */
    void PrintFinalising();

    /**
     * Print initialising message to file.
     */
    void PrintInitialising();
};

#endif /*PROGRESSREPORTER_HPP_*/
