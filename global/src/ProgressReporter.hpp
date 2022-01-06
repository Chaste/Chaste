/*

Copyright (c) 2005-2021, University of Oxford.
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
#ifndef PROGRESSREPORTER_HPP_
#define PROGRESSREPORTER_HPP_

#include <iostream>
#include <ctime>
#include <limits.h>

#include "Timer.hpp"
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

    /** Allow testing of private member functions */
    friend class TestProgressReporter;

    /** Start time of the simulation */
    double mStartTime;

    /** End time of the simulation */
    double mEndTime;

    /** Time step length */
    double mDt;

    /** Percentage increment after which to report change */
    unsigned mPercentageIncrement;

    /** Time step increment after which to report change */
    unsigned mTimestepIncrement;

    /** Time increment (in seconds) after which to report change */
    unsigned mSecondsIncrement;

    /** The most recent percentage that was written */
    unsigned mLastPercentage;

    /** Duration of the simulation in number of time steps */
    unsigned mNumTimesteps;

    /** Max number of digits needed to represent the number of time steps */
    unsigned mNumDigits;

    /** The wall set in the main constructor */
    double mWallTimeAtStart;

    /** Whether to output to console */
    bool mOutputToConsole;

    /** Whether to output to file */
    bool mOutputToFile;

    /** The timer */
    Timer mTimer;

    /** The file to write to */
    out_stream mpFile;

    /**
     * Send message either to the file stream or the console (or both), depending on the values of
     * mOutputToConsole and mOutputToFile
     *
     * @param message the message to print
     */
    void SendMessage(std::string message);

    /**
     * Given a strictly non-negative elapsed time, return a string formatted as (dd:hh:mm:ss)
     *
     * @param timeElapsed a non-negative elapsed time
     * @return a string representing the elapsed time
     */
    std::string GetTimeString(double timeElapsed);

public:

    /**
     * Constructor saves times and opens output file ('progress_status.txt').
     *
     * @param outputDirectory where to open the output file
     * @param startTime the start time
     * @param endTime the end time
     */
    ProgressReporter(std::string outputDirectory, double startTime, double endTime, double dt);

    /** Default constructor */
    ProgressReporter() = default;

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

    /** Print finalising message to file */
    void PrintFinalising();

    /** Print initialising message to file */
    void PrintInitialising();

    /** Set whether to output to console */
    void SetOutputToConsole(bool outputToConsole=true);
};

#endif /*PROGRESSREPORTER_HPP_*/
