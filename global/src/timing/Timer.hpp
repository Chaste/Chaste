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

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <string>

/**
 *  A very simple lightweight benchmarking tool. Call Timer::Reset() to reset the timer
 *  and Timer::Print() to print the time elapsed to stdout.
 *
 *  Usage:
 *
 *  Timer::Reset();
 *  //do something
 *  Timer::PrintAndReset("First thing");
 *  //do something else
 *  Timer::Print("Other thing");
 *
 *  which outputs (for example):
 *
 *  First thing time: 10s
 *  Other thing time: 2s
 */
class Timer
{
private:

    /** The start time. */
    static double msStartTime;

public:

    /**
     * Reset the timer.
     */
    static void Reset();

    /**
     * Print the elapsed wall-clock time (to std::cout and the Log file (under logging-level 2))
     * preceded by the message provided.
     *
     * @param message
     */
    static void Print(std::string message);

    /**
     * Get the elapsed wall-clock time since since midnight of January 1, 1970,
     * or since Reset() was called.
     *
     * @return The elapsed time, in seconds.
     */
    static double GetElapsedTime();

    /**
     * @return The time since since midnight of January 1, 1970 in seconds.
     */
    static double GetWallTime();

    /**
     * Print the elapsed wall-clock time (to std::cout and the Log file (under logging-level 2))
     * preceded by the message provided, and also reset the timer.
     *
     * @param message
     */
    static void PrintAndReset(std::string message);
};

#endif /*TIMER_HPP_*/
