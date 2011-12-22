/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include <ctime>
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
    static time_t StartTime;

public:

    /**
     * Reset the timer.
     */
    static void Reset();

    /**
     * Print the elapsed time (to std::cout and the Log file (under logging-level 2)
     * preceded by the message provided.
     *
     * @param message
     */
    static void Print(std::string message);

    /**
     * Print the elapsed time (to std::cout and the Log file (under logging-level 2)
     * preceded by the message provided, and also reset the timer.
     *
     * @param message
     */
    static void PrintAndReset(std::string message);
};

#endif /*TIMER_HPP_*/
