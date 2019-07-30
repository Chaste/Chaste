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

#ifndef GENERICEVENTHANDLER_HPP_
#define GENERICEVENTHANDLER_HPP_

#include <cassert>
#include <iostream>

#include "Exception.hpp"
#include "PetscTools.hpp"
#include "Timer.hpp"

/**
 * A generic base class providing the functionality for timing various events.
 * Subclasses provide the event codes and names; see HeartEventHandler for an example.
 *
 * Note: this class assume that, for any given concrete class, the last event
 * represents the total time, and thus wraps all other events.
 *
 * The methods in this class are not implemented separately as then they would not be
 * inline, which could impact performance; we generally want timing routines to be very
 * lightweight.
 */
template <unsigned NUM_EVENTS, class CONCRETE>
class GenericEventHandler
{
    friend class TestGenericEventHandler;
    friend class TestCellBasedEventHandler;
    friend class TestHeartEventHandler;

private:

    std::vector<double> mWallTime; /**< Wall time assigned to each event */
    std::vector<bool> mHasBegun; /**< Whether each event is in progress */
    bool mEnabled; /**< Whether the event handler is recording event times */
    bool mInUse; /**< Determines if any of the event have begun */

    /**
     * Sleep for a specified number of milliseconds.
     * Used in testing.
     * Ought to be more portable than sleep() or usleep().
     *
     * @param milliseconds  minimum number of milliseconds for which to sleep (ought to be a multiple of 10)
     */
    static inline void MilliSleep(unsigned milliseconds)
    {
        double min_Wtime = milliseconds/1000.0 + Timer::GetElapsedTime();
        while (Timer::GetElapsedTime() < min_Wtime)
        {
            //pause;
        }
    }

    /**
     * Convert a wall clock time to milliseconds.
     *
     * @param wallTime  the wall time
     * @return time in milliseconds.
     */
    inline double ConvertWallTimeToMilliseconds(double wallTime)
    {
        return wallTime*1000.0;
    }

    /**
     * Convert a wall clock time to seconds.  (Identity function in this implementation.)
     *
     * @param wallTime  the wall time
     * @return Time in seconds
     */
    inline double ConvertWallTimeToSeconds(double wallTime)
    {
        return wallTime;
    }

public:

    /**
     * @return Get the singleton instance of the event handler.
     */
    static GenericEventHandler<NUM_EVENTS, CONCRETE>* Instance()
    {
        static CONCRETE inst;
        return &inst;
    }

    /**
     * Reset the event handler - set all event durations to zero.
     */
    static void Reset()
    {
        Instance()->ResetImpl();
    }

    /**
     * Record the start of an event.
     *
     * @param event  the index of an event (this must be less than NUM_EVENTS)
     */
    static void BeginEvent(unsigned event)
    {
        Instance()->BeginEventImpl(event);
    }

    /**
     * Record the ending of an event.
     *
     * @param event  the index of an event (this must be less than NUM_EVENTS)
     */
    static void EndEvent(unsigned event)
    {
        Instance()->EndEventImpl(event);
    }

    /**
     * @return The time (in milliseconds) accounted so far to the given event.
     *
     * Will automatically determine if the event is currently ongoing or not.
     *
     * @param event  the index of an event (this must be less than NUM_EVENTS)
     */
    static double GetElapsedTime(unsigned event)
    {
        return Instance()->GetElapsedTimeImpl(event);
    }

    /**
     * Print a report on the timed events and reset the handler.
     *
     * Assumes all events have ended.
     *
     * If there is a collection of processes then the report will include an
     * average and maximum over all CPUs.
     */
    static void Report()
    {
        Instance()->ReportImpl();
    }

    /**
     * Output the headings for a report.
     */
    static void Headings()
    {
        Instance()->HeadingsImpl();
    }

    /**
     * Enable the event handler so that it will record event durations.
     */
    static void Enable()
    {
        Instance()->EnableImpl();
    }

    /** Disable the event handler, so that event durations are no longer recorded. */
    static void Disable()
    {
        Instance()->DisableImpl();
    }

    /** Check whether the event handler is enabled.
     * @return true if enabled
     */
    static bool IsEnabled()
    {
        return Instance()->IsEnabledImpl();
    }

protected:

    /**
     * Default constructor.
     */
    GenericEventHandler()
    {
        mEnabled = true;
        mInUse = false;
        mWallTime.resize(NUM_EVENTS, 0.0);
        mHasBegun.resize(NUM_EVENTS, false);
    }

private:

    /**
     * Reset the event handler - set all event durations to zero.
     */
    void ResetImpl()
    {
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            mWallTime[event] = 0.0;
            mHasBegun[event] = false;
        }
        Enable();
        mInUse = false;
    }

    /**
     * Record the start of an event.
     *
     * @param event  the index of an event (this must be less than NUM_EVENTS)
     */
    void BeginEventImpl(unsigned event)
    {
        if (!mEnabled)
        {
            return;
        }
#ifdef CHASTE_EVENT_BARRIERS
        PetscTools::Barrier("BeginEvent");
#endif
        mInUse = true;
        assert(event<NUM_EVENTS);
        // Check that we are recording the total
        if (event != NUM_EVENTS-1) // If use <, Intel complains when NUM_EVENTS==1
        {
            if (!mHasBegun[NUM_EVENTS-1])
            {
                // Silently open the "total" event
                BeginEvent(NUM_EVENTS-1);
            }
        }
        if (mHasBegun[event])
        {
            std::string msg;
            msg += "The event associated with the counter for '";
            msg += CONCRETE::EventName[event];
            msg += "' had already begun when BeginEvent was called.";
            std::cerr << msg << std::endl << std::flush;
            Disable();
            return;
        }
        mWallTime[event] -= Timer::GetWallTime();
        mHasBegun[event] = true;
        //std::cout << PetscTools::GetMyRank()<<": Beginning " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }

    /**
     * Record the ending of an event.
     *
     * @param event  the index of an event (this must be less than NUM_EVENTS)
     */
    void EndEventImpl(unsigned event)
    {
        assert(event<NUM_EVENTS);
        if (!mEnabled)
        {
            return;
        }
#ifdef CHASTE_EVENT_BARRIERS
        PetscTools::Barrier("EndEvent");
#endif
        if (!mHasBegun[event])
        {
            std::string msg;
            msg += "Error: The event associated with the counter for '";
            msg += CONCRETE::EventName[event];
            msg += "' had not begun when EndEvent was called.";
            EXCEPTION(msg);
        }
        mWallTime[event] += Timer::GetWallTime();
        mHasBegun[event] = false;
        //std::cout << PetscTools::GetMyRank()<<": Ending " << EVENT_NAME[event] << " @ " << (clock()/1000) << std::endl;
    }

    /**
     * @return The time (in milliseconds) accounted so far to the given event.
     *
     * Will automatically determine if the event is currently ongoing or not.
     *
     * @param event  the index of an event (this must be less than NUM_EVENTS)
     */
    double GetElapsedTimeImpl(unsigned event)
    {
        assert(event<NUM_EVENTS);
        if (!mEnabled)
        {
            EXCEPTION("Asked to report on a disabled event handler.  Check for contributory errors above.");
        }
        double time;
        if (mHasBegun[event])
        {
            time =  mWallTime[event] + Timer::GetWallTime();
        }
        else
        {
            time = mWallTime[event];
        }
        return ConvertWallTimeToMilliseconds(time);
    }

    /**
     * Print a report on the timed events and reset the handler.
     *
     * Assumes all events have ended.
     *
     * If there is a collection of processes then the report will include an
     * average and maximum over all CPUs.
     */
    void ReportImpl()
    {
        if (!mEnabled)
        {
            EXCEPTION("Asked to report on a disabled event handler.  Check for contributory errors above.");
        }
        if (!mInUse)
        {
            EXCEPTION("Asked to report on an event handler which is set to zero.");
        }
        // Check that all events are finished
        for (unsigned event=0; event<NUM_EVENTS; event++)
        {
            if (mHasBegun[event])
            {
                // Silently close event
                EndEvent(event);
            }
        }
        const unsigned top_event = NUM_EVENTS-1;
        double total = ConvertWallTimeToSeconds(mWallTime[top_event]);

        // Make the output precision depend on total run time
        const char* format;
        if (total > 999999.0)      // 11.5 days
        {
            format = "%8.0f ";     // will allow up to 115 days before columns unaligned
        }
        else if (total > 9999.0)   // 2.7 hours
        {
            format = "%8.1f ";
        }
        else
        {
            format = "%8.3f ";
        }

        PetscTools::BeginRoundRobin();
        {
            std::cout.flush();
            if (PetscTools::IsParallel())
            {
                // Report the process number at the beginning of the line
                printf("%3u: ", PetscTools::GetMyRank()); //5 chars
            }
            for (unsigned event=0; event<NUM_EVENTS; event++)
            {
                const double secs = ConvertWallTimeToSeconds(mWallTime[event]);
                printf(format, secs);
                printf("(%3.0f%%)  ", total == 0.0 ? 0.0 : (secs/total*100.0));
            }
            std::cout << "(seconds) \n";
        }
        PetscTools::EndRoundRobin();

        // If there is a collection of processes then report an average
        if (PetscTools::IsParallel() && !PetscTools::IsIsolated())
        {
            double total_cpu_time[NUM_EVENTS];
            MPI_Reduce(&mWallTime[0], total_cpu_time, NUM_EVENTS, MPI_DOUBLE, MPI_SUM, 0, PetscTools::GetWorld());
            if (PetscTools::AmMaster())
            {
                total = ConvertWallTimeToSeconds(total_cpu_time[top_event]);
                printf("avg: "); //5 chars
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    const double secs = ConvertWallTimeToSeconds(total_cpu_time[event]);
                    printf(format, secs/PetscTools::GetNumProcs());
                    printf("(%3.0f%%)  ", total == 0.0 ? 0.0 : (secs/total*100.0));
                }
                std::cout << "(seconds) \n";
            }

            double max_cpu_time[NUM_EVENTS];
            MPI_Reduce(&mWallTime[0], max_cpu_time, NUM_EVENTS, MPI_DOUBLE, MPI_MAX, 0, PetscTools::GetWorld());
            if (PetscTools::AmMaster())
            {
                total = ConvertWallTimeToSeconds(max_cpu_time[top_event]);
                printf("max: "); //5 chars
                for (unsigned event=0; event<NUM_EVENTS; event++)
                {
                    const double secs = ConvertWallTimeToSeconds(max_cpu_time[event]);
                    printf(format, secs);
                    printf("(%3.0f%%)  ", total == 0.0 ? 0.0 : (secs/total*100.0));
                }
                std::cout << "(seconds) \n";
            }
        }
        std::cout.flush();
        PetscTools::Barrier();
        std::cout.flush();

        Reset();
    }

    /**
     * Output the headings for a report.
     */
    void HeadingsImpl()
    {
        // Make sure that all output (on all processes) is flushed
        std::cout.flush();
        PetscTools::Barrier();
        std::cout.flush();
        if (PetscTools::AmMaster())
        {
            if (PetscTools::IsParallel())
            {
                // Report the process number at the beginning of the line
                printf("Proc "); //5 chars
            }
            for (unsigned event=0; event<NUM_EVENTS; event++)
            {
                printf("%15s%2s", CONCRETE::EventName[event], "");
            }
            std::cout << "\n";
            std::cout.flush();
        }
    }

    /** Enable the event handler so that it will record event durations. */
    void EnableImpl()
    {
        mEnabled = true;
    }

    /** Disable the event handler, so that event durations are no longer recorded. */
    void DisableImpl()
    {
        mEnabled = false;
    }

    /** Check whether the event handler is enabled.
     * @return true if enabled
     */
    bool IsEnabledImpl()
    {
        return mEnabled;
    }
};

#endif /*GENERICEVENTHANDLER_HPP_*/
