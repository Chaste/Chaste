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

#include "ProgressReporter.hpp"

#include <cassert>
#include <cmath>
#include <iomanip>
#include "Exception.hpp"
#include "PetscTools.hpp"

ProgressReporter::ProgressReporter(std::string outputDirectory, double startTime, double endTime, double dt)
        : mStartTime(startTime),
          mEndTime(endTime),
          mDt(dt),
          mPercentageIncrement(5),
          mTimestepIncrement(UINT_MAX),
          mSecondsIncrement(UINT_MAX),
          mLastPercentage(UINT_MAX),
          mOutputToConsole(false),
          mOutputToFile(true)
{
    if (PetscTools::AmMaster())
    {
        assert(startTime < endTime);

        if (outputDirectory.empty())
        {
            mOutputToFile = false;
        }
        else
        {
            OutputFileHandler handler(outputDirectory, false);
            mpFile = handler.OpenOutputFile("progress_status.txt");
        }

        mTimer.Reset();
        mWallTimeAtStart = mTimer.GetWallTime();

        mNumTimesteps = static_cast<unsigned>(std::lround((mEndTime - mStartTime) / dt));
        mNumDigits = 1u + static_cast<unsigned>(std::floor(std::log10(mNumTimesteps)));
    }
}

ProgressReporter::~ProgressReporter()
{
    if (PetscTools::AmMaster())
    {
        std::stringstream message;

        if (mLastPercentage != 100)
        {
            message << "100% completed" << std::endl;
        }

        message << "..done!" << std::endl;

        SendMessage(message.str());

        if (mOutputToFile)
        {
            mpFile->close();
        }

        PrintFinalising();
    }
}

void ProgressReporter::Update(double currentTime)
{
    if (PetscTools::AmMaster())
    {
        if (mOutputToFile || mOutputToConsole)
        {
            auto percentage = std::lround((currentTime - mStartTime) / (mEndTime - mStartTime) * 100);
            auto timesteps_elapsed = std::lround((currentTime - mStartTime) / mDt);
            auto seconds_elapsed = std::lround(mTimer.GetElapsedTime());

            bool percentage_condition = mLastPercentage == UINT_MAX || percentage - mLastPercentage >= mPercentageIncrement;
            bool timesteps_condition = timesteps_elapsed % mTimestepIncrement == 0;
            bool elapsed_time_condition = seconds_elapsed >= mSecondsIncrement;

            if (percentage_condition || timesteps_condition || elapsed_time_condition)
            {
                std::stringstream message;
                message << std::setfill(' ') << std::setw(3) << percentage << "% complete: "
                        << std::setw(mNumDigits) << timesteps_elapsed << "/" << mNumTimesteps
                        << " " << GetTimeString(mTimer.GetWallTime() - mWallTimeAtStart) << std::endl;

                SendMessage(message.str());

                mLastPercentage = static_cast<unsigned>(percentage);

                if (elapsed_time_condition)
                {
                    mTimer.Reset();
                }
            }
        }
    }
}

void ProgressReporter::PrintFinalising()
{
    if (PetscTools::AmMaster())
    {
        std::stringstream message;
        message << "Finalising.." << std::endl;

        SendMessage(message.str());
    }
}

void ProgressReporter::PrintInitialising()
{
    if (PetscTools::AmMaster())
    {
        std::stringstream message;
        message << "Initialising.." << std::endl;

        SendMessage(message.str());
    }
}

void ProgressReporter::SetOutputToConsole(bool outputToConsole)
{
    mOutputToConsole = outputToConsole;
}

void ProgressReporter::SendMessage(std::string message)
{
    if (mOutputToFile)
    {
        *mpFile << message;
    }
    if (mOutputToConsole)
    {
        std::cout << message;
    }
}
std::string ProgressReporter::GetTimeString(double timeElapsed)
{
    if (timeElapsed < 0.0)
    {
        EXCEPTION("Invalid time elapsed.");
    }
    if (timeElapsed > static_cast<double>(LONG_MAX))
    {
        EXCEPTION("Long overflow: something has gone wrong.");
    }

    // Get long int representing the elaspsed time
    auto elapsed = std::lround(timeElapsed);

    // Convert this into seconds, minutes, hours, and days
    auto sec = elapsed % 60;

    elapsed = (elapsed - sec) / 60;
    auto min = elapsed % 60;

    elapsed = (elapsed - min) / 60;
    auto hrs = elapsed % 24;

    elapsed = elapsed - hrs;
    auto day = elapsed / 24;

    // Represent this information as a string
    std::stringstream time_string;
    time_string << "(" << std::setfill('0')
                << std::setw(2) << day << ":"
                << std::setw(2) << hrs << ":"
                << std::setw(2) << min << ":"
                << std::setw(2) << sec << "s)";

    return time_string.str();
}
