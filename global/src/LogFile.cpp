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

#include "LogFile.hpp"
#include "Exception.hpp"

#include <cmath>
#include <sstream>

LogFile* LogFile::mpInstance = nullptr;

LogFile::LogFile()
    : mFileSet(false),
      mInitTime(time(nullptr)),
      mLevel(0),
      mPrecision(6)
{
}

LogFile* LogFile::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new LogFile; // default construtor which doesn't write
    }
    return mpInstance;
}

unsigned LogFile::Level()
{
    if (mpInstance == nullptr)
    {
        return 0;
    }
    else
    {
        return mpInstance->mLevel;
    }
}

void LogFile::Set(unsigned level, const std::string& rDirectory, const std::string& rFileName)
{
    if (level > mMaxLoggingLevel)
    {
        EXCEPTION("Requested level " << level
                  << " should have been less than or equal to "
                  << mMaxLoggingLevel);
    }
    mLevel = level;

//// force log files to be written to the desktop (eg for use on machines which tend to die
//// and need rebooting (/tmp gets wiped on a reboot))
//    std::string file = "/home/chaste/Desktop/" + fileName;
//    out_stream p_file(new std::ofstream(file.c_str()));
//    mpOutStream = p_file;

    OutputFileHandler handler(rDirectory, false);
    mpOutStream = handler.OpenOutputFile(rFileName);
    mFileSet = true;

    // Write header in the log file..?
}

unsigned LogFile::MaxLoggingLevel()
{
    return mMaxLoggingLevel;
}

void LogFile::Close()
{
    if (mpInstance)
    {
        mpInstance->mpOutStream->close();
        delete mpInstance;
        mpInstance = nullptr;
    }
}

void LogFile::WriteHeader(std::string simulationType)
{
    *this << "\nChaste: " << simulationType << " simulation, on " << ctime(&mInitTime) << "\n";
}

void LogFile::WriteElapsedTime(std::string pre)
{
    double fsecs = difftime(time(nullptr),mInitTime);
    long total_secs = static_cast<long>(floor(fsecs+0.5));
    int total_mins = total_secs/60;

    int secs = total_secs%60;
    int mins = total_mins%60;
    int hrs = total_mins/60;

    *this << pre << "Elapsed time is: " <<  hrs << "h " << mins << "m " << secs << "s\n";
}

bool LogFile::IsFileSet()
{
    return mFileSet;
}

void LogFile::SetPrecision(unsigned precision)
{
    assert(precision > 0);
    mPrecision = precision;
}
