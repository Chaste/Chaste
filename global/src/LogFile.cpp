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

#include "LogFile.hpp"
#include "Exception.hpp"

#include <cmath>
#include <sstream>

LogFile* LogFile::mpInstance = NULL;

LogFile::LogFile()
    : mFileSet(false),
      mInitTime(time(NULL)),
      mLevel(0),
      mPrecision(6)
{
}

LogFile* LogFile::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new LogFile; // default construtor which doesn't write
    }
    return mpInstance;
}

unsigned LogFile::Level()
{
    if (mpInstance == NULL)
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
        mpInstance = NULL;
    }
}

void LogFile::WriteHeader(std::string simulationType)
{
    *this << "\nChaste: " << simulationType << " simulation, on " << ctime(&mInitTime) << "\n";
}

void LogFile::WriteElapsedTime(std::string pre)
{
    double fsecs = difftime(time(NULL),mInitTime);
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
