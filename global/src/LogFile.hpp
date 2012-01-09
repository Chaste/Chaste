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

#ifndef LOGFILE_HPP_
#define LOGFILE_HPP_

#include <string>
#include "OutputFileHandler.hpp"
#include <iostream>
#include <cassert>
#include <iomanip>
#include <ctime>

/**
 * A singleton log file class. Allows the user to define log file in the test, which
 * can be written to (without being passed around in the code).
 *
 * Usage (in test):
 *
 * // begining of test
 * LogFile* p_log = LogFile::Instance();
 * p_log->Set(level, "dir","file");
 * p_log->WriteHeader("type_of_sim"); // optional
 *
 * // at end of simulation
 * LogFile::Close();
 *
 * Here 'level' is a number between 0 and LogFile::MaxLoggingLevel, with zero
 * meaning no logging and MaxLoggingLevel meaning full logging.
 *
 * Usage (in source) - use the macro 'LOG'
 * LOG(1, "Info to be written to the log file\n" << "More info\n");
 * LogFile::Instance()->WriteElapsedTime(); // optional
 *
 * This checks whether the given level (here '1') is greater or equal to the given
 * logging level, in which case it writes to the current log file. If there is
 * no log file set up it does nothing.
 *
 * Note the log file can be written to directly, without any level-checking, using
 * (*LogFile::Instance()) << "Info to be written to the log file\n";
 */
class LogFile
{
private:

    /** The static single instance. */
    static LogFile* mpInstance;

    /** Whether a directory and filename has been set. */
    bool mFileSet;

    /** The file to be written to. */
    out_stream mpOutStream;

    /** The current current calendar time. */
    time_t mInitTime;

    /** The level of logging required for this particular log file. */
    unsigned mLevel;

    /** The max level of logging. */
    static const unsigned mMaxLoggingLevel = 2;

    /** The precision with which to output data. */
    unsigned mPrecision;

    /**
     * Constructor. Should never be called directly, call LogFile::Instance() instead.
     */
    LogFile();

public:

    /**
     * Get the single instance of the LogFile object.
     */
    static LogFile* Instance();

    /**
     * Get the logging level.
     */
    static unsigned Level();

    /**
     * Set the logging level, the directory (relative to TEST_OUTPUT) and the file
     * the log should be written to (file defaults to "log.txt").
     *
     * The level should be a number between 0 and LogFile::MaxLoggingLevel() (which is the
     * same as LogFile::mMaxLoggingLevel)
     *
     * Note: we intentionally do NOT check or throw an exception if a file has already
     * been set (i.e. Close() wasn't called the last time a log was used).
     *
     * The directory is never cleaned.
     *
     * @param level  the logging level
     * @param rDirectory  the directory in which to write the data to file
     * @param rFileName  the name of the file to write to, relative to the output directory
     */
    void Set(unsigned level, const std::string& rDirectory, const std::string& rFileName="log.txt");

    /**
     * Get the maximum allowed logging level.
     */
    static unsigned MaxLoggingLevel();

    /**
     * Set the precision to write data (the 'decimal precision', look up
     * documentation for std::setprecision()).
     *
     * @param precision  the precision
     */
    void SetPrecision(unsigned precision);

    /**
     * Write a header in the log file, stating the (given) type of simulation and the
     * date and time.
     *
     * @param simulationType The type of simulation, eg "Bidomain" or "Crypt" or
     * "Cardiac Electromechanics". Defaults to empty.
     */
    void WriteHeader(std::string simulationType="");

    /**
     * Write the elapsed time since the simulation began (since the log file was created).
     *
     * @param pre a string (eg spacings) to write before the elapsed time line.
     */
    void WriteElapsedTime(std::string pre="");

    /**
     * Close the currently open file, and delete the single LogFile instance.
     */
    static void Close();

    /**
     * Whether Set() has been called.
     */
    bool IsFileSet();

    /**
     * Overloaded << operator, to write to the log file, if one has been set, and
     * does nothing if not.
     *
     * @param message the message to write to the log file
     */
    template <class T>
    LogFile& operator<<(T message)
    {
        if (mFileSet)
        {
            (*mpOutStream) << std::setprecision((int)mPrecision) << message << std::flush;
        }

        return *this;
    }
};

// Define the log macro
#define LOG(level, message) assert(level>0); if(level <= LogFile::Level()) { (*LogFile::Instance()) << message << "\n"; }
#define LOG_AND_COUT(level, message) {std::cout << message << std::endl << std::flush; LOG(level, message); }

#endif /*LOGFILE_HPP_*/
