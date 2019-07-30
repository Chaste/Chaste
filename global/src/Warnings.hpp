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

#ifndef _WARNINGS_HPP_
#define _WARNINGS_HPP_

#include <string>
#include <deque>
#include <sstream>

/**
 * The Warnings singleton class collects warnings via the AddWarning() method
 * or the WARNING macro.  This is to provide a mechanism for informing the user of
 * error which are less severe than those demanding exceptions or assertions. (Errors
 * which can be repaired immediately.)
 *
 * Warnings can be polled with GetNumWarnings() and GetNextWarningMessage().
 * Warnings can be ignored and destroyed with QuietDestroy().
 * All warnings left at the close of the test suite (not the individual test), or
 * the close of the executable, will be printed to the screen.
 */
class Warnings
{
private:

    static Warnings* mpInstance; /**<  Pointer to the single instance. For use as singleton */

    /** Container type for warnings */
    typedef std::deque<std::pair<std::string, std::string> > WarningsContainerType;

    WarningsContainerType mWarningMessages; /**< Warnings messages.  First in pair is the context (line number etc.).  Second in pair is the actual warning, */

protected:

    /**
     * Protected constructor.
     * Use Instance() to access the Warnings singleton.
     */
    Warnings();

public:

    /**
     * Default destroyer (takes the place of a destructor, using the std::atexit directive).
     * Prints all warning messages to stdout.
     */
    static void NoisyDestroy();

    /**
     * Public destroyer (takes the place of a destructor but is called by the user to make sure
     * that no warnings appear).
     */
    static void QuietDestroy();

    /**
     * Prints all warning messages to stdout without changing the state of the Warnings.
     */
    static void PrintWarnings();

    /**
     * Make a warning with a message string and add it to the list.
     *
     * @param rMessage  the message
     * @param rFilename  which source file threw the exception
     * @param lineNumber  which line number of the source file threw the exception
     * @param onlyOnce  whether to only log the first warning thrown from this location
     */
    void AddWarning(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber, bool onlyOnce=false);

    /**
     * Get the message associated with the exception with file and line number.
     *
     * @return The message set when the exception was thrown including file and line number information
     */
    std::string PopWarning() const;

    /**
     * @return a pointer to the Warnings object.
     * The object is created the first time this method is called.
     */
    static Warnings* Instance();

    /**
     * @return how many warnings are in the queue.
     */
    unsigned GetNumWarnings();

    /**
     * @return next warning.  Remove and inspect a warning.
     */
    std::string GetNextWarningMessage();
};

#define WARNING(message)                                                    \
{                                                                           \
    std::stringstream msg_stream;                                           \
    msg_stream << message;                                                  \
    Warnings::Instance()->AddWarning(msg_stream.str(), __FILE__, __LINE__); \
}

/**
 * Warn only the first time line is reached. Note: this does not base whether to not warn
 * again on the message content, just on whether the line of code where this macro is placed
 * has been reached. In other words:
 *
 * for (unsigned i=0; i<10; i++)
 * {
 *    WARN_ONCE_ONLY("Don't make me angry");
 * }
 *
 * will print once, whereas
 *
 * WARN_ONCE_ONLY("You wouldn't like me when I'm angry");
 * WARN_ONCE_ONLY("You wouldn't like me when I'm angry");
 *
 * will warn twice.
 */
#define WARN_ONCE_ONLY(message)                                                     \
{                                                                                   \
    std::stringstream msg_stream;                                                   \
    msg_stream << message;                                                          \
    Warnings::Instance()->AddWarning(msg_stream.str(), __FILE__, __LINE__, true);   \
}
#endif // _WARNINGS_HPP_
