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
     * Return a pointer to the Warnings object.
     * The object is created the first time this method is called.
     */
    static Warnings* Instance();

    /**
     * How many warnings are in the queue.
     */
    unsigned GetNumWarnings();

    /**
     * Remove and inspect a warning.
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
