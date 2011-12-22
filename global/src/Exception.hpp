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


#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_

/**
 * @file
 * Contains the Exception class, along with some macros that are widely
 * used throughout the code.
 */
#include <string>
#include <sstream>
#include <cfloat>  // For DBL_MAX
#include <climits> // For UINT_MAX & INT_MAX, necessary in gcc-4.3
#include <cstdlib> // For system() etc., necessary in gcc-4.3

/** Use when initialising an unsigned variable that doesn't have a sensible default value. */
const unsigned UNSIGNED_UNSET = UINT_MAX;
/** Use when initialising an int variable that doesn't have a sensible default value. */
const int INT_UNSET = INT_MAX;
/** Use when initialising a double variable that doesn't have a sensible default value. */
const double DOUBLE_UNSET = DBL_MAX;

/**
 * Exception class.
 * All exceptions thrown by this code are currently instances of this class.
 *
 */
class Exception
{
public:
    /**
     * Construct an exception with a message string.
     *
     * @param rMessage  the message
     * @param rFilename  which source file threw the exception
     * @param lineNumber  which line number of the source file threw the exception
     */
    Exception(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber);

    /**
     * Get the message associated with the exception with file and line number
     *
     * @return The message set when the exception was thrown including file and line number information
     **/
    std::string GetMessage() const;

    /**
     * Get the message associated with the exception
     *
     * @return The message text set when the exception was thrown.
     **/
    std::string GetShortMessage() const;

    /**
     * Helper method for checking we have the right exception.
     *
     * Checks that #mShortMessage matches that given, and returns
     * a suitable error message string if not.  If they do match,
     * returns the empty string.
     *
     * @param expected  the expected value of #mShortMessage
     */
    std::string CheckShortMessage(std::string expected) const;

    /**
     * Helper method for checking we have the right exception.
     *
     * Checks that #mShortMessage contains the given string, and
     * returns a suitable error message string if not.  If it does,
     * returns the empty string.
     *
     * @param expected  some expected substring of #mShortMessage
     */
    std::string CheckShortMessageContains(std::string expected) const;

    /**
     * Level 4 error (Termination).  Execution cannot continue from this point and hence
     * should be terminated (even when running with NDEBUG or in parallel).
     *
     * @param rMessage An error message to appear on the screen
     * @param rFilename  which source file produced the termination error
     * @param lineNumber  which line number of the source file produced the termination error
     */
    static void Terminate(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber);

protected:
    /**
     * Allow subclasses to reset the exception message after construction of the base class,
     * if desired.
     *
     * @param rMessage  the message
     * @param rFilename  which source file threw the exception
     * @param lineNumber  which line number of the source file threw the exception
     */
    void SetMessage(const std::string& rMessage,
                    const std::string& rFilename, unsigned lineNumber);

private:
    std::string mMessage; /**< Full exception message - includes file and line number. */
    std::string mShortMessage; /**< Short exception message - just text of the exception. */
};

/**
 * Convenience macro for throwing an exception, in order to add file and line info.
 *
 * @param message  the error message to use, as a streamed expression
 */
#define EXCEPTION(message)                           \
    do {                                             \
        std::stringstream msg_stream;                \
        msg_stream << message;                       \
        throw Exception(msg_stream.str(), __FILE__, __LINE__); \
    } while (false)

#include <boost/preprocessor/stringize.hpp>

/**
 * Convenience macro for changing an assert into an exception - has the same
 * calling semantics, but throws.
 *
 * @param test  the test that must always be true.
 */
#define EXCEPT_IF_NOT(test) \
    if (!(test)) EXCEPTION("Assertion tripped: " BOOST_PP_STRINGIZE(test))


/**
 * Terminate execution safely, even when running in parallel.  Use for level 4 errors:
 * execution cannot continue from this point and hence should be terminated (even when running with NDEBUG).
 *
 * @param message  explanatory message
 */
#define TERMINATE(message) Exception::Terminate(message, __FILE__, __LINE__)

/**
 * Use for control paths that will never be executed, just to make sure they aren't.
 *
 * The exit statement at the end of NEVER_REACHED is not really needed but prevents g++ from complaining about
 * uninitialised variables when you have code that looks like:
 *
 * \code
 *   RelativeTo::Value relative_to;
 *   switch (rPath.relative_to())
 *   {
 *       case cp::relative_to_type::cwd:
 *           relative_to = RelativeTo::CWD;
 *           break;
 *       case cp::relative_to_type::chaste_test_output:
 *           relative_to = RelativeTo::ChasteTestOutput;
 *           break;
 *       case cp::relative_to_type::chaste_source_root:
 *           relative_to = RelativeTo::ChasteSourceRoot;
 *           break;
 *       case cp::relative_to_type::absolute:
 *           relative_to = RelativeTo::Absolute;
 *           break;
 *       default:
 *           NEVER_REACHED;
 *           break;
 *   }
 * \endcode
 *
 * relative_to is considered potentially uninitialised in the default branch unless the compiler finds a exit,
 * assert or throw statement.
 */
#define NEVER_REACHED  TERMINATE("Should have been impossible to reach this line of code"); exit(EXIT_FAILURE)

/**
 * This is to cope with NDEBUG causing variables to not be used, when they are only
 * used in assert()s.
 * @param var  the "unused" variable
 */
#ifdef NDEBUG
#define UNUSED_OPT(var) var=var
#else
#define UNUSED_OPT(var)
#endif

/**
 * Handy for calling functions like system which return non-zero on error.
 * Throws if an error occurs.
 *
 * @note DO NOT use this macro within an if (PetscTools::AmMaster) block, as then you'll
 * get deadlock if an exception is thrown when running in parallel!
 * (Unless the block is wrapped in a try-catch and exception replication handler.)
 * Instead, use ABORT_IF_NON0.
 *
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define EXPECT0(cmd, arg) {      \
    std::string _arg(arg);       \
    int ret = cmd(_arg.c_str()); \
    if (ret != 0) {              \
        EXCEPTION("Error executing command: " #cmd "(" + _arg + ")"); \
    } }


/**
 * Handy for calling functions like system which return non-zero on error.
 * Terminate if the return code is non-zero, printing a suitable message.
 * @param retcode  command return code
 * @param msg  error message to display
 */
#define ABORT_IF_NON0_WITH_MSG(retcode, msg)     \
    if (retcode != 0) {                          \
        TERMINATE(msg);                          \
    }

/**
 * Handy for calling functions like system which return non-zero on error.
 * Terminate if an error occurs.
 *
 * This macro should be used instead of EXPECT0 within blocks that are only
 * executed by one process, but need to kill all processes if an error occurs.
 *
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define ABORT_IF_NON0(cmd, arg) { \
    std::string _arg(arg);        \
    int ret = cmd(_arg.c_str());  \
    ABORT_IF_NON0_WITH_MSG(ret, "Error executing command: " #cmd "(" + _arg + ")") \
    }

/**
 * Handy for calling functions like system which return non-zero on error.
 * This time we expect failure; throws if the command succeeds.
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define EXPECT_NON0(cmd, arg) {   \
    std::string _arg = (arg);     \
    int ret = cmd(_arg.c_str());  \
    if (ret == 0) {               \
        EXCEPTION("Command: " #cmd "(" + _arg + ") succeeded and it shouldn't have"); \
    } }

/**
 * Handy for calling functions like system which return non-zero on error.
 * This version ignores the return code, in case you don't care about errors for some reason...
 * @param cmd  command to call
 * @param arg  its argument (will be converted to std::string)
 */
#define IGNORE_RET(cmd, arg) {   \
    std::string _arg = (arg);    \
    int ret = cmd(_arg.c_str()); \
    ret = ret;                   \
    }

#endif // _EXCEPTION_HPP_
