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
#include <cstdlib> // Necessary in gcc-4.3 and later which don't include stdlib by default

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
     * @return an empty string when the expected message matches.
     * Checks that #mShortMessage matches that given, and
     * a suitable error message string if not.
     *
     * @param expected  the expected value of #mShortMessage
     */
    std::string CheckShortMessage(std::string expected) const;

    /**
     * Helper method for checking we have the right exception.
     *
     * @return an empty string when the message contains the expected string.
     * Checks that #mShortMessage contains the given string, and
     * returns a suitable error message string if not.
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
#define TERMINATE(message)                           \
    do {                                             \
        std::stringstream msg_stream;                \
        msg_stream << message;                       \
        Exception::Terminate(msg_stream.str(), __FILE__, __LINE__); \
    } while (false)

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
 * Convenience function to convert an exception thrown by a single process into
 * termination of the entire program.
 *
 * @param block  the block of code to execute
 */
#define ABORT_IF_THROWS(block)          \
    try {                               \
        block;                          \
    } catch (const Exception& e) {      \
        TERMINATE(e.GetMessage());      \
    } catch (const std::exception &e) { \
        TERMINATE(e.what());            \
    } catch (...) {                     \
        TERMINATE("Unexpected exception thrown."); \
    }


// The macros below are deprecated.  In most cases, FileFinder routines should be used
// in preference to system() calls.

/**
 * @note This macro is deprecated.
 *
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
    int retcode = cmd(_arg.c_str()); \
    if (retcode != 0) {              \
        EXCEPTION("Error executing command: " #cmd "(" + _arg + ")"); \
    } }


/**
 * @note This macro is deprecated, and no longer used by core code.
 *
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
 * @note This macro is deprecated, and no longer used by core code.
 *
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
 * @note This macro is deprecated, and no longer used by core code.
 *
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
 * @note This macro is deprecated.
 *
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
