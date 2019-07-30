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

#ifndef COMMANDLINEARGUMENTS_HPP_
#define COMMANDLINEARGUMENTS_HPP_

#include "Exception.hpp"
#include <vector>
#include <boost/utility.hpp>

/**
 * A convenient holder for the command line arguments, with helper
 * methods for checking whether an option has been given or getting the
 * value(s) corresponding to a given option.
 *
 * Options are considered to be any argument starting with a '-' followed by a non-numeric character,
 * in order to distinguish between options and negative numbers.  Each option may have zero or more
 * values following it on the command line.  Examples:
 *
 *  \li -option1 on
 *  \li --my_param value
 *  \li --verbose
 *  \li -numbers 1 8 -4.5
 *  \li --models model1 model2
 *
 * The cxxtest harness will fill in the member variables when a test is
 * started.  They can then be read by PETSc when it is initialised, and
 * by other code as required.
 */
class CommandLineArguments : private boost::noncopyable
{
private:

    /** Default constructor. Should never be called directly, call CommandLineArguments::Instance() instead.*/
    CommandLineArguments();

    /** The single instance of the class. */
    static CommandLineArguments* mpInstance;

    /**
     * @return the index for the given argument. Returns -1 if the argument is not found.
     *
     * @param rOption  the argument name as a string.
     * @return the position of the argument in the list (indexed from 1)
     */
    int GetIndexForArgument(std::string rOption);

    /**
     * Throw an exception if the option is not of the required form: '-' followed by a non-numeric character.
     *
     * @param rOption  the option name to check the format of.
     */
    void TestOptionFormat(const std::string& rOption);

public:

    /** The number of command line arguments. */
    int* p_argc;

    /** The arguments themselves. */
    char*** p_argv;

    /** @return the single instance of this class. */
    static CommandLineArguments* Instance();

    /**
     * @return true if a given option exists in the command line arguments.
     *
     * @param rOption  the option name as a string.
     */
    bool OptionExists(const std::string& rOption);

    /**
     * @return the number of arguments for a given option. Throws an Exception if the
     * option does not have any following arguments, or the option does not exist.
     * So you can use OptionExists() to avoid Exceptions in your code when using the other public methods.
     *
     * @param rOption  the option name as a string.
     * @param throwIfNone  Whether to give an exception if there are no arguments for this rOption (defaults to false).
     * @return the number of arguments following this option.
     */
    int GetNumberOfArgumentsForOption(const std::string& rOption, bool throwIfNone=false);

    /**
     * @return the value for a given option, i.e. the argument after the option name in
     * the list of command line arguments. For example, if the following arguments
     * were given
     *  ./heart/build/debug/TestMyClassRunner -timestep 0.04
     * Then calling
     *   CommandLineArguments::Instance()->GetValueCorrespondingToOption("-timestep");
     * will return 0.04 (as a char*).
     * Use atoi or atof to convert the char* to an int or a double(float) respectively,
     * or call one of our convenience methods for common types.
     *
     * @param rOption  the option name as a string.
     * @param valueNumber  the number of the argument following the option definition (defaults to 1, for 1st argument).
     */
    char* GetValueCorrespondingToOption(const std::string& rOption, int valueNumber=1);

    /**
     * @return the double for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char* to a double.
     *
     * @param rOption  the option name as a string.
     * @param valueNumber  the number of the argument following the option definition (defaults to 1, for 1st argument).
     */
    double GetDoubleCorrespondingToOption(const std::string& rOption, int valueNumber=1);

    /**
     * @return the int for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char* to an int.
     *
     * @param rOption  the option name as a string.
     * @param valueNumber  the number of the argument following the option definition (defaults to 1, for 1st argument).
     */
    int GetIntCorrespondingToOption(const std::string& rOption, int valueNumber=1);

    /**
     * @return the unsigned for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char* to an unsigned.
     * Throws an exception if the option converts to a negative integer.
     *
     * @param rOption  the option name as a string.
     * @param valueNumber  the number of the argument following the option definition (defaults to 1, for 1st argument).
     */
    unsigned GetUnsignedCorrespondingToOption(const std::string& rOption, int valueNumber=1);

    /**
     * @return the string for a given option.
     *
     * This uses GetValueCorrespondingToOption and converts the char* to a std::string.
     *
     * @param rOption  the option name as a string.
     * @param valueNumber  the number of the argument following the option definition (defaults to 1, for 1st argument).
     */
    std::string GetStringCorrespondingToOption(const std::string& rOption, int valueNumber=1);

    /**
     * @return a collection of strings for a given option (useful for inputting a list of files for example).
     *
     * This uses GetStringCorrespondingToOption repeatedly.
     *
     * @param rOption  the option name as a string.
     */
    std::vector<std::string> GetStringsCorrespondingToOption(const std::string& rOption);

    /**
     * @return a collection of doubles for a given option.
     *
     * This uses GetDoubleCorrespondingToOption repeatedly.
     *
     * @param rOption  the option name as a string.
     */
    std::vector<double> GetDoublesCorrespondingToOption(const std::string& rOption);

    /**
     * @return a collection of ints for a given option.
     *
     * This uses GetIntCorrespondingToOption repeatedly.
     *
     * @param rOption  the option name as a string.
     */
    std::vector<int> GetIntsCorrespondingToOption(const std::string& rOption);

    /**
     * @return a collection of unsigneds for a given option.
     *
     * This uses GetUnsignedCorrespondingToOption repeatedly.
     *
     * @param rOption  the option name as a string.
     */
    std::vector<unsigned> GetUnsignedsCorrespondingToOption(const std::string& rOption);

    /**
     * Read bools from the command line. This accepts 0 or 1, and
     * also true/True/TRUE and false/False/FALSE.
     *
     * @param rOption  The option name as a string.
     * @return Whether the option was set to true or false.
     */
    bool GetBoolCorrespondingToOption(const std::string& rOption);
};

#endif // COMMANDLINEARGUMENTS_HPP_
