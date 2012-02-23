/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "CommandLineArguments.hpp"

#include <cassert>
#include <cstddef>

CommandLineArguments::CommandLineArguments()
    : p_argc(NULL),
      p_argv(NULL)
{
    // Make doubly sure there's only one instance
    assert(mpInstance == NULL);
}

CommandLineArguments* CommandLineArguments::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new CommandLineArguments;
    }
    return mpInstance;
}

CommandLineArguments* CommandLineArguments::mpInstance = NULL;

bool CommandLineArguments::OptionExists(std::string option)
{
    int index = GetIndexForArgument(option);
    assert(index!=0);

    return (index > 0);
}

char* CommandLineArguments::GetValueCorrespondingToOption(std::string option, int valueNumber)
{
    EXCEPT_IF_NOT(valueNumber>0);
    TestOptionFormat(option);

    int num_args = GetNumberOfArgumentsForOption(option);
    int index = GetIndexForArgument(option);
    EXCEPT_IF_NOT(index>0);
    EXCEPT_IF_NOT(num_args>=valueNumber);
    return (*p_argv)[index+valueNumber];
}

double CommandLineArguments::GetDoubleCorrespondingToOption(std::string option, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(option, valueNumber);
    return atof(val);
}

int CommandLineArguments::GetIntCorrespondingToOption(std::string option, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(option, valueNumber);
    return atoi(val);
}

unsigned CommandLineArguments::GetUnsignedCorrespondingToOption(std::string option, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(option, valueNumber);
    int i = atoi(val);
    if (i < 0)
    {
        EXCEPTION("Option is a negative number and cannot be converted to unsigned.");
    }
    return (unsigned)(i);
}

int CommandLineArguments::GetIndexForArgument(std::string argument)
{
    TestOptionFormat(argument);

    for (int i=1; i<*p_argc; i++)
    {
        if (argument==std::string((*p_argv)[i]))
        {
            return i;
        }
    }
    return -1;
}

int CommandLineArguments::GetNumberOfArgumentsForOption(std::string option)
{
    int start_idx = GetIndexForArgument(option);
    if (start_idx < 0)
    {
        EXCEPTION("Command line option '" + option + "' does not exist");
    }

    int end_idx = start_idx;
    for (int i=start_idx+1; i<*p_argc; i++)
    {
        std::string argument = std::string((*p_argv)[i]);
        if (argument.substr(0,1)=="-" && argument.substr(1,1).find_first_of("0123456789")==std::string::npos)
        {
            break;
        }
        end_idx = i;
    }

    if (end_idx == start_idx)
    {
        EXCEPTION("No value(s) given after command line option '" + option + "'");
    }

    return end_idx - start_idx;
}

std::string CommandLineArguments::GetStringCorrespondingToOption(std::string option, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(option, valueNumber);
    std::string string_arg(val);
    return string_arg;
}

std::vector<std::string> CommandLineArguments::GetStringsCorrespondingToOption(std::string option)
{
    std::vector<std::string> strings;
    int num_args = GetNumberOfArgumentsForOption(option);
    for(int i=1; i<=num_args; ++i)
    {
        strings.push_back(GetStringCorrespondingToOption(option, i));
    }
    return strings;
}

std::vector<double> CommandLineArguments::GetDoublesCorrespondingToOption(std::string option)
{
    std::vector<double> doubles;
    int num_args = GetNumberOfArgumentsForOption(option);
    for(int i=1; i<=num_args; ++i)
    {
        doubles.push_back(GetDoubleCorrespondingToOption(option, i));
    }
    return doubles;
}

std::vector<unsigned> CommandLineArguments::GetUnsignedsCorrespondingToOption(std::string option)
{
    std::vector<unsigned> unsigneds;
    int num_args = GetNumberOfArgumentsForOption(option);
    for(int i=1; i<=num_args; ++i)
    {
        unsigneds.push_back(GetUnsignedCorrespondingToOption(option, i));
    }
    return unsigneds;
}

std::vector<int> CommandLineArguments::GetIntsCorrespondingToOption(std::string option)
{
    std::vector<int> ints;
    int num_args = GetNumberOfArgumentsForOption(option);
    for(int i=1; i<=num_args; ++i)
    {
        ints.push_back(GetIntCorrespondingToOption(option, i));
    }
    return ints;
}

void CommandLineArguments::TestOptionFormat(std::string option)
{
    if ( !( option.substr(0,1) == "-" && option.substr(1,1).find_first_of("0123456789")==std::string::npos ) )
    {
        EXCEPTION("A command line option must begin with '-' followed by a non-numeric character.");
    }
}

