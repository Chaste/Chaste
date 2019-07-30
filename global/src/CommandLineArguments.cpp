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

#include "CommandLineArguments.hpp"

#include <cassert>
#include <cstddef>
#include <algorithm>

CommandLineArguments::CommandLineArguments()
    : p_argc(nullptr),
      p_argv(nullptr)
{
    // Make doubly sure there's only one instance
    assert(mpInstance == nullptr);
}

CommandLineArguments* CommandLineArguments::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new CommandLineArguments;
    }
    return mpInstance;
}

CommandLineArguments* CommandLineArguments::mpInstance = nullptr;

bool CommandLineArguments::OptionExists(const std::string& rOption)
{
    TestOptionFormat(rOption);
    try
    {
        return(GetIndexForArgument(rOption));
    }
    catch (const Exception&)
    {
        return false;
    }
}

char* CommandLineArguments::GetValueCorrespondingToOption(const std::string& rOption, int valueNumber)
{
    EXCEPT_IF_NOT(valueNumber>0);
    TestOptionFormat(rOption);

    int index = GetIndexForArgument(rOption);
    int num_args = GetNumberOfArgumentsForOption(rOption, true);
    if (num_args < valueNumber)
    {
        std::stringstream ss;
        ss<<"Index="<<valueNumber<<" requested for '"<<rOption<<"', but only "<<num_args<<" given.";
        EXCEPTION(ss.str());
    }

    return (*p_argv)[index+valueNumber];
}

double CommandLineArguments::GetDoubleCorrespondingToOption(const std::string& rOption, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(rOption, valueNumber);
    return atof(val);
}

int CommandLineArguments::GetIntCorrespondingToOption(const std::string& rOption, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(rOption, valueNumber);
    return atoi(val);
}

unsigned CommandLineArguments::GetUnsignedCorrespondingToOption(const std::string& rOption, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(rOption, valueNumber);
    int i = atoi(val);
    if (i < 0)
    {
        EXCEPTION("Option is a negative number and cannot be converted to unsigned.");
    }
    return (unsigned)(i);
}

int CommandLineArguments::GetIndexForArgument(std::string rOption)
{
    TestOptionFormat(rOption);

    for (int i=1; i<*p_argc; i++)
    {
        if (rOption==std::string((*p_argv)[i]))
        {
            return i;
        }
    }
    EXCEPTION("Command line option '"+rOption+"' does not exist");
}

int CommandLineArguments::GetNumberOfArgumentsForOption(const std::string& rOption, bool throwIfNone)
{
    int start_idx = GetIndexForArgument(rOption);

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
        if (throwIfNone)
        {
            EXCEPTION("No value(s) given after command line option '" + rOption + "'");
        }
        else
        {
            return 0;
        }
    }

    return end_idx - start_idx;
}

std::string CommandLineArguments::GetStringCorrespondingToOption(const std::string& rOption, int valueNumber)
{
    char* val = GetValueCorrespondingToOption(rOption, valueNumber);
    std::string string_arg(val);
    return string_arg;
}

std::vector<std::string> CommandLineArguments::GetStringsCorrespondingToOption(const std::string& rOption)
{
    std::vector<std::string> strings;
    int num_args = GetNumberOfArgumentsForOption(rOption, true);
    for (int i=1; i<=num_args; ++i)
    {
        strings.push_back(GetStringCorrespondingToOption(rOption, i));
    }
    return strings;
}

std::vector<double> CommandLineArguments::GetDoublesCorrespondingToOption(const std::string& rOption)
{
    std::vector<double> doubles;
    int num_args = GetNumberOfArgumentsForOption(rOption, true);
    for (int i=1; i<=num_args; ++i)
    {
        doubles.push_back(GetDoubleCorrespondingToOption(rOption, i));
    }
    return doubles;
}

std::vector<unsigned> CommandLineArguments::GetUnsignedsCorrespondingToOption(const std::string& rOption)
{
    std::vector<unsigned> unsigneds;
    int num_args = GetNumberOfArgumentsForOption(rOption, true);
    for (int i=1; i<=num_args; ++i)
    {
        unsigneds.push_back(GetUnsignedCorrespondingToOption(rOption, i));
    }
    return unsigneds;
}

std::vector<int> CommandLineArguments::GetIntsCorrespondingToOption(const std::string& rOption)
{
    std::vector<int> ints;
    int num_args = GetNumberOfArgumentsForOption(rOption, true);
    for (int i=1; i<=num_args; ++i)
    {
        ints.push_back(GetIntCorrespondingToOption(rOption, i));
    }
    return ints;
}

bool CommandLineArguments::GetBoolCorrespondingToOption(const std::string& rOption)
{
    bool result;
    std::string string_result = GetStringCorrespondingToOption(rOption);

    // Convert to lowercase
    std::transform(string_result.begin(), string_result.end(), string_result.begin(), ::tolower);

    if (string_result=="0" || string_result=="false")
    {
        result = false;
    }
    else if (string_result=="1" || string_result=="true")
    {
        result = true;
    }
    else
    {
        EXCEPTION("The option '" << rOption << "' argument '" <<
            GetStringCorrespondingToOption(rOption) << "' cannot be converted to a bool.");
    }

    return result;
}

void CommandLineArguments::TestOptionFormat(const std::string& rOption)
{
    if (!( rOption.substr(0,1) == "-" && rOption.substr(1,1).find_first_of("0123456789")==std::string::npos))
    {
        EXCEPTION("A command line option must begin with '-' followed by a non-numeric character.");
    }
}
