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

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "Warnings.hpp"
#include "Exception.hpp"
#include "LogFile.hpp"
#include "BoostFilesystem.hpp"
#include "FileFinder.hpp"
#include "PosixPathFixer.hpp"
#include "GetCurrentWorkingDirectory.hpp"

Warnings* Warnings::mpInstance = nullptr;

Warnings::Warnings()
{
}

void Warnings::NoisyDestroy(void)
{
    PrintWarnings();
    QuietDestroy();
}

void Warnings::QuietDestroy(void)
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

void Warnings::PrintWarnings(void)
{
    if (mpInstance)
    {
        for (WarningsContainerType::iterator it = mpInstance->mWarningMessages.begin();
             it != mpInstance->mWarningMessages.end();
             ++it)
        {
            /*
             * Look at my warnings please.
             * First in pair is the context.
             * Second in pair is that actual warning.
             */
            std::cout << it->first << it->second << std::endl;
        }
    }
}

Warnings* Warnings::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new Warnings();
        std::atexit(NoisyDestroy);
    }
    return mpInstance;
}

void Warnings::AddWarning(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber, bool onlyOnce)
{
    std::string posix_filename(ChastePosixPathFixer::ToPosix(fs::path(rFilename)));
    std::stringstream line_number_stream;
    line_number_stream << lineNumber;
    std::string context("Chaste warning: in file " + posix_filename + " at line "  + line_number_stream.str()  + ": ");
    std::pair<std::string, std::string> item(context, rMessage);

    if (onlyOnce)
    {
        WarningsContainerType::iterator it = find(mWarningMessages.begin(), mWarningMessages.end(), item);
        if (it != mWarningMessages.end())
        {
            return;
        }
    }

    mWarningMessages.push_back(item);
    LOG(1, context + rMessage);
}

unsigned Warnings::GetNumWarnings()
{
    return mWarningMessages.size();
}

std::string Warnings::GetNextWarningMessage()
{
    if (mWarningMessages.empty())
    {
        EXCEPTION("There are no warnings");
    }
    std::string message = mWarningMessages.front().second; // Second in pair is the actual warning
    mWarningMessages.pop_front();

    return message;
}
