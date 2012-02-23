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

#include "GetCurrentWorkingDirectory.hpp"

#include <cstdlib>
#include <cerrno>

#include "Exception.hpp"
#include <unistd.h> // For getcwd()

std::string GetCurrentWorkingDirectory()
{
    size_t bufsize = 1000;
    char* p_buffer = NULL;
    while (true)
    {
        p_buffer = (char*) malloc(bufsize);
        if (!p_buffer)
        {
#define COVERAGE_IGNORE
            // Rather a tricky one to cover...
            EXCEPTION("Run out of memory to allocate CWD buffer");
#undef COVERAGE_IGNORE
        }
        if (getcwd(p_buffer, bufsize) == p_buffer)
        {
            break;
        }
#define COVERAGE_IGNORE
        // Also rather a tricky one to cover...
        else
        {
            free(p_buffer);
            if (errno != ERANGE)
            {
                EXCEPTION("Unable to determine current working directory");
            }
            bufsize *= 2;
        }
#undef COVERAGE_IGNORE
    }

    std::string cwd(p_buffer);
    free(p_buffer);
    return cwd;
}
