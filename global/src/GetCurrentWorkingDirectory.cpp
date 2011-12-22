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
