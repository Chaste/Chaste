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

#include "Debug.hpp"
// If you ever need backtrace in Windows then the RTFM begins at "CaptureStackBackTrace"
#ifndef _MSC_VER
#include <execinfo.h>     //For backtrace
#include <unistd.h>       // For memory profiling
#include <sys/resource.h> // For memory profiling
#endif//_MSC_VER
//#include <cxxabi.h> // For demangling C++ C-style names

double eMemoryAtMarker = 0.0;

std::string FormDebugHead()
{
    std::stringstream header;
    header << "DEBUG: ";
    if (PetscTools::IsParallel())
    {
        header << "proc " << PetscTools::GetMyRank() << ": ";
    }
    return header.str();
}


void PrintTheStack()
{
    // If you ever need backtrace in Windows then the RTFM begins at "CaptureStackBackTrace"
    TRACE("Stack information");
#ifndef _MSC_VER
    // storage array for stack trace address data
    void* address_list[20u]; //20 is about the number of stack symbols we are going to print

    // retrieve current stack addresses
    unsigned num_addresses = backtrace(address_list, sizeof(address_list) / sizeof(void*));

    char** symbol_list = backtrace_symbols(address_list, num_addresses);

    // iterate over the returned symbol lines.
    for (unsigned i = 0; i < num_addresses; i++)
    {
        // We could demangle, but this may not be portable on GNU Linux versus Mac OSX
        // char* ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
        TRACE("Level " << i << ": " << symbol_list[i]);
    }
    free(symbol_list);
#endif //_MSC_VER
}

void MarkMemory()
{
    TRACE("PRINT_MEMORY will now be relative to here");

#ifndef _MSC_VER
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );

    // Store current memory footprint in Mb to allow calculations relative to this point
    eMemoryAtMarker = double(rusage.ru_maxrss)/1024.0;
#endif //_MSC_VER
}

void UnmarkMemory()
{
    TRACE("PRINT_MEMORY will now be absolute footprint");

#ifndef _MSC_VER
    eMemoryAtMarker = 0.0;
#endif //_MSC_VER
}

void PrintMemory()
{
#ifndef _MSC_VER
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );

    // Memory as difference between now and the marker
    double memory = double(rusage.ru_maxrss)/1024.0 - eMemoryAtMarker;

    if (fabs(eMemoryAtMarker) < 1e-6) // If marker is zero, this is the total memory footprint
    {
        std::cout << FormDebugHead() << "Total memory footprint: " << memory << " Mb" << std::endl << std::flush;
    }
    else // If marker non-zero, this is a change in memory
    {
        std::cout << FormDebugHead() << "Memory change from marker: " << memory << " Mb" << std::endl << std::flush;
    }
#endif //_MSC_VER
}
