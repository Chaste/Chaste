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

#ifndef CHASTESYSCALLS_HPP_
#define CHASTESYSCALLS_HPP_

/**
 * @file
 * This file provides access to some 'system calls' in a cross-platform manner.
 * It provides the normal Linux names, even on Windows.
 * Functions provided: chdir, getpid, chmod, setenv.
 */

#ifdef _MSC_VER

#include <direct.h>
#define chdir _chdir

#include <process.h>
#define getpid _getpid

#include <io.h>
#define chmod _chmod
#include <sys/stat.h>
#define CHASTE_READONLY _S_IREAD
#define CHASTE_READ_EXECUTE _S_IREAD | _S_IEXEC
#define CHASTE_READ_WRITE _S_IREAD | _S_IWRITE
#define CHASTE_READ_WRITE_EXECUTE _S_IREAD | _S_IWRITE | _S_IEXEC

/**
 * Windows version of setenv call.  Note that under Linux we always pass 1 (overwrite) for the third arg.
 * @param name  environment variable to set
 * @param value  value to give the variable
 * @param mode  not used on Windows
 */
#define setenv(name, value, mode) _putenv_s(name, value)

#else

#include <unistd.h> // For chdir() and getpid()
#include <sys/stat.h> // For chmod()
/** Mode for chmod() to set readonly permissions for everyone. */
#define CHASTE_READONLY 0444
/** Mode for chmod() to set read & execute permissions for everyone. */
#define CHASTE_READ_EXECUTE 0555
/** Mode for chmod() to set read-write permissions for owner, and readonly for group. */
#define CHASTE_READ_WRITE 0640
/** Mode for chmod() to set full permissions for owner, and read-execute for everyone else. */
#define CHASTE_READ_WRITE_EXECUTE 0755

#endif // _MSC_VER

#endif // CHASTESYSCALLS_HPP_
