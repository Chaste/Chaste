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

#include <cassert>
#include <iostream>
#include <petsc.h>
#include <string>     // std::char_traits

#include "Exception.hpp"
#include "BoostFilesystem.hpp"
#include "FileFinder.hpp"
#include "PosixPathFixer.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "ChasteBuildRoot.hpp"

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 2 || PETSC_VERSION_MAJOR<3 ) // Before PETSc 3.2
typedef PetscTruth PetscBool;
#endif

Exception::Exception(const std::string& rMessage,
                     const std::string& rFilename, unsigned lineNumber)
{
    SetMessage(rMessage, rFilename, lineNumber);
    ///// The following would write the error message to the log file, if one exists.
    ///// It's commented out because you end up with 100s of errors in the log from
    ///// element nodes being swapped around when the mesh is read in.
    //// No way of saying here whether this will be a fatal error, but write
    //// it to the log file (if one exists) in case it is.
    // std::string log_file_message = "Exception occurred (although possibly handled), error message:\n" + message;
    // LOG(1, log_file_message);
}

void Exception::SetMessage(const std::string& rMessage,
                           const std::string& rFilename, unsigned lineNumber)
{
    // Strip off source root dir if exists
    std::string filename(rFilename);
    const size_t root_dir_length = std::char_traits<char>::length(ChasteSourceRootDir());
    if (filename.compare(0,root_dir_length,ChasteSourceRootDir()) == 0)
    {
        filename.replace(0,root_dir_length,"./");
    }

    std::string posix_filename(ChastePosixPathFixer::ToPosix(fs::path(filename)));
    mShortMessage = rMessage;
    std::stringstream line_number_stream;
    line_number_stream << lineNumber;
    mMessage = std::string("\nChaste error: ") + posix_filename + ":"  + line_number_stream.str()  + ": " + mShortMessage;
}

std::string Exception::GetMessage() const
{
    return mMessage;
}

std::string Exception::GetShortMessage() const
{
    return mShortMessage;
}

std::string Exception::CheckShortMessage(std::string expected) const
{
    std::string error;
    if (mShortMessage != expected && mShortMessage != "Another process threw an exception; bailing out.")
    {
        error = "Incorrect exception message thrown: expected (" + expected + "); got (" + mShortMessage + ").";
    }
    return error;
}

std::string Exception::CheckShortMessageContains(std::string expected) const
{
    std::string error;
    if (mShortMessage.find(expected) == std::string::npos && mShortMessage != "Another process threw an exception; bailing out.")
    {
        error = "Incorrect exception message thrown: expected it to contain (" + expected + "); got (" + mShortMessage + ").";
    }
    return error;
}

// LCOV_EXCL_START //Termination NEVER EVER happens under normal testing conditions.
void Exception::Terminate(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber)
{
    std::stringstream error_message;

    error_message << "\nChaste termination: " << rFilename << ":" << lineNumber  << ": " << rMessage<<"\n";
    std::cerr << error_message.str();

    /*
     * Check if we're running in parallel.
     */

    PetscBool is_there;
    PetscInitialized(&is_there);
    if (is_there)
    {
        MPI_Abort(PETSC_COMM_WORLD, EXIT_FAILURE);
    }
    else
    {
        exit(EXIT_FAILURE);
    }
}

// LCOV_EXCL_STOP // Termination NEVER EVER happens under normal testing conditions
