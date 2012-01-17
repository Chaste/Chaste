/*

Copyright (C) University of Oxford, 2005-2012

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

#include <iostream>
#include <petsc.h>

#include "Exception.hpp"

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
typedef PetscBool PetscTruth;
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
    mShortMessage = rMessage;
    std::stringstream line_number_stream;
    line_number_stream << lineNumber;
    mMessage = std::string("\nChaste error: ") + rFilename + ":"  + line_number_stream.str()  + ": " + mShortMessage;
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

#define COVERAGE_IGNORE //Termination NEVER EVER happens under normal testing conditions.
void Exception::Terminate(const std::string& rMessage, const std::string& rFilename, unsigned lineNumber)
{
    std::stringstream error_message;

    error_message << "\nChaste termination: " << rFilename << ":" << lineNumber  << ": " << rMessage<<"\n";
    std::cerr << error_message.str();

    /*
     * Check if we're running in parallel.
     */

    PetscTruth is_there;
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

#undef COVERAGE_IGNORE // Termination NEVER EVER happens under normal testing conditions
