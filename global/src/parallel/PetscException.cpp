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

#include <sstream>

#include "PetscException.hpp"
#include "Exception.hpp"
#include "Warnings.hpp"
#include "ChasteBuildRoot.hpp"


/*
 * Positive codes mean that there's an error.
 * Zero means success.
 * Negative codes should never happen, but we'll throw anyway.
 */
void PetscException(PetscInt petscError,
                    unsigned line,
                    const char* funct,
                    const char* file)
{
    if (petscError != 0)
    {
        const char*  p_text;
        char default_message[30]="Unknown PETSc error code";

        /*
         * PetscErrorMessage will swing p_text to point to the error code's message
         * ...but only if it's a valid code.
         */
        PetscErrorMessage(petscError, &p_text, nullptr);
        if (p_text == nullptr)
        {
            p_text=default_message;
        }

// /todo #2656 - remove ifdef after cmake transition
#ifdef CHASTE_CMAKE
        const size_t root_dir_length = std::char_traits<char>::length(ChasteSourceRootDir());
        EXCEPTION(p_text << " in function '" << funct  << "' on line "
                  << line << " of file ./" << file+root_dir_length);
#else
        EXCEPTION(p_text << " in function '" << funct  << "' on line "
                  << line << " of file " << file);
#endif
    }
}

/*
 * Positive codes mean that the KSP converged.
 * Negative codes mean that the KSP diverged, i.e. there's a problem.
 */
std::string GetKspErrorMessage(PetscInt kspError)
{
    std::string err_string;

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    switch (kspError)
    {
        case KSP_DIVERGED_ITS:
            err_string = "KSP_DIVERGED_ITS";
            break;
        case KSP_DIVERGED_DTOL:
            err_string = "KSP_DIVERGED_DTOL";
            break;
        case KSP_DIVERGED_BREAKDOWN:
            err_string = "KSP_DIVERGED_BREAKDOWN";
            break;
        case KSP_DIVERGED_BREAKDOWN_BICG:
            err_string = "KSP_DIVERGED_BREAKDOWN_BICG";
            break;
        case KSP_DIVERGED_NONSYMMETRIC:
            err_string = "KSP_DIVERGED_NONSYMMETRIC";
            break;
        case KSP_DIVERGED_INDEFINITE_PC:
            err_string = "KSP_DIVERGED_INDEFINITE_PC";
            break;
        default:
            err_string = "Unknown KSP error code";
    }
#else
    // This array contains the strings describing KSP convergence/divergence reasons.
    // It is exported by libpetscksp.a
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
    //External declaration of KSPConvergedReasons already happens
    extern const char * const *KSPConvergedReasons;
#else
    extern const char **KSPConvergedReasons;
#endif

    // The code for the last known error (-10) is hardcoded in PETSc, in future releases it might change.
    // It is defined in src/ksp/ksp/interface/dlregisksp.c
    if (kspError >= -10)
    {
        err_string = KSPConvergedReasons[kspError];
    }
    else
    {
        err_string = "Unknown KSP error code";
    }
#endif

    return err_string;
}

/*
 * Positive codes mean that the KSP converged.
 * Negative codes mean that the KSP diverged, i.e. there's a problem.
 */
void KspException(PetscInt kspError,
                  unsigned line,
                  const char* funct,
                  const char* file)
{
    if (kspError < 0)
    {
        std::stringstream err_string_stream;
        // /todo #2656 - remove ifdef after cmake transition
#ifdef CHASTE_CMAKE
        const size_t root_dir_length = std::char_traits<char>::length(ChasteSourceRootDir());
        err_string_stream << GetKspErrorMessage(kspError) << " in function '" << funct << "' on line "
                      << line << " of file ./" << file + root_dir_length;
#else
        err_string_stream << GetKspErrorMessage(kspError) << " in function '" << funct << "' on line "
                      << line << " of file " << file;
#endif
        EXCEPTION(err_string_stream.str());
    }
}

void KspWarnIfFailed(PetscInt kspError)
{
    if (kspError < 0)
    {
        std::string message = "Linear solve failed: " + GetKspErrorMessage(kspError);
        WARNING(message);
    }
}
