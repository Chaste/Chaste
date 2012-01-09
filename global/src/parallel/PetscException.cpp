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

#include <sstream>

#include "PetscException.hpp"
#include "Exception.hpp"
#include "Warnings.hpp"


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
        PetscErrorMessage(petscError, &p_text, NULL);
        if (p_text == 0)
        {
            p_text=default_message;
        }
        EXCEPTION(p_text << " in function '" << funct  << "' on line "
                  << line << " of file " << file);
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
    extern const char **KSPConvergedReasons;

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
        std::string err_string = GetKspErrorMessage(kspError);

        err_string += " in function '";
        err_string += funct;
        err_string += "' on line ";
        err_string += line;
        err_string += " of file ";
        err_string += file;

        EXCEPTION(err_string);
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
