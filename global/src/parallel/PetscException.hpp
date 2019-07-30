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

#ifndef PETSCEXCEPTION_HPP_
#define PETSCEXCEPTION_HPP_

#include <petsc.h>
#include <petscksp.h>
/**
 * Throw an exception if the PETSc error number (what is generally denoted 'ierr' in
 * PETSc code) is non-zero - see PETSCEXCEPT
 *
 * @param petscError PETSc error number
 * @param line
 * @param funct
 * @param file
 */
extern void PetscException(PetscInt petscError, unsigned line,
                           const char* funct, const char* file);

/**
 * Throw an exception if the KSP error indicates linear solve failure - see
 * KSPEXCEPT.
 *
 * @param kspError KSP error number (obtained using KSPGetConvergedReason)
 * @param line
 * @param funct
 * @param file
 */
extern void KspException(PetscInt kspError, unsigned line,
                         const char* funct, const char* file);

/**
 * Throw a warning if the KSP error indicates linear solve failure - see
 * KSPWARNIFFAILED.
 *
 * @param kspError KSP error number (obtained using KSPGetConvergedReason)
 */
extern void KspWarnIfFailed(PetscInt kspError);


/** Helper function for above functions - convert a KSP error number into an error message */
std::string GetKspErrorMessage(PetscInt kspError);


/**
 * Positive codes mean that there's an error.
 * Zero means success.
 * Negative codes should never happen, but we'll throw anyway.
 */
#define PETSCEXCEPT(n) PetscException(n, __LINE__, __func__, __FILE__)

/**
 * Positive codes mean that the KSP converged.
 * Negative codes mean that the KSP diverged, i.e. there's a problem.
 *
 * Throw an Exception if KSP failed to solve.
 */
#define KSPEXCEPT(n)  KspException(n, __LINE__, __func__, __FILE__)

/**
 * Positive codes mean that the KSP converged.
 * Negative codes mean that the KSP diverged, i.e. there's a problem.
 *
 * Throw a Warning if KSP failed to solve.
 */
#define KSPWARNIFFAILED(n)  KspWarnIfFailed(n)


#endif /*PETSCEXCEPTION_HPP_*/
