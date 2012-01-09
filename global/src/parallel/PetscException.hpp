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
#define PETSCEXCEPT(n) PetscException(n, __LINE__, __FUNCT__,__FILE__)

/**
 * Positive codes mean that the KSP converged.
 * Negative codes mean that the KSP diverged, i.e. there's a problem.
 *
 * Throw an Exception if KSP failed to solve.
 */
#define KSPEXCEPT(n)  KspException(n, __LINE__, __FUNCT__,__FILE__)

/**
 * Positive codes mean that the KSP converged.
 * Negative codes mean that the KSP diverged, i.e. there's a problem.
 *
 * Throw a Warning if KSP failed to solve.
 */
#define KSPWARNIFFAILED(n)  KspWarnIfFailed(n)


#endif /*PETSCEXCEPTION_HPP_*/
