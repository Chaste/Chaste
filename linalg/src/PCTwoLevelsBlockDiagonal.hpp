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

#ifndef PCTWOLEVELSBLOCKDIAGONAL_HPP_
#define PCTWOLEVELSBLOCKDIAGONAL_HPP_

#include <cassert>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscpc.h>
#include "PetscTools.hpp"

/**
 * PETSc will return the control to this function everytime it needs to precondition a vector (i.e. y = inv(M)*x)
 *
 * This function needs to be declared global, since I (migb) haven't found a way of defining it inside a class and
 * be able of passing it by reference.
 *
 * @param pc_context preconditioner context struct. Stores preconditioner state (i.e. PC, Mat, and Vec objects used)
 * @param x unpreconditioned residual.
 * @param y preconditioned residual. y = inv(M)*x
 */
#if (PETSC_VERSION_MAJOR == 3 && (PETSC_VERSION_MINOR == 1 || PETSC_VERSION_MINOR ==  2) ) //PETSc 3.1
PetscErrorCode PCTwoLevelsBlockDiagonalApply(PC pc_context, Vec x, Vec y);
#else
PetscErrorCode PCTwoLevelsBlockDiagonalApply(void* pc_context, Vec x, Vec y);
#endif

/// \todo: #1082 update this description
/**
 *  This class defines a PETSc-compliant purpose-build preconditioner.
 *
 *  Let A be a matrix arising in the FEM discretisation of the bidomain equations with the following block structure:
 *
 *                 A = (A11  B')
 *                     (B   A22)
 *
 *  By creating an instance of this class, one will define the following preconditioner:
 *
 *                 inv(M) = inv( (A11   0)   = (inv(A11)        0)
 *                               (0   A22) )   (0        inv(A22))
 *
 *  The inverses are approximate with one cycle of AMG.
 *
 *  Note: This class requires PETSc to be build including HYPRE library. If it's not available, it will throw the following error:
 *
 *     [0]PETSC ERROR: --------------------- Error Message ------------------------------------
 *     [0]PETSC ERROR: Unknown type. Check for miss-spelling or missing external package needed for type!
 *     [0]PETSC ERROR: Unable to find requested PC type hypre!
 *
 *        and will approximate the inverse of the subblocks with PETSc's default preconditioner (bjacobi at the time of writing this).
 */
class PCTwoLevelsBlockDiagonal
{
public:

    /**
     * This struct defines the state of the preconditioner (initialised data and objects to be reused).
     */
    typedef struct{
        Mat A11_matrix_subblock; /**< Mat object that stores the A11 subblock*/
        Mat A22_B1_matrix_subblock; /**< Mat object that stores the tissue part of the A22 subblock*/
        Mat A22_B2_matrix_subblock; /**< Mat object that stores the bath part of the A22 subblock*/
        PC  PC_amg_A11; /**<  inv(A11) is approximated with ILU*/
        PC  PC_amg_A22_B1; /**<  inv(A22_B1) is approximated with ILU*/
        PC  PC_amg_A22_B2; /**<  inv(A22_B2) is approximated by an AMG cycle. We compute it with HYPRE via a PC object*/
        Vec x1_subvector;/**<  Used to store the first half of the vector to be preconditioned*/
        Vec x21_subvector;/**<  Used to store the tissue part of the second half of the vector to be preconditioned*/
        Vec x22_subvector;/**<  Used to store the bath part of the second half of the vector to be preconditioned*/
        Vec y1_subvector;/**<  Used to store the first half of the preconditioned vector*/
        Vec y21_subvector;/**<  Used to store the tissue part of the second half of the preconditioned vector*/
        Vec y22_subvector;/**<  Used to store the bath part second half of the preconditioned vector*/
        VecScatter A11_scatter_ctx;/**< Scattering context: gather x1 from x and scatter y1 back into y*/
        VecScatter A22_B1_scatter_ctx;/**< Scattering context: gather x21 from x and scatter y21 back into y*/
        VecScatter A22_B2_scatter_ctx;/**< Scattering context: gather x22 from x and scatter y22 back into y*/

    } PCTwoLevelsBlockDiagonalContext;

    PCTwoLevelsBlockDiagonalContext mPCContext; /**< PC context, this will be passed to PCTwoLevelsBlockDiagonalApply when PETSc returns control to our preconditioner subroutine.  See PCShellSetContext().*/
    PC mPetscPCObject;/**< Generic PETSc preconditioner object */

public:

    /**
     * Constructor.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     * @param rBathNodes a list of nodes defining the bath
     */
    PCTwoLevelsBlockDiagonal(KSP& rKspObject, std::vector<PetscInt>& rBathNodes);

    /**
     * Destructor.
     */
    ~PCTwoLevelsBlockDiagonal();

private:

    /**
     * Creates all the state data required by the preconditioner.
     *
     * @param rKspObject KSP object where we want to install the block diagonal preconditioner.
     * @param rBathNodes a list of nodes defining the bath
     */
    void PCTwoLevelsBlockDiagonalCreate(KSP& rKspObject, std::vector<PetscInt>& rBathNodes);

    /**
     * Setups preconditioner.
     */
    void PCTwoLevelsBlockDiagonalSetUp();
};

#endif /*PCTWOLEVELSBLOCKDIAGONAL_HPP_*/
