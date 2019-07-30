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

#ifndef MONODOMAINPURKINJESOLVER_HPP_
#define MONODOMAINPURKINJESOLVER_HPP_

#include "AbstractDynamicLinearPdeSolver.hpp"
#include "MassMatrixAssembler.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "MonodomainCorrectionTermAssembler.hpp"
#include "MonodomainTissue.hpp"
#include "MonodomainPurkinjeVolumeAssembler.hpp"
#include "MonodomainPurkinjeCableAssembler.hpp"


/**
 *  Solver class for Monodomain problems on tissues containing Purkinje fibres.
 *
 *  In such problems, there are a subset of nodes of the mesh which are labelled as
 *  Purkinje nodes (and 1D elements connecting them). There are two variables to be
 *  computed, the (normal, myocardium) transmembrane voltage (V), defined at ALL nodes
 *  of the mesh, and the Purkinje voltage (Vp), defined at the purkinje nodes. Hence, the
 *  Purkinje nodes have two variables defined on them.
 *
 *  For parallelisation/implementation reasons, we choose to have two variables to be
 *  defined at ALL nodes, introducing dummy variables Vp for non-Purkinje nodes. We set
 *  Vp = 0 at non-Purkinje nodes. Hence we solve for {V,Vp} at every node in the mesh,
 *  and PROBLEM_DIM=2. The linear system to be assembled, written as usual in block form
 *  but actually in striped form in the code, is:
 *
 *  [ A1 0 ][V ] = [b1]
 *  [ 0 A2 ][Vp] = [b2]
 *
 *  where each block is of size num_nodes.
 *
 *  A1 and b1 are obtained by integrating  over 3D myocardium elements and are therefore
 *  exactly the same as in a normal monodomain problem.
 *
 *  Suppose the nodes are ordered such that all the Purkinje nodes are last. Then the matrix A2
 *  needs to be of the form
 *  A2 = [I 0 ]
 *       [0 Ap]
 *  where the first block is of size num_non_purkinje_nodes, and the second is of size num_purkinje_nodes.
 *  Ap is obtained by assembling over 1D purkinje elements. After assembly we add in the identity block,
 *  which is just represents the equations Vp=0 for the dummy variables (non-purkinje nodes). Finally,
 *  we similarly have
 *  b2 = [0 ]
 *       [bp]
 *  where bp involves a loop over 1D purkinje elements.
 *
 *  This class implements the above, and is based on (but doesn't inherit from, as the PROBLEM_DIMENSION
 *  is different) MonodomainSolver.
 *
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonodomainPurkinjeSolver
  : public AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>
{
private:
    /** Saved pointer to the mesh in this class, as the pointer saved in the
     *  parent class (AbstractDynamicLinearPdeSolver::mpMesh) is not declared to
     *  be a pointer to a mixed mesh
     */
    MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* mpMixedMesh;

    /** Monodomain tissue class (collection of cells, and conductivities) */
    MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* mpMonodomainTissue;

    /** Boundary conditions */
    BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* mpBoundaryConditions;

    /**
     *  The volume assembler, used to set up volume integral parts of the
     *  LHS matrix
     */
    MonodomainPurkinjeVolumeAssembler<ELEMENT_DIM,SPACE_DIM>* mpVolumeAssembler;
    /**
     *  The cable element assembler, used to set up cable integral parts of the
     *  LHS matrix
     */
    MonodomainPurkinjeCableAssembler<ELEMENT_DIM,SPACE_DIM>* mpCableAssembler;

    /** Assembler for surface integrals coming from any non-zero Neumann boundary conditions */
    NaturalNeumannSurfaceTermAssembler<ELEMENT_DIM,SPACE_DIM,2>* mpNeumannSurfaceTermsAssembler;

    // SVI and Purkinje not yet implemented:
    // MonodomainCorrectionTermAssembler<ELEMENT_DIM,SPACE_DIM>* mpMonodomainCorrectionTermAssembler;

    /** The mass matrix, used to computing the RHS vector */
    Mat mMassMatrix;

    /** The vector multiplied by the mass matrix. Ie, if the linear system to
     *  be solved is Ax=b (excluding surface integrals), this vector is z where b=Mz.
     */
    Vec mVecForConstructingRhs;

    /**
     *  Implementation of SetupLinearSystem() which uses the assembler to compute the
     *  LHS matrix, but sets up the RHS vector using the mass-matrix (constructed
     *  using a separate assembler) multiplied by a vector
     *
     *  @param currentSolution  Solution at current time
     *  @param computeMatrix  Whether to compute the matrix of the linear system
     */
    void SetupLinearSystem(Vec currentSolution, bool computeMatrix);

    /**
     *  For the block in the LHS matrix corresponding to the Purkinje voltage at myocardium nodes:
     *  this block is zero after all elements are assembled so, so we set it to be the
     *  identity block. This is done by just checking which rows of the matrix has zero
     *  diagonal values.
     */
    void SetIdentityBlockToLhsMatrix();

public:
    /**
     *  Overloaded PrepareForSetupLinearSystem() methods which
     *  gets the cell models to solve themselves
     *
     *  @param currentSolution solution at current time
     */
    void PrepareForSetupLinearSystem(Vec currentSolution);

    /**
     *  Overloaded InitialiseForSolve
     *
     *  @param initialSolution initial solution
     */
    virtual void InitialiseForSolve(Vec initialSolution);

    /**
     * Constructor
     *
     * @param pMesh pointer to the mesh
     * @param pTissue pointer to the tissue
     * @param pBoundaryConditions pointer to the boundary conditions
     */
    MonodomainPurkinjeSolver(MixedDimensionMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                             MonodomainTissue<ELEMENT_DIM,SPACE_DIM>* pTissue,
                             BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions);

    /**
     *  Destructor
     */
    virtual ~MonodomainPurkinjeSolver();
};



#endif // MONODOMAINPURKINJESOLVER_HPP_
