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

#ifndef ABSTRACTFEASSEMBLERCOMMON_HPP_
#define ABSTRACTFEASSEMBLERCOMMON_HPP_

#include "AbstractFeAssemblerInterface.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "HeartEventHandler.hpp"
#include "LinearBasisFunction.hpp"
#include "PetscTools.hpp"
#include "AbstractTetrahedralMesh.hpp"

/**
 * Enumeration for defining how much interpolation (onto quadrature points) is
 * required by the concrete class.
 *
 * CARDIAC: only interpolates the first component of the unknown (ie the voltage)
 * NORMAL: interpolates the position X and all components of the unknown u
 * NONLINEAR: interpolates X, u and grad(u). Also computes the gradient of the
 *   basis functions when assembling vectors.
 */
typedef enum InterpolationLevel_
{
    CARDIAC = 0,
    NORMAL,
    NONLINEAR
} InterpolationLevel;

/**
 *   A base class for AbstractFeVolumeIntegralAssembler (the main abstract assembler class), AbstractSurfaceFeObjectAssembler, and
 *   AbstractCableFeObjectAssembler.
 *
 *   The base class of this, AbstractFeAssemblerInterface, defines the interface for these assembler classes. This class
 *   just defines a few pde-folder-specific (ie not continuum-mechanics-related) extra methods.
 *
 *   See AbstractFeVolumeIntegralAssembler documentation for info on these assembler classes.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
class AbstractFeAssemblerCommon : public AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>
{
protected:
    /**
     * If the matrix or vector will be dependent on a current solution, say,
     * this is where that information is put.
     */
    ReplicatableVector mCurrentSolutionOrGuessReplicated;

    /**
     * @return an entry from the solution vector.
     *
     * @param nodeIndex node index
     * @param indexOfUnknown index of unknown
     */
    virtual double GetCurrentSolutionOrGuessValue(unsigned nodeIndex, unsigned indexOfUnknown)
    {
        return mCurrentSolutionOrGuessReplicated[ PROBLEM_DIM*nodeIndex + indexOfUnknown];
    }

    /**
     * The concrete subclass can overload this and IncrementInterpolatedQuantities()
     * if there are some quantities which need to be computed at each Gauss point.
     * They are called in AssembleOnElement().
     */
    virtual void ResetInterpolatedQuantities()
    {}

    /**
     * The concrete subclass can overload this and ResetInterpolatedQuantities()
     * if there are some quantities which need to be computed at each Gauss point.
     * They are called in AssembleOnElement().
     *
     * Note that this method is called over assembly of elements, surface elements and cables.
     *
     * @param phiI
     * @param pNode pointer to a node
     */
    virtual void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
    {}

    /**
     * The concrete subclass can overload this and ResetInterpolatedQuantities()
     * if there are some gradient dependent quantities which need to be computed at each Gauss point.
     * A matrix of all the basis function gradients at the quad point is passed for efficiency reasons.
     * To access the gradient vector use of the current basis function use rGradPhi(:, phi_index);
     * They are called in AssembleOnElement().
     *
     * Note that this method is ONLY called during assembly of elements. NOT during assembly of surface elements or cables.
     *
     * Further, it is ONLY called in the cases where rGradPhi has been computed.  Currently these cases are
     *     - When mAssembleMatrix is set
     *  or - When INTERPOLATION_LEVEL==NONLINEAR
     * If there are use-cases, then allow other cases where interpolated gradients are needed in righthand-side assembly
     * (see #2075)
     *
     * @param rGradPhi A matrix containing the gradient of all the basis functions at this Gauss point.
     * @param phiIndex The index of the current basis function in the rGradPhi matrix.
     * @param pNode pointer to the node associated with the current basis function
     */
    virtual void IncrementInterpolatedGradientQuantities(const c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi, unsigned phiIndex, const Node<SPACE_DIM>* pNode)
    {}
public:

    /**
     * Constructor.
     */
    AbstractFeAssemblerCommon();

    /**
     * Set a current solution vector that will be used in AssembleOnElement and can passed
     * up to ComputeMatrixTerm() or ComputeVectorTerm().
     *
     * @param currentSolution Current solution vector.
     */
    void SetCurrentSolution(Vec currentSolution);

    /**
     * Destructor.
     */
    virtual ~AbstractFeAssemblerCommon()
    {
    }
};

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::AbstractFeAssemblerCommon()
    : AbstractFeAssemblerInterface<CAN_ASSEMBLE_VECTOR,CAN_ASSEMBLE_MATRIX>()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM, bool CAN_ASSEMBLE_VECTOR, bool CAN_ASSEMBLE_MATRIX, InterpolationLevel INTERPOLATION_LEVEL>
void AbstractFeAssemblerCommon<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM, CAN_ASSEMBLE_VECTOR, CAN_ASSEMBLE_MATRIX, INTERPOLATION_LEVEL>::SetCurrentSolution(Vec currentSolution)
{
    assert(currentSolution != nullptr);

    // Replicate the current solution and store so can be used in AssembleOnElement
    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    mCurrentSolutionOrGuessReplicated.ReplicatePetscVector(currentSolution);
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

    // The AssembleOnElement type methods will determine if a current solution or
    // current guess exists by looking at the size of the replicated vector, so
    // check the size is zero if there isn't a current solution.
    assert(mCurrentSolutionOrGuessReplicated.GetSize() > 0);
}


#endif /* ABSTRACTFEASSEMBLERCOMMON_HPP_ */
