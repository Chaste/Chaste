/*

Copyright (C) University of Oxford, 2005-2011

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
     * Useful inline function for getting an entry from the current solution vector.
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
     * @param phiI
     * @param pNode pointer to a node
     */
    virtual void IncrementInterpolatedQuantities(double phiI, const Node<SPACE_DIM>* pNode)
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
    assert(currentSolution != NULL);

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
