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

#ifndef _BOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_
#define _BOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_

#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "DistributedVector.hpp"
#include "Exception.hpp"
#include "HeartEventHandler.hpp"
#include "PetscMatTools.hpp"
#include "PetscVecTools.hpp"
#include "ReplicatableVector.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::BoundaryConditionsContainer(bool deleteConditions)
            : AbstractBoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>(deleteConditions)
{
    mLoadedFromArchive = false;

    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        mpNeumannMap[index_of_unknown] = new std::map< const BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *, const AbstractBoundaryCondition<SPACE_DIM>*>;

        mAnyNonZeroNeumannConditionsForUnknown[index_of_unknown] = false;
        mLastNeumannCondition[index_of_unknown] = mpNeumannMap[index_of_unknown]->begin();

        mpPeriodicBcMap[index_of_unknown] = new std::map< const Node<SPACE_DIM> *, const Node<SPACE_DIM> * >;
    }

    // This zero boundary condition is only used in AddNeumannBoundaryCondition
    mpZeroBoundaryCondition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~BoundaryConditionsContainer()
{
    // Keep track of what boundary condition objects we've deleted
    std::set<const AbstractBoundaryCondition<SPACE_DIM>*> deleted_conditions;
    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        NeumannMapIterator neumann_iterator = mpNeumannMap[i]->begin();
        while (neumann_iterator != mpNeumannMap[i]->end() )
        {

            if (deleted_conditions.count(neumann_iterator->second) == 0)
            {
                deleted_conditions.insert(neumann_iterator->second);
                //Leave the zero boundary condition until last
                if (neumann_iterator->second != mpZeroBoundaryCondition)
                {
                    if (this->mDeleteConditions)
                    {
                        delete neumann_iterator->second;
                    }
                }
            }
            neumann_iterator++;
        }
        delete(mpNeumannMap[i]);
        delete(mpPeriodicBcMap[i]);
    }

    delete mpZeroBoundaryCondition;

    if (this->mDeleteConditions)
    {
        this->DeleteDirichletBoundaryConditions(deleted_conditions);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AddDirichletBoundaryCondition(const Node<SPACE_DIM>* pBoundaryNode,
                                        const AbstractBoundaryCondition<SPACE_DIM>* pBoundaryCondition,
                                        unsigned indexOfUnknown,
                                        bool checkIfBoundaryNode)
{
    assert(indexOfUnknown < PROBLEM_DIM);
    if (checkIfBoundaryNode)
    {
        assert(pBoundaryNode->IsBoundaryNode());
    }

    (*(this->mpDirichletMap[indexOfUnknown]))[pBoundaryNode] = pBoundaryCondition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AddPeriodicBoundaryCondition(const Node<SPACE_DIM>* pNode1,
                                                                                                  const Node<SPACE_DIM>* pNode2)
{
    assert(pNode1->IsBoundaryNode());
    assert(pNode2->IsBoundaryNode());

    // will assume the periodic BC is to be applied to ALL unknowns, can't really imagine a
    // situation where this isn't going to be true. If necessary can easily change this method
    // to take in the index of the unknown
    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        (*(this->mpPeriodicBcMap[i]))[pNode1] = pNode2;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AddNeumannBoundaryCondition( const BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> * pBoundaryElement,
                                      const AbstractBoundaryCondition<SPACE_DIM> * pBoundaryCondition,
                                      unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    /*
     * If this condition is constant, we can test whether it is zero.
     * Otherwise we assume that this could be a non-zero boundary condition.
     */
    const ConstBoundaryCondition<SPACE_DIM>* p_const_cond = dynamic_cast<const ConstBoundaryCondition<SPACE_DIM>*>(pBoundaryCondition);
    if (p_const_cond)
    {
        if (p_const_cond->GetValue(pBoundaryElement->GetNode(0)->GetPoint()) != 0.0)
        {
            mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = true;
        }
    }
    else
    {
        mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = true;
    }

    for (unsigned unknown=0; unknown<PROBLEM_DIM; unknown++)
    {
        if (unknown == indexOfUnknown)
        {
            (*(mpNeumannMap[indexOfUnknown]))[pBoundaryElement] = pBoundaryCondition;
        }
        else
        {
            // If can't find pBoundaryElement in map[unknown]
            if (mpNeumannMap[unknown]->find(pBoundaryElement)==mpNeumannMap[unknown]->end())
            {
                // Add zero bc to other unknowns (so all maps are in sync)
                (*(mpNeumannMap[unknown]))[pBoundaryElement] = mpZeroBoundaryCondition;
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineZeroDirichletOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                           unsigned indexOfUnknown)
{
    this->DefineConstantDirichletOnMeshBoundary(pMesh, 0.0, indexOfUnknown);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineConstantDirichletOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                               double value,
                                               unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);
    // In applying a condition to the boundary, we need to be sure that the boundary exists
    assert(PetscTools::ReplicateBool( pMesh->GetNumBoundaryNodes() > 0 ) );

    ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(value);

    typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryNodeIterator iter;
    iter = pMesh->GetBoundaryNodeIteratorBegin();
    while (iter != pMesh->GetBoundaryNodeIteratorEnd())
    {
        AddDirichletBoundaryCondition(*iter, p_boundary_condition, indexOfUnknown);
        iter++;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineZeroNeumannOnMeshBoundary(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                         unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    // In applying a condition to the boundary, we need to be sure that the boundary exists
    assert(pMesh->GetNumBoundaryElements() > 0);
    ConstBoundaryCondition<SPACE_DIM>* p_zero_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>( 0.0 );

    typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::BoundaryElementIterator iter;
    iter = pMesh->GetBoundaryElementIteratorBegin();
    while (iter != pMesh->GetBoundaryElementIteratorEnd())
    {
        AddNeumannBoundaryCondition(*iter, p_zero_boundary_condition, indexOfUnknown);
        iter++;
    }

    mAnyNonZeroNeumannConditionsForUnknown[indexOfUnknown] = false;
}

/**
 * Modifies a linear system to incorporate Dirichlet boundary conditions
 *
 * The BCs are imposed in such a way as to ensure that a symmetric linear system remains symmetric.
 * For each node with a boundary condition applied, both the corresponding row and column are zero'd
 * and the RHS vector modified to take into account the zero'd column. See #577.
 *
 * Suppose we have a matrix
 * [a b c] [x] = [ b1 ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and we want to apply the boundary condition x=v without losing symmetry if the matrix is
 * symmetric. We apply the boundary condition
 * [1 0 0] [x] = [ v  ]
 * [d e f] [y]   [ b2 ]
 * [g h i] [z]   [ b3 ]
 * and then zero the column as well, adding a term to the RHS to take account for the
 * zero-matrix components
 * [1 0 0] [x] = [ v  ] - v[ 0 ]
 * [0 e f] [y]   [ b2 ]    [ d ]
 * [0 h i] [z]   [ b3 ]    [ g ]
 * Note the last term is the first column of the matrix, with one component zeroed, and
 * multiplied by the boundary condition value. This last term is then stored in
 * rLinearSystem.rGetDirichletBoundaryConditionsVector(), and in general form is the
 * SUM_{d=1..D} v_d a'_d
 * where v_d is the boundary value of boundary condition d (d an index into the matrix),
 * and a'_d is the dth-column of the matrix but with the d-th component zeroed, and where
 * there are D boundary conditions
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ApplyDirichletToLinearProblem(
        LinearSystem& rLinearSystem,
        bool applyToMatrix,
        bool applyToRhsVector)
{
    HeartEventHandler::BeginEvent(HeartEventHandler::DIRICHLET_BCS);

    if (applyToMatrix)
    {
        if (!this->HasDirichletBoundaryConditions())
        {
            // Short-circuit the replication if there are no conditions
            HeartEventHandler::EndEvent(HeartEventHandler::DIRICHLET_BCS);
            return;
        }

        bool matrix_is_symmetric = rLinearSystem.IsMatrixSymmetric();

        if (matrix_is_symmetric)
        {
            /*
             * Modifications to the RHS are stored in the Dirichlet boundary
             * conditions vector. This is done so that they can be reapplied
             * at each time step.
             * Make a new vector to store the Dirichlet offsets in.
             */
            Vec& r_bcs_vec = rLinearSystem.rGetDirichletBoundaryConditionsVector();
            if (!r_bcs_vec)
            {
                VecDuplicate(rLinearSystem.rGetRhsVector(), &r_bcs_vec);
            }
            PetscVecTools::Zero(r_bcs_vec);
            /*
             * If the matrix is symmetric, calls to GetMatrixRowDistributed()
             * require the matrix to be in assembled state. Otherwise we can
             * defer it.
             */
            rLinearSystem.AssembleFinalLinearSystem();
        }

        // Work out where we're setting Dirichlet boundary conditions *everywhere*, not just those locally known
        ReplicatableVector dirichlet_conditions(rLinearSystem.GetSize());
        unsigned lo, hi;
        {
            PetscInt lo_s, hi_s;
            rLinearSystem.GetOwnershipRange(lo_s, hi_s);
            lo = lo_s; hi = hi_s;
        }
        // Initialise all local entries to DBL_MAX, i.e. don't know if there's a condition
        for (unsigned i=lo; i<hi; i++)
        {
            dirichlet_conditions[i] = DBL_MAX;
        }
        // Now fill in the ones we know
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

            while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
            {
                unsigned node_index = this->mDirichIterator->first->GetIndex();
                double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());
                assert(value != DBL_MAX);

                unsigned row = PROBLEM_DIM*node_index + index_of_unknown;
                dirichlet_conditions[row] = value;

                this->mDirichIterator++;
            }
        }

        // And replicate
        dirichlet_conditions.Replicate(lo, hi);

        // Which rows have conditions?
        std::vector<unsigned> rows_to_zero;
        for (unsigned i=0; i<dirichlet_conditions.GetSize(); i++)
        {
            if (dirichlet_conditions[i] != DBL_MAX)
            {
                rows_to_zero.push_back(i);
            }
        }

        if (matrix_is_symmetric)
        {
            // Modify the matrix columns
            for (unsigned i=0; i<rows_to_zero.size(); i++)
            {
                unsigned col = rows_to_zero[i];
                double minus_value = -dirichlet_conditions[col];

                /*
                 * Get a vector which will store the column of the matrix (column d,
                 * where d is the index of the row (and column) to be altered for the
                 * boundary condition. Since the matrix is symmetric when get row
                 * number "col" and treat it as a column. PETSc uses compressed row
                 * format and therefore getting rows is far more efficient than getting
                 * columns.
                 */
                Vec matrix_col = rLinearSystem.GetMatrixRowDistributed(col);

                // Zero the correct entry of the column
                PetscVecTools::SetElement(matrix_col, col, 0.0);

                /*
                 * Set up the RHS Dirichlet boundary conditions vector.
                 * Assuming one boundary at the zeroth node (x_0 = value), this is equal to
                 *   -value*[0 a_21 a_31 .. a_N1]
                 * and will be added to the RHS.
                 */
                PetscVecTools::AddScaledVector(rLinearSystem.rGetDirichletBoundaryConditionsVector(), matrix_col, minus_value);
                PetscTools::Destroy(matrix_col);
            }
        }

        /*
         * Now zero the appropriate rows and columns of the matrix. If the matrix
         * is symmetric we apply the boundary conditions in a way the symmetry isn't
         * lost (rows and columns). If not only the row is zeroed.
         */
        if (matrix_is_symmetric)
        {
            rLinearSystem.ZeroMatrixRowsAndColumnsWithValueOnDiagonal(rows_to_zero, 1.0);
        }
        else
        {
            rLinearSystem.ZeroMatrixRowsWithValueOnDiagonal(rows_to_zero, 1.0);
        }
    }

    if (applyToRhsVector)
    {
        // Apply the RHS boundary conditions modification if required.
        if (rLinearSystem.rGetDirichletBoundaryConditionsVector())
        {
            PetscVecTools::AddScaledVector(rLinearSystem.rGetRhsVector(), rLinearSystem.rGetDirichletBoundaryConditionsVector(), 1.0);
        }

        /*
         * Apply the actual boundary condition to the RHS, note this must be
         * done after the modification to the RHS vector.
         */
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

            while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
            {
                unsigned node_index = this->mDirichIterator->first->GetIndex();
                double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());

                unsigned row = PROBLEM_DIM*node_index + index_of_unknown;

                rLinearSystem.SetRhsVectorElement(row, value);

                this->mDirichIterator++;
            }
        }
    }

    HeartEventHandler::EndEvent(HeartEventHandler::DIRICHLET_BCS);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ApplyPeriodicBcsToLinearProblem(LinearSystem& rLinearSystem,
                                                                                                     bool applyToMatrix,
                                                                                                     bool applyToRhsVector)
{
    bool has_periodic_bcs = false;
    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        if (!mpPeriodicBcMap[i]->empty())
        {
            has_periodic_bcs = true;
            break;
        }
    }

    EXCEPT_IF_NOT(has_periodic_bcs);

    if (applyToMatrix)
    {
        std::vector<unsigned> rows_to_zero;
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            for (typename std::map< const Node<SPACE_DIM> *, const Node<SPACE_DIM> * >::const_iterator iter = mpPeriodicBcMap[index_of_unknown]->begin();
                iter != mpPeriodicBcMap[index_of_unknown]->end();
                ++iter)
            {
                unsigned node_index_1 = iter->first->GetIndex();
                unsigned row_index_1 = PROBLEM_DIM*node_index_1 + index_of_unknown;
                rows_to_zero.push_back(row_index_1);
            }
        }

        rLinearSystem.ZeroMatrixRowsWithValueOnDiagonal(rows_to_zero, 1.0);

        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            for (typename std::map< const Node<SPACE_DIM> *, const Node<SPACE_DIM> * >::const_iterator iter = mpPeriodicBcMap[index_of_unknown]->begin();
                iter != mpPeriodicBcMap[index_of_unknown]->end();
                ++iter)
            {
                unsigned node_index_1 = iter->first->GetIndex();
                unsigned node_index_2 = iter->second->GetIndex();

                unsigned mat_index1 = PROBLEM_DIM*node_index_1 + index_of_unknown;
                unsigned mat_index2 = PROBLEM_DIM*node_index_2 + index_of_unknown;
                PetscMatTools::SetElement(rLinearSystem.rGetLhsMatrix(), mat_index1, mat_index2, -1.0);
            }
        }
    }

    if (applyToRhsVector)
    {
        for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
        {
            for (typename std::map< const Node<SPACE_DIM> *, const Node<SPACE_DIM> * >::const_iterator iter = mpPeriodicBcMap[index_of_unknown]->begin();
                 iter != mpPeriodicBcMap[index_of_unknown]->end();
                 ++iter)
            {
                unsigned node_index = iter->first->GetIndex();
                unsigned row_index = PROBLEM_DIM*node_index + index_of_unknown;
                rLinearSystem.SetRhsVectorElement(row_index, 0.0);
            }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ApplyDirichletToNonlinearResidual(
        const Vec currentSolution,
        Vec residual,
        DistributedVectorFactory& rFactory)
{
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

        DistributedVector solution_distributed = rFactory.CreateDistributedVector(currentSolution, true /*Read-only*/);
        DistributedVector residual_distributed = rFactory.CreateDistributedVector(residual);


        while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
        {
            DistributedVector::Stripe solution_stripe(solution_distributed, index_of_unknown);
            DistributedVector::Stripe residual_stripe(residual_distributed, index_of_unknown);

            unsigned node_index = this->mDirichIterator->first->GetIndex();

            double value = this->mDirichIterator->second->GetValue(this->mDirichIterator->first->GetPoint());

            if (solution_distributed.IsGlobalIndexLocal(node_index))
            {
                residual_stripe[node_index]=solution_stripe[node_index] - value;
            }
            this->mDirichIterator++;
        }
        // Don't restore the read-only one: solution_distributed.Restore();
        residual_distributed.Restore();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::ApplyDirichletToNonlinearJacobian(Mat jacobian)
{
    unsigned num_boundary_conditions = 0;
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        num_boundary_conditions += this->mpDirichletMap[index_of_unknown]->size();
    }

    std::vector<unsigned> rows_to_zero(num_boundary_conditions);

    unsigned counter=0;
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        this->mDirichIterator = this->mpDirichletMap[index_of_unknown]->begin();

        while (this->mDirichIterator != this->mpDirichletMap[index_of_unknown]->end() )
        {
            unsigned node_index = this->mDirichIterator->first->GetIndex();
            rows_to_zero[counter++] = PROBLEM_DIM*node_index + index_of_unknown;
            this->mDirichIterator++;
        }
    }
    PetscMatTools::Finalise(jacobian);
    PetscMatTools::ZeroRowsWithValueOnDiagonal(jacobian, rows_to_zero, 1.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::Validate(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    bool valid = true;

    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        // Iterate over surface elements
        typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::BoundaryElementIterator elt_iter
        = pMesh->GetBoundaryElementIteratorBegin();
        while (valid && elt_iter != pMesh->GetBoundaryElementIteratorEnd())
        {
            if (!HasNeumannBoundaryCondition(*elt_iter, index_of_unknown))
            {
                // Check for Dirichlet conditions on this element's nodes
                for (unsigned i=0; i<(*elt_iter)->GetNumNodes(); i++)
                {
                    if (!this->HasDirichletBoundaryCondition((*elt_iter)->GetNode(i)))
                    {
                        valid = false;
                    }
                }
            }
            elt_iter++;
        }
    }
    return valid;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetNeumannBCValue(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* pSurfaceElement,
                             const ChastePoint<SPACE_DIM>& rX,
                             unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    // Did we see this condition on the last search we did?
    if (mLastNeumannCondition[indexOfUnknown] == mpNeumannMap[indexOfUnknown]->end() ||
        mLastNeumannCondition[indexOfUnknown]->first != pSurfaceElement)
    {
        mLastNeumannCondition[indexOfUnknown] = mpNeumannMap[indexOfUnknown]->find(pSurfaceElement);
    }
    if (mLastNeumannCondition[indexOfUnknown] == mpNeumannMap[indexOfUnknown]->end())
    {
        // No Neumann condition is equivalent to a zero Neumann condition
        return 0.0;
    }
    else
    {
        return mLastNeumannCondition[indexOfUnknown]->second->GetValue(rX);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::HasNeumannBoundaryCondition(const BoundaryElement<ELEMENT_DIM-1,SPACE_DIM>* pSurfaceElement,
                                     unsigned indexOfUnknown)
{
    assert(indexOfUnknown < PROBLEM_DIM);

    mLastNeumannCondition[indexOfUnknown] = mpNeumannMap[indexOfUnknown]->find(pSurfaceElement);

    return (mLastNeumannCondition[indexOfUnknown] != mpNeumannMap[indexOfUnknown]->end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AnyNonZeroNeumannConditions()
{
    bool ret = false;
    for (unsigned index_of_unknown=0; index_of_unknown<PROBLEM_DIM; index_of_unknown++)
    {
        if (mAnyNonZeroNeumannConditionsForUnknown[index_of_unknown] == true)
        {
            ret = true;
        }
    }
    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
typename BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::NeumannMapIterator BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::BeginNeumann()
{
    // [0] is ok as all maps will be in sync due to the way ApplyNeumannBoundaryCondition works
    return mpNeumannMap[0]->begin();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
typename BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::NeumannMapIterator BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::EndNeumann()
{
    // [0] is ok as all maps will be in sync due to the way ApplyNeumannBoundaryCondition works
    return mpNeumannMap[0]->end();
}

#endif // _BOUNDARYCONDITIONSCONTAINERIMPLEMENTATION_HPP_
