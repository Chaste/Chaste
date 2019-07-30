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


#include "AbstractBidomainSolver.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscMatTools.hpp"
#include "PetscVecTools.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    // The base class method that calls this function will only call it with a null linear system
    assert(this->mpLinearSystem == NULL);

    // linear system created here
    AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 2>::InitialiseForSolve(initialSolution);

    if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
#ifdef TRACE_KSP
        if (PetscTools::AmMaster())
        {
            std::cout << "Using absolute tolerance: " << mpConfig->GetAbsoluteTolerance() << "\n";
        }
#endif
        this->mpLinearSystem->SetAbsoluteTolerance(mpConfig->GetAbsoluteTolerance());
    }
    else
    {
#ifdef TRACE_KSP
        if (PetscTools::AmMaster())
        {
            std::cout << "Using relative tolerance: " << mpConfig->GetRelativeTolerance() << "\n";
        }
#endif
        this->mpLinearSystem->SetRelativeTolerance(mpConfig->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());

    // Note that twolevelblockdiagonal was never finished.  It was a preconditioner specific to the Parabolic-Parabolic formulation of Bidomain
    // Two levels block diagonal only worked in serial with TetrahedralMesh.
    assert(std::string(HeartConfig::Instance()->GetKSPPreconditioner()) != std::string("twolevelsblockdiagonal"));
    this->mpLinearSystem->SetPcType(HeartConfig::Instance()->GetKSPPreconditioner());

    if (mRowForAverageOfPhiZeroed == INT_MAX)
    {
        // not applying average(phi)=0 constraint, so matrix is symmetric
        this->mpLinearSystem->SetMatrixIsSymmetric(true);
    }
    else
    {
        //Turn off preallocation so that we can have one dense row in the matrix.
        PetscMatTools::TurnOffVariableAllocationError(this->mpLinearSystem->rGetLhsMatrix());

        // applying average(phi)=0 constraint, so matrix is not symmetric
        this->mpLinearSystem->SetMatrixIsSymmetric(false);
    }

    this->mpLinearSystem->SetUseFixedNumberIterations(
        HeartConfig::Instance()->GetUseFixedNumberIterationsLinearSolver(),
        HeartConfig::Instance()->GetEvaluateNumItsEveryNSolves());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec existingSolution)
{
    mpBidomainTissue->SolveCellSystems(existingSolution, PdeSimulationTime::GetTime(), PdeSimulationTime::GetNextTime());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Vec AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::GenerateNullBasis() const
{
    double sqrt_num_nodes = sqrt((double) this->mpMesh->GetNumNodes());

    Vec null_basis;
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    null_basis = p_factory->CreateVec(2);

    DistributedVector dist_null_basis = p_factory->CreateDistributedVector(null_basis);
    DistributedVector::Stripe null_basis_stripe_0(dist_null_basis,0);
    DistributedVector::Stripe null_basis_stripe_1(dist_null_basis,1);
    for (DistributedVector::Iterator index = dist_null_basis.Begin();
         index != dist_null_basis.End();
         ++index)
    {
        null_basis_stripe_0[index] = 0.0;
        null_basis_stripe_1[index] = 1.0/sqrt_num_nodes; // normalised vector
    }
    dist_null_basis.Restore();

    return null_basis;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::FinaliseLinearSystem(Vec existingSolution)
{
    // If no dirichlet boundary conditions
    //  (i) Check compatibility condition to check we are solving
    //      a linear system that can be solved
    // Then either
    //  (a) If not setting average(phi)=0, we are solving a singular system,
    //      so set up a null space
    //  (b) Apply average(phi)=0 constraint by altering the last row, to
    //      get a non-singular system
    //
    if (!(GetBoundaryConditions()->HasDirichletBoundaryConditions()))
    {
        // first check compatibility condition
        CheckCompatibilityCondition();

        // Check whether applying average(phi_e)=0 constraint
        if (mRowForAverageOfPhiZeroed==INT_MAX)
        {
            // We're not using the constraint, hence use a null space
            if (!mNullSpaceCreated)
            {
                // No null space set up, so create one and pass it to the linear system
                Vec null_basis[] = {GenerateNullBasis()};

                this->mpLinearSystem->SetNullBasis(null_basis, 1);

                PetscTools::Destroy(null_basis[0]);
                mNullSpaceCreated = true;
            }
        }
        else
        {
            // Using the average(phi_e)=0 constraint

            // CG (default solver) won't work since the system isn't symmetric anymore. Switch to GMRES
            this->mpLinearSystem->SetKspType("gmres"); // Switches the solver
            mpConfig->SetKSPSolver("gmres", true); // Makes sure this change will be reflected in the XML file written to disk at the end of the simulation.
            //(If the user doesn't have gmres then the "true" warns the user about the switch)

            // Set average phi_e to zero
            unsigned matrix_size = this->mpLinearSystem->GetSize();
            if (!this->mMatrixIsAssembled)
            {

                // Set the mRowForAverageOfPhiZeroed-th matrix row to 0 1 0 1 ...
                std::vector<unsigned> row_for_average;
                row_for_average.push_back(mRowForAverageOfPhiZeroed);
                this->mpLinearSystem->ZeroMatrixRowsWithValueOnDiagonal(row_for_average, 0.0);
                for (unsigned col_index=0; col_index<matrix_size; col_index++)
                {
                    if (col_index%2 == 1)
                    {
                        this->mpLinearSystem->SetMatrixElement(mRowForAverageOfPhiZeroed, col_index, 1);
                    }

                }
                this->mpLinearSystem->FinaliseLhsMatrix();

            }
            // Set the mRowForAverageOfPhiZeroed-th rhs vector row to 0
            this->mpLinearSystem->SetRhsVectorElement(mRowForAverageOfPhiZeroed, 0);

            this->mpLinearSystem->FinaliseRhsVector();
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::CheckCompatibilityCondition()
{
    if (GetBoundaryConditions()->HasDirichletBoundaryConditions() || this->mRowForAverageOfPhiZeroed!=INT_MAX)
    {
        // not a singular system, no compatibility condition
        return;
    }

#ifndef NDEBUG
    DistributedVector distributed_rhs=this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(this->mpLinearSystem->rGetRhsVector());
    DistributedVector::Stripe distributed_rhs_phi_entries(distributed_rhs,1);

    double local_sum=0;
    for (DistributedVector::Iterator index = distributed_rhs.Begin();
         index!= distributed_rhs.End();
         ++index)
    {
        local_sum += distributed_rhs_phi_entries[index];
    }

    double global_sum;
    int mpi_ret = MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    assert(mpi_ret == MPI_SUCCESS);

    if (fabs(global_sum)>1e-6) // magic number! sum should really be a sum of zeros and exactly zero though anyway (or a-a+b-b+c-c.. etc in the case of electrodes)
    {
        // LCOV_EXCL_START
        // shouldn't ever reach this line but useful to have the error printed out if you do
        //std::cout << "Sum of b_{2i+1} = " << global_sum << " (should be zero for compatibility)\n";
        EXCEPTION("Linear system does not satisfy compatibility constraint!");
        // LCOV_EXCL_STOP
    }
#endif
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::AbstractBidomainSolver(
            bool bathSimulation,
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
            BidomainTissue<SPACE_DIM>* pTissue,
            BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,2>* pBoundaryConditions)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,2>(pMesh),
      mBathSimulation(bathSimulation),
      mpBidomainTissue(pTissue),
      mpBoundaryConditions(pBoundaryConditions)
{
    assert(pTissue != NULL);
    assert(pBoundaryConditions != NULL);

    mNullSpaceCreated = false;

    // important!
    this->mMatrixIsConstant = true;

    mRowForAverageOfPhiZeroed = INT_MAX; //this->mpLinearSystem->GetSize() - 1;
    mpConfig = HeartConfig::Instance();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::~AbstractBidomainSolver()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::SetFixedExtracellularPotentialNodes(
            std::vector<unsigned> fixedExtracellularPotentialNodes)
{
    for (unsigned i=0; i<fixedExtracellularPotentialNodes.size(); i++)
    {
        if (fixedExtracellularPotentialNodes[i] >= this->mpMesh->GetNumNodes() )
        {
            EXCEPTION("Fixed node number must be less than total number nodes");
        }
    }

    mFixedExtracellularPotentialNodes = fixedExtracellularPotentialNodes;

    // We will need to recalculate this when HasDirichletBoundaryConditions() is called.
    GetBoundaryConditions()->ResetDirichletCommunication();

    for (unsigned i=0; i<mFixedExtracellularPotentialNodes.size(); i++)
    {
        if (this->mpMesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(mFixedExtracellularPotentialNodes[i]))
        {
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition
                 = new ConstBoundaryCondition<SPACE_DIM>(0.0);

            //Throws if node is not owned locally
            Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(mFixedExtracellularPotentialNodes[i]);

            GetBoundaryConditions()->AddDirichletBoundaryCondition(p_node, p_boundary_condition, 1);

        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::SetRowForAverageOfPhiZeroed(unsigned row)
{
    // Row should be odd in C++-like indexing
    if (row%2 == 0)
    {
        EXCEPTION("Row for applying the constraint 'Average of phi_e = zero' should be odd in C++ like indexing");
    }

    mRowForAverageOfPhiZeroed = row;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractBidomainSolver<ELEMENT_DIM,SPACE_DIM>::FinaliseForBath(bool computeMatrix, bool computeVector)
{
    assert(mBathSimulation);
    PetscBool is_matrix_assembled;
    MatAssembled(this->mpLinearSystem->GetLhsMatrix(), &is_matrix_assembled);
    assert(is_matrix_assembled == PETSC_TRUE);

    /*
     *  Before revision 6516, we used to zero out i-th row and column here. It seems to be redundant because they are already zero after assembly.
     *  When assembling a bath element you get a matrix subblock that looks like (2D example):
     *
     *     Vm   0 0 0 0 0 0
     *     Vm   0 0 0 0 0 0
     *     Vm   0 0 0 0 0 0
     *     Phie 0 0 0 x x x
     *     Phie 0 0 0 x x x  -> the x subblock is assembled from div(grad_phi) = 0
     *     Phie 0 0 0 x x x
     *
     *  Therefore, all the Vm entries of this node are already 0.
     *
     *  Explicitly checking it in non-production builds.
     */
#ifndef NDEBUG
    if (computeMatrix)
    {
        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter=this->mpMesh->GetNodeIteratorBegin();
             iter != this->mpMesh->GetNodeIteratorEnd();
             ++iter)
        {
            if (HeartRegionCode::IsRegionBath( (*iter).GetRegion() ))
            {
                int num_equation = 2*iter->GetIndex(); // assumes Vm and Phie are interleaved

                PetscInt local_lo, local_hi;
                this->mpLinearSystem->GetOwnershipRange(local_lo, local_hi);

                // If this processor owns i-th row, check it.
                if ((local_lo <= (int)num_equation) && ((int)num_equation < local_hi))
                {
                    for (unsigned column=0; column < this->mpLinearSystem->GetSize(); column++)
                    {
                        assert(this->mpLinearSystem->GetMatrixElement(num_equation, column)==0.0);
                    }
                }

                // Check the local entries of the i-th column
                for (int row=local_lo; row<local_hi; row++)
                {
                    assert(this->mpLinearSystem->GetMatrixElement(row, num_equation)==0);
                }
            }
        }
    }
#endif

    /*
     *  These two loops are decoupled because interleaving calls to GetMatrixElement and MatSetValue
     *  require reassembling the matrix before GetMatrixElement which generates massive communication
     *  overhead for large models and/or large core counts.
     */

    for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::NodeIterator iter=this->mpMesh->GetNodeIteratorBegin();
         iter != this->mpMesh->GetNodeIteratorEnd();
         ++iter)
    {
        if (HeartRegionCode::IsRegionBath((*iter).GetRegion() ))
        {
            PetscInt index = 2*iter->GetIndex(); // assumes Vm and Phie are interleaved

            if (computeMatrix)
            {
                // put 1.0 on the diagonal
                PetscMatTools::SetElement(this->mpLinearSystem->rGetLhsMatrix(), index, index, 1.0);
            }

            if (computeVector)
            {
                // zero rhs vector entry
                PetscVecTools::SetElement(this->mpLinearSystem->rGetRhsVector(), index, 0.0);
            }
        }
    }
}

// Explicit instantiation
template class AbstractBidomainSolver<1,1>;
template class AbstractBidomainSolver<2,2>;
template class AbstractBidomainSolver<3,3>;
