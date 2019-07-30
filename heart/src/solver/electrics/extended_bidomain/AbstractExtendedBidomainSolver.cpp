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

#include "AbstractExtendedBidomainSolver.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractExtendedBidomainSolver<ELEMENT_DIM, SPACE_DIM>::InitialiseForSolve(Vec initialSolution)
{
    // The base class method that calls this function will only call it with a null linear system
    assert(this->mpLinearSystem == NULL);

    // linear system created here
    AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, 3>::InitialiseForSolve(initialSolution);

    if (HeartConfig::Instance()->GetUseAbsoluteTolerance())
    {
#ifdef TRACE_KSP
        std::cout << "Using absolute tolerance: " << mpConfig->GetAbsoluteTolerance() << "\n";
#endif
        this->mpLinearSystem->SetAbsoluteTolerance(mpConfig->GetAbsoluteTolerance());
    }
    else
    {
#ifdef TRACE_KSP
        std::cout << "Using relative tolerance: " << mpConfig->GetRelativeTolerance() << "\n";
#endif
        this->mpLinearSystem->SetRelativeTolerance(mpConfig->GetRelativeTolerance());
    }

    this->mpLinearSystem->SetKspType(HeartConfig::Instance()->GetKSPSolver());
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
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::PrepareForSetupLinearSystem(Vec existingSolution)
{
    mpExtendedBidomainTissue->SolveCellSystems(existingSolution, PdeSimulationTime::GetTime(), PdeSimulationTime::GetNextTime());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Vec AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::GenerateNullBasis() const
{
    double sqrrt_num_nodes = sqrt(  this->mpMesh->GetNumNodes());
    double normalised_basis_value = 1.0 / sqrrt_num_nodes;

    Vec nullbasis;
    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    nullbasis=p_factory->CreateVec(3);
    DistributedVector dist_null_basis = p_factory->CreateDistributedVector(nullbasis);
    DistributedVector::Stripe null_basis_stripe_0(dist_null_basis,0);
    DistributedVector::Stripe null_basis_stripe_1(dist_null_basis,1);
    DistributedVector::Stripe null_basis_stripe_2(dist_null_basis,2);
    for (DistributedVector::Iterator index = dist_null_basis.Begin();
         index != dist_null_basis.End();
         ++index)
    {
        null_basis_stripe_0[index] = 0.0;
        null_basis_stripe_1[index] = 0.0;
        null_basis_stripe_2[index] = normalised_basis_value;
    }
    dist_null_basis.Restore();
    return nullbasis;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::FinaliseLinearSystem(Vec existingSolution)
{
    if (!(this->mpBoundaryConditions->HasDirichletBoundaryConditions()))
    {
        // We're not pinning any nodes.
        if (mRowForAverageOfPhiZeroed==INT_MAX)
        {
            // We're not using the 'Average phi_e = 0' method, hence use a null space
            if (!mNullSpaceCreated)
            {
                // No null space set up, so create one and pass it to the linear system
                Vec nullbasis[] = {GenerateNullBasis()};

                this->mpLinearSystem->SetNullBasis(nullbasis, 1);

                PetscTools::Destroy(nullbasis[0]);
                mNullSpaceCreated = true;
            }
        }
        else  // mRowForAverageOfPhiZeroed!=INT_MAX, i.e. we're using the 'Average phi_e = 0' method
        {
            // CG (default solver) won't work since the system isn't symmetric anymore. Switch to GMRES
            this->mpLinearSystem->SetKspType("gmres"); // Switches the solver
            mpConfig->SetKSPSolver("gmres", true); // Makes sure this change will be reflected in the XML file written to disk at the end of the simulation.
            //(If the user doesn't have gmres then the "true" warns the user about the switch)

            // Set average phi_e to zero
            unsigned matrix_size = this->mpLinearSystem->GetSize();
            if (!this->mMatrixIsAssembled)
            {

                // Set the mRowForAverageOfPhiZeroed-th matrix row to 0 0 1 0 0 1 ...
                std::vector<unsigned> row_for_average;
                row_for_average.push_back(mRowForAverageOfPhiZeroed);
                this->mpLinearSystem->ZeroMatrixRowsWithValueOnDiagonal(row_for_average, 0.0);
                for (unsigned col_index=0; col_index<matrix_size; col_index++)
                {
                    if (((col_index-2)%3 == 0) && (col_index>1))//phi_e column indices are 2,5,8,11 etc....
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
    CheckCompatibilityCondition();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::CheckCompatibilityCondition()
{
    if (this->mpBoundaryConditions->HasDirichletBoundaryConditions() || mRowForAverageOfPhiZeroed!=INT_MAX )
    {
        // not a singular system, no compability condition
        return;
    }

#ifndef NDEBUG
    DistributedVector distributed_rhs=this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(this->mpLinearSystem->rGetRhsVector());
    DistributedVector::Stripe distributed_rhs_phi_entries(distributed_rhs,2); // stripe number 2 -> phi_e

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
        std::cout << "Sum of b_{every 3 items} = " << global_sum << " (should be zero for compatibility)\n";
        EXCEPTION("Linear system does not satisfy compatibility constraint!");
        // LCOV_EXCL_STOP
    }
#endif
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::AbstractExtendedBidomainSolver(
        bool bathSimulation,
        AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
        ExtendedBidomainTissue<SPACE_DIM>* pTissue,
        BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 3>* pBcc)
    : AbstractDynamicLinearPdeSolver<ELEMENT_DIM,SPACE_DIM,3>(pMesh),
              mBathSimulation(bathSimulation),
              mpExtendedBidomainTissue(pTissue),
              mpBoundaryConditions(pBcc)
{
    assert(pTissue != NULL);
    assert(pBcc != NULL);

    mNullSpaceCreated = false;

    // important!
    this->mMatrixIsConstant = true;

    mRowForAverageOfPhiZeroed = INT_MAX; //this->mpLinearSystem->GetSize() - 1;
    mpConfig = HeartConfig::Instance();

    mpExtendedBidomainAssembler = NULL; // can't initialise until know what dt is
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::~AbstractExtendedBidomainSolver()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::SetFixedExtracellularPotentialNodes(
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
    this->mpBoundaryConditions->ResetDirichletCommunication();

    for (unsigned i=0; i<mFixedExtracellularPotentialNodes.size(); i++)
    {
        if (this->mpMesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(mFixedExtracellularPotentialNodes[i]))
        {
            ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(0.0);

            //Throws if node is not owned locally
            Node<SPACE_DIM>* p_node = this->mpMesh->GetNode(mFixedExtracellularPotentialNodes[i]);

            //the "false" passed in tells not to check that it is a boundary node (since we might have grounded electrodes within the tissue)
            GetBoundaryConditions()->AddDirichletBoundaryCondition(p_node, p_boundary_condition, 2, false);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractExtendedBidomainSolver<ELEMENT_DIM,SPACE_DIM>::SetRowForAverageOfPhiZeroed(unsigned row)
{
    // Row should be every 3 entries, starting from zero...
    if (((row-2)%3 != 0) && row!=INT_MAX)
    {
        EXCEPTION("Row for applying the constraint 'Average of phi_e = zero' should be every 3 rows");
    }

    mRowForAverageOfPhiZeroed = row;
}

// Explicit instantiation
template class AbstractExtendedBidomainSolver<1,1>;
template class AbstractExtendedBidomainSolver<2,2>;
template class AbstractExtendedBidomainSolver<3,3>;
