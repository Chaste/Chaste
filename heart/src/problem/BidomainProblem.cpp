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


#include "BidomainProblem.hpp"
#include "BidomainSolver.hpp"
#include "HeartConfig.hpp"
#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"

template <unsigned DIM>
void BidomainProblem<DIM>::AnalyseMeshForBath()
{
    // Annotate bath notes with correct region code
    if (mHasBath)
    {
        // Initialize all nodes to be bath nodes
        for (typename AbstractTetrahedralMesh<DIM,DIM>::NodeIterator iter=this->mpMesh->GetNodeIteratorBegin();
             iter != this->mpMesh->GetNodeIteratorEnd();
            ++iter)
        {
            (*iter).SetRegion(HeartRegionCode::GetValidBathId());
        }

        bool any_bath_element_found = false;

        // Set nodes that are part of a heart element to be heart nodes
        //for (unsigned i=0; i<this->mpMesh->GetNumElements(); i++)
        for (typename AbstractTetrahedralMesh<DIM,DIM>::ElementIterator it = this->mpMesh->GetElementIteratorBegin();
             it != this->mpMesh->GetElementIteratorEnd();
             ++it)
        {
            Element<DIM, DIM>& r_element = *it;

            if (HeartRegionCode::IsRegionTissue( r_element.GetUnsignedAttribute() ))
            {
                for (unsigned j=0; j<r_element.GetNumNodes(); j++)
                {
                    r_element.GetNode(j)->SetRegion(HeartRegionCode::GetValidTissueId());
                }
            }
            else
            {
                assert(HeartRegionCode::IsRegionBath( r_element.GetUnsignedAttribute() ));
                any_bath_element_found = true;
            }
        }

        if (!PetscTools::ReplicateBool(any_bath_element_found))
        {
           EXCEPTION("No bath element found");
        }
    }
    else
    {
        // The problem hasn't been set up with a bath, so check that the user hasn't set any options
        // which would suggest they're expecting there to be a bath!
        std::set<unsigned> bath_identifiers = HeartConfig::Instance()->rGetBathIdentifiers();
        if (!(bath_identifiers.size()==1 && bath_identifiers.find(1)==bath_identifiers.begin())) // IF NOT only 1 in the bath identifiers set
        {
            EXCEPTION("User has set bath identifiers, but the BidomainProblem isn't expecting a bath. Did you mean to use BidomainProblem(..., true)? Or alternatively, BidomainWithBathProblem(...)?");
        }
    }
}

template<unsigned DIM>
Vec BidomainProblem<DIM>::CreateInitialCondition()
{
    Vec init_cond = AbstractCardiacProblem<DIM,DIM,2>::CreateInitialCondition();
    if (mHasBath)
    {
        // get the voltage stripe
        DistributedVector ic = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(init_cond);
        DistributedVector::Stripe voltage_stripe = DistributedVector::Stripe(ic,0);

        for (DistributedVector::Iterator index = ic.Begin();
             index!= ic.End();
             ++index)
        {
            if (HeartRegionCode::IsRegionBath( this->mpMesh->GetNode( index.Global )->GetRegion() ))
            {
                voltage_stripe[index] = 0.0;
            }
        }
        ic.Restore();
    }

    return init_cond;
}

template<unsigned DIM>
AbstractCardiacTissue<DIM> * BidomainProblem<DIM>::CreateCardiacTissue()
{
    AnalyseMeshForBath();
    mpBidomainTissue = new BidomainTissue<DIM>(this->mpCellFactory, HeartConfig::Instance()->GetUseStateVariableInterpolation());
    return mpBidomainTissue;
}

template<unsigned DIM>
AbstractDynamicLinearPdeSolver<DIM, DIM, 2>* BidomainProblem<DIM>::CreateSolver()
{
    /*
     * NOTE: The this->mpBoundaryConditionsContainer.get() lines below convert a
     * boost::shared_ptr to a normal pointer, as this is what the solvers are
     * expecting. We have to be a bit careful though as boost could decide to delete
     * them whenever it feels like as it won't count the assembers as using them.
     *
     * As long as they are kept as member variables here for as long as they are
     * required in the solvers it should all work OK.
     */
    mpSolver = new BidomainSolver<DIM,DIM>(mHasBath,
                                           this->mpMesh,
                                           mpBidomainTissue,
                                           this->mpBoundaryConditionsContainer.get());

    try
    {
        mpSolver->SetFixedExtracellularPotentialNodes(mFixedExtracellularPotentialNodes);
        mpSolver->SetRowForAverageOfPhiZeroed(mRowForAverageOfPhiZeroed);
    }
    catch (const Exception& e)
    {
        delete mpSolver;
        throw e;
    }

    return mpSolver;
}

template<unsigned DIM>
BidomainProblem<DIM>::BidomainProblem(
            AbstractCardiacCellFactory<DIM>* pCellFactory, bool hasBath)
    : AbstractCardiacProblem<DIM,DIM, 2>(pCellFactory),
      mpBidomainTissue(NULL),
      mRowForAverageOfPhiZeroed(INT_MAX),
      mHasBath(hasBath)
{
    mFixedExtracellularPotentialNodes.resize(0);
}

template<unsigned DIM>
BidomainProblem<DIM>::BidomainProblem()
    : AbstractCardiacProblem<DIM, DIM, 2>(),
      mpBidomainTissue(NULL),
      mRowForAverageOfPhiZeroed(INT_MAX)
{
    mFixedExtracellularPotentialNodes.resize(0);
}

template<unsigned DIM>
void BidomainProblem<DIM>::SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes)
{
    mFixedExtracellularPotentialNodes.resize(nodes.size());
    for (unsigned i=0; i<nodes.size(); i++)
    {
        // the solver checks that the nodes[i] is less than
        // the number of nodes in the mesh so this is not done here
        mFixedExtracellularPotentialNodes[i] = nodes[i];
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::SetNodeForAverageOfPhiZeroed(unsigned node)
{
    mRowForAverageOfPhiZeroed = 2*node+1;
}

template<unsigned DIM>
BidomainTissue<DIM>* BidomainProblem<DIM>::GetBidomainTissue()
{
    assert(mpBidomainTissue!=NULL);
    return mpBidomainTissue;
}

template<unsigned DIM>
void BidomainProblem<DIM>::WriteInfo(double time)
{
    if (PetscTools::AmMaster())
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
    }

    double v_max, v_min, phi_max, phi_min;

    VecStrideMax( this->mSolution, 0, PETSC_NULL, &v_max );
    VecStrideMin( this->mSolution, 0, PETSC_NULL, &v_min );

    VecStrideMax( this->mSolution, 1, PETSC_NULL, &phi_max );
    VecStrideMin( this->mSolution, 1, PETSC_NULL, &phi_min );

    if (PetscTools::AmMaster())
    {
        std::cout << " V; phi_e = " << "[" <<v_min << ", " << v_max << "]" << ";\t"
                  << "[" <<phi_min << ", " << phi_max << "]" << "\n"
                  << std::flush;
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::DefineWriterColumns(bool extending)
{
    AbstractCardiacProblem<DIM,DIM,2>::DefineWriterColumns(extending);
    if (extending)
    {
        mExtracelluarColumnId = this->mpWriter->GetVariableByName("Phi_e");
    }
    else
    {
        mExtracelluarColumnId = this->mpWriter->DefineVariable("Phi_e","mV");
    }
    AbstractCardiacProblem<DIM,DIM,2>::DefineExtraVariablesWriterColumns(extending);
}

template<unsigned DIM>
void BidomainProblem<DIM>::WriteOneStep(double time, Vec voltageVec)
{
    this->mpWriter->PutUnlimitedVariable(time);
    std::vector<int> variable_ids;
    variable_ids.push_back(this->mVoltageColumnId);
    variable_ids.push_back(mExtracelluarColumnId);
    this->mpWriter->PutStripedVector(variable_ids, voltageVec);
    AbstractCardiacProblem<DIM,DIM,2>::WriteExtraVariablesOneStep();
}

template<unsigned DIM>
void BidomainProblem<DIM>::PreSolveChecks()
{
    AbstractCardiacProblem<DIM,DIM, 2>::PreSolveChecks();
    if (mFixedExtracellularPotentialNodes.empty())
    {
        // We're not pinning any nodes.
        if (mRowForAverageOfPhiZeroed==INT_MAX)
        {
            // We're not using the constrain Average phi_e to 0 method, hence use a null space
            // Check that the KSP solver isn't going to do anything stupid:
            // phi_e is not bounded, so it'd be wrong to use a relative tolerance
            if (HeartConfig::Instance()->GetUseRelativeTolerance())
            {
                EXCEPTION("Bidomain external voltage is not bounded in this simulation - use KSP *absolute* tolerance");
            }
        }
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::SetElectrodes()
{
    if (!mHasBath)
    {
        //Cannot set electrodes when problem has been defined to not have a bath
        return;
    }

    assert(this->mpMesh!=NULL);

    if (HeartConfig::Instance()->IsElectrodesPresent())
    {
        mpElectrodes = (boost::shared_ptr<Electrodes<DIM> >) new Electrodes<DIM>(*(this->mpMesh));
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::AtBeginningOfTimestep(double time)
{
    if (mpElectrodes && mpElectrodes->SwitchOn(time))
    {
        // At the moment mpBcc and mpDefaultBcc point to a set default BC
        assert(this->mpBoundaryConditionsContainer);
        //assert(this->mpDefaultBoundaryConditionsContainer);

        // Note, no point calling this->SetBoundaryConditionsContainer() as the
        // solver has already been created..
        mpSolver->ResetBoundaryConditionsContainer(mpElectrodes->GetBoundaryConditionsContainer().get());

        // ..but we set mpBcc anyway, so the local mpBcc is
        // the same as the one being used in the solver...
        this->mpBoundaryConditionsContainer = mpElectrodes->GetBoundaryConditionsContainer();

        /// \todo #1159 #1324 heart/src/problem/AbstractCardiacProblem.hpp:657 expects both pointing at the same place when unarchiving
        this->mpDefaultBoundaryConditionsContainer = this->mpBoundaryConditionsContainer;

        // At t==0 or after checkpointing we won't have a system assembled at this stage: BCs will be applied once the matrix
        // is assembled. Dirichlet BCs will be present at the time of assembly and no null space will be created either.
        if (mpSolver->GetLinearSystem() != NULL)
        {
            // System matrix is assembled once at the beginning of the simulation. After that, nobody will take care
            // of applying new BC to the system matrix. Must be triggered explicitly.
            if (mpElectrodes->HasGroundedElectrode())
            {
                this->mpBoundaryConditionsContainer->ApplyDirichletToLinearProblem( *(mpSolver->GetLinearSystem()),
                                                                                   true, false);
            }

            // If a grounded electrode is switched on, the linear system is not singular anymore. Remove the null space.
            if (mpElectrodes->HasGroundedElectrode())
            {
                mpSolver->GetLinearSystem()->RemoveNullSpace();
            }
        }
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::OnEndOfTimestep(double time)
{
    if (mpElectrodes && mpElectrodes->SwitchOff(time))
    {
        // At the moment mpBcc should exist and therefore
        // mpDefaultBcc should be empty (not if electrodes switched on after 0ms)
        assert(this->mpBoundaryConditionsContainer);
        //assert(! this->mpDefaultBoundaryConditionsContainer);

        // Set up default boundary conditions container - no Neumann fluxes
        // or Dirichlet fixed nodes
        this->mpDefaultBoundaryConditionsContainer.reset(new BoundaryConditionsContainer<DIM,DIM,2>);
        for (unsigned problem_index=0; problem_index<2; problem_index++)
        {
            this->mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(this->mpMesh, problem_index);
        }

        // If there's a grounded electrode, we must remove BC from linear system. At the moment, we don't
        // have a sensible way of doing this, therefore we reassemble the system.
        if (mpElectrodes->HasGroundedElectrode())
        {
            delete mpSolver;
            AbstractCardiacProblem<DIM,DIM,2>::mpSolver = CreateSolver();
            mpSolver->SetTimeStep(HeartConfig::Instance()->GetPdeTimeStep());
        }

        // Note, no point calling this->SetBoundaryConditionsContainer() as the
        // solver has already been created..
        mpSolver->ResetBoundaryConditionsContainer(this->mpDefaultBoundaryConditionsContainer.get());
        // ..but we set mpBcc to be mpDefaultBcc anyway, so the local mpBcc is
        // the same as the one being used in the solver...
        this->mpBoundaryConditionsContainer = this->mpDefaultBoundaryConditionsContainer;
    }
}

template<unsigned DIM>
void BidomainProblem<DIM>::SetUpAdditionalStoppingTimes(std::vector<double>& rAdditionalStoppingTimes)
{
    if (mpElectrodes)
    {
        rAdditionalStoppingTimes.push_back(mpElectrodes->GetSwitchOnTime());
        rAdditionalStoppingTimes.push_back(mpElectrodes->GetSwitchOffTime());
    }
}

template<unsigned DIM>
bool BidomainProblem<DIM>::GetHasBath()
{
    return mHasBath;
}

// Explicit instantiation
template class BidomainProblem<1>;
template class BidomainProblem<2>;
template class BidomainProblem<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BidomainProblem)
