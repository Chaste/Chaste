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


#include "ExtendedBidomainProblem.hpp"
#include "ExtendedBidomainSolver.hpp"
#include "AbstractDynamicLinearPdeSolver.hpp"
#include "HeartConfig.hpp"
#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
ExtendedBidomainProblem<DIM>::ExtendedBidomainProblem(
            AbstractCardiacCellFactory<DIM>* pCellFactory, AbstractCardiacCellFactory<DIM>* pSecondCellFactory, bool hasBath)
    : AbstractCardiacProblem<DIM,DIM, 3>(pCellFactory),
      mpSecondCellFactory(pSecondCellFactory),
      mpExtendedBidomainTissue(NULL),
      mUserSpecifiedSecondCellConductivities(false),
      mUserHasSetBidomainValuesExplicitly(false),
      mpExtracellularStimulusFactory(NULL),
      mRowForAverageOfPhiZeroed(INT_MAX),
      mApplyAveragePhieZeroConstraintAfterSolving(false),
      mUserSuppliedExtracellularStimulus(false),
      mHasBath(hasBath)
{
    mFixedExtracellularPotentialNodes.resize(0);
}

template<unsigned DIM>
ExtendedBidomainProblem<DIM>::ExtendedBidomainProblem()
    : AbstractCardiacProblem<DIM,DIM, 3>(),
      mpSecondCellFactory(NULL),
      mpExtendedBidomainTissue(NULL),
      mUserSpecifiedSecondCellConductivities(false),
      mUserHasSetBidomainValuesExplicitly(false),
      mpExtracellularStimulusFactory(NULL),
      mRowForAverageOfPhiZeroed(INT_MAX),
      mApplyAveragePhieZeroConstraintAfterSolving(false),
      mUserSuppliedExtracellularStimulus(false)
{
    mFixedExtracellularPotentialNodes.resize(0);
}

template<unsigned DIM>
Vec ExtendedBidomainProblem<DIM>::CreateInitialCondition()
{

    DistributedVectorFactory* p_factory = this->mpMesh->GetDistributedVectorFactory();
    Vec initial_condition = p_factory->CreateVec(3);
    DistributedVector ic = p_factory->CreateDistributedVector(initial_condition);
    std::vector<DistributedVector::Stripe> stripe;
    stripe.reserve(3);

    stripe.push_back(DistributedVector::Stripe(ic, 0));
    stripe.push_back(DistributedVector::Stripe(ic, 1));
    stripe.push_back(DistributedVector::Stripe(ic, 2));

    for (DistributedVector::Iterator index = ic.Begin();
         index != ic.End();
         ++index)
    {
        stripe[0][index] = mpExtendedBidomainTissue->GetCardiacCell(index.Global)->GetVoltage();//phi_i of frrst cell = Vm first cell at the start
        stripe[1][index] = mpExtendedBidomainTissue->GetCardiacSecondCell(index.Global)->GetVoltage();//phi_i of second cell = Vm second cell at the start
        stripe[2][index] = 0.0;//
    }
    ic.Restore();

    return initial_condition;
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::ProcessExtracellularStimulus()
{
    if (mpExtracellularStimulusFactory == NULL) // user has not passed in any extracellular stimulus in any form
    {
        mpExtracellularStimulusFactory = new AbstractStimulusFactory<DIM>();
        // Create one (with default implementation to zero stimulus everywhere)
    }

    assert(mpExtracellularStimulusFactory); // should be created by now, either above or by the user...
    mpExtracellularStimulusFactory->SetMesh(this->mpMesh);//so, set the mesh into it.
    mpExtracellularStimulusFactory->SetCompatibleExtracellularStimulus();//make sure compatibility condition will be valid

    std::vector<AbstractChasteRegion<DIM>* > grounded_regions = mpExtracellularStimulusFactory->GetRegionsToBeGrounded();

    if ((mUserSuppliedExtracellularStimulus) && grounded_regions.size() > 0) //we check for grunded nodes here
    {
        std::vector<unsigned> grounded_indices;
        for (unsigned global_node_index = 0; global_node_index < this->mpMesh->GetNumNodes(); global_node_index++)
        {
            if (this->mpMesh->GetDistributedVectorFactory()->IsGlobalIndexLocal(global_node_index))
            {
                for (unsigned region_index = 0; region_index <grounded_regions.size(); region_index++)
                {
                    if (grounded_regions[region_index]->DoesContain(this->mpMesh->GetNode(global_node_index)->GetPoint()))
                    {
                        grounded_indices.push_back( this->mpMesh->GetNode(global_node_index)->GetIndex() );
                    }
                }
            }
        }
        PetscTools::Barrier();
        SetFixedExtracellularPotentialNodes(grounded_indices);
    }
}

template<unsigned DIM>
AbstractCardiacTissue<DIM> * ExtendedBidomainProblem<DIM>::CreateCardiacTissue()
{
    //set the mesh into the second cell factory as well.
    mpSecondCellFactory->SetMesh(this->mpMesh);

    //deal with extracellular stimulus, if any
    ProcessExtracellularStimulus();

    //Now create the tissue object
    mpExtendedBidomainTissue = new ExtendedBidomainTissue<DIM>(this->mpCellFactory, mpSecondCellFactory,mpExtracellularStimulusFactory);

    //Let the Tissue know if the user wants an extracellular stimulus (or if we had to create a default zero one).
    mpExtendedBidomainTissue->SetUserSuppliedExtracellularStimulus(mUserSuppliedExtracellularStimulus);

    //if the user remembered to set a different value for the sigma of the second cell...
    if (mUserSpecifiedSecondCellConductivities)
    {
        mpExtendedBidomainTissue->SetIntracellularConductivitiesSecondCell(mIntracellularConductivitiesSecondCell);
    }
    else //..otherwise it gets the same as the first cell (according to heartconfig...)
    {
        c_vector<double, DIM> intra_conductivities;
        HeartConfig::Instance()->GetIntracellularConductivities(intra_conductivities);
        mpExtendedBidomainTissue->SetIntracellularConductivitiesSecondCell(intra_conductivities);
    }

    //the conductivities for the first cell are created within the tissue constructor in the abstract class
    //here we create the ones for the second cell
    mpExtendedBidomainTissue->CreateIntracellularConductivityTensorSecondCell();

    if (mUserHasSetBidomainValuesExplicitly)
    {
        mpExtendedBidomainTissue->SetAmFirstCell(mAmFirstCell);
        mpExtendedBidomainTissue->SetAmSecondCell(mAmSecondCell);
        mpExtendedBidomainTissue->SetAmGap(mAmGap);
        mpExtendedBidomainTissue->SetGGap(mGGap);
        mpExtendedBidomainTissue->SetCmFirstCell(mCmFirstCell);
        mpExtendedBidomainTissue->SetCmSecondCell(mCmSecondCell);
    }
    else//we set all the Am and Cm to the values set by the heartconfig (only one value for all Am and one value for all Cms)
    {
        mpExtendedBidomainTissue->SetAmFirstCell(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        mpExtendedBidomainTissue->SetAmSecondCell(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        mpExtendedBidomainTissue->SetAmGap(HeartConfig::Instance()->GetSurfaceAreaToVolumeRatio());
        mpExtendedBidomainTissue->SetGGap(0.0);
        mpExtendedBidomainTissue->SetCmFirstCell(HeartConfig::Instance()->GetCapacitance());
        mpExtendedBidomainTissue->SetCmSecondCell(HeartConfig::Instance()->GetCapacitance());
    }

    mpExtendedBidomainTissue->SetGgapHeterogeneities(mGgapHeterogeneityRegions, mGgapHeterogenousValues);//set user input into the tissue class
    mpExtendedBidomainTissue->CreateGGapConductivities();//if vectors are empty, mGgap will be put everywhere by this method.

    return mpExtendedBidomainTissue;
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::SetExtendedBidomainParameters(double Am1, double Am2, double AmGap, double Cm1, double Cm2, double Ggap)
{
     mAmFirstCell = Am1;
     mAmSecondCell = Am2;
     mAmGap = AmGap;
     mCmFirstCell = Cm1;
     mCmSecondCell = Cm2;
     mGGap = Ggap;

     mUserHasSetBidomainValuesExplicitly = true;
}

template <unsigned DIM>
void ExtendedBidomainProblem<DIM>::SetGgapHeterogeneities ( std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rGgapHeterogeneityRegions, std::vector<double>& rGgapValues)
{
    if (rGgapHeterogeneityRegions.size() != rGgapValues.size() )
    {
        EXCEPTION  ("Gap junction heterogeneity areas must be of the same number as the heterogeneity values");
    }
    mGgapHeterogeneityRegions = rGgapHeterogeneityRegions;
    mGgapHeterogenousValues =rGgapValues;
}

template <unsigned DIM>
void ExtendedBidomainProblem<DIM>::SetExtracellularStimulusFactory( AbstractStimulusFactory<DIM>* pFactory)
{
    mpExtracellularStimulusFactory = pFactory;
    mUserSuppliedExtracellularStimulus = true;
}

template<unsigned DIM>
AbstractDynamicLinearPdeSolver<DIM, DIM, 3>* ExtendedBidomainProblem<DIM>::CreateSolver()
{
    /*
     * NOTE: The this->mpBoundaryConditionsContainer.get() lines below convert a
     * boost::shared_ptr to a normal pointer, as this is what the assemblers are
     * expecting. We have to be a bit careful though as boost could decide to delete
     * them whenever it feels like as it won't count the assemblers as using them.
     *
     * As long as they are kept as member variables here for as long as they are
     * required in the assemblers it should all work OK.
     */

    mpSolver = new ExtendedBidomainSolver<DIM,DIM>( mHasBath,
                                                    this->mpMesh,
                                                    mpExtendedBidomainTissue,
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
ExtendedBidomainProblem<DIM>::~ExtendedBidomainProblem()
{
    if (!mUserSuppliedExtracellularStimulus)
    {
        delete mpExtracellularStimulusFactory;
    }
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::SetIntracellularConductivitiesForSecondCell(c_vector<double, DIM> conductivities)
{
    for (unsigned i = 0; i < DIM; i++)
    {
        mIntracellularConductivitiesSecondCell[i] = conductivities[i];
    }
    mUserSpecifiedSecondCellConductivities = true;
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::SetFixedExtracellularPotentialNodes(std::vector<unsigned> nodes)
{
    assert(mFixedExtracellularPotentialNodes.size() == 0); ///\todo turn this into an exception if the user calls this twice...
    mFixedExtracellularPotentialNodes.resize(nodes.size());
    for (unsigned i=0; i<nodes.size(); i++)
    {
        // the assembler checks that the nodes[i] is less than
        // the number of nodes in the mesh so this is not done here
        mFixedExtracellularPotentialNodes[i] = nodes[i];
    }
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::SetNodeForAverageOfPhiZeroed(unsigned node)
{
    if (node==0)
    {
        mRowForAverageOfPhiZeroed = 2;
    }
    else
    {
        //Phie is every three lines, starting from zero...
        mRowForAverageOfPhiZeroed = 3*node  - 1;
    }
}

template<unsigned DIM>
ExtendedBidomainTissue<DIM>* ExtendedBidomainProblem<DIM>::GetExtendedBidomainTissue()
{
    assert(mpExtendedBidomainTissue!=NULL);
    return mpExtendedBidomainTissue;
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::WriteInfo(double time)
{
    if (PetscTools::AmMaster())
    {
        std::cout << "Solved to time " << time << "\n" << std::flush;
    }

    double V_max_first_cell, V_min_first_cell, V_max_second_cell, V_min_second_cell, phi_e_min, phi_e_max;

    VecStrideMax( this->mSolution, 0, PETSC_NULL, &V_max_first_cell );
    VecStrideMin( this->mSolution, 0, PETSC_NULL, &V_min_first_cell );

    VecStrideMax( this->mSolution, 1, PETSC_NULL, &V_max_second_cell );
    VecStrideMin( this->mSolution, 1, PETSC_NULL, &V_min_second_cell );

    VecStrideMax( this->mSolution, 2, PETSC_NULL, &phi_e_max );
    VecStrideMin( this->mSolution, 2, PETSC_NULL, &phi_e_min );

    if (PetscTools::AmMaster())
    {
        std::cout << " V first cell = "  << "[" <<V_min_first_cell << ",   " << V_max_first_cell << "]" << ";\n"
                  << " V second cell = " << "[" <<V_min_second_cell << ", " << V_max_second_cell << "]" << ";\n"
                  << " Phi_e = "         << "[" <<phi_e_min << ", "         << phi_e_max << "]" << ";\n"
                  << std::flush;
    }
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::DefineWriterColumns(bool extending)
{
    if (!extending)
    {
        if (this->mNodesToOutput.empty())
        {
            // Set writer to output all nodes
            this->mpWriter->DefineFixedDimension(this->mpMesh->GetNumNodes());
        }
//        else
//        {
//            // Output only the nodes indicted
//            this->mpWriter->DefineFixedDimension( this->mNodesToOutput, this->mpMesh->GetNumNodes() );
//        }
        //mNodeColumnId = mpWriter->DefineVariable("Node", "dimensionless");
        mVoltageColumnId_Vm1 = this->mpWriter->DefineVariable("V","mV");
        mVoltageColumnId_Vm2 = this->mpWriter->DefineVariable("V_2","mV");
        mVoltageColumnId_Phie = this->mpWriter->DefineVariable("Phi_e","mV");
        mVariablesIDs.push_back(mVoltageColumnId_Vm1);
        mVariablesIDs.push_back(mVoltageColumnId_Vm2);
        mVariablesIDs.push_back(mVoltageColumnId_Phie);

        // Only used to get an estimate of the # of timesteps below (copied from Abstract class)
        TimeStepper stepper(AbstractCardiacProblem<DIM,DIM,3>::mCurrentTime,
                            HeartConfig::Instance()->GetSimulationDuration(),
                            HeartConfig::Instance()->GetPrintingTimeStep());
        this->mpWriter->DefineUnlimitedDimension("Time", "msecs", stepper.EstimateTimeSteps()+1); // +1 for start and end
    }
    else
    {
        mVoltageColumnId_Vm1 = this->mpWriter->GetVariableByName("V");
        mVoltageColumnId_Vm2 = this->mpWriter->GetVariableByName("V_2");
        mVoltageColumnId_Phie = this->mpWriter->GetVariableByName("Phi_e");
    }
    //define any extra variable. NOTE: it must be in the first cell (not the second)
    AbstractCardiacProblem<DIM,DIM,3>::DefineExtraVariablesWriterColumns(extending);

}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::WriteOneStep(double time, Vec voltageVec)
{
    this->mpWriter->PutUnlimitedVariable(time);

   // Create a striped vector
    Vec ordered_voltages =  this->mpMesh->GetDistributedVectorFactory()->CreateVec(3);
    DistributedVector wrapped_ordered_solution = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(ordered_voltages);
    DistributedVector wrapped_solution = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(voltageVec);

    DistributedVector::Stripe V_first_cell_stripe(wrapped_solution,0);
    DistributedVector::Stripe V_second_cell_stripe(wrapped_solution,1);
    DistributedVector::Stripe phi_e_stripe(wrapped_solution,2);

    DistributedVector::Stripe wrapped_ordered_solution_first_stripe(wrapped_ordered_solution,0);
    DistributedVector::Stripe wrapped_ordered_solution_second_stripe(wrapped_ordered_solution,1);
    DistributedVector::Stripe wrapped_ordered_solution_third_stripe(wrapped_ordered_solution,2);

    for (DistributedVector::Iterator index = wrapped_solution.Begin();
         index != wrapped_solution.End();
         ++index)
    {
        wrapped_ordered_solution_first_stripe[index] = V_first_cell_stripe[index];
        wrapped_ordered_solution_second_stripe[index] = V_second_cell_stripe[index];
        wrapped_ordered_solution_third_stripe[index] = phi_e_stripe[index];
    }
    wrapped_solution.Restore();
    wrapped_ordered_solution.Restore();

    this->mpWriter->PutStripedVector(mVariablesIDs, ordered_voltages);
    PetscTools::Destroy(ordered_voltages);
    //write any extra variable. Note that this method in the parent class will
    //take the extra variable only from the first cell.
    ///\todo write a specific method for this class
    AbstractCardiacProblem<DIM,DIM,3>::WriteExtraVariablesOneStep();
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::PreSolveChecks()
{
    AbstractCardiacProblem<DIM,DIM, 3>::PreSolveChecks();
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
bool ExtendedBidomainProblem<DIM>::GetHasBath()
{
    return mHasBath;
}

template<unsigned DIM>
void ExtendedBidomainProblem<DIM>::SetHasBath(bool hasBath)
{
    mHasBath = hasBath;
}

// Explicit instantiation
template class ExtendedBidomainProblem<1>;
template class ExtendedBidomainProblem<2>;
template class ExtendedBidomainProblem<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtendedBidomainProblem)
