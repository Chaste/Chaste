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

#include "AbstractCardiacTissue.hpp"

#include <boost/scoped_array.hpp>

#include "DistributedVector.hpp"
#include "AxisymmetricConductivityTensors.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "Exception.hpp"
#include "ChastePoint.hpp"
#include "AbstractChasteRegion.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"
#include "PetscVecTools.hpp"
#include "AbstractCvodeCell.hpp"
#include "Warnings.hpp"

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacTissue(
            AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory,
            bool exchangeHalos)
    : mpMesh(pCellFactory->GetMesh()),
      mpDistributedVectorFactory(mpMesh->GetDistributedVectorFactory()),
      mpConductivityModifier(NULL),
      mHasPurkinje(false),
      mDoCacheReplication(true),
      mMeshUnarchived(false),
      mExchangeHalos(exchangeHalos)
{
    //This constructor is called from the Initialise() method of the CardiacProblem class
    assert(pCellFactory != NULL);
    assert(pCellFactory->GetMesh() != NULL);

    if (PetscTools::IsSequential())
    {
        //Remove the request for a halo exchange
        mExchangeHalos = false;
    }

    unsigned num_local_nodes = mpDistributedVectorFactory->GetLocalOwnership();
    bool process_has_no_nodes = (num_local_nodes == 0u);
    if (PetscTools::ReplicateBool(process_has_no_nodes))
    {
        /* If there are no nodes on a process then there is a potential for an error in preconditioning.
         * This is dangerous because the process without nodes may segfault and not propagate the error
         * to the other processes.  Therefore, to avoid deadlock, we share this potential for error between
         * processes and throw an exception.
         */
// LCOV_EXCL_START
        // This problem normally occurs on 3 or more processes, so we can't cover it - coverage only runs with 1 and 2 processes.
        EXCEPTION("No cells were assigned some process in AbstractCardiacTissue constructor. Advice: Make total number of processors no greater than number of nodes in the mesh");
// LCOV_EXCL_STOP
    }
    unsigned ownership_range_low = mpDistributedVectorFactory->GetLow();
    mCellsDistributed.resize(num_local_nodes);

    // Figure out if we're dealing with Purkinje
    AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>* p_purkinje_cell_factory
       = dynamic_cast<AbstractPurkinjeCellFactory<ELEMENT_DIM,SPACE_DIM>*>(pCellFactory);
    if (p_purkinje_cell_factory)
    {
        mHasPurkinje = true;
        mPurkinjeCellsDistributed.resize(num_local_nodes);
    }

    /////////////////////////////////////////////////////
    // Set up cells
    /////////////////////////////////////////////////////
    try
    {
        for (unsigned local_index = 0; local_index < num_local_nodes; local_index++)
        {
            unsigned global_index = ownership_range_low + local_index;
            Node<SPACE_DIM>* p_node = mpMesh->GetNode(global_index);
            mCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(p_node);
            mCellsDistributed[local_index]->SetUsedInTissueSimulation();

            if (mHasPurkinje)
            {
                mPurkinjeCellsDistributed[local_index]
                    = p_purkinje_cell_factory->CreatePurkinjeCellForNode(p_node, mCellsDistributed[local_index]);
                mPurkinjeCellsDistributed[local_index]->SetUsedInTissueSimulation();
            }
        }

        pCellFactory->FinaliseCellCreation(&mCellsDistributed,
                                           mpDistributedVectorFactory->GetLow(),
                                           mpDistributedVectorFactory->GetHigh());
        if (mHasPurkinje)
        {
            p_purkinje_cell_factory->FinalisePurkinjeCellCreation(&mPurkinjeCellsDistributed,
                                                                  mpDistributedVectorFactory->GetLow(),
                                                                  mpDistributedVectorFactory->GetHigh());
        }
    }
    catch (const Exception& e)
    {
        // This catch statement is quite tricky to cover, but it is actually done
        // in TestCardiacSimulation::TestMono1dSodiumBlockBySettingNamedParameter()

        // Errors thrown creating cells will often be process-specific
        PetscTools::ReplicateException(true);

        // Delete cells
        // Should really do this for other processes too, but this is all we need
        // to get memory testing to pass, and leaking when we're about to die isn't
        // that bad!
        for (std::vector<AbstractCardiacCellInterface*>::iterator cell_iterator = mCellsDistributed.begin();
             cell_iterator != mCellsDistributed.end();
             ++cell_iterator)
        {
            delete (*cell_iterator);
        }

        throw e;
    }
    PetscTools::ReplicateException(false);

    // Halo nodes (if required)
    SetUpHaloCells(pCellFactory);

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    mIionicCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );
    mIntracellularStimulusCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );

    if (mHasPurkinje)
    {
        mPurkinjeIionicCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );
        mPurkinjeIntracellularStimulusCacheReplicated.Resize( pCellFactory->GetNumberOfCells() );
    }
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

    if (HeartConfig::Instance()->IsMeshProvided() && HeartConfig::Instance()->GetLoadMesh())
    {
        mFibreFilePathNoExtension = HeartConfig::Instance()->GetMeshName();
    }
    else
    {
        // As of r10671 fibre orientation can only be defined when loading a mesh from disc.
        mFibreFilePathNoExtension = "";
    }
    CreateIntracellularConductivityTensor();
}

// Constructor used for archiving
template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::AbstractCardiacTissue(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : mpMesh(pMesh),
      mpDistributedVectorFactory(mpMesh->GetDistributedVectorFactory()),
      mHasPurkinje(false),
      mDoCacheReplication(true),
      mMeshUnarchived(true),
      mExchangeHalos(false)
{
    mIionicCacheReplicated.Resize(mpDistributedVectorFactory->GetProblemSize());
    mIntracellularStimulusCacheReplicated.Resize(mpDistributedVectorFactory->GetProblemSize());

    mFibreFilePathNoExtension = ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename();
    CreateIntracellularConductivityTensor();
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::~AbstractCardiacTissue()
{
    // Delete cells
    for (std::vector<AbstractCardiacCellInterface*>::iterator iter = mCellsDistributed.begin();
         iter != mCellsDistributed.end();
         ++iter)
    {
        delete (*iter);
    }

    // Delete cells for halo nodes
    for (std::vector<AbstractCardiacCellInterface*>::iterator iter = mHaloCellsDistributed.begin();
         iter != mHaloCellsDistributed.end();
         ++iter)
    {
        delete (*iter);
    }

    delete mpIntracellularConductivityTensors;

    // Delete Purkinje cells
    for (std::vector<AbstractCardiacCellInterface*>::iterator iter = mPurkinjeCellsDistributed.begin();
         iter != mPurkinjeCellsDistributed.end();
         ++iter)
    {
        delete (*iter);
    }

    // If the mesh was unarchived we need to free it explicitly.
    if (mMeshUnarchived)
    {
        delete mpMesh;
    }
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
bool AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::HasPurkinje()
{
    return mHasPurkinje;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::CreateIntracellularConductivityTensor()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
    mpConfig = HeartConfig::Instance();

    if (mpConfig->IsMeshProvided() && mpConfig->GetLoadMesh())
    {
        assert(mFibreFilePathNoExtension != "");

        switch (mpConfig->GetConductivityMedia())
        {
            case cp::media_type::Orthotropic:
            {
                mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
                FileFinder ortho_file(mFibreFilePathNoExtension + ".ortho", RelativeTo::AbsoluteOrCwd);
                assert(ortho_file.Exists());
                mpIntracellularConductivityTensors->SetFibreOrientationFile(ortho_file);
                break;
            }

            case cp::media_type::Axisymmetric:
            {
                mpIntracellularConductivityTensors = new AxisymmetricConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
                FileFinder axi_file(mFibreFilePathNoExtension + ".axi", RelativeTo::AbsoluteOrCwd);
                assert(axi_file.Exists());
                mpIntracellularConductivityTensors->SetFibreOrientationFile(axi_file);
                break;
            }

            case cp::media_type::NoFibreOrientation:
                /// \todo #1316 Create a class defining constant tensors to be used when no fibre orientation is provided.
                mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
                break;

            default :
                NEVER_REACHED;
        }
    }
    else // Slab defined in config file or SetMesh() called; no fibre orientation assumed
    {
        /// \todo #1316 Create a class defining constant tensors to be used when no fibre orientation is provided.
        mpIntracellularConductivityTensors = new OrthotropicConductivityTensors<ELEMENT_DIM,SPACE_DIM>;
    }

    c_vector<double, SPACE_DIM> intra_conductivities;
    mpConfig->GetIntracellularConductivities(intra_conductivities);

    // this definition must be here (and not inside the if statement) because SetNonConstantConductivities() will keep
    // a pointer to it and we don't want it to go out of scope before Init() is called
    unsigned num_local_elements = mpMesh->GetNumLocalElements();
    std::vector<c_vector<double, SPACE_DIM> > hetero_intra_conductivities;

    if (mpConfig->GetConductivityHeterogeneitiesProvided())
    {
        try
        {
            assert(hetero_intra_conductivities.size()==0);
            hetero_intra_conductivities.resize(num_local_elements, intra_conductivities);
        }
        // LCOV_EXCL_START
        catch(std::bad_alloc &r_bad_alloc)
        {
            std::cout << "Failed to allocate std::vector of size " << num_local_elements << std::endl;
            PetscTools::ReplicateException(true);
            throw r_bad_alloc;
        }
        // LCOV_EXCL_STOP

        PetscTools::ReplicateException(false);

        std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > conductivities_heterogeneity_areas;
        std::vector< c_vector<double,3> > intra_h_conductivities;
        std::vector< c_vector<double,3> > extra_h_conductivities;
        HeartConfig::Instance()->GetConductivityHeterogeneities(conductivities_heterogeneity_areas,
                                                                intra_h_conductivities,
                                                                extra_h_conductivities);

        unsigned local_element_index = 0;

        for (typename AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>::ElementIterator it = mpMesh->GetElementIteratorBegin();
             it != mpMesh->GetElementIteratorEnd();
             ++it)
        {
//            unsigned element_index = it->GetIndex();
            // if element centroid is contained in the region
            ChastePoint<SPACE_DIM> element_centroid(it->CalculateCentroid());
            for (unsigned region_index=0; region_index< conductivities_heterogeneity_areas.size(); region_index++)
            {
                if (conductivities_heterogeneity_areas[region_index]->DoesContain(element_centroid))
                {
                    //We don't use ublas vector assignment here, because we might be getting a subvector of a 3-vector
                    for (unsigned i=0; i<SPACE_DIM; i++)
                    {
                        hetero_intra_conductivities[local_element_index][i] = intra_h_conductivities[region_index][i];
                    }
                }
            }
            local_element_index++;
        }

        mpIntracellularConductivityTensors->SetNonConstantConductivities(&hetero_intra_conductivities);
    }
    else
    {
        mpIntracellularConductivityTensors->SetConstantConductivities(intra_conductivities);
    }

    mpIntracellularConductivityTensors->Init(this->mpMesh);
    HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::SetCacheReplication(bool doCacheReplication)
{
    mDoCacheReplication = doCacheReplication;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
bool AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::GetDoCacheReplication()
{
    return mDoCacheReplication;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetIntracellularConductivityTensor(unsigned elementIndex)
{
    assert( mpIntracellularConductivityTensors);
    if (mpConductivityModifier==NULL)
    {
        return (*mpIntracellularConductivityTensors)[elementIndex];
    }
    else
    {
        return mpConductivityModifier->rGetModifiedConductivityTensor(elementIndex, (*mpIntracellularConductivityTensors)[elementIndex], 0u);
    }
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const c_matrix<double, SPACE_DIM, SPACE_DIM>& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetExtracellularConductivityTensor(unsigned elementIndex)
{
     EXCEPTION("Monodomain tissues do not have extracellular conductivity tensors.");
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacCellInterface* AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::GetCardiacCell( unsigned globalIndex )
{
    assert(mpDistributedVectorFactory->GetLow() <= globalIndex &&
           globalIndex < mpDistributedVectorFactory->GetHigh());
    return mCellsDistributed[globalIndex - mpDistributedVectorFactory->GetLow()];
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacCellInterface* AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::GetPurkinjeCell( unsigned globalIndex )
{
    assert(mpDistributedVectorFactory->GetLow() <= globalIndex &&
           globalIndex < mpDistributedVectorFactory->GetHigh());
    EXCEPT_IF_NOT(mHasPurkinje);
    return mPurkinjeCellsDistributed[globalIndex - mpDistributedVectorFactory->GetLow()];
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
AbstractCardiacCellInterface* AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::GetCardiacCellOrHaloCell( unsigned globalIndex )
{
    std::map<unsigned, unsigned>::const_iterator node_position;
    // First search the halo
    if ((node_position=mHaloGlobalToLocalIndexMap.find(globalIndex)) != mHaloGlobalToLocalIndexMap.end())
    {
        // Found a halo node
        return mHaloCellsDistributed[node_position->second];
    }
    // Then search the owned node
    if (mpDistributedVectorFactory->IsGlobalIndexLocal(globalIndex))
    {
        // Found an owned node
        return mCellsDistributed[globalIndex - mpDistributedVectorFactory->GetLow()];
    }
    // Not here
    EXCEPTION("Requested node/halo " << globalIndex << " does not belong to processor " << PetscTools::GetMyRank());
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::CalculateHaloNodesFromNodeExchange()
{
    std::set<unsigned> halos_as_set;
    for (unsigned proc=0; proc<PetscTools::GetNumProcs(); proc++)
    {
        halos_as_set.insert(mNodesToReceivePerProcess[proc].begin(), mNodesToReceivePerProcess[proc].end());
    }
    mHaloNodes = std::vector<unsigned>(halos_as_set.begin(), halos_as_set.end());
    //PRINT_VECTOR(mHaloNodes);
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::SetUpHaloCells(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory)
{
    if (mExchangeHalos)
    {
        mpMesh->CalculateNodeExchange(mNodesToSendPerProcess, mNodesToReceivePerProcess);
        // Note that the following call will not work for a TetrahedralMesh which has
        // no concept of halo nodes.
        //mpMesh->GetHaloNodeIndices( mHaloNodes );
        CalculateHaloNodesFromNodeExchange();
        unsigned num_halo_nodes = mHaloNodes.size();
        mHaloCellsDistributed.resize( num_halo_nodes );
        for (unsigned local_index = 0; local_index < num_halo_nodes; local_index++)
        {
            unsigned global_index = mHaloNodes[local_index];
            // These are all halo nodes, so we use the "GetNodeOrHaloNode" variety of GetNode
            Node<SPACE_DIM>* p_node = mpMesh->GetNodeOrHaloNode(global_index);
            mHaloCellsDistributed[local_index] = pCellFactory->CreateCardiacCellForNode(p_node);
            mHaloCellsDistributed[local_index]->SetUsedInTissueSimulation();
            mHaloGlobalToLocalIndexMap[global_index] = local_index;
        }
        // No need to call FinaliseCellCreation() as halo node cardiac cells will
        // never be stimulated (their values are communicated from the process that
        // owns them).
    }
}


template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::SolveCellSystems(Vec existingSolution, double time, double nextTime, bool updateVoltage)
{
    if (mHasPurkinje)
    {
        // can't do Purkinje and operator splitting
        assert(!updateVoltage);
        // The code below assumes Purkinje is are monodomain, so the vector has two stripes.
        // The assert will fail the first time bidomain purkinje is coded - need to decide what
        // ordering the three stripes (V, V_purk, phi_e) are in
        assert(PetscVecTools::GetSize(existingSolution)==2*mpMesh->GetNumNodes());
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::SOLVE_ODES);

    DistributedVector dist_solution = mpDistributedVectorFactory->CreateDistributedVector(existingSolution);

    /////////////////////////////////////////////////////////////
    // Solve cell models (except purkinje cell models)
    /////////////////////////////////////////////////////////////
    DistributedVector::Stripe voltage(dist_solution, 0);
    try
    {
        double voltage_before_update;
        for (DistributedVector::Iterator index = dist_solution.Begin();
             index != dist_solution.End();
             ++index)
        {
            voltage_before_update = voltage[index];
            mCellsDistributed[index.Local]->SetVoltage( voltage_before_update );

            // Added a try-catch here to provide more output to screen when an error occurs.
            /// \todo This may want to go to std::cerr ??
            try
            {
                if (!updateVoltage)
                {
                    // solve ODE system at this node.
                    // Note: Voltage is not being updated. The voltage is updated in the PDE solve.
#ifndef CHASTE_CVODE
                    mCellsDistributed[index.Local]->ComputeExceptVoltage(time, nextTime);
#else
                    // If CVODE is enabled, and this is a CVODE cell
                    // there's a chance we can recover this by doing a reset so put the above call in a try...catch.
                    try
                    {
                        mCellsDistributed[index.Local]->ComputeExceptVoltage(time, nextTime);
                    }
                    catch (Exception &e)
                    {
                        // Try an 'emergency' reset if this is a CVODE cell.
                        // See #2594 for why we think this may be necessary.
                        if (dynamic_cast<AbstractCvodeCell*>(mCellsDistributed[index.Local]))
                        {
                            // Reset the CVODE cell, this leads to a call to CVodeReInit.
                            static_cast<AbstractCvodeCell*>(mCellsDistributed[index.Local])->ResetSolver();
                            mCellsDistributed[index.Local]->ComputeExceptVoltage(time, nextTime);
                            WARNING("Global node " << index.Global << " had an ODE solving problem in t = [" << time <<
                                    ", " << nextTime << "] ms. This was fixed by a reset of CVODE, but may suggest PDE time"
                                    " step should be reduced, or CVODE tolerances relaxed.");
                        }
                        else
                        {
                            throw e;
                        }
                    }
#endif // CHASTE_CVODE
                }
                else
                {
                    // solve, including updating the voltage (for the operator-splitting implementation of the monodomain solver)
                    mCellsDistributed[index.Local]->SolveAndUpdateState(time, nextTime);
                    voltage[index] = mCellsDistributed[index.Local]->GetVoltage();
                }
            }
            catch (Exception &e)
            {
                std::cout << std::setprecision(16);
                std::cout << "Global node " << index.Global << " had problems with ODE solve between "
                        "t = " << time << " and " << nextTime << "ms.\n";

                std::cout << "Voltage at this node before solve was " << voltage_before_update << "mV\n"
                        "(this SHOULD NOT necessarily be the same as the one in the state variables,\n"
                        "which can be ignored and stay at the initial condition - the voltage is dictated by PDE instead of state variable.)\n";

                std::cout << "Stimulus current (NB converted to micro-Amps per cm^3) applied here is equal to:\n\t"
                    << mCellsDistributed[index.Local]->GetIntracellularStimulus(time) << " at t = " << time     << "ms,\n\t"
                    << mCellsDistributed[index.Local]->GetIntracellularStimulus(nextTime) << " at t = " << nextTime << "ms.\n";

                std::cout << "Cell model: " << dynamic_cast<AbstractUntemplatedParameterisedSystem*>(mCellsDistributed[index.Local])->GetSystemName() << "\n";

                std::cout << "All state variables are now:\n";
                std::vector<double> state_vars = mCellsDistributed[index.Local]->GetStdVecStateVariables();
                std::vector<std::string> state_var_names = mCellsDistributed[index.Local]->rGetStateVariableNames();
                for (unsigned i=0; i<state_vars.size(); i++)
                {
                    std::cout << "\t" << state_var_names[i] << "\t:\t" << state_vars[i] << "\n";
                }
                std::cout << std::flush;

                throw e;
            }
            // update the Iionic and stimulus caches
            UpdateCaches(index.Global, index.Local, nextTime);
        }

        if (updateVoltage)
        {
            dist_solution.Restore();
        }
    }
    catch (Exception &e)
    {
        PetscTools::ReplicateException(true);
        throw e;
    }

    /////////////////////////////////////////////////////////////
    // Solve purkinje cell models
    /////////////////////////////////////////////////////////////
    if (mHasPurkinje)
    {
        DistributedVector::Stripe purkinje_voltage(dist_solution, 1);
        try
        {
            for (DistributedVector::Iterator index = dist_solution.Begin();
                 index != dist_solution.End();
                 ++index)
            {
                // overwrite the voltage with the input value
                mPurkinjeCellsDistributed[index.Local]->SetVoltage( purkinje_voltage[index] );

                // solve
                // Note: Voltage is not being updated. The voltage is updated in the PDE solve.
                mPurkinjeCellsDistributed[index.Local]->ComputeExceptVoltage(time, nextTime);

                // update the Iionic and stimulus caches
                UpdatePurkinjeCaches(index.Global, index.Local, nextTime);
            }
        }
        // LCOV_EXCL_START
        catch (Exception& e)
        {
            ///\todo #2017 This code will be needed in the future, not being covered now.
            /// Add a test which covers this, e.g.  by doing a simulation with a bad stimulus for one of the
            /// Purkinje cells - stimulating with the wrong choice of sign say.

            NEVER_REACHED;
            PetscTools::ReplicateException(true);
            throw e;
        }
        // LCOV_EXCL_STOP
    }

    PetscTools::ReplicateException(false);
    HeartEventHandler::EndEvent(HeartEventHandler::SOLVE_ODES);

    // Communicate new state variable values to halo nodes
    if (mExchangeHalos)
    {
        assert(!mHasPurkinje);

        for ( unsigned rank_offset = 1; rank_offset < PetscTools::GetNumProcs(); rank_offset++ )
        {
            unsigned send_to      = (PetscTools::GetMyRank() + rank_offset) % (PetscTools::GetNumProcs());
            unsigned receive_from = (PetscTools::GetMyRank() + PetscTools::GetNumProcs()- rank_offset ) % (PetscTools::GetNumProcs());

            unsigned number_of_cells_to_send    = mNodesToSendPerProcess[send_to].size();
            unsigned number_of_cells_to_receive = mNodesToReceivePerProcess[receive_from].size();

            // Pack send buffer
            unsigned send_size = 0;
            for (unsigned i=0; i<number_of_cells_to_send; i++)
            {
                unsigned global_cell_index = mNodesToSendPerProcess[send_to][i];
                send_size += mCellsDistributed[global_cell_index - mpDistributedVectorFactory->GetLow()]->GetNumberOfStateVariables();
            }

            boost::scoped_array<double> send_data(new double[send_size]);

            unsigned send_index = 0;
            for (unsigned cell = 0; cell < number_of_cells_to_send; cell++)
            {
                unsigned global_cell_index = mNodesToSendPerProcess[send_to][cell];
                AbstractCardiacCellInterface* p_cell = mCellsDistributed[global_cell_index - mpDistributedVectorFactory->GetLow()];
                std::vector<double> cell_data = p_cell->GetStdVecStateVariables();
                const unsigned num_state_vars = p_cell->GetNumberOfStateVariables();
                for (unsigned state_variable = 0; state_variable < num_state_vars; state_variable++)
                {
                    send_data[send_index++] = cell_data[state_variable];
                }
            }
            // Receive buffer
            unsigned receive_size = 0;
            for (unsigned i=0; i<number_of_cells_to_receive; i++)
            {
                unsigned halo_cell_index = mHaloGlobalToLocalIndexMap[mNodesToReceivePerProcess[receive_from][i]];
                receive_size += mHaloCellsDistributed[halo_cell_index]->GetNumberOfStateVariables();
            }

            boost::scoped_array<double> receive_data(new double[receive_size]);

            // Send and receive
            int ret;
            MPI_Status status;
            ret = MPI_Sendrecv(send_data.get(), send_size,
                               MPI_DOUBLE,
                               send_to, 0,
                               receive_data.get(), receive_size,
                               MPI_DOUBLE,
                               receive_from, 0,
                               PETSC_COMM_WORLD, &status);
            UNUSED_OPT(ret);
            assert ( ret == MPI_SUCCESS);

            // Unpack
            unsigned receive_index = 0;
            for ( unsigned cell = 0; cell < number_of_cells_to_receive; cell++ )
            {
                AbstractCardiacCellInterface* p_cell = mHaloCellsDistributed[mHaloGlobalToLocalIndexMap[mNodesToReceivePerProcess[receive_from][cell]]];
                const unsigned number_of_state_variables = p_cell->GetNumberOfStateVariables();

                std::vector<double> cell_data(number_of_state_variables);
                for (unsigned state_variable = 0; state_variable < number_of_state_variables; state_variable++)
                {
                    cell_data[state_variable] = receive_data[receive_index++];
                }
                p_cell->SetStateVariables(cell_data);
            }
        }
    }

    HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
    if (mDoCacheReplication)
    {
        ReplicateCaches();
    }
    HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetIionicCacheReplicated()
{
    return mIionicCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetIntracellularStimulusCacheReplicated()
{
    return mIntracellularStimulusCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetPurkinjeIionicCacheReplicated()
{
    EXCEPT_IF_NOT(mHasPurkinje);
    return mPurkinjeIionicCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
ReplicatableVector& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetPurkinjeIntracellularStimulusCacheReplicated()
{
    EXCEPT_IF_NOT(mHasPurkinje);
    return mPurkinjeIntracellularStimulusCacheReplicated;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
{
    mIionicCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIIonic();
    mIntracellularStimulusCacheReplicated[globalIndex] = mCellsDistributed[localIndex]->GetIntracellularStimulus(nextTime);
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::UpdatePurkinjeCaches(unsigned globalIndex, unsigned localIndex, double nextTime)
{
    assert(mHasPurkinje);
    mPurkinjeIionicCacheReplicated[globalIndex] = mPurkinjeCellsDistributed[localIndex]->GetIIonic();
    mPurkinjeIntracellularStimulusCacheReplicated[globalIndex] = mPurkinjeCellsDistributed[localIndex]->GetIntracellularStimulus(nextTime);
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::ReplicateCaches()
{
    // ReplicateCaches only needed for SVI (and non-matrix based assembly which is no longer in code)
    // which is not implemented with Purkinje. See commented code below if introducing this.
    assert(!mHasPurkinje);

    mIionicCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());
    mIntracellularStimulusCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());

    //if (mHasPurkinje)
    //{
    //    mPurkinjeIionicCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());
    //    mPurkinjeIntracellularStimulusCacheReplicated.Replicate(mpDistributedVectorFactory->GetLow(), mpDistributedVectorFactory->GetHigh());
    //}
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const std::vector<AbstractCardiacCellInterface*>& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetCellsDistributed() const
{
    return mCellsDistributed;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const std::vector<AbstractCardiacCellInterface*>& AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::rGetPurkinjeCellsDistributed() const
{
    EXCEPT_IF_NOT(mHasPurkinje);
    return mPurkinjeCellsDistributed;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
const AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::pGetMesh() const
{
    return mpMesh;
}

template <unsigned ELEMENT_DIM,unsigned SPACE_DIM>
void AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>::SetConductivityModifier(AbstractConductivityModifier<ELEMENT_DIM,SPACE_DIM>* pModifier)
{
    assert(pModifier!=NULL);
    assert(mpConductivityModifier==NULL); // shouldn't be called twice for example, or with two different modifiers (remove this assert
                                          // if for whatever reason want to be able to overwrite modifiers)
    mpConductivityModifier = pModifier;
}

// Explicit instantiation
template class AbstractCardiacTissue<1,1>;
template class AbstractCardiacTissue<1,2>;
template class AbstractCardiacTissue<1,3>;
template class AbstractCardiacTissue<2,2>;
template class AbstractCardiacTissue<3,3>;
