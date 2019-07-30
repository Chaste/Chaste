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
#ifndef ABSTRACTCARDIACTISSUE_HPP_
#define ABSTRACTCARDIACTISSUE_HPP_

#include <set>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "UblasMatrixInclude.hpp"

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "ChasteSerializationVersion.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>

#include "AbstractCardiacCellInterface.hpp"
#include "FakeBathCell.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractConductivityTensors.hpp"
#include "AbstractPurkinjeCellFactory.hpp"
#include "ReplicatableVector.hpp"
#include "HeartConfig.hpp"
#include "ArchiveLocationInfo.hpp"
#include "AbstractDynamicallyLoadableEntity.hpp"
#include "DynamicModelLoaderRegistry.hpp"
#include "AbstractConductivityModifier.hpp"

/**
 * Class containing "tissue-like" functionality used in monodomain and bidomain
 * problems.
 *
 * Contains the cardiac cells (ODE systems for each node of the mesh) and
 * conductivity tensors (dependent on fibre directions).
 *
 * Also contains knowledge of parallelisation in the form of the
 * distributed vector factory. This class deals with created a distributed
 * vector of cells, and getting the ionic current and stimuli from these
 * cells and putting them in replicated arrays for the PDE solvers to call.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AbstractCardiacTissue : private boost::noncopyable
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    friend class TestMonodomainTissue;

    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        if (version >= 3)
        {
            archive & mHasPurkinje;
        }
        if (version >= 2)
        {
            archive & mExchangeHalos;
        }
        // Don't use the std::vector serialization for cardiac cells, so that we can load them
        // more cleverly when migrating checkpoints.
        SaveCardiacCells(*ProcessSpecificArchive<Archive>::Get(), version);

        // archive & mpMesh; Archived in save/load_constructs at the bottom of Mono/BidomainTissue.hpp
        // archive & mpIntracellularConductivityTensors; Loaded from HeartConfig every time constructor is called
        if (HeartConfig::Instance()->IsMeshProvided() && HeartConfig::Instance()->GetLoadMesh())
        {
            switch (HeartConfig::Instance()->GetConductivityMedia())
            {
                case cp::media_type::Orthotropic:
                {
                    FileFinder source_file(mFibreFilePathNoExtension + ".ortho", RelativeTo::AbsoluteOrCwd);
                    assert(source_file.Exists());
                    FileFinder dest_file(ArchiveLocationInfo::GetArchiveRelativePath() + ArchiveLocationInfo::GetMeshFilename() + ".ortho", RelativeTo::ChasteTestOutput);

                    TRY_IF_MASTER(source_file.CopyTo(dest_file));
                    break;
                }

                case cp::media_type::Axisymmetric:
                {
                    FileFinder source_file(mFibreFilePathNoExtension + ".axi", RelativeTo::AbsoluteOrCwd);
                    assert(source_file.Exists());
                    FileFinder dest_file(ArchiveLocationInfo::GetArchiveRelativePath()
                                       + ArchiveLocationInfo::GetMeshFilename() + ".axi", RelativeTo::ChasteTestOutput);

                    TRY_IF_MASTER(source_file.CopyTo(dest_file));
                    break;
                }

                case cp::media_type::NoFibreOrientation:
                    break;

                default :
                    NEVER_REACHED;
            }
        }

        // archive & mIionicCacheReplicated; // will be regenerated
        // archive & mIntracellularStimulusCacheReplicated; // will be regenerated
        archive & mDoCacheReplication;
        // archive & mMeshUnarchived; Not archived since set to true when archiving constructor is called.

        (*ProcessSpecificArchive<Archive>::Get()) & mpDistributedVectorFactory;

        // Paranoia: check we agree with the mesh on who owns what
        assert(mpDistributedVectorFactory == mpMesh->GetDistributedVectorFactory());
        assert(mpDistributedVectorFactory->GetLow()==mpMesh->GetDistributedVectorFactory()->GetLow());
        assert(mpDistributedVectorFactory->GetLocalOwnership()==mpMesh->GetDistributedVectorFactory()->GetLocalOwnership());
    }

    /**
     * Unarchive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        // archive & mpMesh; Archived in save/load_constructs at the bottom of Mono/BidomainTissue.hpp
        // archive & mpIntracellularConductivityTensors; Loaded from HeartConfig every time constructor is called

        if (version >= 3)
        {
            archive & mHasPurkinje;
            if (mHasPurkinje)
            {
                mPurkinjeIionicCacheReplicated.Resize(mpDistributedVectorFactory->GetProblemSize());
            }
        }
        if (version >= 2)
        {
            archive & mExchangeHalos;
            if (mExchangeHalos)
            {
                mpMesh->CalculateNodeExchange(mNodesToSendPerProcess, mNodesToReceivePerProcess);
                CalculateHaloNodesFromNodeExchange();
                unsigned num_halo_nodes = mHaloNodes.size();
                mHaloCellsDistributed.resize( num_halo_nodes );
                for (unsigned local_index = 0; local_index < num_halo_nodes; local_index++)
                {
                    unsigned global_index = mHaloNodes[local_index];
                    mHaloGlobalToLocalIndexMap[global_index] = local_index;
                }
            }
        }

        // mCellsDistributed & mHaloCellsDistributed:
        LoadCardiacCells(*ProcessSpecificArchive<Archive>::Get(), version);

        // archive & mIionicCacheReplicated; // will be regenerated
        // archive & mIntracellularStimulusCacheReplicated; // will be regenerated
        archive & mDoCacheReplication;

        // we no longer have a bool mDoOneCacheReplication, but to maintain backwards compatibility
        // we archive something if version==0
        if (version==0)
        {
            bool do_one_cache_replication = true;
            archive & do_one_cache_replication;
        }

        (*ProcessSpecificArchive<Archive>::Get()) & mpDistributedVectorFactory;

        // Paranoia: check we agree with the mesh on who owns what
        assert(mpDistributedVectorFactory == mpMesh->GetDistributedVectorFactory());
        assert(mpDistributedVectorFactory->GetLow()==mpMesh->GetDistributedVectorFactory()->GetLow());
        assert(mpDistributedVectorFactory->GetLocalOwnership()==mpMesh->GetDistributedVectorFactory()->GetLocalOwnership());
        // archive & mMeshUnarchived; Not archived since set to true when archiving constructor is called.

        // not archiving mpConductivityModifier for the time being (mechanics simulations are only use-case at the moment, and they
        // do not get archived...). mpConductivityModifier has to be reset to NULL upon load.
        mpConductivityModifier = NULL;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * Convenience method for intracellular conductivity tensor creation
     */
    void CreateIntracellularConductivityTensor();

protected:

    /** It's handy to keep a pointer to the mesh object*/
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* mpMesh;

    /** Intracellular conductivity tensors. Not archived, since it's loaded from the
     *  HeartConfig singleton. */
    AbstractConductivityTensors<ELEMENT_DIM,SPACE_DIM>* mpIntracellularConductivityTensors;

    /** The vector of cells. Distributed. */
    std::vector< AbstractCardiacCellInterface* > mCellsDistributed;

    /** The vector of the purkinje cells. Distributed.
     * Empty unless a AbstractPurkinjeCellFactory is given to the constructor. */
    std::vector< AbstractCardiacCellInterface* > mPurkinjeCellsDistributed;

    /**
     *  Cache containing all the ionic currents for each node,
     *  replicated over all processes.
     */
    ReplicatableVector mIionicCacheReplicated;

    /**
     *  Cache containing all the ionic currents for each purkinje node,
     *  replicated over all processes.
     */
    ReplicatableVector mPurkinjeIionicCacheReplicated;

    /**
     *  Cache containing all the stimulus currents for each node,
     *  replicated over all processes.
     */
    ReplicatableVector mIntracellularStimulusCacheReplicated;

    /**
     *  Cache containing all the stimulus currents for each Purkinje node,
     *  replicated over all processes.
     */
    ReplicatableVector mPurkinjeIntracellularStimulusCacheReplicated;

    /** Local pointer to the HeartConfig singleton instance, for convenience. */
    HeartConfig* mpConfig;

    /**
     * Local pointer to the distributed vector factory associated with the mesh object used.
     *
     * Used to retrieve node ownership range when needed.
     *
     * NB: This is set from mpMesh->GetDistributedVectorFactory() and thus always equal to
     * that.  We never assume ownership of the object.
     */
    DistributedVectorFactory* mpDistributedVectorFactory;

    /**
     * Path to the location of the fibre file without extension.
     */
    std::string mFibreFilePathNoExtension;

    /**
     * This class, if not NULL, will be used to modify the conductivity that is obtained from
     * mpIntracellularConductivityTensors when rGetIntracellularConductivityTensor() is called.
     * For example, it is required when conductivities become deformation dependent.
     */
    AbstractConductivityModifier<ELEMENT_DIM,SPACE_DIM>* mpConductivityModifier;

    /** Whether this tissue has any Purkinje cells. */
    bool mHasPurkinje;

    /**
     * Whether we need to replicate the caches.
     *
     * When doing matrix-based RHS assembly, we only actually need information from
     * cells/nodes local to the processor, so replicating the caches is an
     * unnecessary communication overhead.
     *
     * Defaults to true.
     */
    bool mDoCacheReplication;

    /**
     * Whether the mesh was unarchived or got from elsewhere.
     */
    bool mMeshUnarchived;

    /**
     * Whether to exchange cell models across the halo boundaries.
     * Used in state variable interpolation.
     */
    bool mExchangeHalos;

    /** Vector of halo node indices for current process */
    std::vector<unsigned> mHaloNodes;

    /** The vector of halo cells. Distributed. */
    std::vector< AbstractCardiacCellInterface* > mHaloCellsDistributed;

    /** Map of global to local indices for halo nodes. */
    std::map<unsigned, unsigned> mHaloGlobalToLocalIndexMap;

    /**
     * A vector which will be of size GetNumProcs() where each internal vector except
     * i=GetMyRank() contains an ordered list of indices of nodes to send to process i
     * during data exchange
     */
    std::vector<std::vector<unsigned> > mNodesToSendPerProcess;

    /**
     * A vector which will be of size GetNumProcs() for information to receive for
     * process i.
     */
    std::vector<std::vector<unsigned> > mNodesToReceivePerProcess;

    /**
     * If the mesh is a tetrahedral mesh then all elements and nodes are known.
     * The halo nodes to the ones which are actually used as cardiac cells
     * must be calculated explicitly.
     */
    void CalculateHaloNodesFromNodeExchange();

    /**
     * If #mExchangeHalos is true, this method calls CalculateHaloNodesFromNodeExchange
     * and sets up the halo cell data structures #mHaloCellsDistributed and #mHaloGlobalToLocalIndexMap.
     *
     * @param pCellFactory  cell factory to use to create halo cells
     */
    void SetUpHaloCells(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory);

public:
    /**
     * This constructor is called from the Initialise() method of the CardiacProblem class.
     * It creates all the cell objects, and sets up the conductivities.
     *
     * Note that pCellFactory contains a pointer to the mesh
     *
     * @param pCellFactory  factory to use to create cardiac cells.
     *     If this is actually an AbstractPurkinjeCellFactory it creates purkinje cells.
     * @param exchangeHalos used in state-variable interpolation.  Defaults to false.
     */
    AbstractCardiacTissue(AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory, bool exchangeHalos=false);

    /**
     * This constructor is called by the archiver only.
     *
     * @param pMesh  a pointer to the AbstractTetrahedral mesh.
     */
    AbstractCardiacTissue(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /** Virtual destructor */
    virtual ~AbstractCardiacTissue();

    /** @return whether this tissue contains Purkinje fibres. */
    bool HasPurkinje();

    /**
     * Set whether or not to replicate the caches across all processors.
     *
     * See also mDoCacheReplication.
     * @param doCacheReplication - true if the cache needs to be replicated
     */
    void SetCacheReplication(bool doCacheReplication);

    /**
     * Get whether or not to replicate the caches across all processors.
     *
     * @return mDoCacheReplication - true if the cache needs to be replicated
     */
    bool GetDoCacheReplication();

    /** @return the intracellular conductivity tensor for the given element
     * @param elementIndex  index of the element of interest
     */
    const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetIntracellularConductivityTensor(unsigned elementIndex);

    /**
     * @return the extracellular conductivity tensor for the given element
     * (this throws an exception in this abstract class since monodomain don't have extracellular ones)
     * it is overridden in the BidomainTissue.
     *
     * @param elementIndex  index of the element of interest
     */
    virtual const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetExtracellularConductivityTensor(unsigned elementIndex);

    /**
     * @return a pointer to a cell, indexed by the global node index.
     *
     * \note Should only called by the process owning the cell -
     * triggers an assertion otherwise.
     *
     * @param globalIndex  global node index for which to retrieve a cell
     */
    AbstractCardiacCellInterface* GetCardiacCell( unsigned globalIndex );

    /**
     * @return a pointer to a Purkinje cell, indexed by the global node index.
     *
     * \note Should only called by the process owning the cell -
     * triggers an assertion otherwise.
     *
     * @param globalIndex  global node index for which to retrieve a cell
     */
    AbstractCardiacCellInterface* GetPurkinjeCell( unsigned globalIndex );

    /**
     * @return a pointer to a halo cell, indexed by the global node index.
     *
     * \note Should only called by the process halo owning the cell -
     * triggers an assertion otherwise.
     *
     * @param globalIndex  global node index for which to retrieve a cell
     */
    AbstractCardiacCellInterface* GetCardiacCellOrHaloCell( unsigned globalIndex );

    /**
     * Integrate the cell ODEs and update ionic current etc for each of the
     * cells, between the two times provided.
     *
     * @param existingSolution  the current voltage solution vector
     * @param time  the current simulation time
     * @param nextTime  when to simulate the cells until
     * @param updateVoltage whether to also solve for the voltage (generally false, true for operator splitting methods). Defaults to false
     */
    virtual void SolveCellSystems(Vec existingSolution, double time, double nextTime, bool updateVoltage=false);

    /** @return the entire ionic current cache */
    ReplicatableVector& rGetIionicCacheReplicated();

    /** @return the entire stimulus current cache */
    ReplicatableVector& rGetIntracellularStimulusCacheReplicated();

    /** @return the entire Purkinje ionic current cache */
    ReplicatableVector& rGetPurkinjeIionicCacheReplicated();

    /** @return the entire Purkinje stimulus current cache */
    ReplicatableVector& rGetPurkinjeIntracellularStimulusCacheReplicated();


    /**
     * Update the Iionic and intracellular stimulus caches.
     *
     * @param globalIndex  global index of the entry to update
     * @param localIndex  local index of the entry to update
     * @param nextTime  the next PDE time point, at which to evaluate the stimulus current
     */
    void UpdateCaches(unsigned globalIndex, unsigned localIndex, double nextTime);

    /**
     * Update the Iionic and intracellular stimulus caches for Purkinje cells.
     *
     * @param globalIndex  global index of the entry to update
     * @param localIndex  local index of the entry to update
     * @param nextTime  the next PDE time point, at which to evaluate the stimulus current
     */
    void UpdatePurkinjeCaches(unsigned globalIndex, unsigned localIndex, double nextTime);

    /**
     *  Replicate the Iionic and intracellular stimulus caches.
     */
    void ReplicateCaches();

    /**
     *  @return a reference to the vector of distributed cells. Needed for archiving.
     */
    const std::vector<AbstractCardiacCellInterface*>& rGetCellsDistributed() const;

    /**
     *  @return a reference to the vector of distributed Purkinje cells. Needed for archiving.
     */
    const std::vector<AbstractCardiacCellInterface*>& rGetPurkinjeCellsDistributed() const;

    /**
     *  @return a pointer to the mesh object
     *
     *  @return pointer to mesh object
     */
    const AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pGetMesh() const;

    /**
     * Set a modifier class which will be used to modifier a conductivity obtained from mpIntracellularConductivityTensors
     * when rGetIntracellularConductivityTensor() is called. For example, it is required when conductivities become deformation-dependent.
     * @param pModifier Pointer to the concrete modifier class
     */
    void SetConductivityModifier(AbstractConductivityModifier<ELEMENT_DIM,SPACE_DIM>* pModifier);

    /**
     * Save our tissue to an archive.
     *
     * Writes:
     *  -# #mpDistributedVectorFactory
     *  -# number of cells on this process
     *  -# each cell pointer in turn, interleaved with Purkinje cells if present
     *
     * @param archive  the process-specific archive to write cells to.
     * @param version
     */
    template<class Archive>
    void SaveCardiacCells(Archive & archive, const unsigned int version) const
    {
        const std::vector<AbstractCardiacCellInterface*> & r_cells_distributed = rGetCellsDistributed();
        assert(mpDistributedVectorFactory == this->mpMesh->GetDistributedVectorFactory());
        archive & mpDistributedVectorFactory; // Needed when loading
        const unsigned num_cells = r_cells_distributed.size();
        archive & num_cells;
        for (unsigned i=0; i<num_cells; i++)
        {
            AbstractDynamicallyLoadableEntity* p_entity = dynamic_cast<AbstractDynamicallyLoadableEntity*>(r_cells_distributed[i]);
            bool is_dynamic = (p_entity != NULL);
            archive & is_dynamic;
            if (is_dynamic)
            {
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
                archive & p_entity->GetLoader()->GetLoadableModulePath();
#else
                // We should have thrown an exception before this point
                NEVER_REACHED;
#endif // CHASTE_CAN_CHECKPOINT_DLLS
            }
            archive & r_cells_distributed[i];
            if (mHasPurkinje)
            {
                archive & rGetPurkinjeCellsDistributed()[i];
            }
        }
    }

    /**
     * Load our tissue from an archive.
     *
     * Handles the checkpoint migration case, deleting loaded cells immediately if they are
     * not local to this process.
     *
     * Also loads halo cells if we're doing halo exchange, by using the non-local cells from the
     * archive.
     *
     * @param archive  the process-specific archive to load from
     * @param version  archive version
     */
    template<class Archive>
    void LoadCardiacCells(Archive & archive, const unsigned int version)
    {
        // Note that p_factory loaded from this archive might not be the same as our mesh's factory,
        // since we're loading a process-specific archive that could have been written by any process.
        // We therefore need to use p_mesh_factory to determine the partitioning in use for the resumed
        // simulation, and p_factory to determine what the original partitioning was when the simulation
        // was saved.
        DistributedVectorFactory* p_factory;
        DistributedVectorFactory* p_mesh_factory = this->mpMesh->GetDistributedVectorFactory();
        archive & p_factory;
        unsigned num_cells;
        archive & num_cells;
        if (mCellsDistributed.empty())
        {
            mCellsDistributed.resize(p_mesh_factory->GetLocalOwnership());
#ifndef NDEBUG
            // Paranoia
            for (unsigned i=0; i<mCellsDistributed.size(); i++)
            {
                assert(mCellsDistributed[i] == NULL);
            }
#endif
        }
        else
        {
            assert(mCellsDistributed.size() == p_mesh_factory->GetLocalOwnership());
        }
        if (mHasPurkinje)
        {
            if (mPurkinjeCellsDistributed.empty())
            {
                mPurkinjeCellsDistributed.resize(p_mesh_factory->GetLocalOwnership());
            }
            else
            {
                assert(mPurkinjeCellsDistributed.size() == p_mesh_factory->GetLocalOwnership());
            }
        }

        // We don't store a cell index in the archive, so need to work out what global index this tissue starts at.
        // If we have an original factory we use the original low index; otherwise we use the current low index.
        unsigned index_low = p_factory->GetOriginalFactory() ? p_factory->GetOriginalFactory()->GetLow() : p_mesh_factory->GetLow();

        // Track fake cells (which might have multiple pointers to the same object) to make sure we only delete non-local ones
        std::set<FakeBathCell*> fake_cells_non_local, fake_cells_local;

        /*
         * Historical note:
         *
         * We always do a dumb partition when we unarchive.
         *
         * When unarchive was first implemented in parallel (migration #1199) it was thought that we might want to repartition the mesh. This would be feasible and would give
         * better partitions when we move to different numbers of processes.  However it would require us to apply a new permutation to entire data structure.
         *
         * In the case where the original mesh was permuted and *copied* into the archive, we need to apply the stored permutation to the mesh but not to the archive (cells).  That
         * is, any permutation stored with the mesh can be safely ignored.  (If we also had to repartition/permute the archive, then we would be applying a double permutation to the
         * mesh and a single permutation to archive.)
         *
         */
        for (unsigned local_index=0; local_index<num_cells; local_index++)
        {
            // Figure out where this cell goes
            unsigned global_index = index_low + local_index;
            bool local = p_mesh_factory->IsGlobalIndexLocal(global_index);

            // Check if this will be a halo cell
            std::map<unsigned, unsigned>::const_iterator halo_position;
            bool halo = ((halo_position=mHaloGlobalToLocalIndexMap.find(global_index)) != mHaloGlobalToLocalIndexMap.end());
            // halo_position->second is local halo index

            bool is_dynamic;
            archive & is_dynamic;
            if (is_dynamic)
            {
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
                // Ensure the shared object file for this cell model is loaded.
                // We need to do this here, rather than in the class' serialization code,
                // because that code won't be available until this is done...
                std::string shared_object_path;
                archive & shared_object_path;
                DynamicModelLoaderRegistry::Instance()->GetLoader(shared_object_path);
#else
                // Could only happen on Mac OS X, and will probably be trapped earlier.
                NEVER_REACHED;
#endif // CHASTE_CAN_CHECKPOINT_DLLS
            }
            AbstractCardiacCellInterface* p_cell;
            archive & p_cell;
            AbstractCardiacCellInterface* p_purkinje_cell = NULL;
            if (mHasPurkinje)
            {
                archive & p_purkinje_cell;
            }
            // Check if it's a fake cell
            FakeBathCell* p_fake = dynamic_cast<FakeBathCell*>(p_cell);
            if (p_fake)
            {
                if (halo || local)
                {
                    fake_cells_local.insert(p_fake);
                }
                else
                {
                    fake_cells_non_local.insert(p_fake);
                }
            }
            FakeBathCell* p_fake_purkinje = dynamic_cast<FakeBathCell*>(p_purkinje_cell);
            if (p_fake_purkinje)
            {
                if (halo || local)
                {
                    fake_cells_local.insert(p_fake_purkinje);
                }
                else
                {
                    fake_cells_non_local.insert(p_fake_purkinje);
                }
            }
            // Add real cells to the local or halo vectors
            if (local)
            {
                // Note that the original local_index was relative to the archived partition (distributed vector)
                // The new_local_index is local relative to the new partition in memory
                unsigned new_local_index = global_index - p_mesh_factory->GetLow();
                assert(mCellsDistributed[new_local_index] == NULL);
                mCellsDistributed[new_local_index] = p_cell;
                if (mHasPurkinje)
                {
                    assert(mPurkinjeCellsDistributed[new_local_index] == NULL);
                    mPurkinjeCellsDistributed[new_local_index] = p_purkinje_cell;
                }
            }
            else if (halo)
            {
                assert(mHaloCellsDistributed[halo_position->second] == NULL);
                mHaloCellsDistributed[halo_position->second] = p_cell;
            }
            else
            {
                if (!p_fake)
                {
                    // Non-local real cell, so free the memory.
                    delete p_cell;
                }
                if (!p_fake_purkinje)
                {
                    // This will be NULL if there's no Purkinje, so a delete is OK.
                    delete p_purkinje_cell;
                }
            }
        }

        // Delete any unused fake cells
        for (std::set<FakeBathCell*>::iterator it = fake_cells_non_local.begin();
             it != fake_cells_non_local.end();
             ++it)
        {
            if (fake_cells_local.find(*it) == fake_cells_local.end())
            {
                delete (*it);
            }
        }
    }
};

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractCardiacTissue)

namespace boost {
namespace serialization {
/**
 * Specify a version number for archive backwards compatibility.
 *
 * This is how to do BOOST_CLASS_VERSION(AbstractCardiacTissue, 1)
 * with a templated class.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct version<AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM> >
{
    ///Macro to set the version number of templated archive in known versions of Boost
    CHASTE_VERSION_CONTENT(3);
};
} // namespace serialization
} // namespace boost

#endif /*ABSTRACTCARDIACTISSUE_HPP_*/

