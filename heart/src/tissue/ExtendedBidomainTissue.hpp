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


#ifndef EXTENDEDBIDOMAINTISSUE_HPP_
#define EXTENDEDBIDOMAINTISSUE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>
#include "UblasMatrixInclude.hpp"

#include "AbstractStimulusFunction.hpp"
#include "AbstractStimulusFactory.hpp"
#include "AbstractConductivityTensors.hpp"

#include "AbstractCardiacTissue.hpp"

/**
 * Class that provides functionalities to specify a tissue within the context of the extended bidomain framework.
 *
 * The extended bidomain equations are of the form:
 *
 *  - div ( sigma_i1 grad Phi_1 ) + Am1*Cm1*d(phi_1)/dt - Am2*Cm2*d(phi_e)/dt + Am1*I_ion1 - Am1*I_stim1 + Amgap*G_gap*(phi_1 - phi_2)
 *  - div ( sigma_i2 grad Phi_2 ) + Am2*Cm2*d(phi_2)/dt - Am2*Cm2*d(phi_e)/dt + Am2*I_ion2 - Am2*I_stim2 - Amgap*G_gap*(phi_1 - phi_2)
 *   div ( sigma_e grad Phi_e ) + div ( sigma_i1 grad Phi_1 ) + div ( sigma_i2 grad Phi_2 )  = I_stim
 *
 *   The unknowns are:
 *
 *   - Phi_1 (intracellular potential of the first cell).
 *   - Phi_2 (intracellular potential of the second cell).
 *   - Phi_e (extracellular potential).
 *
 *   Am1, Am2 and Amgap are surface-to-volume ratios for first cell, second cell and gap junction. User can set their values.
 *   Cm1 and cm2 are capaciatnce values of first and second cell respectively
 *   sigma_i1 and sigma_i2 are intracellular conductivity tensors of first and second cell respectively
 *   sigma_e is the conductivity tensor for the extracellular space
 *   G_gap is the conductance (in ms/cm2) of the gap junction channel.
 *
 *
 *
 */
template <unsigned SPACE_DIM>
class ExtendedBidomainTissue : public virtual AbstractCardiacTissue<SPACE_DIM>
{
private:
    friend class TestExtendedBidomainTissue; // for testing.

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacTissue<SPACE_DIM> >(*this);
        // Conductivity tensors are dealt with by HeartConfig, and the caches get regenerated.

        archive & mAmFirstCell;
        archive & mAmSecondCell;
        archive & mAmGap;
        archive & mCmFirstCell;
        archive & mCmSecondCell;
        archive & mGGap;
        archive & mUserSuppliedExtracellularStimulus;
    }

    /** Intracellular conductivity tensors for the second cell.*/
    AbstractConductivityTensors<SPACE_DIM, SPACE_DIM> *mpIntracellularConductivityTensorsSecondCell;

    /**
     * Stores the values of the conductivities for the second cell. Accessible via get and set methods. The problem class will set it
     * This variable is a convenient interface for other classes. It is used to fill in mpIntracellularConductivityTensorsSecondCell.
     */
    c_vector<double, SPACE_DIM>  mIntracellularConductivitiesSecondCell;


    /** Extracellular conductivity tensors. */
    AbstractConductivityTensors<SPACE_DIM, SPACE_DIM> *mpExtracellularConductivityTensors;

    /**
     *  Cache containing all the stimulus currents for each node,
     *  replicated over all processes.
     */
    ReplicatableVector mExtracellularStimulusCacheReplicated;

    /**
     *  Cache containing all the stimulus currents for each node,
     *  replicated over all processes.
     */
    ReplicatableVector mGgapCacheReplicated;

    /**
     *  Cache containing all the ionic currents for each node for the seconed cell,
     *  replicated over all processes.
     */
    ReplicatableVector mIionicCacheReplicatedSecondCell;

    /**
     *  Cache containing all the stimulus currents for each node for the second cell,
     *  replicated over all processes.
     */
    ReplicatableVector mIntracellularStimulusCacheReplicatedSecondCell;

    /** The vector of cells (the second one). Distributed. */
    std::vector< AbstractCardiacCellInterface* > mCellsDistributedSecondCell;

    /** The vector of stimuli for the extracellular stimulus. Distributed. */
    std::vector<boost::shared_ptr<AbstractStimulusFunction> > mExtracellularStimuliDistributed;

    /** The vector of gap junction conductances. Distributed*/
    std::vector<double> mGgapDistributed;

    /**the Am for the first cell, set by the problem class and picked up by the assembler*/
    double mAmFirstCell;
    /**the Am for the second cell, set by the problem class and picked up by the assembler*/
    double mAmSecondCell;
    /**the Am for the gap junction, set by the problem class and picked up by the assembler*/
    double mAmGap;
    /**the Cm for the first cell, set by the problem class and picked up by the assembler*/
    double mCmFirstCell;
    /**the Cm for the second cell, set by the problem class and picked up by the assembler*/
    double mCmSecondCell;
    /**the conductance of the gap junction, in mS/cm2. Set by the problem class and picked up by the assembler*/
    double mGGap;

    /**
     * Whether the extracellular stimulus that is passed in was supplied by the user or not
     * (it could be the default zero implementation). Initialise to false (user did not pass in anything).
     */
    bool mUserSuppliedExtracellularStimulus;

    /**
     * Convenience method for extracellular conductivity tensors creation
     */
    void CreateExtracellularConductivityTensors();

    /**
     * The parent class AbstractCardiacTissue has a method UpdateCaches that updates some caches of general use.
     * This method updates more caches that are specific to extended bidomain problems, namely:
     *
     * - Iionic and intracellular stimulus for the second cell
     * - Extracellular stimulus
     * - Gap junction conductivities (Ggap)
     *
     * It is typically called right after the UpdateCaches method in the parent class.
     *
     * @param globalIndex  global index of the entry to update
     * @param localIndex  local index of the entry to update
     * @param nextTime  the next PDE time point, at which to evaluate the stimulus current
     */
    void UpdateAdditionalCaches(unsigned globalIndex, unsigned localIndex, double nextTime);

    /**
     * The parent class AbstractCardiacTissue has a method ReplicateCaches that replicates some caches of general use.
     * This method replicates more caches that are specific to extended bidomain problems, namely:
     *
     * - Iionic and intracellular stimulus for the second cell
     * - Extracellular stimulus
     * - Gap junction conductivities (Ggap)
     *
     * It is typically called right after the ReplicateCaches method in the parent class.
     */
    void ReplicateAdditionalCaches();

    /** vector of regions for Ggap heterogeneities*/
    std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > mGgapHeterogeneityRegions;
    /**values of heterogeneous Ggaps corresponding to mGgapHeterogeneityRegions. This has the same size as mGgapHeterogeneityRegions*/
    std::vector<double> mGgapValues;

public:

    /**
     * Constructor sets up extracellular conductivity tensors.
     * @param pCellFactory factory to pass on to the base class constructor
     * @param pCellFactorySecondCell factory to pass on to the base class constructor for the second cell
     * @param pExtracellularStimulusFactory factory for creating extracellular stimuli
     */
    ExtendedBidomainTissue(AbstractCardiacCellFactory<SPACE_DIM>* pCellFactory, AbstractCardiacCellFactory<SPACE_DIM>* pCellFactorySecondCell, AbstractStimulusFactory<SPACE_DIM>* pExtracellularStimulusFactory);

    /**
     *  Archiving constructor
     * @param rCellsDistributed  local cell models (recovered from archive)
     * @param rSecondCellsDistributed  local cell models for second cells (recovered from archive)
     * @param rExtraStimuliDistributed local extracellular stimuli (recovered from archive)
     * @param rGgapsDistributed distributed Ggaps (recovered from archive)
     * @param pMesh  a pointer to the AbstractTetrahedral mesh (recovered from archive).
     * @param intracellularConductivitiesSecondCell a vector with the orthotropic conductivities for the second cell (this is needed because the second cell values may not be taken from HeartConfig as the the ones for the first cell are).
     */
    ExtendedBidomainTissue(std::vector<AbstractCardiacCellInterface*> & rCellsDistributed,
                           std::vector<AbstractCardiacCellInterface*> & rSecondCellsDistributed,
                           std::vector<boost::shared_ptr<AbstractStimulusFunction> > & rExtraStimuliDistributed,
                           std::vector<double>& rGgapsDistributed,
                           AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh,
                           c_vector<double, SPACE_DIM>  intracellularConductivitiesSecondCell);

    /**
     * Destructor
     */
    virtual ~ExtendedBidomainTissue();

    /**
     * Sets the value of the conductivities for the second cell.
     *
     * @param conductivities the conductivities to be set.
     */
    void SetIntracellularConductivitiesSecondCell(c_vector<double, SPACE_DIM> conductivities);

    /**
     * @return a pointer to the second cell
     *
     * @param globalIndex the global index in the mesh
     */
    AbstractCardiacCellInterface* GetCardiacSecondCell( unsigned globalIndex );


    /**
     * @return a pointer to the extracellular stimulus. Useful for testing
     *
     * @param globalIndex the global index in the mesh
     */
    boost::shared_ptr<AbstractStimulusFunction> GetExtracellularStimulus( unsigned globalIndex );

    /**
     *  @return a reference to the vector of distributed cells (second cell). Needed for archiving.
     */
    const std::vector<AbstractCardiacCellInterface*>& rGetSecondCellsDistributed() const;

    /**
     *  @return a reference to the vector of distributed values of Ggaps. Needed for archiving.
     */
    const std::vector<double>& rGetGapsDistributed() const;


    /**
     *  @return a reference to the vector of distributed extracellular stimuli. Needed for archiving.
     */
    const std::vector<boost::shared_ptr<AbstractStimulusFunction> >& rGetExtracellularStimulusDistributed() const;


    /**
     * @return the intracellular conductivities of the second cell
     */
    c_vector<double, SPACE_DIM> GetIntracellularConductivitiesSecondCell() const;

    /**
     * Integrate the cell ODEs and update ionic current etc for each of the
     * cells, between the two times provided. This is a re-implementation from the version in the base class.
     *
     * @param existingSolution  the current voltage solution vector
     * @param time  the current simulation time
     * @param nextTime  when to simulate the cells until
     * @param updateVoltage (unused here)
     */
    virtual void SolveCellSystems(Vec existingSolution, double time, double nextTime, bool updateVoltage = false);

    /**
     * Convenience method for intracellular conductivity tensors creation for the second cell
     */
    void CreateIntracellularConductivityTensorSecondCell();

    /**
     *  Set the values of mCellHeterogeneityRegions and mGgapValues for the heterogeneities of Ggap.
     *
     *  @param rGgapHeterogeneityRegions a vector of (pointers to) heterogeneity regions for gap junctions
     *  @param rGgapValues a vector (of the same size as rGgapHeterogeneityRegions) with the respective values of Ggap for every region.
     */
    void SetGgapHeterogeneities ( std::vector<boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > > & rGgapHeterogeneityRegions, std::vector<double> rGgapValues);

    /**
     * Create the pattern of Ggap across the mesh based upon mCellHeterogeneityRegions, mGgapValues and mGgap. This will fill in mGgapDistributed.
     * It will set mGgap everywhere except in  the areas mCellHeterogeneityRegions[i] where it will put mGgapValues[i] instead.
     * If mCellHeterogeneityRegions (and mGgapValues) are empty, mGgap will be set everywhere.
     */
    void CreateGGapConductivities();

    /**
     * @return the extracellular conductivity tensor for the given element
     * @param elementIndex  index of the element of interest
     */
     const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetExtracellularConductivityTensor(unsigned elementIndex);

     /**
      * @return the intracellular conductivity tensor for the given element for tehs econd cell
      * @param elementIndex  index of the element of interest
      */
      const c_matrix<double, SPACE_DIM, SPACE_DIM>& rGetIntracellularConductivityTensorSecondCell(unsigned elementIndex);


     /** @return the entire ionic current cache for the second cell*/
     ReplicatableVector& rGetIionicCacheReplicatedSecondCell();

     /** @return the entire stimulus current cache for the second cell*/
     ReplicatableVector& rGetIntracellularStimulusCacheReplicatedSecondCell();

     /** @return the extracellular stimulus*/
     ReplicatableVector& rGetExtracellularStimulusCacheReplicated();

     /** @return the values of ggap*/
     ReplicatableVector& rGetGgapCacheReplicated();

     /**
      * @return Am for the first cell
      */
     double GetAmFirstCell();

     /**
      * @return Am for the second cell
      */
     double GetAmSecondCell();

     /**
      * @return Am for the gap junction
      */
     double GetAmGap();

     /**
      * @return Cm for the first cell
      */
     double GetCmFirstCell();

     /**
      * @return Cm for the second cell
      */
     double GetCmSecondCell();

     /**
      *  @return the conducatnce of the gap junction (mGGap)
      */
     double GetGGap();

     /**
      * @param value Am for the first cell
      */
     void SetAmFirstCell(double value);

     /**
      * @param value Am for the second cell
      */
     void SetAmSecondCell(double value);

     /**
      * @param value Am for the gap junction
      */
     void SetAmGap(double value);

     /**
      * @param value Cm for the first cell
      */
     void SetCmFirstCell(double value);

     /**
      * @param value Cm for the first cell
      */
     void SetCmSecondCell(double value);

     /**
      * @param value conductance, in mS of the gap junction
      */
     void SetGGap(double value);

     /**
      * This method gives access to the member variable mUserSuppliedExtracellularStimulus,
      * which is false by default but turned true if the user supplies an extracellular stimulus
      * in any form.
      *
      * @return true if the user supplied an extracellular stimulus.
      */
     bool HasTheUserSuppliedExtracellularStimulus();

     /**
      * This method allows modifications of the  mUserSuppliedExtracellularStimulus flag (false by default).
      * Other classes (e.g., Problem classes) can use this method to tell the Tissue that the user
      * specified an extracellular stimulus.
      *
      * @param flag ; true if you want to tell the Tissue object that the user supplied an extracellular stimulus explicitly
      */
     void SetUserSuppliedExtracellularStimulus(bool flag);


     /**
      * This method is the equivalent of SaveCardiacCells in the abstract class but save both cells of the extended bidomain tissue
      *
      * @param archive  the master archive; cells will actually be written to the process-specific archive.
      * @param version
      */
     template<class Archive>
     void SaveExtendedBidomainCells(Archive & archive, const unsigned int version) const
     {
         Archive& r_archive = *ProcessSpecificArchive<Archive>::Get();
         const std::vector<AbstractCardiacCellInterface*> & r_cells_distributed = this->rGetCellsDistributed();
         const std::vector<AbstractCardiacCellInterface*> & r_cells_distributed_second_cell = rGetSecondCellsDistributed();
         const std::vector<double> & r_ggaps_distributed = rGetGapsDistributed();

         r_archive & this->mpDistributedVectorFactory; // Needed when loading
         const unsigned num_cells = r_cells_distributed.size();
         r_archive & num_cells;
         for (unsigned i=0; i<num_cells; i++)
         {
             AbstractDynamicallyLoadableEntity* p_entity = dynamic_cast<AbstractDynamicallyLoadableEntity*>(r_cells_distributed[i]);
             bool is_dynamic = (p_entity != NULL);
             r_archive & is_dynamic;
             if (is_dynamic)
             {
 #ifdef CHASTE_CAN_CHECKPOINT_DLLS
                 ///\todo Dynamically loaded cell models aren't saved to archive in extended Bidomain
                 NEVER_REACHED;
                 //r_archive & p_entity->GetLoader()->GetLoadableModulePath();
 #else
                 // We should have thrown an exception before this point
                 NEVER_REACHED;
 #endif // CHASTE_CAN_CHECKPOINT_DLLS
             }
             r_archive & r_cells_distributed[i];
             r_archive & r_cells_distributed_second_cell[i];
             r_archive & r_ggaps_distributed[i];
         }
     }

     /**
      * Load our tissue from an archive. This is the equivalent of LoadCardiacCells in the abstract class.
      * it loads the two cells instead of only one.
      *
      * Handles the checkpoint migration case, deleting loaded cells immediately if they are
      * not local to this process.
      *
      * @param archive  the process-specific archive to load from
      * @param version  archive version
      * @param rCells  vector to fill in with pointers to local cells
      * @param rSecondCells vector to fill in with pointers to the second cells
      * @param rGgaps vector of values of gap junctions
      * @param pMesh  the mesh, so we can get at the node permutation, if any
      */
     template<class Archive>
     static void LoadExtendedBidomainCells(Archive & archive, const unsigned int version,
                                  std::vector<AbstractCardiacCellInterface*>& rCells,
                                  std::vector<AbstractCardiacCellInterface*>& rSecondCells,
                                  std::vector<double>& rGgaps,
                                  AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
     {
         assert(pMesh!=NULL);
         DistributedVectorFactory* p_factory;
         archive & p_factory;
         unsigned num_cells;
         archive & num_cells;
         rCells.resize(p_factory->GetLocalOwnership());
         rSecondCells.resize(p_factory->GetLocalOwnership());
         rGgaps.resize(p_factory->GetLocalOwnership());
 #ifndef NDEBUG
         // Paranoia
         assert(rCells.size() == rSecondCells.size());
         for (unsigned i=0; i<rCells.size(); i++)
         {
             assert(rCells[i] == NULL);
             assert(rSecondCells[i] == NULL);
         }
 #endif

         // We don't store a cell index in the archive, so need to work out what global
         // index this tissue starts up.  If we're migrating (so have an
         // original factory) we use the original low index; otherwise we use the current
         // low index.
         unsigned index_low = p_factory->GetOriginalFactory() ? p_factory->GetOriginalFactory()->GetLow() : p_factory->GetLow();

         for (unsigned local_index=0; local_index<num_cells; local_index++)
         {
             unsigned global_index = index_low + local_index;
             unsigned new_local_index = global_index - p_factory->GetLow();
             bool local = p_factory->IsGlobalIndexLocal(global_index);

             bool is_dynamic;
             archive & is_dynamic;

             if (is_dynamic)
             {
 #ifdef CHASTE_CAN_CHECKPOINT_DLLS
                 ///\todo Dynamically loaded cell models aren't loaded from archive in extended Bidomain
                 NEVER_REACHED;
                 // Ensure the shared object file for this cell model is loaded.
                 // We need to do this here, rather than in the class' serialization code,
                 // because that code won't be available until this is done...
//                 std::string shared_object_path;
//                 archive & shared_object_path;
//                 DynamicModelLoaderRegistry::Instance()->GetLoader(shared_object_path);
 #else
                 // Could only happen on Mac OS X, and will probably be trapped earlier.
                 NEVER_REACHED;
 #endif // CHASTE_CAN_CHECKPOINT_DLLS
             }

             AbstractCardiacCellInterface* p_cell;
             AbstractCardiacCellInterface* p_second_cell;
             double g_gap;
             archive & p_cell;
             archive & p_second_cell;
             archive & g_gap;
             if (local)
             {
                 rCells[new_local_index] = p_cell; // Add to local cells
                 rSecondCells[new_local_index] = p_second_cell;
                 rGgaps[new_local_index] = g_gap;
             }
             else
             {
                 //not sure how to cover this, we are already looping over local cells...
                 NEVER_REACHED;
                // Non-local real cell, so free the memory.
                // delete p_cell;
                // delete p_second_cell;
             }
         }
     }

     /**
      * This method is the equivalent of SaveCardiacCells but Saves the extracellular stimulus instead
      *
      * @param archive  the master archive; cells will actually be written to the process-specific archive.
      * @param version
      */
     template<class Archive>
     void SaveExtracellularStimulus(Archive & archive, const unsigned int version) const
     {
         Archive& r_archive = *ProcessSpecificArchive<Archive>::Get();
         const std::vector<boost::shared_ptr<AbstractStimulusFunction> > & r_stimulus_distributed = rGetExtracellularStimulusDistributed();
         r_archive & this->mpDistributedVectorFactory; // Needed when loading
         const unsigned num_cells = r_stimulus_distributed.size();
         r_archive & num_cells;
         for (unsigned i=0; i<num_cells; i++)
         {
             r_archive & r_stimulus_distributed[i];
         }
     }

     /**
      * This method is the equivalent of LoadCardiacCells but Load the extracellular stimulus instead
      *
      * @param archive  the master archive; cells will actually be written to the process-specific archive.
      * @param version
      * @param rStimuli the extracellular stimuli (will be filled from the archive).
      * @param pMesh the mesh (needed to work out number of nodes). Here it is assumed we have already unarchived the mesh somewhere and the pointer passed in is not NULL.
      */
     template<class Archive>
     void LoadExtracellularStimulus(Archive & archive, const unsigned int version,
                                               std::vector<boost::shared_ptr<AbstractStimulusFunction> >& rStimuli,
                                               AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* pMesh)
     {

        DistributedVectorFactory* p_factory;
        archive & p_factory;
        unsigned num_cells;
        archive & num_cells;
        rStimuli.resize(p_factory->GetLocalOwnership());
#ifndef NDEBUG
          // Paranoia
          for (unsigned i=0; i<rStimuli.size(); i++)
          {
              assert(rStimuli[i] == NULL);
          }
#endif

        // We don't store a cell index in the archive, so need to work out what global
        // index this tissue starts up.  If we're migrating (so have an
        // original factory) we use the original low index; otherwise we use the current
        // low index.
        unsigned index_low = p_factory->GetOriginalFactory() ? p_factory->GetOriginalFactory()->GetLow() : p_factory->GetLow();

        assert(pMesh!=NULL);
        //unsigned num_cells = pMesh->GetNumNodes();
        for (unsigned local_index=0; local_index<num_cells; local_index++)
        {
          unsigned global_index = index_low + local_index;

          unsigned new_local_index = global_index - p_factory->GetLow();
          bool local = p_factory->IsGlobalIndexLocal(global_index);

          boost::shared_ptr<AbstractStimulusFunction> p_stim;
          archive & p_stim;//get from archive

          if (local)
          {
              rStimuli[new_local_index] = p_stim; // Add stimulus to local cells
          }
          //otherwise we should delete, but I think shared pointers delete themselves?
        }
    }
};

 // Declare identifier for the serializer
 #include "SerializationExportWrapper.hpp"
 EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExtendedBidomainTissue)

 namespace boost
 {
 namespace serialization
 {

 template<class Archive, unsigned SPACE_DIM>
 inline void save_construct_data(
     Archive & ar, const ExtendedBidomainTissue<SPACE_DIM> * t, const unsigned int file_version)
 {
     //archive the conductivity tensor of the second cell (which may not be dealt with by heartconfig)
     c_vector<double, SPACE_DIM>  intracellular_conductivities_second_cell = t->GetIntracellularConductivitiesSecondCell();
     //note that simple: ar & intracellular_conductivities_second_cell may not be liked by some boost versions
     for (unsigned i = 0; i < SPACE_DIM; i++)
     {
         ar & intracellular_conductivities_second_cell(i);
     }

     const AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_mesh = t->pGetMesh();
     ar & p_mesh;

     // Don't use the std::vector serialization for cardiac cells, so that we can load them
     // more cleverly when migrating checkpoints.
     t->SaveExtendedBidomainCells(ar, file_version);
     t->SaveExtracellularStimulus(ar, file_version);

     // Creation of conductivity tensors are called by constructor and uses HeartConfig. So make sure that it is
     // archived too (needs doing before construction so appears here instead of usual archive location).
     HeartConfig* p_config = HeartConfig::Instance();
     ar & *p_config;
     ar & p_config;
 }

 /**
  * Allow us to not need a default constructor, by specifying how Boost should
  * instantiate an instance (using existing constructor)
  */
 template<class Archive, unsigned SPACE_DIM>
 inline void load_construct_data(
     Archive & ar, ExtendedBidomainTissue<SPACE_DIM> * t, const unsigned int file_version)
 {
     //Load conductivities of the conductivity of the second cell.
     c_vector<double, SPACE_DIM>  intra_cond_second_cell;
     //note that simple: ar & intra_cond_second_cell may not be liked by some boost versions
     for (unsigned i = 0; i < SPACE_DIM; i++)
     {
         double cond;
         ar & cond;
         intra_cond_second_cell(i) = cond;
     }


     std::vector<AbstractCardiacCellInterface*> cells_distributed;
     std::vector<AbstractCardiacCellInterface*> cells_distributed_second_cell;
     std::vector<boost::shared_ptr<AbstractStimulusFunction> > extra_stim;
     std::vector<double> g_gaps;
     AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>* p_mesh;
     ar & p_mesh;

     // Load only the cells we actually own
     t->LoadExtendedBidomainCells(
             *ProcessSpecificArchive<Archive>::Get(), file_version, cells_distributed, cells_distributed_second_cell, g_gaps, p_mesh);

     t->LoadExtracellularStimulus(
             *ProcessSpecificArchive<Archive>::Get(), file_version, extra_stim, p_mesh);

     // CreateIntracellularConductivityTensor() is called by AbstractCardiacTissue constructor and uses HeartConfig.
     // (as does CreateExtracellularConductivityTensor). So make sure that it is
     // archived too (needs doing before construction so appears here instead of usual archive location).
     HeartConfig* p_config = HeartConfig::Instance();
     ar & *p_config;
     ar & p_config;

     ::new(t)ExtendedBidomainTissue<SPACE_DIM>(cells_distributed, cells_distributed_second_cell, extra_stim, g_gaps, p_mesh, intra_cond_second_cell);
 }
 }
 } // namespace ...

#endif /*EXTENDEDBIDOMAINTISSUE_HPP_*/
