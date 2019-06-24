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


#ifndef HEARTCONFIG_HPP_
#define HEARTCONFIG_HPP_

#include <string>
#include <vector>
#include <set>
#include <map>
#include <boost/shared_ptr.hpp>

#include "UblasVectorInclude.hpp"

#include "ChasteParameters_2017_1.hpp"

#include "AbstractStimulusFunction.hpp"
#include "AbstractChasteRegion.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "ChasteEllipsoid.hpp"
#include "DistributedTetrahedralMeshPartitionType.hpp"
#include "PetscTools.hpp"
#include "FileFinder.hpp"

#include "ChasteSerialization.hpp"
#include "ChasteSerializationVersion.hpp"
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

namespace cp = chaste::parameters::v2017_1;

// Forward declaration to avoid circular includes
class HeartFileFinder;


/**
 * A singleton class containing configuration parameters for heart simulations.
 *
 * This class wraps the settings from the XML configuration file in a more friendly
 * interface, providing methods to read and write all the settings, and round-trip
 * them to/from XML format.  It also deals with the complexities of supporting
 * multiple versions of CodeSynthesis XSD.
 *
 * chaste_parameters_type is a convenience class created by CodeSynthesis XSD
 */
class HeartConfig
{
private:
    /**
     * Throws if the time steps don't obey constraints (within machine precision)
     * ode_step > 0.0
     * pde_step = n1 * ode_step (where n1 is a positive integer)
     * printing_step = n2 * pde_step (where n2 is a positive integer)
     */
    void CheckTimeSteps() const;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        //Only the Master should be writing the configuration file
        if (PetscTools::AmMaster())
        {
            mpInstance->Write( true );
        }

        // Archive other member variables that don't appear in the XML
        if (version > 1)
        {
            archive & mEpiFraction;
            archive & mEndoFraction;
            archive & mMidFraction;
            archive & mIndexMid;
            archive & mIndexEpi;
            archive & mIndexEndo;
            archive & mUserAskedForCellularTransmuralHeterogeneities;
            archive & mUseMassLumping;
            archive & mUseMassLumpingForPrecond;
            archive & mUseReactionDiffusionOperatorSplitting;
            archive & mBathConductivities;
            archive & mTissueIdentifiers;
            archive & mBathIdentifiers;
            archive & mUseFixedNumberIterations;
            archive & mEvaluateNumItsEveryNSolves;
        }

        PetscTools::Barrier("HeartConfig::save");
    }

    /**
     * Un-archive the object.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        LoadFromCheckpoint();

        // Load other member variables
        if (version > 1)
        {
            archive & mEpiFraction;
            archive & mEndoFraction;
            archive & mMidFraction;
            archive & mIndexMid;
            archive & mIndexEpi;
            archive & mIndexEndo;
            archive & mUserAskedForCellularTransmuralHeterogeneities;
            archive & mUseMassLumping;
            archive & mUseMassLumpingForPrecond;
            archive & mUseReactionDiffusionOperatorSplitting;
            archive & mBathConductivities;
            archive & mTissueIdentifiers;
            archive & mBathIdentifiers;
            archive & mUseFixedNumberIterations;
            archive & mEvaluateNumItsEveryNSolves;
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * This method is called by load() to do the actual work - we don't need the Boost archives
     * since we load from our serialised XML.
     */
    void LoadFromCheckpoint();

    /**
     * When loading a simulation from archive, some parameters can get overridden by the content of the ResumeSimulation
     * element.  This method does that.
     *
     * @param pResumeParameters  the parameters containing the ResumeSimulation element.
     */
    void UpdateParametersFromResumeSimulation(boost::shared_ptr<cp::chaste_parameters_type> pResumeParameters);

public:
    /**
     * Our type for specifying schema location properties: a map from namespace URI
     * to schema URI.  The default namespace is specified by an empty namespace URI.
     */
    typedef std::map<std::string, std::string> SchemaLocationsMap;

private:
    /**
     * Fixed location of schema files for the different Chaste parameters namespaces.
     */
    SchemaLocationsMap mSchemaLocations;

    /**
     * Set default schema locations in the Chaste source tree.
     */
    void SetDefaultSchemaLocations();

public:
    /**
     * Call this method to access the global parameters holder.
     *
     * @return a single instance of the class
     */
    static HeartConfig* Instance();

    /**
     * @param useFixedSchemaLocation  whether to read the schema location from the XML
     *    file (false) or use the schema located at heart/src/io/ChasteParameters.xsd
     *    in the Chaste source tree (or specified with SetFixedSchemaLocations()) (true).
     */
    void SetUseFixedSchemaLocation(bool useFixedSchemaLocation);

    /**
     * Set the schema files to use.
     * Also calls SetUseFixedSchemaLocation(true).
     *
     * @param rSchemaLocations  map from namespace URI to schema URI
     */
    void SetFixedSchemaLocations(const SchemaLocationsMap& rSchemaLocations);

    /**
     * #mpParameters  is set to a new context associated with a parameters file
     * @param rFileName The name of the parameters file
     */
    void SetParametersFile(const std::string& rFileName);

    /**
     * Write out the complete configuration set (ChasteParameters
     * and ChasteDefaults) as an XML file.
     * Note that the location of ChasteParameters.xsd (schema definition)
     * will be hard-coded in the XML file.
     * @param useArchiveLocationInfo  if false, then use self's GetOutputDirectory() and open in *named* subfolder
     *                                if true, then use ArchiveLocationInfo
     * @param subfolderName -- where to store with respect to GetOutputDirectory()
     *
     * @note This method is collective if useArchiveLocationInfo is false
     */
    void Write(bool useArchiveLocationInfo=false, std::string subfolderName="output");

    /**
     * Try to copy the latest version of the schema to the given directory.
     * If we can't find the latest version of the schema, generate a warning.
     *
     * @note Must be called by the master only, preferably in a TRY_IF_MASTER()
     * in case of failure, to handle exceptions nicely.
     *
     * @param rToDirectory  directory to copy to
     */
    void CopySchema(const std::string& rToDirectory);

    /**
     * Utility method to parse an XML parameters file.
     * @param rFileName  Name of XML file
     * @return a pointer to the parameters in a convenience class created by CodeSynthesis XSD
     */
    boost::shared_ptr<cp::chaste_parameters_type> ReadFile(const std::string& rFileName);

    /**
     * Throw away the current instance by resetting unique_ptr #mpInstance to NULL.
     * "New" another #mpInstance
     */
    static void Reset();

    ~HeartConfig(); /**< Destructor*/

    /**
     * @return the Chaste version of a parameters file, given its namespace URI.
     * The version will be encoded as major*1000+minor.
     *
     * @param rNamespaceUri  the namespace URI of the parameters file
     */
    unsigned GetVersionFromNamespace(const std::string& rNamespaceUri);

    ///////////////////////////////////////////////////////////////
    //
    //  Get methods
    //
    ///////////////////////////////////////////////////////////////

    /**
     * @return where the user parameters were read from.
     * The result is undefined if no parameters file has been read.
     */
    FileFinder GetParametersFilePath();

    // Methods for asking the configuration file about which sections are defined.

    /**
     *  Returns whether the configuration file defines a new simulation.
     *
     *  @return is a new simulation?
     */
    bool IsSimulationDefined() const;

    /**
     *  Returns whether the configuration file resumes an archived simulation.
     *
     *  @return is a resumed simulation?
     */
    bool IsSimulationResumed() const;

    // Simulation
    unsigned GetSpaceDimension() const; /**< @return space dimension 1, 2 or 3.*/
    double GetSimulationDuration() const; /**< @return duration of the simulation (ms)*/

    /**
     * cp::domain_type is an xsd convenience class type
     *
     * @return domain type of simulation: bi- or mono-domain
     */
    cp::domain_type GetDomain() const;

    /**
     * Default cardiac cell model to use at all mesh nodes
     * (unless otherwise specified by GetIonicModelRegions).
     * cp::ionic_model_selection_type is generated automatically from the XML Schema.
     *
     * @return  type of model
     */
    cp::ionic_model_selection_type GetDefaultIonicModel() const;

    /**
     * Regions where we need to use a different cell model (think infarction).
     * cp::ionic_model_selection_type is generated automatically from the XML Schema.
     *
     * The supplied vectors are first cleared, then filled in with the information from the
     * parameters files.  On return, both vectors will be the same length (one entry per region).
     *
     * @param rDefinedRegions vector of axis-aligned box regions (one per cellular heterogeneity)
     * @param rIonicModels vector of models (one per cellular heterogeneity)
     */
    template<unsigned DIM>
    void GetIonicModelRegions(std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rDefinedRegions,
                              std::vector<cp::ionic_model_selection_type>& rIonicModels) const;

    /**
     * Set the regions where we need to use a different cell model (think infarction).
     * Unlike the get method, this is currently only supported in 3d.
     * cp::ionic_model_selection_type is generated automatically from the XML Schema.
     *
     * The input standard vectors must be of the same length (one entry per region)
     * otherwise the method throws.
     *
     * @param rDefinedRegions vector of axis-aligned box regions (one per cellular heterogeneity)
     * @param rIonicModels vector of models (one per cellular heterogeneity)
     */
    void SetIonicModelRegions(std::vector<ChasteCuboid<3> >& rDefinedRegions,
                              std::vector<cp::ionic_model_selection_type>& rIonicModels) const;

    bool IsMeshProvided() const; /**< @return true if a mesh file name is given.  (Otherwise it's assumed that this is a cuboid simulation.)*/
    bool GetCreateMesh() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
    bool GetCreateSlab() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
    bool GetCreateSheet() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
    bool GetCreateFibre() const; /**< @return true if it's a cuboid simulation (no mesh on disk)*/
    bool GetLoadMesh() const; /**< @return true if a mesh file name is given and we are expecting to load a mesh from file*/

    /**
     * @param slabDimensions  return vector for the (cuboid) mesh dimensions (cm)
     */
    void GetSlabDimensions(c_vector<double, 3>& slabDimensions) const;
    /**
     * @param sheetDimensions  return vector for the (cuboid) mesh dimensions (cm)
     */
    void GetSheetDimensions(c_vector<double, 2>& sheetDimensions) const;
    /**
     * @param fibreLength  return vector for the (cuboid) mesh dimensions (cm)
     */
    void GetFibreLength(c_vector<double, 1>& fibreLength) const;
    double GetInterNodeSpace() const; /**< @return internode space of cuboid mesh (cm)*/

    std::string GetMeshName() const;/**< @return path/basename of mesh files*/

    cp::media_type GetConductivityMedia() const;/**< @return media (Orthotropic/Axisymmetric/NoFibreOrientation) so that we know whether to read a .ortho/.axi file*/

    /**
     * Return a number of stimulated regions (Axis-aligned boxes)
     * The returned std::vectors are all of the same length
     * @param rStimuliApplied  rStimuliApplied[0] is stimulus for the first region
     * @param rStimulatedAreas  rStimulatedAreas[0] is the first region to be stimulated
     * \todo - do we assume the vectors are initially empty?
     *
     * \todo There is no set method
     */
    template<unsigned DIM>
    void GetStimuli(std::vector<boost::shared_ptr<AbstractStimulusFunction> >& rStimuliApplied,
                    std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rStimulatedAreas) const;

    /**
     * Reads from the XML file the cellular hetrogeneities. It fugures out whether the user specified a cuboid
     * or a transmural-type of hetrogeneities. In the latter case, it stores the percentage values of Epi and Endo layers
     * in two member variables, accessible via get methods. It also checks if the user-supplied numbers are consistent (i.e., positive and add up to less than 1)
     * Return a number of heterogeneous regions for special gating variable changes
     * The returned std::vectors are all of the same length
     * @param rCellHeterogeneityRegions  cellHeterogeneityAreas[0] is the first region
     * @param rScaleFactorGks  scaleFactorGks[0] is a scaling factor for the first region
     * @param rScaleFactorIto  scaleFactorIto[0] is a scaling factor for the first region
     * @param rScaleFactorGkr  scaleFactorGkr[0] is a scaling factor for the first region
     * @param pParameterSettings  specification of named parameters to set on the cell models; each entry is a map
     *     from parameter name to value.
     * \todo - do we assume the vectors are initially empty?
     * \todo There is no set method
     */
    template<unsigned DIM>
    void GetCellHeterogeneities(std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rCellHeterogeneityRegions,
                                std::vector<double>& rScaleFactorGks,
                                std::vector<double>& rScaleFactorIto,
                                std::vector<double>& rScaleFactorGkr,
                                std::vector<std::map<std::string, double> >* pParameterSettings);

    bool GetConductivityHeterogeneitiesProvided() const; /**< @return  true if there are conductivity heterogeneities for GetConductivityHeterogeneities to return*/

    /**
     * @return the value of the flag that tells whether the user asked for cellular transmural heterogeneities
     */
    bool AreCellularTransmuralHeterogeneitiesRequested();

    /**
     * @return the fraction of epicardial layer
     */
    double GetEpiLayerFraction();

    /**
     * @return the fraction of endocardial layer
     */
    double GetEndoLayerFraction();

    /**
     * @return the fraction of endocardial layer
     */
    double GetMidLayerFraction();

    /**
     * @return the index with which the epicardial layer is supplied (i.e., the order it comes in the XML file)
     */
    unsigned GetEpiLayerIndex();

    /**
     * @return the index with which the endocardial layer is supplied (i.e., the order it comes in the XML file)
     */
    unsigned GetEndoLayerIndex();

    /**
     * @return the index with which the midmyocardial layer is supplied (i.e., the order it comes in the XML file)
     */
    unsigned GetMidLayerIndex();


    /**
     * Return a number of heterogeneous regions (Axis-aligned boxes)
     * The returned std::vectors are all of the same length
     * @param conductivitiesHeterogeneityAreas  conductivitiesHeterogeneityAreas[0] is the first region
     * @param intraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param extraConductivities  extraConductivities[0] is conductivity vector for the first region
     * \todo - do we assume the vectors are initially empty?
     */
    template<unsigned DIM>
    void GetConductivityHeterogeneities(std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& conductivitiesHeterogeneityAreas,
                                        std::vector< c_vector<double,3> >& intraConductivities,
                                        std::vector< c_vector<double,3> >& extraConductivities) const;
    std::string GetOutputDirectory() const; /**< @return output directory path name*/

    /**
     * @return  Prefix for files
     * If set to "res" this produces
     * [path]/res.h5
     * [path]/output/res_mesh.pts
     * [path]/output/res_mesh.tri
     * [path]/output/res_parameters.xml  (a copy of this configuration at the end of the simulation)
     * [path]/output/res_times.info
     * [path]/output/res_V.dat
     */
    std::string GetOutputFilenamePrefix() const;

    /**
     * @return true iff any extra output variables have been requested
     */
    bool GetOutputVariablesProvided() const;

    /**
     * @return the extra output variables from the xml file.
     *
     * @param rOutputVariables reference to std::vector to contain the output variables requested.
     *    Note: will be cleared before being filled.
     */
    void GetOutputVariables(std::vector<std::string>& rOutputVariables) const;

    /**
     * @return whether to write output HDF5 file using the original
     * mesh permutation (in situations where a parallel partition may have
     * permuted the node).  The default is to use the new, not original permutation,
     */
    bool GetOutputUsingOriginalNodeOrdering();

    /**
     * Get whether simulation should be checkpointed or not
     *
     * @return archive simulation
     */
    bool GetCheckpointSimulation() const;

    /**
     * Get checkpointing timestep
     *
     * @return checkpointing timestep
     */
    double GetCheckpointTimestep() const;

    /**
     * Get number of checkpoints to keep on disk
     *
     * @return checkpointing timestep
     */
    unsigned GetMaxCheckpointsOnDisk() const;

    // ResumeSimulation
    /**
     * @return directory where the archived simulation to resume is defined
     */
    HeartFileFinder GetArchivedSimulationDir() const;


    // Physiological
    /**
     * 3D version
     * @param rIntraConductivities  DIM-vector for returning intracellular conductivities (mS/cm)
     */
    void GetIntracellularConductivities(c_vector<double, 3>& rIntraConductivities) const;

    /**
     * 2D version
     * @param rIntraConductivities  DIM-vector for returning intracellular conductivities (mS/cm)
     */
    void GetIntracellularConductivities(c_vector<double, 2>& rIntraConductivities) const;

    /**
     * 1D version
     * @param rIntraConductivities  DIM-vector for returning intracellular conductivities (mS/cm)
     */
    void GetIntracellularConductivities(c_vector<double, 1>& rIntraConductivities) const;

    /**
     * 3D version
     * @param rExtraConductivities  DIM-vector for returning extracellular conductivities (mS/cm)
     */
    void GetExtracellularConductivities(c_vector<double, 3>& rExtraConductivities) const;

    /**
     * 2D version
     * @param rExtraConductivities  DIM-vector for returning extracellular conductivities (mS/cm)
     */
    void GetExtracellularConductivities(c_vector<double, 2>& rExtraConductivities) const;

    /**
     * 1D version
     * @param rExtraConductivities  DIM-vector for returning extracellular conductivities (mS/cm)
     */
    void GetExtracellularConductivities(c_vector<double, 1>& rExtraConductivities) const;

    /**
     *  Returns bath conductivities for different regions of the bath. When called without a
     * region identifier, it will return whatever has been defined as BathConductivity
     *
     *  @param bathRegion region identifier
     *  @return bath conductivity (mS/cm)
     */
    double GetBathConductivity(unsigned bathRegion=UINT_MAX) const;

    /**
     *  Gets  region identifiers that have to be considered as cardiac tissue.
     *
     *  @return set of identifiers
     */
    const std::set<unsigned>& rGetTissueIdentifiers();

    /**
     *  Gets region identifiers that have to be considered as bath.
     *
     *  @return set of identifiers
     */
    const std::set<unsigned>& rGetBathIdentifiers();


    double GetSurfaceAreaToVolumeRatio() const; /**< @return surface area to volume ratio chi a.k.a Am for PDE (1/cm)*/

    double GetCapacitance() const; /**< @return surface capacitance Cm for PDE (uF/cm^2)*/

    // Numerical
    double GetOdeTimeStep() const; /**< @return ODE time-step (ms)*/
    double GetPdeTimeStep() const; /**< @return PDE time-step (ms)*/
    double GetPrintingTimeStep() const; /**< @return priting time-step (ms)*/

    bool GetUseAbsoluteTolerance() const; /**< @return true if we are using KSP absolute tolerance*/
    double GetAbsoluteTolerance() const; /**< @return KSP absolute tolerance (or throw if we are using relative)*/

    bool GetUseRelativeTolerance() const; /**< @return true if we are using KSP relative tolerance*/
    double GetRelativeTolerance() const;  /**< @return KSP relative tolerance (or throw if we are using absolute)*/

    const char* GetKSPSolver() const; /**< @return name of -ksp_type from {"gmres", "cg", "symmlq"}*/
    const char* GetKSPPreconditioner() const; /**< @return name of -pc_type from {"jacobi", "bjacobi", "hypre", "ml", "spai", "blockdiagonal", "ldufactorisation", "none"}*/

    DistributedTetrahedralMeshPartitionType::type GetMeshPartitioning() const; /**< @return the mesh partitioning method to use */

    // Adaptivity
    /**
     * Adaptivity is now deprecated.  This method now gives a warning before returning true.
     * @return true if there is an adaptivity section
     */
    bool IsAdaptivityParametersPresent() const;

    // Post processing
    /**
     * @return true if there is a post-processing section
     */
    bool IsPostProcessingSectionPresent() const;

    /**
     * Create a PostProcessing section in the user parameters if one doesn't exist.
     */
    void EnsurePostProcessingSectionPresent();

    /**
     * @return true if any post-processing information has been requested
     */
    bool IsPostProcessingRequested() const;

    /**
     * @return true if APD maps have been requested
     */
    bool IsApdMapsRequested() const;

    /**
     * @param rApdMaps  each entry is a request for a map with
     *  - a percentage in the range [1, 100)
     *  - a threshold (in mV)
     */
    void GetApdMaps(std::vector<std::pair<double,double> >& rApdMaps) const;

    /**
     * @return true if upstroke time maps have been requested
     */
    bool IsUpstrokeTimeMapsRequested() const;
    /**
     * @param rUpstrokeTimeMaps  each entry is a request for a map with
     *  - a threshold (in mV)
     */
    void GetUpstrokeTimeMaps (std::vector<double>& rUpstrokeTimeMaps) const;

    /**
     * @return true maximum upstroke velocity maps have been requested
     */
    bool IsMaxUpstrokeVelocityMapRequested() const;

    /**
     * @param rUpstrokeVelocityMaps  each entry is a request for a map with
     *  - a threshold (in mV, defaulted to -30 mV)
     */
    void GetMaxUpstrokeVelocityMaps(std::vector<double>& rUpstrokeVelocityMaps) const;

    /**
     * @return true if conduction velocity maps have been requested
     */
    bool IsConductionVelocityMapsRequested() const;

    /**
     * @param rConductionVelocityMaps  each entry is a request for a map with
     *  - an index to treat as the source for wave propagation
     */
    void GetConductionVelocityMaps(std::vector<unsigned>& rConductionVelocityMaps) const;

    /**
     * @return true if any nodal time trace is requested
     */
    bool IsAnyNodalTimeTraceRequested() const;

    /**
     * @param rRequestedNodes vector of indices of requested nodes that will be filled in by this function
     */
    void GetNodalTimeTraceRequested(std::vector<unsigned>& rRequestedNodes) const;

    /**
     * @return true iff pseudo-ECG calculation has been requested
     */
    bool IsPseudoEcgCalculationRequested() const;

    /**
     * @param rPseudoEcgElectrodePositions  will be filled in with the positions of electrodes
     * to use in calculating pseudo-ECGs (if any)
     */
    template<unsigned SPACE_DIM>
    void GetPseudoEcgElectrodePositions(std::vector<ChastePoint<SPACE_DIM> >& rPseudoEcgElectrodePositions) const;

    /**
     * @return true if state variable interpolation is used
     */
    bool GetUseStateVariableInterpolation() const;


    // Output visualization

    /** @return whether there is an OutputVisualizer element present. */
    bool IsOutputVisualizerPresent() const;

    /** @return whether to convert the output from HDF5 to meshalyzer readable format */
    bool GetVisualizeWithMeshalyzer() const;

    /** @return whether to convert the output from HDF5 to Cmgui readable format */
    bool GetVisualizeWithCmgui() const;

    /** @return whether to convert the output from HDF5 to Vtk readable format */
    bool GetVisualizeWithVtk() const;

    /** @return whether to convert the output from HDF5 to parallel Vtk readable format */
    bool GetVisualizeWithParallelVtk() const;

    /** @return the number of digits to be output in the converted solution data files. */
    unsigned GetVisualizerOutputPrecision();

    /**
     * @return true if there is an electrodes section
     */
    bool IsElectrodesPresent() const;


    /**
     * Get electrode parameters.
     *
     *  @param rGroundSecondElectrode Whether to ground the second electrode (see class documentation)
     *  @param rIndex The value i when applying the electrodes to x_i=a and x_i=b (a<b)
     *  @param rMagnitude Magnitude of the stimulus
     *  @param rStartTime Switch on time
     *  @param rDuration Duration of the stimulus.
     */
    void GetElectrodeParameters(bool& rGroundSecondElectrode,
                                unsigned& rIndex, double& rMagnitude,
                                double& rStartTime, double& rDuration );

    /**
     * @return whether to use mass lumping in the FE solver or not.
     */
    bool GetUseMassLumping();

    /**
     * @return whether to use mass lumping in the construction of the preconditioner of the FE solver or not.
     */
    bool GetUseMassLumpingForPrecond();

    /**
     *  @return whether to use Strang operator splitting of the reaction and diffusion terms (see
     *  Set method documentation).
     */
    bool GetUseReactionDiffusionOperatorSplitting();

    /**
     *  @return whether to use a fixed number of iterations in the linear solver
     */
    bool GetUseFixedNumberIterationsLinearSolver();

    /**
     *  @return how often perform a solve with residual-based stop criteria in order to decide how many iterations to perform in following linear solves.
     */
    unsigned GetEvaluateNumItsEveryNSolves();


    ///////////////////////////////////////////////////////////////
    //
    //  Set methods
    //
    ///////////////////////////////////////////////////////////////

    // Simulation
    /** Set the configuration dimension
     * @param spaceDimension 1, 2 or 3.
     */
    void SetSpaceDimension(unsigned spaceDimension);

    /**
     * Set the configuration simulation end time.
     * @param simulationDuration end time for the next call to Solve() (in ms).
     */
    void SetSimulationDuration(double simulationDuration);

    /**
     * Set the configuration to run mono or bidomain
     * cp::domain_type is an xsd convenience class type
     *
     * @param rDomain type of simulation bi- mono-domain
     */
    void SetDomain(const cp::domain_type& rDomain);

    /**
     * Set the configuration to place the given cardiac cell models at all mesh nodes
     * (unless otherwise specified by SetIonicModelRegions).
     * cp::ionic_models_available_type is generated automatically from the XML Schema.
     *
     * @param rIonicModel  type of model
     */
    void SetDefaultIonicModel(const cp::ionic_models_available_type& rIonicModel);

    /**
     * Set dimensions of simulation for use with a cuboid mesh generated on the fly.  3-D.
     * @param x  length in 1st dimension (cm)
     * @param y  length in 2nd dimension (cm)
     * @param z  length in 3rd dimension (cm)
     * @param inter_node_space  Spacing in cartesian direction (cm). Diagonals will be longer.
     */
    void SetSlabDimensions(double x, double y, double z, double inter_node_space);

    /**
     * Set dimensions of simulation for use with a cuboid mesh generated on the fly.  2-D.
     * @param x  length in 1st dimension (cm)
     * @param y  length in 2nd dimension (cm)
     * @param inter_node_space  Spacing in cartesian direction (cm). Diagonals will be longer.
     */
    void SetSheetDimensions(double x, double y, double inter_node_space);

    /**
     * Set dimensions of simulation for use with a cuboid mesh generated on the fly.  1-D.
     * @param x  length in 1st dimension (cm)
     * @param inter_node_space  Spacing in cartesian direction (cm).
     */
    void SetFibreLength(double x, double inter_node_space);

    /**
     * Sets the name of a mesh to be read from disk for this simulation
     * @param meshPrefix  path and basename of a set of mesh files (.nodes .ele etc) in triangle/tetget format
     * @param fibreDefinition  if set (Orthotropic/Axisymmetric) then a (.ortho/.axi) file should also be read
     * \todo There is no Get method
     */
    void SetMeshFileName(std::string meshPrefix, cp::media_type fibreDefinition=cp::media_type::NoFibreOrientation);

    /**
     * Set a number of heterogeneous regions (Axis-aligned boxes)
     * It is assumed that the std::vectors are all of the same length
     * @param rConductivityAreas conductivityAreas[0] is the first region
     * @param rIntraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param rExtraConductivities  extraConductivities[0] is conductivity vector for the first region
     */
    void SetConductivityHeterogeneities(std::vector<ChasteCuboid<3> >& rConductivityAreas,
                                        std::vector< c_vector<double,3> >& rIntraConductivities,
                                        std::vector< c_vector<double,3> >& rExtraConductivities);
    /**
     * Set a number of heterogeneous regions (Axis-aligned ellipsoids)
     * It is assumed that the std::vectors are all of the same length
     * @param rConductivityAreas conductivityAreas[0] is the first region
     * @param rIntraConductivities  intraConductivities[0] is conductivity vector for the first region
     * @param rExtraConductivities  extraConductivities[0] is conductivity vector for the first region
     */
    void SetConductivityHeterogeneitiesEllipsoid(std::vector<ChasteEllipsoid<3> >& rConductivityAreas,
                                                 std::vector< c_vector<double,3> >& rIntraConductivities,
                                                 std::vector< c_vector<double,3> >& rExtraConductivities);
    /**
     * @param rOutputDirectory  Full path to output directory (will be created if necessary)
     */
    void SetOutputDirectory(const std::string& rOutputDirectory);

    /**
     * @param rOutputFilenamePrefix  Prefix for files
     * If set to "res" this will produce
     * [path]/res.h5
     * [path]/output/res_mesh.pts
     * [path]/output/res_mesh.tri
     * [path]/output/res_parameters.xml  (a copy of this configuration at the end of the simulation)
     * [path]/output/res_times.info
     * [path]/output/res_V.dat
     */
    void SetOutputFilenamePrefix(const std::string& rOutputFilenamePrefix);

    /**
     * @param rOutputVariables  a vector of std::strings of the names
     * of each variable that should be outputted at each time step.
     *
     * USING THIS METHOD WILL OVERRIDE ANY OUTPUT VARIABLES SET IN THE XML FILE
     */
    void SetOutputVariables(const std::vector<std::string>& rOutputVariables);

    /**
     * This method may set the output HDF5 file to be written using the original
     * mesh permutation (in situations where a parallel partition may have
     * permuted the node).  The default is to use the new, not original permutation,
     * i.e.  useOriginal=false
     *
     * @param useOriginal  whether to use the original permutation
     */
    void SetOutputUsingOriginalNodeOrdering(bool useOriginal);

    /**
     * Set whether the simulation should be checkpointed or not.
     *
     * @param checkpointSimulation whether to do checkpointing
     * @param checkpointTimestep checkpointing timestep
     * @param maxCheckpointsOnDisk maximum number of checkpoint archives to keep on disk
     */
     void SetCheckpointSimulation(bool checkpointSimulation, double checkpointTimestep=-1.0, unsigned maxCheckpointsOnDisk=UINT_MAX);

    // Physiological
    /**
     * 3D version
     * @param rIntraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 3>& rIntraConductivities);

    /**
     * 2D version
     * @param rIntraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 2>& rIntraConductivities);

    /**
     * 1D version
     * @param rIntraConductivities  DIM-vector of intracellular conductivities (mS/cm)
     */
    void SetIntracellularConductivities(const c_vector<double, 1>& rIntraConductivities);

    /**
     * 3D version
     * @param rExtraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 3>& rExtraConductivities);

    /**
     * 2D version
     * @param rExtraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 2>& rExtraConductivities);

    /**
     * 1D version
     * @param rExtraConductivities  DIM-vector of extracellular conductivities (mS/cm)
     */
    void SetExtracellularConductivities(const c_vector<double, 1>& rExtraConductivities);

    /**
     * Set bath default conductivity
     * @param bathConductivity default conductivity for perfusing bath (mS/cm)
     * \todo Is this used anywhere?
     */
    void SetBathConductivity(double bathConductivity);

    /**
     * Set multiple bath conductivities based on element region label (mS/cm)
     *
     * @param bathConductivities map between different bath region identifier and their conductivity (if different from default)
     */
    void SetBathMultipleConductivities(std::map<unsigned, double> bathConductivities);

    /**
     *  Sets which region identifiers have to be considered cardiac tissue and bath.
     *
     *  @param rTissueIds set of identifiers
     *  @param rBathIds set of identifiers
     */
    void SetTissueAndBathIdentifiers(const std::set<unsigned>& rTissueIds, const std::set<unsigned>& rBathIds);

    /**
     *  Sets which region identifiers have to be considered cardiac tissue.
     *
     *  param tissueIds set of identifiers
     * \todo #1703 Think about adding this convenience method either copying the existing BathIds, resetting them out of the way, or making them empty...
     */
    //void SetTissueIdentifiers(const std::set<unsigned>& tissueIds);

    /**
     * Set surface area to volume ratio Am (for PDE)
     * @param ratio (1/cm)
     */
    void SetSurfaceAreaToVolumeRatio(double ratio);

    /**
     * Set surface capacitance Cm (for PDE)
     * @param capacitance (uF/cm^2)
     */
    void SetCapacitance(double capacitance);

    // Numerical
    /**
     * Set the configuration to use ode, pde and printing times of given values
     *
     * The ODE step is used by explicit solver schemes such as ForwardEuler
     * to evolve the ODE system at each node in the mesh.
     * AbstractCardiacCells set their internal timestep to this value
     * in their constructors. (AbstractCvodeCells will ignore this setting
     * and use an adaptive time step scheme between PDE or sampling times.)
     *
     * The PDE time step dictates how long a PDE solve should run before re-evaluating
     * the ODE and cell model states and recalculating current contributions.
     * The ODE time step should be a subdivision of this PDE timestep.
     *
     * The sampling timestep should be a multiple of the PDE timestep, and dictates
     * how frequently the output to file of results should occur.
     *
     * This method calls CheckTimeSteps() to ensure the above compatibility
     * conditions are met.
     *
     * @param odeTimeStep  ode value to use
     * @param pdeTimeStep  pde value to use
     * @param printingTimeStep  printing value to use
     */
    void SetOdePdeAndPrintingTimeSteps(double odeTimeStep, double pdeTimeStep, double printingTimeStep);

    /**
     * Set the configuration to use ode time step of given value,
     * for explicit solver schemes such as ForwardEuler.
     * AbstractCardiacCells set their internal timestep to this
     * in their constructors.
     *
     * AbstractCvodeCells will ignore this setting and use an
     * adaptive time step scheme between PDE or sampling times.
     *
     * Calls CheckTimeSteps via SetOdePdeAndPrintingTimeSteps
     * @param odeTimeStep  the value to use
     */
    void SetOdeTimeStep(double odeTimeStep);

    /** Set the configuration to use pde time of given value
     * Calls CheckTimeSteps via SetOdePdeAndPrintingTimeSteps
     * @param pdeTimeStep  the value to use
     */
    void SetPdeTimeStep(double pdeTimeStep);

    /**
     * Set the configuration to use printing time of given value.
     * The printing time step is how long between timesteps that
     * are written to the HDF5 file.
     *
     * Calls CheckTimeSteps via SetOdePdeAndPrintingTimeSteps
     *
     * @param printingTimeStep  the value to use
     */
     void SetPrintingTimeStep(double printingTimeStep);

    /** Set the configuration to use KSP relative tolerance of given value
     * @param relativeTolerance  the value to use
     */
    void SetUseRelativeTolerance(double relativeTolerance);

    /** Set the configuration to use KSP absolute tolerance of given value
     * @param absoluteTolerance  the value to use
     */
    void SetUseAbsoluteTolerance(double absoluteTolerance);

    /** Set the type of KSP solver as with the flag "-ksp_type"
     * @param kspSolver  a string from {"gmres", "cg", "symmlq"}
     * @param warnOfChange  Warn if this set is changing the current value because the calling
     * code may be (silently) overwriting a user setting
     */
    void SetKSPSolver(const char* kspSolver, bool warnOfChange=false);

    /** Set the type of preconditioner as with the flag "-pc_type"
     * @param kspPreconditioner  a string from {"jacobi", "bjacobi", "hypre", "ml", "spai", "blockdiagonal", "ldufactorisation", "none"}
     */
    void SetKSPPreconditioner(const char* kspPreconditioner);

    /** Set the type of mesh partitioning method
     * @param meshPartioningMethod  a string from {"dumb", "metis", "parmetis", "petsc"}
     */
    void SetMeshPartitioning(const char* meshPartioningMethod);

    /** Set the parameters of the apd map requested
     *
     *  @param rApdMaps  each entry is a request for a map with
     *  - a percentage in the range [1, 100) (ranges are not checked by this method, but during the calculation)
     *  - a threshold (in mV)
     */
    void SetApdMaps(const std::vector<std::pair<double,double> >& rApdMaps);

    /** Set the parameters of the upstroke time map requested
     *
     * @param rUpstrokeTimeMaps  is the list of thresholds (mV) with respect to which the upstroke time maps are calculated.
     *     The threshold is used for determining when an action potential occurs.
     */
    void SetUpstrokeTimeMaps(std::vector<double>& rUpstrokeTimeMaps);

    /** Set the parameters of the maximal upstroke velocity map requested
     *
     *  @param rMaxUpstrokeVelocityMaps is the list of thresholds (mV) with respect to which the upstroke velocity maps are calculated.
     *     The threshold is used for determining when an action potential occurs.
     */
    void SetMaxUpstrokeVelocityMaps(std::vector<double>& rMaxUpstrokeVelocityMaps);

    /** Set the parameters of the conduction velocity map requested
     *
     *  @param rConductionVelocityMaps is a list of origin node indices. One map is created for each origin node.
     */
    void SetConductionVelocityMaps(std::vector<unsigned>& rConductionVelocityMaps);

    /**
     * Sets some requested nodes for printing of their variables over time in separate files.
     * The node numbering is referred to the original numbering (unpermuted).
     *
     * @param requestedNodes the node indices (in the unpermuted mesh) where we want the plot over time
     */
    void SetRequestedNodalTimeTraces(std::vector<unsigned>& requestedNodes);

    /** Set the parameters for pseudo-ECG calculation.
     *
     * @param rPseudoEcgElectrodePositions  should contan the positions of electrodes
     * to use in calculating pseudo-ECGs (if any)
     */
    template<unsigned SPACE_DIM>
    void SetPseudoEcgElectrodePositions(const std::vector<ChastePoint<SPACE_DIM> >& rPseudoEcgElectrodePositions);


    // Output visualization

    /** Create the OutputVisualizer element if it doesn't exist */
    void EnsureOutputVisualizerExists(void);

    /** Set whether to convert the output from HDF5 to meshalyzer readable format
     *
     * @param useMeshalyzer
     */
    void SetVisualizeWithMeshalyzer(bool useMeshalyzer=true);

    /** Set whether to convert the output from HDF5 to Cmgui readable format
     *
     * @param useCmgui
     */
    void SetVisualizeWithCmgui(bool useCmgui=true);

    /** Set whether to convert the output from HDF5 to Vtk readable format
     *
     * @param useVtk
     */
    void SetVisualizeWithVtk(bool useVtk=true);

    /** Set whether to convert the output from HDF5 to parallel Vtk readable format
     *
     * @param useParallelVtk
     */
    void SetVisualizeWithParallelVtk(bool useParallelVtk=true);

    /**
     * Set the precision with which to output textual visualizer formats (e.g. meshalyzer).
     * Use '0' for the implementation-defined default precision.
     *
     * @param numberOfDigits  how many digits of precision to use
     */
    void SetVisualizerOutputPrecision(unsigned numberOfDigits);

    /**
     * Setup electrode parameters.
     *
     *  @param groundSecondElectrode Whether to ground the second electrode (see class documentation)
     *  @param index The value i when applying the electrodes to x_i=a and x_i=b (a<b)
     *  @param magnitude Magnitude of the stimulus
     *  @param startTime Switch on time
     *  @param duration Duration of the stimulus.
     */
    void SetElectrodeParameters(bool groundSecondElectrode,
                                unsigned index, double magnitude,
                                double startTime, double duration);

   /**
     * Set the use of State Variable Interpolation in the computation of ionic currents.
     * See documentation page ChasteGuides/StateVariableInterpolation.
     *
     * @param  useStateVariableInterpolation Whether to use it.
     */
    void SetUseStateVariableInterpolation(bool useStateVariableInterpolation = true);

    /**
     * Set the use of mass lumping in the FE solver.
     *
     * @param useMassLumping Whether to use it
     */
    void SetUseMassLumping(bool useMassLumping = true);

    /**
     * Set the use of mass lumping in the construction of the preconditioner in the FE solver.
     *
     * @param useMassLumping Whether to use it
     */
    void SetUseMassLumpingForPrecond(bool useMassLumping = true);

    /**
     * Use Strang operator splitting of the diffusion (conductivity) term and the reaction (ionic current) term,
     * instead of solving the full reaction-diffusion PDE. This does NOT refer to operator splitting of the
     * two PDEs in the bidomain equations. For details see for example Sundnes et al "Computing the Electrical
     * Activity of the Heart".
     *
     * @param useOperatorSplitting Whether to use operator splitting (defaults to true).
     */
    void SetUseReactionDiffusionOperatorSplitting(bool useOperatorSplitting = true);

    /**
     * Set the use of fixed number of iterations in the linear solver
     *
     * @param useFixedNumberIterations Whether to use a fixed number of iterations for the linear solver
     * @param evaluateNumItsEveryNSolves Perform a solve with convergence-based stop criteria every n solves to decide how many iterations perform for the next n-1 solves. Default is perfoming a single evaluation at the beginning of the simulation.
     */
    void SetUseFixedNumberIterationsLinearSolver(bool useFixedNumberIterations = true, unsigned evaluateNumItsEveryNSolves=UINT_MAX);

    /**
     * @return whether HeartConfig has a drug concentration and any IC50s set up
     */
    bool HasDrugDose() const;

    /**
     * @return the dose of the drug (in the same units as the IC50s)
     */
    double GetDrugDose() const;

    /**
     * @param drugDose  The dose of the drug to use (should be in units consistent with the IC50s).
     */
    void SetDrugDose(double drugDose);

    /**
     * Add a new conductance block model for a particular channel.
     *
     * @param rCurrentName  The Oxford metadata name of the current (e.g. membrane_fast_sodium_current)
     * @param ic50  The IC50 value for this channel (should be in consistent units with drug dose)
     * @param hill  The hill coefficient to use (usually default to 1)
     */
    void SetIc50Value(const std::string& rCurrentName, double ic50, double hill=1.0);

    /**
     * Get the parameters for the model of "conductance-block" drug action on a set of ion channels.
     *
     * @return  A map between the current/channel name, and a pair giving IC50 value and hill coefficient.
     */
    std::map<std::string, std::pair<double, double> > GetIc50Values();

    //
    // Purkinje-related methods
    //

    /**
     * @return whether this simulation contains a Purkinje system.
     */
    bool HasPurkinje();

    /**
     * @return the surface capacitance for Purkinje myocytes.
     */
    double GetPurkinjeCapacitance();

    /**
     * Set the surface capacitance for Purkinje myocytes.
     * @param capacitance  Purkinje capacitance (Cm) (units uF/cm^2)
     */
    void SetPurkinjeCapacitance(double capacitance);

    /**
     * @return the surface area to volume ratio for Purkinje fibres.
     */
    double GetPurkinjeSurfaceAreaToVolumeRatio();

    /**
     * Set the surface area to volume ratio for Purkinje fibres.
     * @param ratio  the ratio (Am) (units 1/cm)
     */
    void SetPurkinjeSurfaceAreaToVolumeRatio(double ratio);

    /**
     * @return the default conductivity for Purkinje fibres.
     */
    double GetPurkinjeConductivity();

    /**
     * Set the default conductivity for Purkinje fibres.
     * @param conductivity  Purkinje conductivity (units mS/cm)
     */
    void SetPurkinjeConductivity(double conductivity);

private:
    // Only to be accessed by the tests
    friend class TestHeartConfig;

    /*Constructor is private, since the class is only accessed by the singleton instance() method*/
    HeartConfig();

    /** Pointer to parameters read from the user's input XML file  */
    boost::shared_ptr<cp::chaste_parameters_type> mpParameters;

    /** The single instance of the class */
    static boost::shared_ptr<HeartConfig> mpInstance;

    /**
     * Where the user parameters were read from.
     */
    FileFinder mParametersFilePath;

    /**
     * @return whether to read the schema location from the XML file (false) or use the schema
     * located at heart/src/io/ChasteParameters.xsd in the Chaste source tree (true).
     */
    bool mUseFixedSchemaLocation;

    /**
     * Fraction of epicardial layer
     */
    double mEpiFraction;

    /**
     * Fraction of endocardial layer
     */
    double mEndoFraction;

    /**
     * Fraction of midmyocardial layer
     */
    double mMidFraction;

    /**
     * Order index in which the midmyocardial heterogeneities are supplied
     */
    unsigned mIndexMid;

    /**
     * Order index in which the epicardial heterogeneities are supplied
     */
    unsigned mIndexEpi;

    /**
     * Order index in which the endocardial heterogeneities are supplied
     */
    unsigned mIndexEndo;

    /**
     * Flag to check whether the user asked for cellular transmural heterogeneities
     */
    bool mUserAskedForCellularTransmuralHeterogeneities;

    /**
     * Flag telling whether to use mass lumping or not.
     */
    bool mUseMassLumping;

    /**
     * Flag telling whether to use mass lumping in the preconditioner or not.
     */
    bool mUseMassLumpingForPrecond;

    /**
     *  @return whether to use Strang operator splitting of the diffusion and reaction terms (see
     *  Set method documentation).
     */
    bool mUseReactionDiffusionOperatorSplitting;

    /**
     *  Map defining bath conductivity for multiple bath regions
     */
    std::map<unsigned, double> mBathConductivities;

    /**
     * Mesh region identifiers to be considered as cardiac tissue
     */
    std::set<unsigned> mTissueIdentifiers;

    /**
     * Mesh region identifiers to be considered as Bath
     */
    std::set<unsigned> mBathIdentifiers;

    /**
     * @return whether to use a fixed number of iterations for the linear solver
     */
    bool mUseFixedNumberIterations;

    /**
     * Perform a solve with convergence-based stop criteria every n solves
     * to decide how many iterations perform for the next n-1 solves. Default
     * is performing a single evaluation at the beginning of the simulation.
     */
    unsigned mEvaluateNumItsEveryNSolves;

    /**
     * CheckSimulationIsDefined is a convenience method for checking if the "<"Simulation">" element
     * has been defined and therefore is safe to use the Simulation().get() pointer to access
     * other data.
     *
     * Throws an exception if not.
     *
     * @param callingMethod string describing the get method performing the check.
     */
    void CheckSimulationIsDefined(std::string callingMethod="") const;

    /**
     * CheckSimulationIsDefined is a convenience method for checking if the "<"ResumeSimulation">" element
     * has been defined and therefore is safe to use the ResumeSimulation().get() pointer to access
     * other data.
     *
     * Throws an exception if not.
     *
     * @param callingMethod string describing the get method performing the check.
     */
    void CheckResumeSimulationIsDefined(std::string callingMethod="") const;
};


BOOST_CLASS_VERSION(HeartConfig, 2)
#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(HeartConfig)

#endif /*HEARTCONFIG_HPP_*/
