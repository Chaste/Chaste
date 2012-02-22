/*

Copyright (C) University of Oxford, 2005-2012

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

// Work-around for newer Boost versions
#include "CheckpointArchiveTypesIfNeeded.hpp"

#include "UblasCustomFunctions.hpp"

#include "HeartConfig.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
#include "ChastePoint.hpp"
#include "Version.hpp"
#include "AbstractChasteRegion.hpp"
#include "HeartFileFinder.hpp"
#include "Warnings.hpp"

#include "HeartRegionCodes.hpp"

#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"

#include <string>
#include <istream>
#include <fstream>
#include <cassert>
#include <map>

#include "XmlTools.hpp"
#include <xsd/cxx/tree/exceptions.hxx>
using namespace xsd::cxx::tree;

// Coping with changes to XSD interface
#if (XSD_INT_VERSION >= 3000000L)
#define XSD_SEQUENCE_TYPE(base) base##_sequence
#define XSD_ITERATOR_TYPE(base) base##_iterator
#define XSD_NESTED_TYPE(t) t##_type
#define XSD_ANON_TYPE(t1, t2) \
    t1::t2##_type
#else
#define XSD_SEQUENCE_TYPE(base) base::container
#define XSD_ITERATOR_TYPE(base) base::iterator
#define XSD_NESTED_TYPE(t) t::type
#define XSD_ANON_TYPE(t1, t2) \
    t1::t2::_xsd_##t2##_::t2
#endif

// These are for convenience
#define XSD_ANON_SEQUENCE_TYPE(t1, t2, t3) \
    XSD_SEQUENCE_TYPE(XSD_ANON_TYPE(t1, t2)::t3)
#define XSD_ANON_ITERATOR_TYPE(t1, t2, t3) \
    XSD_ITERATOR_TYPE(XSD_ANON_TYPE(t1, t2)::t3)

// Newer versions don't allow you to set fixed attributes
#if (XSD_INT_VERSION >= 3020000L)
#define XSD_CREATE_WITH_FIXED_ATTR(type, name, attr) \
    type name
#define XSD_CREATE_WITH_FIXED_ATTR1(type, name, arg1, attr) \
    type name(arg1)
#define XSD_CREATE_WITH_FIXED_ATTR2(type, name, arg1, arg2, attr) \
    type name(arg1, arg2)
#define XSD_CREATE_WITH_FIXED_ATTR3(type, name, arg1, arg2, arg3, attr) \
    type name(arg1, arg2, arg3)
#else
#define XSD_CREATE_WITH_FIXED_ATTR(type, name, attr) \
    type name(attr)
#define XSD_CREATE_WITH_FIXED_ATTR1(type, name, arg1, attr) \
    type name(arg1, attr)
#define XSD_CREATE_WITH_FIXED_ATTR2(type, name, arg1, arg2, attr) \
    type name(arg1, arg2, attr)
#define XSD_CREATE_WITH_FIXED_ATTR3(type, name, arg1, arg2, arg3, attr) \
    type name(arg1, arg2, arg3, attr)
#endif

/**
 * This gets used by various set methods.
 */
#define ENSURE_SECTION_PRESENT(location, type) \
    if (!location.present())                   \
    {                                          \
        type empty_item;                       \
        location.set(empty_item);              \
    }


/**
 * A class of utility methods for transforming parameters files from previous versions
 * of the Chaste schema to the latest version.
 */
class XmlTransforms
{
public:
    /**
     * Edits the DOM tree to wrap ionic model definitions from old (release 1 or 1.1)
     * configuration files in a 'Hardcoded' element.
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pRootElement  the root of the tree to be transformed
     */
    static void TransformIonicModelDefinitions(xercesc::DOMDocument* pDocument,
                                               xercesc::DOMElement* pRootElement);

    /**
     * Edits the DOM tree to change the 'ArchiveDirectory' element from a simple string
     * to a cp::path_type.  This is used for 2.0 -> 2.1 migration.  We assume that the
     * path is relative to CHASTE_TEST_OUTPUT.
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pRootElement  the root of the tree to be transformed
     */
    static void TransformArchiveDirectory(xercesc::DOMDocument* pDocument,
                                          xercesc::DOMElement* pRootElement);

    /**
     * Release 2.1 removes the ilu preconditioner as an option, so throw an exception
     * if it is used.
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pRootElement  the root of the tree to be transformed
     */
    static void CheckForIluPreconditioner(xercesc::DOMDocument* pDocument,
                                          xercesc::DOMElement* pRootElement);

    /**
     * Release 3.1 moved the ConductivityHeterogeneities element from Simulation to
     * Physiological, to be next to the default conductivity definitions.
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pRootElement  the root of the tree to be transformed
     */
    static void MoveConductivityHeterogeneities(xercesc::DOMDocument* pDocument,
                                                xercesc::DOMElement* pRootElement);
};

//
// Default settings
//
#include "HeartConfigDefaults.hpp"

//
// Definition of static member variables
//
std::auto_ptr<HeartConfig> HeartConfig::mpInstance;

//
// Methods
//

HeartConfig* HeartConfig::Instance()
{
    if (mpInstance.get() == NULL)
    {
        mpInstance.reset(new HeartConfig);
    }
    return mpInstance.get();
}

HeartConfig::HeartConfig()
    : mUseMassLumping(false),
      mUseMassLumpingForPrecond(false),
      mUseFixedNumberIterations(false),
      mEvaluateNumItsEveryNSolves(UINT_MAX)
{
    assert(mpInstance.get() == NULL);
    mUseFixedSchemaLocation = true;
    SetDefaultSchemaLocations();

    mpDefaultParameters = CreateDefaultParameters();
    mpUserParameters = mpDefaultParameters;
    //CheckTimeSteps(); // necessity of this line of code is not tested -- remove with caution!

    //initialise the member variable of the layers
    mEpiFraction = -1.0;
    mEndoFraction =  -1.0;
    mMidFraction = -1.0;
    mUserAskedForCellularTransmuralHeterogeneities = false;
    // initialise to senseless values (these should be only 0, 1 and 2)
    // note: the 'minus 3' is for checking purposes as we need to add 0, 1 or 2 to this initial value
    // and UINT_MAX+1 seems to be 0
    mIndexMid = UINT_MAX-3u;
    mIndexEpi = UINT_MAX-3u;
    mIndexEndo = UINT_MAX-3u;

    mUseReactionDiffusionOperatorSplitting = false;

    /// \todo #1703 This defaults should be set in HeartConfigDefaults.hpp
    mTissueIdentifiers.insert(0);
    mBathIdentifiers.insert(1);
}

HeartConfig::~HeartConfig()
{
}

void HeartConfig::SetDefaultsFile(const std::string& rFileName)
{
    bool same_target = (mpUserParameters == mpDefaultParameters);

    mpDefaultParameters = ReadFile(rFileName);

    if (same_target)
    {
        mpUserParameters = mpDefaultParameters;
    }
    CheckTimeSteps();
}

void HeartConfig::Write(bool useArchiveLocationInfo, std::string subfolderName)
{
    //Output file
    std::string output_dirname;
    if (useArchiveLocationInfo)
    {
        output_dirname = ArchiveLocationInfo::GetArchiveDirectory();
    }
    else
    {
        OutputFileHandler handler(GetOutputDirectory() + "/" + subfolderName, false);
        output_dirname = handler.GetOutputDirectoryFullPath();
    }
    if (!PetscTools::AmMaster())
    {
        //Only the master process is writing the configuration files
        return;
    }
    out_stream p_defaults_file( new std::ofstream( (output_dirname+"ChasteDefaults.xml").c_str() ) );
    out_stream p_parameters_file( new std::ofstream( (output_dirname+"ChasteParameters.xml").c_str() ) );

    if (!p_defaults_file->is_open() || !p_parameters_file->is_open())
    {
        EXCEPTION("Could not open XML file in HeartConfig");
    }

    //Schema map
    //Note - this location is relative to where we are storing the xml
    ::xml_schema::namespace_infomap map;
    // Release 1.1 (and earlier) didn't use a namespace
    map[""].schema = "ChasteParameters_1_1.xsd";
    // Later releases use namespaces of the form https://chaste.comlab.ox.ac.uk/nss/parameters/N_M
    map["cp20"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/2_0";
    map["cp20"].schema = "ChasteParameters_2_0.xsd";
    map["cp21"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/2_1";
    map["cp21"].schema = "ChasteParameters_2_1.xsd";
    map["cp22"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/2_2";
    map["cp22"].schema = "ChasteParameters_2_2.xsd";
    map["cp23"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/2_3";
    map["cp23"].schema = "ChasteParameters_2_3.xsd";
    map["cp30"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/3_0";
    map["cp30"].schema = "ChasteParameters_3_0.xsd";
    // We use 'cp' as prefix for the latest version to avoid having to change saved
    // versions for comparison at every release.
    map["cp"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/3_1";
    map["cp"].schema = "ChasteParameters_3_1.xsd";

    cp::ChasteParameters(*p_parameters_file, *mpUserParameters, map);
    cp::ChasteParameters(*p_defaults_file, *mpDefaultParameters, map);

    // If we're archiving, try to save a copy of the latest schema too
    if (useArchiveLocationInfo)
    {
        CopySchema(output_dirname);
    }
}

void HeartConfig::CopySchema(const std::string& rToDirectory)
{
    if (PetscTools::AmMaster())
    {
        std::string schema_name("ChasteParameters_3_1.xsd");
        FileFinder schema_location("heart/src/io/" + schema_name, RelativeTo::ChasteSourceRoot);
        if (!schema_location.Exists())
        {
            // Try a relative path instead
            schema_location.SetPath(schema_name, RelativeTo::CWD);
            if (!schema_location.Exists())
            {
                // Warn the user
                std::string message("Unable to locate schema file " + schema_name +
                                    ". You will need to ensure it is available when resuming from the checkpoint.");
                WARN_ONCE_ONLY(message);
            }
        }
        if (schema_location.Exists())
        {
            // Called by master only so use this instead of EXPECT0
            ABORT_IF_NON0(system,("cp " + schema_location.GetAbsolutePath() + " " + rToDirectory).c_str());
        }
    }
}

void HeartConfig::SetDefaultSchemaLocations()
{
    mSchemaLocations.clear();
    // Location of schemas in the source tree
    std::string root_dir = std::string(ChasteBuildInfo::GetRootDir()) + "/heart/src/io/";
    // Release 1.1 (and earlier) didn't use a namespace
    mSchemaLocations[""] = root_dir + "ChasteParameters_1_1.xsd";
    // Later releases use namespaces of the form https://chaste.comlab.ox.ac.uk/nss/parameters/N_M
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_0"] = root_dir + "ChasteParameters_2_0.xsd";
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_1"] = root_dir + "ChasteParameters_2_1.xsd";
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_2"] = root_dir + "ChasteParameters_2_2.xsd";
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_3"] = root_dir + "ChasteParameters_2_3.xsd";
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_0"] = root_dir + "ChasteParameters_3_0.xsd";
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_1"] = root_dir + "ChasteParameters_3_1.xsd";
}

unsigned HeartConfig::GetVersionFromNamespace(const std::string& rNamespaceUri)
{
    unsigned version_major = 0;
    unsigned version_minor = 0;
    if (rNamespaceUri == "")
    {
        version_major = 1;
        version_minor = 1;
    }
    else
    {
        std::string uri_base("https://chaste.comlab.ox.ac.uk/nss/parameters/");
        if (rNamespaceUri.substr(0, uri_base.length()) == uri_base)
        {
            std::istringstream version_string(rNamespaceUri.substr(uri_base.length()));
            version_string >> version_major;
            version_string.ignore(1);
            version_string >> version_minor;
            if (version_string.fail())
            {
                version_major = 0;
                version_minor = 0;
            }
        }
    }

    unsigned version = version_major * 1000 + version_minor;
    if (version == 0)
    {
        EXCEPTION(rNamespaceUri + " is not a recognised Chaste parameters namespace.");
    }
    return version;
}

void HeartConfig::SetFixedSchemaLocations(const SchemaLocationsMap& rSchemaLocations)
{
    mSchemaLocations = rSchemaLocations;
    SetUseFixedSchemaLocation(true);
}

void HeartConfig::SetUseFixedSchemaLocation(bool useFixedSchemaLocation)
{
    mUseFixedSchemaLocation = useFixedSchemaLocation;
}

boost::shared_ptr<cp::chaste_parameters_type> HeartConfig::ReadFile(const std::string& rFileName)
{
    // Determine whether to use the schema path given in the input XML, or our own schema
    ::xml_schema::properties props;
    if (mUseFixedSchemaLocation)
    {
        for (SchemaLocationsMap::iterator it = mSchemaLocations.begin();
             it != mSchemaLocations.end();
             ++it)
        {
            if (it->first == "")
            {
                props.no_namespace_schema_location(XmlTools::EscapeSpaces(it->second));
            }
            else
            {
                props.schema_location(it->first, XmlTools::EscapeSpaces(it->second));
            }
        }
    }

    // Get the parameters using the method 'ChasteParameters(rFileName)',
    // which returns a std::auto_ptr. We convert to a shared_ptr for easier semantics.
    try
    {
        // Make sure Xerces finalization happens
        XmlTools::Finalizer finalizer(false);
        // Parse XML to DOM
        xsd::cxx::xml::dom::auto_ptr<xercesc::DOMDocument> p_doc = XmlTools::ReadXmlFile(rFileName, props);
        // Test the namespace on the root element
        xercesc::DOMElement* p_root_elt = p_doc->getDocumentElement();
        std::string namespace_uri(X2C(p_root_elt->getNamespaceURI()));
        const unsigned version = GetVersionFromNamespace(namespace_uri);
        if (version < 2000) // Changes made in release 2.0
        {
            XmlTransforms::TransformIonicModelDefinitions(p_doc.get(), p_root_elt);
        }
        if (version < 2001) // Changes made in release 2.1
        {
            XmlTransforms::TransformArchiveDirectory(p_doc.get(), p_root_elt);
            XmlTransforms::CheckForIluPreconditioner(p_doc.get(), p_root_elt);
        }
        if (version < 3001) // Changes made in release 3.1
        {
            XmlTransforms::MoveConductivityHeterogeneities(p_doc.get(), p_root_elt);
        }
        if (version < 3001) // Not the latest release
        {
            XmlTools::SetNamespace(p_doc.get(), p_root_elt, "https://chaste.comlab.ox.ac.uk/nss/parameters/3_1");
        }
        // Parse DOM to object model
        std::auto_ptr<cp::chaste_parameters_type> p_params(cp::ChasteParameters(*p_doc, ::xml_schema::flags::dont_initialize, props));
        // Get rid of the DOM stuff
        p_doc.reset();

        return boost::shared_ptr<cp::chaste_parameters_type>(p_params);
    }
    catch (const xml_schema::exception& e)
    {
        std::cerr << e << std::endl;
        // Make sure we don't store invalid parameters
        mpUserParameters.reset();
        mpDefaultParameters.reset();
        EXCEPTION("XML parsing error in configuration file: " + rFileName);
    }
    catch (...)
    {
        // Make sure we don't store invalid parameters
        mpUserParameters.reset();
        mpDefaultParameters.reset();
        throw;
    }
}

void HeartConfig::SetParametersFile(const std::string& rFileName)
{
    mpUserParameters = ReadFile(rFileName);
    mParametersFilePath.SetPath(rFileName, RelativeTo::AbsoluteOrCwd);

    CheckTimeSteps(); // For consistency with SetDefaultsFile
}

FileFinder HeartConfig::GetParametersFilePath()
{
    return mParametersFilePath;
}

void HeartConfig::UpdateParametersFromResumeSimulation(boost::shared_ptr<cp::chaste_parameters_type> pResumeParameters)
{
    // Check for user foolishness
    if ( (pResumeParameters->ResumeSimulation()->SpaceDimension() != HeartConfig::Instance()->GetSpaceDimension())
         ||(pResumeParameters->ResumeSimulation()->Domain() != HeartConfig::Instance()->GetDomain()))
    {
        EXCEPTION("Problem type and space dimension should match when restarting a simulation.");
    }

    // New simulation duration
    HeartConfig::Instance()->SetSimulationDuration(pResumeParameters->ResumeSimulation()->SimulationDuration());

    // Stimulus definition.  For these we always replace any previous definitions (at least for now...)
    if (pResumeParameters->ResumeSimulation()->Stimuli().present())
    {
        mpUserParameters->Simulation()->Stimuli().set(pResumeParameters->ResumeSimulation()->Stimuli().get());
    }

    // Cell heterogeneities.  Note that while we copy the elements here, other code in CardiacSimulation actually updates
    // the loaded simulation to take account of the new settings.
    if (pResumeParameters->ResumeSimulation()->CellHeterogeneities().present())
    {
        if (!mpUserParameters->Simulation()->CellHeterogeneities().present())
        {
            // Original parameters had no heterogeneities, so just copy the whole element
            mpUserParameters->Simulation()->CellHeterogeneities().set(pResumeParameters->ResumeSimulation()->CellHeterogeneities().get());
        }
        else
        {
            // Need to append the new heterogeneity defitions to the original sequence
            XSD_SEQUENCE_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity)&
                new_seq = pResumeParameters->ResumeSimulation()->CellHeterogeneities()->CellHeterogeneity();
            XSD_SEQUENCE_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity)&
                orig_seq = mpUserParameters->Simulation()->CellHeterogeneities()->CellHeterogeneity();
            for (XSD_ITERATOR_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity) i = new_seq.begin();
                 i != new_seq.end();
                 ++i)
            {
                orig_seq.push_back(*i);
            }
        }
    }

    // Whether to checkpoint the resumed simulation
    if (pResumeParameters->ResumeSimulation()->CheckpointSimulation().present())
    {
        HeartConfig::Instance()->SetCheckpointSimulation(true,
                                                         pResumeParameters->ResumeSimulation()->CheckpointSimulation()->timestep(),
                                                         pResumeParameters->ResumeSimulation()->CheckpointSimulation()->max_checkpoints_on_disk());
    }

    //Visualization parameters are compulsory
    HeartConfig::Instance()->SetVisualizeWithParallelVtk(pResumeParameters->ResumeSimulation()->OutputVisualizer().parallel_vtk() == cp::yesno_type::yes);
    HeartConfig::Instance()->SetVisualizeWithVtk(pResumeParameters->ResumeSimulation()->OutputVisualizer().vtk() == cp::yesno_type::yes);
    HeartConfig::Instance()->SetVisualizeWithCmgui(pResumeParameters->ResumeSimulation()->OutputVisualizer().cmgui() == cp::yesno_type::yes);
    HeartConfig::Instance()->SetVisualizeWithMeshalyzer(pResumeParameters->ResumeSimulation()->OutputVisualizer().meshalyzer() == cp::yesno_type::yes);

    // Numerical parameters may be overridden
    {
        cp::numerical_type& r_resume = pResumeParameters->Numerical();
        cp::numerical_type& r_user = mpUserParameters->Numerical();
        if (r_resume.TimeSteps().present())
        {
            r_user.TimeSteps().set(r_resume.TimeSteps().get());
        }
        if (r_resume.KSPTolerances().present())
        {
            r_user.KSPTolerances().set(r_resume.KSPTolerances().get());
        }
        if (r_resume.KSPSolver().present())
        {
            r_user.KSPSolver().set(r_resume.KSPSolver().get());
        }
        if (r_resume.KSPPreconditioner().present())
        {
            r_user.KSPPreconditioner().set(r_resume.KSPPreconditioner().get());
        }
        if (r_resume.AdaptivityParameters().present())
        {
            r_user.AdaptivityParameters().set(r_resume.AdaptivityParameters().get());
        }
    }

    // Post-processing parameters may be overridden
    if (pResumeParameters->PostProcessing().present())
    {
        EnsurePostProcessingSectionPresent();
        cp::postprocessing_type& r_resume = pResumeParameters->PostProcessing().get();
        cp::postprocessing_type& r_user = mpUserParameters->PostProcessing().get();
        if (!r_resume.ActionPotentialDurationMap().empty())
        {
            r_user.ActionPotentialDurationMap() = r_resume.ActionPotentialDurationMap();
        }
        if (!r_resume.UpstrokeTimeMap().empty())
        {
            r_user.UpstrokeTimeMap() = r_resume.UpstrokeTimeMap();
        }
        if (!r_resume.MaxUpstrokeVelocityMap().empty())
        {
            r_user.MaxUpstrokeVelocityMap() = r_resume.MaxUpstrokeVelocityMap();
        }
        if (!r_resume.ConductionVelocityMap().empty())
        {
            r_user.ConductionVelocityMap() = r_resume.ConductionVelocityMap();
        }
    }
}


void HeartConfig::Reset()
{
    // Throw it away first, so that mpInstance is NULL when we...
    mpInstance.reset();
    // ...make a new one
    mpInstance.reset(new HeartConfig);
}

bool HeartConfig::IsSimulationDefined() const
{
    return mpUserParameters->Simulation().present();
}

bool HeartConfig::IsSimulationResumed() const
{
    return mpUserParameters->ResumeSimulation().present();
}


template<class TYPE>
TYPE* HeartConfig::DecideLocation(TYPE* ptr1, TYPE* ptr2, const std::string& nameParameter) const
{
    if (ptr1->present())
    {
        return ptr1;
    }
    if (ptr2->present())
    {
        return ptr2;
    }

    EXCEPTION("No " + nameParameter + " provided (neither default nor user defined)");
}

void HeartConfig::CheckSimulationIsDefined(std::string callingMethod) const
{
    if (IsSimulationResumed())
    {
        EXCEPTION(callingMethod + " information is not available in a resumed simulation.");
    }
}

void HeartConfig::CheckResumeSimulationIsDefined(std::string callingMethod) const
{
    if (IsSimulationDefined())
    {
        EXCEPTION(callingMethod + " information is not available in a standard (non-resumed) simulation.");
    }
}

unsigned HeartConfig::GetSpaceDimension() const
{
    if (IsSimulationDefined())
    {
        return DecideLocation( & mpUserParameters->Simulation()->SpaceDimension(),
                               & mpDefaultParameters->Simulation()->SpaceDimension(),
                               "SpaceDimension")->get();
    }
    else
    {
        return mpUserParameters->ResumeSimulation()->SpaceDimension();
    }
}

double HeartConfig::GetSimulationDuration() const
{
    if (IsSimulationDefined())
    {
        return DecideLocation( & mpUserParameters->Simulation()->SimulationDuration(),
                               & mpDefaultParameters->Simulation()->SimulationDuration(),
                               "Simulation/SimulationDuration")->get();
    }
    else // IsSimulationResumed
    {
        return mpUserParameters->ResumeSimulation()->SimulationDuration();
    }
}

cp::domain_type HeartConfig::GetDomain() const
{
    if (IsSimulationDefined())
    {
        return DecideLocation( & mpUserParameters->Simulation()->Domain(),
                               & mpDefaultParameters->Simulation()->Domain(),
                               "Domain")->get();
    }
    else
    {
        return mpUserParameters->ResumeSimulation()->Domain();
    }
}

cp::ionic_model_selection_type HeartConfig::GetDefaultIonicModel() const
{
    CheckSimulationIsDefined("DefaultIonicModel");

    return DecideLocation( & mpUserParameters->Simulation()->IonicModels(),
                           & mpDefaultParameters->Simulation()->IonicModels(),
                           "IonicModel")->get().Default();
}

template<unsigned DIM>
void HeartConfig::GetIonicModelRegions(std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& definedRegions,
                                       std::vector<cp::ionic_model_selection_type>& ionicModels) const
{
    CheckSimulationIsDefined("IonicModelRegions");
    definedRegions.clear();
    ionicModels.clear();

    XSD_SEQUENCE_TYPE(cp::ionic_models_type::Region)&
         regions = DecideLocation( & mpUserParameters->Simulation()->IonicModels(),
                                   & mpDefaultParameters->Simulation()->IonicModels(),
                                   "IonicModel")->get().Region();

    for (XSD_ITERATOR_TYPE(cp::ionic_models_type::Region) i = regions.begin();
         i != regions.end();
         ++i)
    {
        cp::ionic_model_region_type ionic_model_region(*i);

        if (ionic_model_region.Location().Cuboid().present() || ionic_model_region.Location().Ellipsoid().present())
        {
            if (ionic_model_region.Location().Cuboid().present())
            {
                cp::point_type point_a = ionic_model_region.Location().Cuboid()->LowerCoordinates();
                cp::point_type point_b = ionic_model_region.Location().Cuboid()->UpperCoordinates();

                switch (DIM)
                {
                    case 1:
                    {
                        ChastePoint<DIM> chaste_point_a ( point_a.x() );
                        ChastePoint<DIM> chaste_point_b ( point_b.x() );
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>( chaste_point_a, chaste_point_b )));
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a ( point_a.x(), point_a.y() );
                        ChastePoint<DIM> chaste_point_b ( point_b.x(), point_b.y() );
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>( chaste_point_a, chaste_point_b )));
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a ( point_a.x(), point_a.y(), point_a.z() );
                        ChastePoint<DIM> chaste_point_b ( point_b.x(), point_b.y(), point_b.z() );
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>( chaste_point_a, chaste_point_b )) );
                        break;
                    }
                    default:
                        NEVER_REACHED;
                        break;
                }
            }
            else if (ionic_model_region.Location().Ellipsoid().present())
            {
                cp::point_type centre = ionic_model_region.Location().Ellipsoid()->Centre();
                cp::point_type radii  = ionic_model_region.Location().Ellipsoid()->Radii();
                switch (DIM)
                {
                    case 1:
                    {
                        ChastePoint<DIM> chaste_point_a ( centre.x() );
                        ChastePoint<DIM> chaste_point_b ( radii.x() );
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b )) );
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a ( centre.x(), centre.y() );
                        ChastePoint<DIM> chaste_point_b ( radii.x(), radii.y() );
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b )) );
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a ( centre.x(), centre.y(), centre.z() );
                        ChastePoint<DIM> chaste_point_b ( radii.x(), radii.y(), radii.z() );
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b )) );
                        break;
                    }
                    default:
                    {
                        NEVER_REACHED;
                        break;
                    }
                }
            }
            else
            {
                NEVER_REACHED;
            }

            ionicModels.push_back(ionic_model_region.IonicModel());
        }
        else if(ionic_model_region.Location().EpiLayer().present() || ionic_model_region.Location().MidLayer().present() || ionic_model_region.Location().EndoLayer().present() )
        {
            ///\todo When this is implemented, then we require an example in ChasteParametersFullFormat.xml
            EXCEPTION("Definition of transmural layers is not yet supported for defining different ionic models, please use cuboids instead");
        }
        else
        {
            EXCEPTION("Invalid region type for ionic model definition");
        }
    }
}


bool HeartConfig::IsMeshProvided() const
{
    CheckSimulationIsDefined("Mesh");

    try
    {
        DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                        & mpDefaultParameters->Simulation()->Mesh(),
                        "Mesh");
        return true;
    }
    catch (Exception& e)
    {
        return false;
    }
}

bool HeartConfig::GetCreateMesh() const
{
    CheckSimulationIsDefined("Mesh");

    cp::mesh_type mesh = DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                         & mpDefaultParameters->Simulation()->Mesh(),
                                         "Mesh")->get();

    return (mesh.Slab().present() || mesh.Sheet().present() || mesh.Fibre().present());
}

bool HeartConfig::GetCreateSlab() const
{
    CheckSimulationIsDefined("Mesh");

    cp::mesh_type mesh = DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                         & mpDefaultParameters->Simulation()->Mesh(),
                                         "Mesh")->get();

    return (mesh.Slab().present());
}

bool HeartConfig::GetCreateSheet() const
{
    CheckSimulationIsDefined("Mesh");

    cp::mesh_type mesh = DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                         & mpDefaultParameters->Simulation()->Mesh(),
                                         "Mesh")->get();

    return (mesh.Sheet().present());
}

bool HeartConfig::GetCreateFibre() const
{
    CheckSimulationIsDefined("Mesh");

    cp::mesh_type mesh = DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                         & mpDefaultParameters->Simulation()->Mesh(),
                                         "Mesh")->get();

    return (mesh.Fibre().present());
}


bool HeartConfig::GetLoadMesh() const
{
    CheckSimulationIsDefined("Mesh");

    return (DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                            & mpDefaultParameters->Simulation()->Mesh(),
                            "Mesh")->get().LoadMesh().present());
}

void HeartConfig::GetSlabDimensions(c_vector<double, 3>& slabDimensions) const
{
    CheckSimulationIsDefined("Slab");

    if (GetSpaceDimension()!=3 || !GetCreateSlab())
    {
        EXCEPTION("Tissue slabs can only be defined in 3D");
    }

    optional<cp::slab_type, false> slab_dimensions = DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                                                     & mpDefaultParameters->Simulation()->Mesh(),
                                                                     "Slab")->get().Slab();

    slabDimensions[0] = slab_dimensions->x();
    slabDimensions[1] = slab_dimensions->y();
    slabDimensions[2] = slab_dimensions->z();
}

void HeartConfig::GetSheetDimensions(c_vector<double, 2>& sheetDimensions) const
{
    CheckSimulationIsDefined("Sheet");

    if (GetSpaceDimension()!=2 || !GetCreateSheet())
    {
        EXCEPTION("Tissue sheets can only be defined in 2D");
    }

    optional<cp::sheet_type, false> sheet_dimensions = DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                                                       & mpDefaultParameters->Simulation()->Mesh(),
                                                                       "Sheet")->get().Sheet();

    sheetDimensions[0] = sheet_dimensions->x();
    sheetDimensions[1] = sheet_dimensions->y();
}

void HeartConfig::GetFibreLength(c_vector<double, 1>& fibreLength) const
{
    CheckSimulationIsDefined("Fibre");

    if (GetSpaceDimension()!=1 || !GetCreateFibre())
    {
        EXCEPTION("Tissue fibres can only be defined in 1D");
    }

    optional<cp::fibre_type, false> fibre_length = DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                                                   & mpDefaultParameters->Simulation()->Mesh(),
                                                                   "Fibre")->get().Fibre();

    fibreLength[0] = fibre_length->x();
}

double HeartConfig::GetInterNodeSpace() const
{
    CheckSimulationIsDefined("InterNodeSpace");
    assert(GetCreateMesh());

    switch(GetSpaceDimension())
    {
        case 3:
            return DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                   & mpDefaultParameters->Simulation()->Mesh(),
                                   "Slab")->get().Slab()->inter_node_space();
            break;
        case 2:
            return DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                   & mpDefaultParameters->Simulation()->Mesh(),
                                   "Sheet")->get().Sheet()->inter_node_space();
            break;
        case 1:
            return DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                                   & mpDefaultParameters->Simulation()->Mesh(),
                                   "Fibre")->get().Fibre()->inter_node_space();
            break;
        default:
            NEVER_REACHED;
#define COVERAGE_IGNORE
            return 0.0; //To fool the compiler
#undef COVERAGE_IGNORE
    }
}

std::string HeartConfig::GetMeshName() const
{
    CheckSimulationIsDefined("LoadMesh");
    assert(GetLoadMesh());

    return DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                           & mpDefaultParameters->Simulation()->Mesh(),
                           "LoadMesh")->get().LoadMesh()->name();
}

cp::media_type HeartConfig::GetConductivityMedia() const
{
    CheckSimulationIsDefined("LoadMesh");
    assert(GetLoadMesh());

    return DecideLocation( & mpUserParameters->Simulation()->Mesh(),
                           & mpDefaultParameters->Simulation()->Mesh(),
                           "LoadMesh")->get().LoadMesh()->conductivity_media();
}

template <unsigned DIM>
void HeartConfig::GetStimuli(std::vector<boost::shared_ptr<AbstractStimulusFunction> >& rStimuliApplied,
                             std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rStimulatedAreas) const
{
    CheckSimulationIsDefined("Stimuli");
    XSD_SEQUENCE_TYPE(cp::stimuli_type::Stimulus) stimuli;


    try
    {
         stimuli = DecideLocation( & mpUserParameters->Simulation()->Stimuli(),
                                   & mpDefaultParameters->Simulation()->Stimuli(),
                                   "Stimuli")->get().Stimulus();
    }
    catch(Exception& e)
    {
        // Finding no stimuli defined is allowed (although HeartConfigRelatedFactory does
        // throw an exception is no stimuli and no electrodes)
        return;
    }

    for (XSD_ITERATOR_TYPE(cp::stimuli_type::Stimulus) i = stimuli.begin();
         i != stimuli.end();
         ++i)
    {
        cp::stimulus_type stimulus(*i);
        if (stimulus.Location().Cuboid().present() || stimulus.Location().Ellipsoid().present())
        {
            boost::shared_ptr<AbstractChasteRegion<DIM> > area_ptr;
            if (stimulus.Location().Cuboid().present() )
            {
                cp::point_type point_a = stimulus.Location().Cuboid()->LowerCoordinates();
                cp::point_type point_b = stimulus.Location().Cuboid()->UpperCoordinates();
                switch (DIM)
                {
                    case 1:
                    {
                        ChastePoint<DIM> chaste_point_a ( point_a.x() );
                        ChastePoint<DIM> chaste_point_b ( point_b.x() );
                        area_ptr.reset(new ChasteCuboid<DIM>( chaste_point_a, chaste_point_b ) );
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a ( point_a.x(), point_a.y() );
                        ChastePoint<DIM> chaste_point_b ( point_b.x(), point_b.y() );
                        area_ptr.reset(new ChasteCuboid<DIM>( chaste_point_a, chaste_point_b ) );
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a ( point_a.x(), point_a.y(), point_a.z() );
                        ChastePoint<DIM> chaste_point_b ( point_b.x(), point_b.y(), point_b.z() );
                        area_ptr.reset(new ChasteCuboid<DIM>( chaste_point_a, chaste_point_b ) );
                        break;
                    }
                    default:
                        NEVER_REACHED;
                        break;
                }
            }
            else if (stimulus.Location().Ellipsoid().present())
            {
                cp::point_type centre = stimulus.Location().Ellipsoid()->Centre();
                cp::point_type radii  = stimulus.Location().Ellipsoid()->Radii();
                switch (DIM)
                {
                    case 1:
                    {
                        ChastePoint<DIM> chaste_point_a ( centre.x() );
                        ChastePoint<DIM> chaste_point_b ( radii.x() );
                        area_ptr.reset( new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b ) );
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a ( centre.x(), centre.y() );
                        ChastePoint<DIM> chaste_point_b ( radii.x(), radii.y() );
                        area_ptr.reset( new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b ) );
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a ( centre.x(), centre.y(), centre.z() );
                        ChastePoint<DIM> chaste_point_b ( radii.x(), radii.y(), radii.z() );
                        area_ptr.reset( new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b ) );
                        break;
                    }
                    default:
                    {
                        NEVER_REACHED;
                        break;
                    }
                }
            }
            rStimulatedAreas.push_back(area_ptr);

            boost::shared_ptr<AbstractStimulusFunction> stim;

            if (stimulus.Period().present())
            {
                if (stimulus.StopTime().present())
                {
                    stim.reset(new RegularStimulus(stimulus.Strength(),
                                                   stimulus.Duration(),
                                                   stimulus.Period().get(),
                                                   stimulus.Delay(),
                                                   stimulus.StopTime().get()));
                }
                else
                {
                    stim.reset(new RegularStimulus(stimulus.Strength(),
                                                   stimulus.Duration(),
                                                   stimulus.Period().get(),
                                                   stimulus.Delay()));
                }

            }
            else
            {
                if (stimulus.StopTime().present())
                {
                    EXCEPTION("Stop time can not be defined for SimpleStimulus. Use Duration instead.");
                }

                stim.reset(new SimpleStimulus(stimulus.Strength(),
                                              stimulus.Duration(),
                                              stimulus.Delay()));
            }
            rStimuliApplied.push_back( stim );
        }
        else if(stimulus.Location().EpiLayer().present() || stimulus.Location().MidLayer().present() || stimulus.Location().EndoLayer().present() )
        {
            EXCEPTION("Definition of transmural layers is not yet supported for specifying stimulated areas, please use cuboids instead");
        }
        else
        {
            EXCEPTION("Invalid region type for stimulus definition");
        }
    }
}


template<unsigned DIM>
void HeartConfig::GetCellHeterogeneities(std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rCellHeterogeneityRegions,
                                         std::vector<double>& rScaleFactorGks,
                                         std::vector<double>& rScaleFactorIto,
                                         std::vector<double>& rScaleFactorGkr,
                                         std::vector<std::map<std::string, double> >* pParameterSettings)
{
    CheckSimulationIsDefined("CellHeterogeneities");
    XSD_SEQUENCE_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity) cell_heterogeneity;

    try
    {
         cell_heterogeneity = DecideLocation( & mpUserParameters->Simulation()->CellHeterogeneities(),
                                              & mpDefaultParameters->Simulation()->CellHeterogeneities(),
                                              "CellHeterogeneities")->get().CellHeterogeneity();
    }
    catch(Exception& e)
    {
        // finding no heterogeneities defined is allowed
        return;
    }

    bool user_supplied_negative_value = false;
    bool user_asking_for_transmural_layers = false;
    bool user_asked_for_cuboids_or_ellipsoids = false;
    unsigned counter_of_heterogeneities = 0;

    for (XSD_ITERATOR_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity) i = cell_heterogeneity.begin();
         i != cell_heterogeneity.end();
         ++i)
    {
        cp::cell_heterogeneity_type ht(*i);


        if (ht.Location().Cuboid().present())
        {
            user_asked_for_cuboids_or_ellipsoids = true;
            cp::point_type point_a = ht.Location().Cuboid()->LowerCoordinates();
            cp::point_type point_b = ht.Location().Cuboid()->UpperCoordinates();

            ChastePoint<DIM> chaste_point_a ( point_a.x(), point_a.y(), point_a.z() );
            ChastePoint<DIM> chaste_point_b ( point_b.x(), point_b.y(), point_b.z() );

            rCellHeterogeneityRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM> ( chaste_point_a, chaste_point_b )) );
        }
        else if (ht.Location().Ellipsoid().present())
        {
            user_asked_for_cuboids_or_ellipsoids = true;
            cp::point_type centre = ht.Location().Ellipsoid()->Centre();
            cp::point_type radii  = ht.Location().Ellipsoid()->Radii();

            ChastePoint<DIM> chaste_point_a ( centre.x(), centre.y(), centre.z() );
            ChastePoint<DIM> chaste_point_b ( radii.x(), radii.y(), radii.z() );
            rCellHeterogeneityRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b )) );
        }
        else if (ht.Location().EpiLayer().present())
        {
            mEpiFraction  =  ht.Location().EpiLayer().get();

            user_asking_for_transmural_layers = true;
            if (mEpiFraction <0)
            {
                user_supplied_negative_value=true;
            }
            mIndexEpi = counter_of_heterogeneities;
        }
        else if (ht.Location().EndoLayer().present())
        {
            mEndoFraction  =  ht.Location().EndoLayer().get();

            user_asking_for_transmural_layers = true;
            if (mEndoFraction <0)
            {
                user_supplied_negative_value=true;
            }
            mIndexEndo = counter_of_heterogeneities;
        }
        else if (ht.Location().MidLayer().present())
        {
            mMidFraction  =  ht.Location().MidLayer().get();

            user_asking_for_transmural_layers = true;
            if (mMidFraction <0)
            {
                user_supplied_negative_value=true;
            }
            mIndexMid =  counter_of_heterogeneities;
        }
        else
        {
            EXCEPTION("Invalid region type for cell heterogeneity definition");
        }

        // Old scale factors
        rScaleFactorGks.push_back(ht.ScaleFactorGks().present() ? (double)ht.ScaleFactorGks().get() : 1.0);
        rScaleFactorIto.push_back(ht.ScaleFactorIto().present() ? (double)ht.ScaleFactorIto().get() : 1.0);
        rScaleFactorGkr.push_back(ht.ScaleFactorGkr().present() ? (double)ht.ScaleFactorGkr().get() : 1.0);

        // Named parameters
        if (pParameterSettings)
        {
            std::map<std::string, double> param_settings;
            XSD_SEQUENCE_TYPE(cp::cell_heterogeneity_type::SetParameter)& params = ht.SetParameter();
            for (XSD_ITERATOR_TYPE(cp::cell_heterogeneity_type::SetParameter) param_it = params.begin();
                 param_it != params.end();
                 ++param_it)
            {
                cp::set_parameter_type param(*param_it);
                param_settings[param.name()] = param.value();
            }
            pParameterSettings->push_back(param_settings);
        }

        counter_of_heterogeneities++;
    }

    //set the flag for request of transmural layers
    mUserAskedForCellularTransmuralHeterogeneities = user_asking_for_transmural_layers;

    if ( mUserAskedForCellularTransmuralHeterogeneities )
    {
        // cuboids/ellipsoids and layers at the same time are not yet supported
        if (user_asked_for_cuboids_or_ellipsoids )
        {
            EXCEPTION ("Specification of cellular heterogeneities by cuboids/ellipsoids and layers at the same time is not yet supported");
        }

        //check that the user supplied all three layers, the indexes should be 0, 1 and 2.
        // As they are initialised to a higher value, if their summation is higher than 3,
        // one (or more) is missing
        if ((mIndexMid+mIndexEndo+mIndexEpi) > 3)
        {
            EXCEPTION ("Three specifications of layers must be supplied");
        }
        if (fabs((mEndoFraction+mMidFraction+mEpiFraction)-1)>1e-2)
        {
            EXCEPTION ("Summation of epicardial, midmyocardial and endocardial fractions should be 1");
        }
        if (user_supplied_negative_value)
        {
           EXCEPTION ("Fractions must be positive");
        }
    }
}

bool HeartConfig::AreCellularTransmuralHeterogeneitiesRequested()
{
    return mUserAskedForCellularTransmuralHeterogeneities;
}

double HeartConfig::GetEpiLayerFraction()
{
    return mEpiFraction;
}

double HeartConfig::GetEndoLayerFraction()
{
    return mEndoFraction;
}

double HeartConfig::GetMidLayerFraction()
{
    return mMidFraction;
}

unsigned HeartConfig::GetEpiLayerIndex()
{
    return mIndexEpi;
}

unsigned HeartConfig::GetEndoLayerIndex()
{
    return mIndexEndo;
}

unsigned HeartConfig::GetMidLayerIndex()
{
    return mIndexMid;
}

bool HeartConfig::GetConductivityHeterogeneitiesProvided() const
{
    CheckSimulationIsDefined("ConductivityHeterogeneities");
    try
    {
        DecideLocation( & mpUserParameters->Physiological().ConductivityHeterogeneities(),
                        & mpDefaultParameters->Physiological().ConductivityHeterogeneities(),
                        "ConductivityHeterogeneities");
        return true;
    }
    catch (Exception& e)
    {
        return false;
    }
}

template<unsigned DIM>
void HeartConfig::GetConductivityHeterogeneities(
        std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& conductivitiesHeterogeneityAreas,
        std::vector< c_vector<double,3> >& intraConductivities,
        std::vector< c_vector<double,3> >& extraConductivities) const
{
    CheckSimulationIsDefined("ConductivityHeterogeneities");
    XSD_ANON_SEQUENCE_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity)&
         conductivity_heterogeneity = DecideLocation( & mpUserParameters->Physiological().ConductivityHeterogeneities(),
                                                      & mpDefaultParameters->Physiological().ConductivityHeterogeneities(),
                                                      "ConductivityHeterogeneities")->get().ConductivityHeterogeneity();

    for (XSD_ANON_ITERATOR_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity) i = conductivity_heterogeneity.begin();
         i != conductivity_heterogeneity.end();
         ++i)
    {
        cp::conductivity_heterogeneity_type ht(*i);

        if (ht.Location().Cuboid().present())
        {
            cp::point_type point_a = ht.Location().Cuboid()->LowerCoordinates();
            cp::point_type point_b = ht.Location().Cuboid()->UpperCoordinates();
            ChastePoint<DIM> chaste_point_a ( point_a.x(), point_a.y(), point_a.z() );
            ChastePoint<DIM> chaste_point_b ( point_b.x(), point_b.y(), point_b.z() );
            conductivitiesHeterogeneityAreas.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >   (  new ChasteCuboid<DIM> ( chaste_point_a, chaste_point_b )) );
        }
        else if (ht.Location().Ellipsoid().present())
        {
            cp::point_type centre = ht.Location().Ellipsoid()->Centre();
            cp::point_type radii  = ht.Location().Ellipsoid()->Radii();
            ChastePoint<DIM> chaste_point_a ( centre.x(), centre.y(), centre.z() );
            ChastePoint<DIM> chaste_point_b ( radii.x(), radii.y(), radii.z() );
            conductivitiesHeterogeneityAreas.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >   (  new ChasteEllipsoid<DIM> ( chaste_point_a, chaste_point_b )) );
        }
        else if (ht.Location().EpiLayer().present() || ht.Location().MidLayer().present() || ht.Location().EndoLayer().present() )
        {
            ///\todo When this is implemented, then we require an example in ChasteParametersFullFormat.xml
            EXCEPTION("Definition of transmural layers is not allowed for conductivities heterogeneities, you may use fibre orientation support instead");
        }
        else
        {
            EXCEPTION("Invalid region type for conductivity definition");
        }

        if (ht.IntracellularConductivities().present())
        {
            double intra_x = ht.IntracellularConductivities()->longi();
            double intra_y = ht.IntracellularConductivities()->trans();
            double intra_z = ht.IntracellularConductivities()->normal();

            intraConductivities.push_back( Create_c_vector(intra_x, intra_y, intra_z) );
        }
        else
        {
            c_vector<double, 3> intra_conductivities;
            GetIntracellularConductivities(intra_conductivities);
            intraConductivities.push_back(intra_conductivities);
        }

        if (ht.ExtracellularConductivities().present())
        {
            double extra_x = ht.ExtracellularConductivities()->longi();
            double extra_y = ht.ExtracellularConductivities()->trans();
            double extra_z = ht.ExtracellularConductivities()->normal();

            extraConductivities.push_back( Create_c_vector(extra_x, extra_y, extra_z) );
        }
        else
        {
            c_vector<double, 3> extra_conductivities;
            GetExtracellularConductivities(extra_conductivities);
            extraConductivities.push_back(extra_conductivities);
        }

    }
}

std::string HeartConfig::GetOutputDirectory() const
{
    CheckSimulationIsDefined("Simulation/OutputDirectory");
    return DecideLocation( & mpUserParameters->Simulation()->OutputDirectory(),
                           & mpDefaultParameters->Simulation()->OutputDirectory(),
                           "Simulation/OutputDirectory")->get();
}

std::string HeartConfig::GetOutputFilenamePrefix() const
{
    CheckSimulationIsDefined("Simulation/OutputFilenamePrefix");
    return DecideLocation( & mpUserParameters->Simulation()->OutputFilenamePrefix(),
                           & mpDefaultParameters->Simulation()->OutputFilenamePrefix(),
                           "Simulation/OutputFilenamePrefix")->get();
}

bool HeartConfig::GetOutputVariablesProvided() const
{
    CheckSimulationIsDefined("OutputVariables");

    try
    {
        DecideLocation( & mpUserParameters->Simulation()->OutputVariables(),
                        & mpDefaultParameters->Simulation()->OutputVariables(),
                        "OutputVariables");
        return true;
    }
    catch (Exception& e)
    {
        return false;
    }
}

void HeartConfig::GetOutputVariables(std::vector<std::string>& rOutputVariables) const
{
    CheckSimulationIsDefined("OutputVariables");
    XSD_SEQUENCE_TYPE(cp::output_variables_type::Var)&
         output_variables = DecideLocation( & mpUserParameters->Simulation()->OutputVariables(),
                                            & mpDefaultParameters->Simulation()->OutputVariables(),
                                            "OutputVariables")->get().Var();
    rOutputVariables.clear();

    for (XSD_ITERATOR_TYPE(cp::output_variables_type::Var) i = output_variables.begin();
         i != output_variables.end();
         ++i)
    {
        cp::var_type& r_var(*i);

        // Add to outputVariables the string returned by var.name()
        rOutputVariables.push_back(r_var.name());
    }
}
bool HeartConfig::GetOutputUsingOriginalNodeOrdering()
{
    try
    {
        return (DecideLocation(& mpUserParameters->Simulation()->OutputUsingOriginalNodeOrdering(),
                               & mpDefaultParameters->Simulation()->OutputUsingOriginalNodeOrdering(),
                              "OutputUsingOriginalNodeOrdering")->get()  == cp::yesno_type::yes);
    }
    catch (Exception &e)
    {
        //If it didn't exist, then we default to false
        return false;
    }
}

bool HeartConfig::GetCheckpointSimulation() const
{
    try
    {
        if (IsSimulationDefined())
        {
            CheckSimulationIsDefined("GetCheckpointSimulation");
            DecideLocation(& mpUserParameters->Simulation()->CheckpointSimulation(),
                           & mpDefaultParameters->Simulation()->CheckpointSimulation(),
                           "CheckpointSimulation");
        }
        else
        {
            CheckResumeSimulationIsDefined("GetCheckpointSimulation");
            DecideLocation(& mpUserParameters->ResumeSimulation()->CheckpointSimulation(),
                           & mpDefaultParameters->Simulation()->CheckpointSimulation(),
                           "CheckpointSimulation");
        }
        return true;
    }
    catch (Exception& e)
    {
        return false;
    }
}

double HeartConfig::GetCheckpointTimestep() const
{
    if (IsSimulationDefined())
    {
        CheckSimulationIsDefined("GetCheckpointTimestep");
        return DecideLocation(& mpUserParameters->Simulation()->CheckpointSimulation(),
                              & mpDefaultParameters->Simulation()->CheckpointSimulation(),
                              "CheckpointSimulation")->get().timestep();
    }
    else
    {
        CheckResumeSimulationIsDefined("GetCheckpointTimestep");
        return DecideLocation(& mpUserParameters->ResumeSimulation()->CheckpointSimulation(),
                              & mpDefaultParameters->Simulation()->CheckpointSimulation(),
                              "CheckpointSimulation")->get().timestep();
    }
}

unsigned HeartConfig::GetMaxCheckpointsOnDisk() const
{
    if (IsSimulationDefined())
    {
        CheckSimulationIsDefined("GetMaxCheckpointsOnDisk");
        return DecideLocation(& mpUserParameters->Simulation()->CheckpointSimulation(),
                              & mpDefaultParameters->Simulation()->CheckpointSimulation(),
                              "CheckpointSimulation")->get().max_checkpoints_on_disk();
    }
    else
    {
        CheckResumeSimulationIsDefined("GetMaxCheckpointsOnDisk");
        return DecideLocation(& mpUserParameters->ResumeSimulation()->CheckpointSimulation(),
                              & mpDefaultParameters->Simulation()->CheckpointSimulation(),
                              "CheckpointSimulation")->get().max_checkpoints_on_disk();
    }
}


HeartFileFinder HeartConfig::GetArchivedSimulationDir() const
{
    CheckResumeSimulationIsDefined("GetArchivedSimulationDir");

    return HeartFileFinder(mpUserParameters->ResumeSimulation()->ArchiveDirectory());
}


void HeartConfig::GetIntracellularConductivities(c_vector<double, 3>& intraConductivities) const
{
    optional<cp::conductivities_type, false>* intra_conductivities
        = DecideLocation( & mpUserParameters->Physiological().IntracellularConductivities(),
                          & mpDefaultParameters->Physiological().IntracellularConductivities(),
                          "IntracellularConductivities");
    double intra_x_cond = intra_conductivities->get().longi();
    double intra_y_cond = intra_conductivities->get().trans();
    double intra_z_cond = intra_conductivities->get().normal();;

    assert(intra_y_cond != DBL_MAX);
    assert(intra_z_cond != DBL_MAX);

    intraConductivities[0] = intra_x_cond;
    intraConductivities[1] = intra_y_cond;
    intraConductivities[2] = intra_z_cond;
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 2>& intraConductivities) const
{
    optional<cp::conductivities_type, false>* intra_conductivities
        = DecideLocation( & mpUserParameters->Physiological().IntracellularConductivities(),
                          & mpDefaultParameters->Physiological().IntracellularConductivities(),
                          "IntracellularConductivities");
    double intra_x_cond = intra_conductivities->get().longi();
    double intra_y_cond = intra_conductivities->get().trans();

    assert(intra_y_cond != DBL_MAX);

    intraConductivities[0] = intra_x_cond;
    intraConductivities[1] = intra_y_cond;
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 1>& intraConductivities) const
{
    optional<cp::conductivities_type, false>* intra_conductivities
        = DecideLocation( & mpUserParameters->Physiological().IntracellularConductivities(),
                          & mpDefaultParameters->Physiological().IntracellularConductivities(),
                          "IntracellularConductivities");
    double intra_x_cond = intra_conductivities->get().longi();

    intraConductivities[0] = intra_x_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 3>& extraConductivities) const
{
    optional<cp::conductivities_type, false>* extra_conductivities
        = DecideLocation( & mpUserParameters->Physiological().ExtracellularConductivities(),
                          & mpDefaultParameters->Physiological().ExtracellularConductivities(),
                          "ExtracellularConductivities");
    double extra_x_cond = extra_conductivities->get().longi();
    double extra_y_cond = extra_conductivities->get().trans();
    double extra_z_cond = extra_conductivities->get().normal();;

    assert(extra_y_cond != DBL_MAX);
    assert(extra_z_cond != DBL_MAX);

    extraConductivities[0] = extra_x_cond;
    extraConductivities[1] = extra_y_cond;
    extraConductivities[2] = extra_z_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 2>& extraConductivities) const
{
    optional<cp::conductivities_type, false>* extra_conductivities
        = DecideLocation( & mpUserParameters->Physiological().ExtracellularConductivities(),
                          & mpDefaultParameters->Physiological().ExtracellularConductivities(),
                          "ExtracellularConductivities");
    double extra_x_cond = extra_conductivities->get().longi();
    double extra_y_cond = extra_conductivities->get().trans();

    assert(extra_y_cond != DBL_MAX);

    extraConductivities[0] = extra_x_cond;
    extraConductivities[1] = extra_y_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 1>& extraConductivities) const
{
    optional<cp::conductivities_type, false>* extra_conductivities
        = DecideLocation( & mpUserParameters->Physiological().ExtracellularConductivities(),
                          & mpDefaultParameters->Physiological().ExtracellularConductivities(),
                          "ExtracellularConductivities");
    double extra_x_cond = extra_conductivities->get().longi();

    extraConductivities[0] = extra_x_cond;
}

double HeartConfig::GetBathConductivity(unsigned bathRegion) const
{
    /*
     *  We have to consider three cases: The user asks for ...
     *    a) ... the default conductivity (bathRegion=UINT_MAX)
     *    b) ... the conductivity of region defined to be heterogeneous
     *    c) ... the conductivity of region NOT defined to be heterogeneous
     *
     *  a) and c) should return the same
     */

    if (bathRegion == UINT_MAX)
    {
        /*bath conductivity mS/cm*/
        return DecideLocation( & mpUserParameters->Physiological().BathConductivity(),
                               & mpDefaultParameters->Physiological().BathConductivity(),
                               "BathConductivity")->get();
    }
    else
    {
        assert(HeartRegionCode::IsRegionBath(bathRegion));

        std::map<unsigned, double>::const_iterator map_entry = mBathConductivities.find(bathRegion);

        if (map_entry != mBathConductivities.end())
        {
            return map_entry->second;
        }
        else
        {
            /*bath conductivity mS/cm*/
            return DecideLocation( & mpUserParameters->Physiological().BathConductivity(),
                                   & mpDefaultParameters->Physiological().BathConductivity(),
                                   "BathConductivity")->get();
        }
    }
}
const std::set<unsigned>&  HeartConfig::rGetTissueIdentifiers()
{
    return mTissueIdentifiers;
}

const std::set<unsigned>&  HeartConfig::rGetBathIdentifiers()
{
    return mBathIdentifiers;
}

double HeartConfig::GetSurfaceAreaToVolumeRatio() const
{
    /*surface area to volume ratio: 1/cm*/
    return DecideLocation( & mpUserParameters->Physiological().SurfaceAreaToVolumeRatio(),
                           & mpDefaultParameters->Physiological().SurfaceAreaToVolumeRatio(),
                           "SurfaceAreaToVolumeRatio")->get();
}

double HeartConfig::GetCapacitance() const
{
    //         capacitance                 : uF/cm^2
    return DecideLocation( & mpUserParameters->Physiological().Capacitance(),
                           & mpDefaultParameters->Physiological().Capacitance(),
                           "Capacitance")->get();
}

double HeartConfig::GetOdeTimeStep() const
{
    return DecideLocation( & mpUserParameters->Numerical().TimeSteps(),
                           & mpDefaultParameters->Numerical().TimeSteps(),
                           "ode TimeStep")->get().ode();
}

double HeartConfig::GetPdeTimeStep() const
{
    return DecideLocation( & mpUserParameters->Numerical().TimeSteps(),
                           & mpDefaultParameters->Numerical().TimeSteps(),
                           "pde TimeStep")->get().pde();
}

double HeartConfig::GetPrintingTimeStep() const
{
    return DecideLocation( & mpUserParameters->Numerical().TimeSteps(),
                           & mpDefaultParameters->Numerical().TimeSteps(),
                           "printing TimeStep")->get().printing();
}

bool HeartConfig::GetUseAbsoluteTolerance() const
{
    return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                            & mpDefaultParameters->Numerical().KSPTolerances(),
                            "KSPTolerances")->get().KSPAbsolute().present();
}

double HeartConfig::GetAbsoluteTolerance() const
{
    if (!GetUseAbsoluteTolerance())
    {
        EXCEPTION("Absolute tolerance is not set in Chaste parameters");
    }
    return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                           & mpDefaultParameters->Numerical().KSPTolerances(),
                           "KSPTolerances")->get().KSPAbsolute().get();
}

bool HeartConfig::GetUseRelativeTolerance() const
{
     return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                            & mpDefaultParameters->Numerical().KSPTolerances(),
                            "KSPTolerances")->get().KSPRelative().present();
}

double HeartConfig::GetRelativeTolerance() const
{
    if (!GetUseRelativeTolerance())
    {
        EXCEPTION("Relative tolerance is not set in Chaste parameters");
    }
    return DecideLocation( & mpUserParameters->Numerical().KSPTolerances(),
                           & mpDefaultParameters->Numerical().KSPTolerances(),
                           "KSPTolerances")->get().KSPRelative().get();
}

const char* HeartConfig::GetKSPSolver() const
{
    switch ( DecideLocation( & mpUserParameters->Numerical().KSPSolver(),
                             & mpDefaultParameters->Numerical().KSPSolver(),
                            "KSPSolver")->get() )
    {
        case cp::ksp_solver_type::gmres :
            return "gmres";
        case cp::ksp_solver_type::cg :
            return "cg";
        case cp::ksp_solver_type::symmlq :
            return "symmlq";
        case cp::ksp_solver_type::chebychev :
            return "chebychev";
    }
#define COVERAGE_IGNORE
    EXCEPTION("Unknown ksp solver");
#undef COVERAGE_IGNORE
}

const char* HeartConfig::GetKSPPreconditioner() const
{
    switch ( DecideLocation( & mpUserParameters->Numerical().KSPPreconditioner(),
                             & mpDefaultParameters->Numerical().KSPPreconditioner(),
                             "KSPPreconditioner")->get() )
    {
        case cp::ksp_preconditioner_type::jacobi :
            return "jacobi";
        case cp::ksp_preconditioner_type::bjacobi :
            return "bjacobi";
        case cp::ksp_preconditioner_type::hypre :
            return "hypre";
        case cp::ksp_preconditioner_type::ml :
            return "ml";
        case cp::ksp_preconditioner_type::spai :
            return "spai";
        case cp::ksp_preconditioner_type::blockdiagonal :
            return "blockdiagonal";
        case cp::ksp_preconditioner_type::ldufactorisation :
            return "ldufactorisation";
        case cp::ksp_preconditioner_type::twolevelsblockdiagonal :
            return "twolevelsblockdiagonal";
        case cp::ksp_preconditioner_type::none :
            return "none";

    }
#define COVERAGE_IGNORE
    EXCEPTION("Unknown ksp preconditioner");
#undef COVERAGE_IGNORE
}

DistributedTetrahedralMeshPartitionType::type HeartConfig::GetMeshPartitioning() const
{
    switch ( DecideLocation( & mpUserParameters->Numerical().MeshPartitioning(),
                             & mpDefaultParameters->Numerical().MeshPartitioning(),
                             "MeshPartitioning")->get() )
    {
        case cp::mesh_partitioning_type::dumb :
            return DistributedTetrahedralMeshPartitionType::DUMB;
        case cp::mesh_partitioning_type::metis :
            return DistributedTetrahedralMeshPartitionType::METIS_LIBRARY;
        case cp::mesh_partitioning_type::parmetis :
            return DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY;
        case cp::mesh_partitioning_type::petsc :
            return DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION;
    }
#define COVERAGE_IGNORE
    EXCEPTION("Unknown mesh partitioning type");
#undef COVERAGE_IGNORE
}

bool HeartConfig::IsAdaptivityParametersPresent() const
{
    try
    {
        DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                        & mpDefaultParameters->Numerical().AdaptivityParameters(),
                        "AdaptivityParameters")->present();
        //If there's a section
        return true;
    }
    catch (Exception &e)
    {
        //No section
        return false;
    }
}

double HeartConfig::GetTargetErrorForAdaptivity() const
{
    if ( IsAdaptivityParametersPresent() )
    {
        return DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                               & mpDefaultParameters->Numerical().AdaptivityParameters(),
                               "TargetError")->get().target_error();
    }
    else
    {
        EXCEPTION("Adaptivity parameters have not been set");
    }
}

double HeartConfig::GetSigmaForAdaptivity() const
{
    if ( IsAdaptivityParametersPresent() )
    {
        return DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                               & mpDefaultParameters->Numerical().AdaptivityParameters(),
                               "TargetError")->get().sigma();
    }
    else
    {
        EXCEPTION("Adaptivity parameters have not been set");
    }
}

double HeartConfig::GetMaxEdgeLengthForAdaptivity() const
{
    if ( IsAdaptivityParametersPresent() )
    {
    return DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                           & mpDefaultParameters->Numerical().AdaptivityParameters(),
                           "TargetError")->get().max_edge_length();
    }
    else
    {
        EXCEPTION("Adaptivity parameters have not been set");
    }
}

double HeartConfig::GetMinEdgeLengthForAdaptivity() const
{
    if ( IsAdaptivityParametersPresent() )
    {
        return DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                               & mpDefaultParameters->Numerical().AdaptivityParameters(),
                               "TargetError")->get().min_edge_length();
    }
    else
    {
        EXCEPTION("Adaptivity parameters have not been set");
    }
}

double HeartConfig::GetGradationForAdaptivity() const
{
    if ( IsAdaptivityParametersPresent() )
    {
        return DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                               & mpDefaultParameters->Numerical().AdaptivityParameters(),
                               "TargetError")->get().gradation();
    }
    else
    {
        EXCEPTION("Adaptivity parameters have not been set");
    }
}

unsigned HeartConfig::GetMaxNodesForAdaptivity() const
{
    if ( IsAdaptivityParametersPresent() )
    {
        return DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                               & mpDefaultParameters->Numerical().AdaptivityParameters(),
                               "TargetError")->get().max_nodes();
    }
    else
    {
        EXCEPTION("Adaptivity parameters have not been set");
    }
}

unsigned HeartConfig::GetNumberOfAdaptiveSweeps() const
{
    if ( IsAdaptivityParametersPresent() )
    {
        return DecideLocation( & mpUserParameters->Numerical().AdaptivityParameters(),
                               & mpDefaultParameters->Numerical().AdaptivityParameters(),
                               "TargetError")->get().num_sweeps();
    }
    else
    {
        EXCEPTION("Adaptivity parameters have not been set");
    }
}

/*
 * PostProcessing
 */

bool HeartConfig::IsPostProcessingSectionPresent() const
{
    try
    {
        DecideLocation( & mpUserParameters->PostProcessing(),
                        & mpDefaultParameters->PostProcessing(),
                        "PostProcessing")->present();
        //If there's a section
        return true;
    }
    catch (Exception &e)
    {
        //No section
        return false;
    }
}

void HeartConfig::EnsurePostProcessingSectionPresent()
{
    ENSURE_SECTION_PRESENT(mpUserParameters->PostProcessing(), cp::postprocessing_type);
}

bool HeartConfig::IsPostProcessingRequested() const
{
    if (IsPostProcessingSectionPresent() == false)
    {
        return false;
    }
    else
    {
        return(IsApdMapsRequested() ||
               IsUpstrokeTimeMapsRequested() ||
               IsMaxUpstrokeVelocityMapRequested() ||
               IsConductionVelocityMapsRequested() ||
               IsAnyNodalTimeTraceRequested()||
               IsPseudoEcgCalculationRequested());
    }
}
bool HeartConfig::IsApdMapsRequested() const
{
    assert(IsPostProcessingSectionPresent());

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ActionPotentialDurationMap)&
        apd_maps = DecideLocation( & mpUserParameters->PostProcessing(),
                                   & mpDefaultParameters->PostProcessing(),
                                   "ActionPotentialDurationMap")->get().ActionPotentialDurationMap();
    return (apd_maps.begin() != apd_maps.end());
}

void HeartConfig::GetApdMaps(std::vector<std::pair<double,double> >& apd_maps) const
{
    assert(IsApdMapsRequested());
    apd_maps.clear();

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ActionPotentialDurationMap)&
        apd_maps_sequence = DecideLocation( & mpUserParameters->PostProcessing(),
                                            & mpDefaultParameters->PostProcessing(),
                                            "ActionPotentialDurationMap")->get().ActionPotentialDurationMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::ActionPotentialDurationMap) i = apd_maps_sequence.begin();
         i != apd_maps_sequence.end();
         ++i)
    {
        std::pair<double,double> map(i->repolarisation_percentage(),i->threshold());

        apd_maps.push_back(map);
    }
}

bool HeartConfig::IsUpstrokeTimeMapsRequested() const
{
    assert(IsPostProcessingSectionPresent());

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::UpstrokeTimeMap)&
        upstroke_map = DecideLocation( & mpUserParameters->PostProcessing(),
                                       & mpDefaultParameters->PostProcessing(),
                                       "UpstrokeTimeMap")->get().UpstrokeTimeMap();
    return (upstroke_map.begin() != upstroke_map.end());
}
void HeartConfig::GetUpstrokeTimeMaps (std::vector<double>& upstroke_time_maps) const
{
    assert(IsUpstrokeTimeMapsRequested());
    assert(upstroke_time_maps.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::UpstrokeTimeMap)&
        upstroke_maps_sequence = DecideLocation( & mpUserParameters->PostProcessing(),
                                                 & mpDefaultParameters->PostProcessing(),
                                                 "UpstrokeTimeMap")->get().UpstrokeTimeMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::UpstrokeTimeMap) i = upstroke_maps_sequence.begin();
         i != upstroke_maps_sequence.end();
         ++i)
    {
        upstroke_time_maps.push_back(i->threshold());
    }
}

bool HeartConfig::IsMaxUpstrokeVelocityMapRequested() const
{
    assert(IsPostProcessingSectionPresent());

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap)&
        max_upstroke_velocity_map = DecideLocation( & mpUserParameters->PostProcessing(),
                                                    & mpDefaultParameters->PostProcessing(),
                                                    "MaxUpstrokeVelocityMap")->get().MaxUpstrokeVelocityMap();

    return (max_upstroke_velocity_map.begin() != max_upstroke_velocity_map.end());
}

void HeartConfig::GetMaxUpstrokeVelocityMaps(std::vector<double>& upstroke_velocity_maps) const
{
    assert(IsMaxUpstrokeVelocityMapRequested());
    assert(upstroke_velocity_maps.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap)&
        max_upstroke_velocity_maps_sequence = DecideLocation( & mpUserParameters->PostProcessing(),
                                                              & mpDefaultParameters->PostProcessing(),
                                                              "MaxUpstrokeVelocityMap")->get().MaxUpstrokeVelocityMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap) i = max_upstroke_velocity_maps_sequence.begin();
         i != max_upstroke_velocity_maps_sequence.end();
         ++i)
    {
        upstroke_velocity_maps.push_back(i->threshold());
    }
}

bool HeartConfig::IsConductionVelocityMapsRequested() const
{
    assert(IsPostProcessingSectionPresent());

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ConductionVelocityMap)&
        cond_vel_maps = DecideLocation( & mpUserParameters->PostProcessing(),
                                        & mpDefaultParameters->PostProcessing(),
                                        "ConductionVelocityMap")->get().ConductionVelocityMap();
    return (cond_vel_maps.begin() != cond_vel_maps.end());
}

void HeartConfig::GetConductionVelocityMaps(std::vector<unsigned>& conduction_velocity_maps) const
{
    assert(IsConductionVelocityMapsRequested());
    assert(conduction_velocity_maps.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ConductionVelocityMap)&
        cond_vel_maps_sequence = DecideLocation( & mpUserParameters->PostProcessing(),
                                                 & mpDefaultParameters->PostProcessing(),
                                                 "ConductionVelocityMap")->get().ConductionVelocityMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::ConductionVelocityMap) i = cond_vel_maps_sequence.begin();
         i != cond_vel_maps_sequence.end();
         ++i)
    {
        conduction_velocity_maps.push_back(i->origin_node());
    }
}

bool HeartConfig::IsAnyNodalTimeTraceRequested() const
{
    assert(IsPostProcessingSectionPresent());

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::TimeTraceAtNode)&
        requested_nodes = DecideLocation( & mpUserParameters->PostProcessing(),
                                        & mpDefaultParameters->PostProcessing(),
                                        "TimeTraceAtNode")->get().TimeTraceAtNode();
    return (requested_nodes.begin() != requested_nodes.end());
}

void HeartConfig::GetNodalTimeTraceRequested(std::vector<unsigned>& rRequestedNodes) const
{
    assert(IsAnyNodalTimeTraceRequested());
    assert(rRequestedNodes.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::TimeTraceAtNode)&
        req_nodes = DecideLocation( & mpUserParameters->PostProcessing(),
                                          & mpDefaultParameters->PostProcessing(),
                                          "TimeTraceAtNode")->get().TimeTraceAtNode();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::TimeTraceAtNode) i = req_nodes.begin();
         i != req_nodes.end();
         ++i)
    {
        rRequestedNodes.push_back(i->node_number());
    }
}


bool HeartConfig::IsPseudoEcgCalculationRequested() const
{
    assert(IsPostProcessingSectionPresent());

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::PseudoEcgElectrodePosition)&
        electrodes = DecideLocation( & mpUserParameters->PostProcessing(),
                                     & mpDefaultParameters->PostProcessing(),
                                     "PseudoEcgElectrodePosition")->get().PseudoEcgElectrodePosition();
    return (electrodes.begin() != electrodes.end());
}

template<unsigned SPACE_DIM>
void HeartConfig::GetPseudoEcgElectrodePositions(std::vector<ChastePoint<SPACE_DIM> >& rPseudoEcgElectrodePositions) const
{
    rPseudoEcgElectrodePositions.clear();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::PseudoEcgElectrodePosition)&
        electrodes = DecideLocation( & mpUserParameters->PostProcessing(),
                                     & mpDefaultParameters->PostProcessing(),
                                     "PseudoEcgElectrodePosition")->get().PseudoEcgElectrodePosition();
    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::PseudoEcgElectrodePosition) i = electrodes.begin();
         i != electrodes.end();
         ++i)
    {
        rPseudoEcgElectrodePositions.push_back(ChastePoint<SPACE_DIM>(i->x(), i->y(), i->z()));
    }
}


/*
 * Output visualization
 */

bool HeartConfig::IsOutputVisualizerPresent() const
{
    CheckSimulationIsDefined("OutputVisualizer");

    try
    {
        DecideLocation( & mpUserParameters->Simulation()->OutputVisualizer(),
                        & mpDefaultParameters->Simulation()->OutputVisualizer(),
                        "OutputVisualizer")->present();
        // If there's an element
        return true;
    }
    catch (Exception &e)
    {
        // No element
        return false;
    }
}

bool HeartConfig::GetVisualizeWithMeshalyzer() const
{
    if (!IsOutputVisualizerPresent())
    {
        return true;
    }
    else
    {
        return DecideLocation( & mpUserParameters->Simulation()->OutputVisualizer(),
                               & mpDefaultParameters->Simulation()->OutputVisualizer(),
                               "OutputVisualizer")->get().meshalyzer() == cp::yesno_type::yes;
    }
}

bool HeartConfig::GetVisualizeWithCmgui() const
{
    if (!IsOutputVisualizerPresent())
    {
        return false;
    }
    else
    {
        return DecideLocation( & mpUserParameters->Simulation()->OutputVisualizer(),
                               & mpDefaultParameters->Simulation()->OutputVisualizer(),
                               "OutputVisualizer")->get().cmgui() == cp::yesno_type::yes;
    }
}

bool HeartConfig::GetVisualizeWithParallelVtk() const
{
    if (!IsOutputVisualizerPresent())
    {
        return false;
    }
    else
    {
        return DecideLocation( & mpUserParameters->Simulation()->OutputVisualizer(),
                               & mpDefaultParameters->Simulation()->OutputVisualizer(),
                               "OutputVisualizer")->get().parallel_vtk() == cp::yesno_type::yes;
    }
}

bool HeartConfig::GetVisualizeWithVtk() const
{
    if (!IsOutputVisualizerPresent())
    {
        return false;
    }
    else
    {
        return DecideLocation( & mpUserParameters->Simulation()->OutputVisualizer(),
                               & mpDefaultParameters->Simulation()->OutputVisualizer(),
                               "OutputVisualizer")->get().vtk() == cp::yesno_type::yes;
    }
}

unsigned HeartConfig::GetVisualizerOutputPrecision()
{
    if (!IsOutputVisualizerPresent())
    {
        return 0u;
    }
    else
    {
        return DecideLocation( & mpUserParameters->Simulation()->OutputVisualizer(),
                               & mpDefaultParameters->Simulation()->OutputVisualizer(),
                               "OutputVisualizer")->get().precision();
    }
}


bool HeartConfig::IsElectrodesPresent() const
{
    try
    {
        DecideLocation( & mpUserParameters->Simulation()->Electrodes(),
                        & mpDefaultParameters->Simulation()->Electrodes(),
                        "Electrodes")->present();
        //If there's a section
        return true;
    }
    catch (Exception &e)
    {
        //No section
        return false;
    }
}

/*
 *  Set methods
 */
void HeartConfig::SetSpaceDimension(unsigned spaceDimension)
{
    mpUserParameters->Simulation()->SpaceDimension().set(spaceDimension);
}

void HeartConfig::SetSimulationDuration(double simulationDuration)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::time_type, time, simulationDuration, "ms");
    mpUserParameters->Simulation()->SimulationDuration().set(time);
}

void HeartConfig::SetDomain(const cp::domain_type& rDomain)
{
    mpUserParameters->Simulation()->Domain().set(rDomain);
}

void HeartConfig::SetDefaultIonicModel(const cp::ionic_models_available_type& rIonicModel)
{
    cp::ionic_model_selection_type ionic_model;
    ionic_model.Hardcoded(rIonicModel);
    cp::ionic_models_type container(ionic_model);
    mpUserParameters->Simulation()->IonicModels().set(container);
}

void HeartConfig::SetSlabDimensions(double x, double y, double z, double inter_node_space)
{
    if ( ! mpUserParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpUserParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    cp::slab_type slab_definition(x, y, z, inter_node_space);
    mpUserParameters->Simulation()->Mesh()->Slab().set(slab_definition);
}

void HeartConfig::SetSheetDimensions(double x, double y, double inter_node_space)
{
    if ( ! mpUserParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpUserParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    cp::sheet_type sheet_definition(x, y, inter_node_space);
    mpUserParameters->Simulation()->Mesh()->Sheet().set(sheet_definition);
}

void HeartConfig::SetFibreLength(double x, double inter_node_space)
{
    if ( ! mpUserParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpUserParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    cp::fibre_type fibre_definition(x, inter_node_space);
    mpUserParameters->Simulation()->Mesh()->Fibre().set(fibre_definition);
}

void HeartConfig::SetMeshFileName(std::string meshPrefix, cp::media_type fibreDefinition)
{
    if ( ! mpUserParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpUserParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    XSD_NESTED_TYPE(cp::mesh_type::LoadMesh) mesh_prefix(meshPrefix, fibreDefinition);
    mpUserParameters->Simulation()->Mesh()->LoadMesh().set(mesh_prefix);
}

void HeartConfig::SetIonicModelRegions(std::vector<ChasteCuboid<3> >& rDefinedRegions,
                                       std::vector<cp::ionic_model_selection_type>& rIonicModels) const
{
    assert(rDefinedRegions.size() == rIonicModels.size());
    // You need to have defined a default model first...
    assert(mpUserParameters->Simulation()->IonicModels().present());
    XSD_SEQUENCE_TYPE(cp::ionic_models_type::Region)&
        regions = mpUserParameters->Simulation()->IonicModels()->Region();
    regions.clear();
    for (unsigned region_index=0; region_index<rDefinedRegions.size(); region_index++)
    {
        cp::point_type point_a(rDefinedRegions[region_index].rGetLowerCorner()[0],
                               rDefinedRegions[region_index].rGetLowerCorner()[1],
                               rDefinedRegions[region_index].rGetLowerCorner()[2]);

        cp::point_type point_b(rDefinedRegions[region_index].rGetUpperCorner()[0],
                               rDefinedRegions[region_index].rGetUpperCorner()[1],
                               rDefinedRegions[region_index].rGetUpperCorner()[2]);

        XSD_CREATE_WITH_FIXED_ATTR(cp::location_type, locn, "cm");
        locn.Cuboid().set(cp::box_type(point_a, point_b));

        cp::ionic_model_region_type region(rIonicModels[region_index], locn);
        regions.push_back(region);
    }
}

void HeartConfig::SetConductivityHeterogeneities(std::vector<ChasteCuboid<3> >& conductivityAreas,
        std::vector< c_vector<double,3> >& intraConductivities,
        std::vector< c_vector<double,3> >& extraConductivities)
{
    assert ( conductivityAreas.size() == intraConductivities.size() );
    assert ( intraConductivities.size() == extraConductivities.size());

    XSD_ANON_SEQUENCE_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity) heterogeneities_container;

    for (unsigned region_index=0; region_index<conductivityAreas.size(); region_index++)
    {
        cp::point_type point_a(conductivityAreas[region_index].rGetLowerCorner()[0],
                           conductivityAreas[region_index].rGetLowerCorner()[1],
                           conductivityAreas[region_index].rGetLowerCorner()[2]);

        cp::point_type point_b(conductivityAreas[region_index].rGetUpperCorner()[0],
                           conductivityAreas[region_index].rGetUpperCorner()[1],
                           conductivityAreas[region_index].rGetUpperCorner()[2]);

        XSD_CREATE_WITH_FIXED_ATTR(cp::location_type, locn, "cm");
        locn.Cuboid().set(cp::box_type(point_a, point_b));
        cp::conductivity_heterogeneity_type ht(locn);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                    intraConductivities[region_index][0],
                                    intraConductivities[region_index][1],
                                    intraConductivities[region_index][2],
                                    "mS/cm");

        ht.IntracellularConductivities(intra);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                    extraConductivities[region_index][0],
                                    extraConductivities[region_index][1],
                                    extraConductivities[region_index][2],
                                    "mS/cm");

        ht.ExtracellularConductivities(extra);

        heterogeneities_container.push_back(ht);
    }

    XSD_ANON_TYPE(cp::physiological_type, ConductivityHeterogeneities) heterogeneities_object;
    heterogeneities_object.ConductivityHeterogeneity(heterogeneities_container);

    mpUserParameters->Physiological().ConductivityHeterogeneities().set(heterogeneities_object);
}

void HeartConfig::SetConductivityHeterogeneitiesEllipsoid(std::vector<ChasteEllipsoid<3> >& conductivityAreas,
        std::vector< c_vector<double,3> >& intraConductivities,
        std::vector< c_vector<double,3> >& extraConductivities)
{
    assert ( conductivityAreas.size() == intraConductivities.size() );
    assert ( intraConductivities.size() == extraConductivities.size());

    XSD_ANON_SEQUENCE_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity) heterogeneities_container;

    for (unsigned region_index=0; region_index<conductivityAreas.size(); region_index++)
    {
        cp::point_type centre(conductivityAreas[region_index].rGetCentre()[0],
                              conductivityAreas[region_index].rGetCentre()[1],
                              conductivityAreas[region_index].rGetCentre()[2]);

        cp::point_type radii(conductivityAreas[region_index].rGetRadii()[0],
                             conductivityAreas[region_index].rGetRadii()[1],
                             conductivityAreas[region_index].rGetRadii()[2]);

        XSD_CREATE_WITH_FIXED_ATTR(cp::location_type, locn, "cm");
        locn.Ellipsoid().set(cp::ellipsoid_type(centre, radii));
        cp::conductivity_heterogeneity_type ht(locn);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                    intraConductivities[region_index][0],
                                    intraConductivities[region_index][1],
                                    intraConductivities[region_index][2],
                                    "mS/cm");

        ht.IntracellularConductivities(intra);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                    extraConductivities[region_index][0],
                                    extraConductivities[region_index][1],
                                    extraConductivities[region_index][2],
                                    "mS/cm");

        ht.ExtracellularConductivities(extra);

        heterogeneities_container.push_back(ht);
    }

    XSD_ANON_TYPE(cp::physiological_type, ConductivityHeterogeneities) heterogeneities_object;
    heterogeneities_object.ConductivityHeterogeneity(heterogeneities_container);

    mpUserParameters->Physiological().ConductivityHeterogeneities().set(heterogeneities_object);
}

void HeartConfig::SetOutputDirectory(const std::string& rOutputDirectory)
{
    mpUserParameters->Simulation()->OutputDirectory().set(rOutputDirectory);
}

void HeartConfig::SetOutputFilenamePrefix(const std::string& rOutputFilenamePrefix)
{
    mpUserParameters->Simulation()->OutputFilenamePrefix().set(rOutputFilenamePrefix);
}

void HeartConfig::SetOutputVariables(const std::vector<std::string>& rOutputVariables)
{
    if ( ! mpUserParameters->Simulation()->OutputVariables().present())
    {
        cp::output_variables_type variables_requested;
        mpUserParameters->Simulation()->OutputVariables().set(variables_requested);
    }

    XSD_SEQUENCE_TYPE(cp::output_variables_type::Var)&
        var_type_sequence = mpUserParameters->Simulation()->OutputVariables()->Var();
    //Erase or create a sequence
    var_type_sequence.clear();

    for (unsigned i=0; i<rOutputVariables.size(); i++)
    {
        cp::var_type temp(rOutputVariables[i]);
        var_type_sequence.push_back(temp);
    }
}

void  HeartConfig::SetOutputUsingOriginalNodeOrdering(bool useOriginal)
{
    //What if it doesn't exist?
    mpUserParameters->Simulation()->OutputUsingOriginalNodeOrdering().set(useOriginal? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetCheckpointSimulation(bool saveSimulation, double checkpointTimestep, unsigned maxCheckpointsOnDisk)
{
    if (saveSimulation)
    {
        // Make sure values for the optional parameters have been provided
        assert(checkpointTimestep!=-1.0 && maxCheckpointsOnDisk!=UINT_MAX);

        XSD_CREATE_WITH_FIXED_ATTR2(cp::simulation_type::XSD_NESTED_TYPE(CheckpointSimulation),
                                    cs,
                                    checkpointTimestep,
                                    maxCheckpointsOnDisk,
                                    "ms");
        mpUserParameters->Simulation()->CheckpointSimulation().set(cs);
    }
    else
    {
        mpUserParameters->Simulation()->CheckpointSimulation().reset();
    }

    CheckTimeSteps();
}

// Physiological

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 3>& intraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                intraConductivities[0],
                                intraConductivities[1],
                                intraConductivities[2],
                                "mS/cm");

    mpUserParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 2>& intraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                intraConductivities[0],
                                intraConductivities[1],
                                0.0, "mS/cm");

    mpUserParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 1>& intraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                intraConductivities[0],
                                0.0, 0.0, "mS/cm");

    mpUserParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 3>& extraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                extraConductivities[0],
                                extraConductivities[1],
                                extraConductivities[2],
                                "mS/cm");

    mpUserParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 2>& extraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                extraConductivities[0],
                                extraConductivities[1],
                                0.0, "mS/cm");

    mpUserParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 1>& extraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                extraConductivities[0],
                                0.0, 0.0, "mS/cm");

    mpUserParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetBathConductivity(double bathConductivity)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::conductivity_type, cond, bathConductivity, "mS/cm");
    mpUserParameters->Physiological().BathConductivity().set(cond);
}

void HeartConfig::SetBathMultipleConductivities(std::map<unsigned, double> bathConductivities)
{
    /// \todo: This implementation is temporary until we incorporate the bath heterogeneities to the XML schema
    mBathConductivities = bathConductivities;
}

//void HeartConfig::SetTissueIdentifiers(const std::set<unsigned>& tissueIds)
//{
//    std::set<unsigned> empty_bath_identifiers;  //Too dangerous (see GetValidBathId)
//    SetTissueAndBathIdentifiers(tissueIds, mBathIdentifiers);
//}

void HeartConfig::SetTissueAndBathIdentifiers(const std::set<unsigned>& tissueIds, const std::set<unsigned>& bathIds)
{
    if (tissueIds.empty() || bathIds.empty() )
    {
        EXCEPTION("Identifying set must be non-empty");
    }
    std::set<unsigned> shared_identifiers;
    std::set_intersection(tissueIds.begin(),
                          tissueIds.end(),
                          bathIds.begin(),
                          bathIds.end(),
                          std::inserter(shared_identifiers, shared_identifiers.begin()));

    if (!shared_identifiers.empty())
    {
        EXCEPTION("Tissue identifiers and bath identifiers overlap");
    }
    mTissueIdentifiers=tissueIds;
    mBathIdentifiers=bathIds;
}

void HeartConfig::SetSurfaceAreaToVolumeRatio(double ratio)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::inverse_length_type, ratio_object, ratio, "1/cm");
    mpUserParameters->Physiological().SurfaceAreaToVolumeRatio().set(ratio_object);
}

void HeartConfig::SetCapacitance(double capacitance)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::capacitance_type, capacitance_object, capacitance, "uF/cm^2");
    mpUserParameters->Physiological().Capacitance().set(capacitance_object);
}


// Numerical
void HeartConfig::SetOdePdeAndPrintingTimeSteps(double odeTimeStep, double pdeTimeStep, double printingTimeStep)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::time_steps_type, time_steps,
                                odeTimeStep, pdeTimeStep, printingTimeStep, "ms");
    mpUserParameters->Numerical().TimeSteps().set(time_steps);
    CheckTimeSteps();
}

void HeartConfig::SetOdeTimeStep(double odeTimeStep)
{
    SetOdePdeAndPrintingTimeSteps(odeTimeStep, GetPdeTimeStep(), GetPrintingTimeStep());
}

void HeartConfig::SetPdeTimeStep(double pdeTimeStep)
{
    SetOdePdeAndPrintingTimeSteps(GetOdeTimeStep(), pdeTimeStep, GetPrintingTimeStep());
}

void HeartConfig::SetPrintingTimeStep(double printingTimeStep)
{
    SetOdePdeAndPrintingTimeSteps(GetOdeTimeStep(), GetPdeTimeStep(), printingTimeStep);
}

void HeartConfig::CheckTimeSteps() const
{
    if (GetOdeTimeStep() <= 0)
    {
        EXCEPTION("Ode time-step should be positive");
    }
    if (GetPdeTimeStep() <= 0)
    {
        EXCEPTION("Pde time-step should be positive");
    }
    if (GetPrintingTimeStep() <= 0.0)
    {
        EXCEPTION("Printing time-step should be positive");
    }

    if (GetPdeTimeStep()>GetPrintingTimeStep())
    {
        EXCEPTION("Printing time-step should not be smaller than PDE time-step");
    }

    if ( !Divides(GetPdeTimeStep(), GetPrintingTimeStep()) )
    {
        EXCEPTION("Printing time-step should be a multiple of PDE time step");
    }

    if ( GetOdeTimeStep() > GetPdeTimeStep() )
    {
        EXCEPTION("Ode time-step should not be greater than PDE time-step");
    }

    if (GetCheckpointSimulation())
    {
        if (GetCheckpointTimestep() <= 0.0)
        {
            EXCEPTION("Checkpoint time-step should be positive");
        }

        if ( !Divides(GetPrintingTimeStep(), GetCheckpointTimestep()) )
        {
            EXCEPTION("Checkpoint time-step should be a multiple of printing time-step");
        }
    }
}


void HeartConfig::SetUseRelativeTolerance(double relativeTolerance)
{
    ENSURE_SECTION_PRESENT(mpUserParameters->Numerical().KSPTolerances(), cp::ksp_tolerances_type);
    //Remove any reference to tolerances is user parameters
    mpUserParameters->Numerical().KSPTolerances()->KSPAbsolute().reset();
    mpUserParameters->Numerical().KSPTolerances()->KSPRelative().set(relativeTolerance);
}

void HeartConfig::SetUseAbsoluteTolerance(double absoluteTolerance)
{
    ENSURE_SECTION_PRESENT(mpUserParameters->Numerical().KSPTolerances(), cp::ksp_tolerances_type);
    //Remove any reference to tolerances is user parameters
    mpUserParameters->Numerical().KSPTolerances()->KSPRelative().reset();
    mpUserParameters->Numerical().KSPTolerances()->KSPAbsolute().set(absoluteTolerance);
}

void HeartConfig::SetKSPSolver(const char* kspSolver)
{
    /* Note that changes in these conditions need to be reflected in the Doxygen*/
    if ( strcmp(kspSolver, "gmres") == 0)
    {
        mpUserParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::gmres);
        return;
    }
    if ( strcmp(kspSolver, "cg") == 0)
    {
        mpUserParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::cg);
        return;
    }
    if ( strcmp(kspSolver, "symmlq") == 0)
    {
        mpUserParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::symmlq);
        return;
    }
    if ( strcmp(kspSolver, "chebychev") == 0)
    {
        mpUserParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::chebychev);
        return;
    }

    EXCEPTION("Unknown solver type provided");
}

void HeartConfig::SetKSPPreconditioner(const char* kspPreconditioner)
{
    /* Note that changes in these conditions need to be reflected in the Doxygen*/
    if ( strcmp(kspPreconditioner, "jacobi") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::jacobi);
        return;
    }
    if ( strcmp(kspPreconditioner, "bjacobi") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::bjacobi);
        return;
    }
    if ( strcmp(kspPreconditioner, "hypre") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::hypre);
        return;
    }
    if ( strcmp(kspPreconditioner, "ml") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::ml);
        return;
    }
    if ( strcmp(kspPreconditioner, "spai") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::spai);
        return;
    }
    if ( strcmp(kspPreconditioner, "twolevelsblockdiagonal") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::twolevelsblockdiagonal);
        return;
    }
    if ( strcmp(kspPreconditioner, "blockdiagonal") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::blockdiagonal);
        return;
    }
    if ( strcmp(kspPreconditioner, "ldufactorisation") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::ldufactorisation);
        return;
    }
    if ( strcmp(kspPreconditioner, "none") == 0)
    {
        mpUserParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::none);
        return;
    }

    EXCEPTION("Unknown preconditioner type provided");
}

void HeartConfig::SetMeshPartitioning(const char* meshPartioningMethod)
{
    /* Note that changes in these conditions need to be reflected in the Doxygen*/
    if ( strcmp(meshPartioningMethod, "dumb") == 0)
    {
        mpUserParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::dumb);
        return;
    }
    if ( strcmp(meshPartioningMethod, "metis") == 0)
    {
        mpUserParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::metis);
        return;
    }
    if ( strcmp(meshPartioningMethod, "parmetis") == 0)
    {
        mpUserParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::parmetis);
        return;
    }
    if ( strcmp(meshPartioningMethod, "petsc") == 0)
    {
        mpUserParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::petsc);
        return;
    }

    EXCEPTION("Unknown mesh partitioning method provided");
}

void HeartConfig::SetAdaptivityParameters(double targetError,
                                          double sigma,
                                          double maxEdgeLength,
                                          double minEdgeLength,
                                          double gradation,
                                          unsigned maxNodes,
                                          unsigned numSweeps )
{
    if ( maxEdgeLength < minEdgeLength )
    {
        EXCEPTION("AdaptivityParameters: maxEdgeLength must be greater than minEdgeLength.");
    }

    cp::adaptivity_parameters_type element(targetError,
                                           sigma,
                                           maxEdgeLength,
                                           minEdgeLength,
                                           gradation,
                                           maxNodes,
                                           numSweeps );
    mpUserParameters->Numerical().AdaptivityParameters().set(element);

    if (IsAdaptivityParametersPresent())
    {
        mpUserParameters->Numerical().AdaptivityParameters()->target_error(targetError);
        mpUserParameters->Numerical().AdaptivityParameters()->sigma(sigma);
        mpUserParameters->Numerical().AdaptivityParameters()->max_edge_length(maxEdgeLength);
        mpUserParameters->Numerical().AdaptivityParameters()->min_edge_length(minEdgeLength);
        mpUserParameters->Numerical().AdaptivityParameters()->gradation(gradation);
        mpUserParameters->Numerical().AdaptivityParameters()->max_nodes(maxNodes);
        mpUserParameters->Numerical().AdaptivityParameters()->num_sweeps(numSweeps);
    }
}

void HeartConfig::SetTargetErrorForAdaptivity(double targetError)
{
    SetAdaptivityParameters( targetError,
                             GetSigmaForAdaptivity(),
                             GetMaxEdgeLengthForAdaptivity(),
                             GetMinEdgeLengthForAdaptivity(),
                             GetGradationForAdaptivity(),
                             GetMaxNodesForAdaptivity(),
                             GetNumberOfAdaptiveSweeps() );
}

void HeartConfig::SetSigmaForAdaptivity(double sigma)
{
    SetAdaptivityParameters( GetTargetErrorForAdaptivity(),
                             sigma,
                             GetMaxEdgeLengthForAdaptivity(),
                             GetMinEdgeLengthForAdaptivity(),
                             GetGradationForAdaptivity(),
                             GetMaxNodesForAdaptivity(),
                             GetNumberOfAdaptiveSweeps() );
}

void HeartConfig::SetMaxEdgeLengthForAdaptivity(double maxEdgeLength)
{
    SetAdaptivityParameters( GetTargetErrorForAdaptivity(),
                             GetSigmaForAdaptivity(),
                             maxEdgeLength,
                             GetMinEdgeLengthForAdaptivity(),
                             GetGradationForAdaptivity(),
                             GetMaxNodesForAdaptivity(),
                             GetNumberOfAdaptiveSweeps() );
}

void HeartConfig::SetMinEdgeLengthForAdaptivity(double minEdgeLength)
{
    SetAdaptivityParameters( GetTargetErrorForAdaptivity(),
                             GetSigmaForAdaptivity(),
                             GetMaxEdgeLengthForAdaptivity(),
                             minEdgeLength,
                             GetGradationForAdaptivity(),
                             GetMaxNodesForAdaptivity(),
                             GetNumberOfAdaptiveSweeps() );
}

void HeartConfig::SetGradationForAdaptivity(double gradation)
{
    SetAdaptivityParameters( GetTargetErrorForAdaptivity(),
                             GetSigmaForAdaptivity(),
                             GetMaxEdgeLengthForAdaptivity(),
                             GetMinEdgeLengthForAdaptivity(),
                             gradation,
                             GetMaxNodesForAdaptivity(),
                             GetNumberOfAdaptiveSweeps() );
}

void HeartConfig::SetMaxNodesForAdaptivity(unsigned maxNodes)
{
    SetAdaptivityParameters( GetTargetErrorForAdaptivity(),
                             GetSigmaForAdaptivity(),
                             GetMaxEdgeLengthForAdaptivity(),
                             GetMinEdgeLengthForAdaptivity(),
                             GetGradationForAdaptivity(),
                             maxNodes,
                             GetNumberOfAdaptiveSweeps() );
}

void HeartConfig::SetNumberOfAdaptiveSweeps(unsigned numSweeps)
{
    SetAdaptivityParameters( GetTargetErrorForAdaptivity(),
                             GetSigmaForAdaptivity(),
                             GetMaxEdgeLengthForAdaptivity(),
                             GetMinEdgeLengthForAdaptivity(),
                             GetGradationForAdaptivity(),
                             GetMaxNodesForAdaptivity(),
                             numSweeps );
}

void HeartConfig::SetApdMaps(const std::vector<std::pair<double,double> >& apdMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ActionPotentialDurationMap)& apd_maps_sequence
        = mpUserParameters->PostProcessing()->ActionPotentialDurationMap();
    //Erase or create a sequence
    apd_maps_sequence.clear();

    for (unsigned i=0; i<apdMaps.size(); i++)
    {
        XSD_CREATE_WITH_FIXED_ATTR2(cp::apd_map_type, temp,
                                    apdMaps[i].first, apdMaps[i].second,
                                    "mV");
        apd_maps_sequence.push_back( temp);
    }
}


void HeartConfig::SetUpstrokeTimeMaps (std::vector<double>& upstrokeTimeMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::UpstrokeTimeMap)& var_type_sequence
        = mpUserParameters->PostProcessing()->UpstrokeTimeMap();

    //Erase or create a sequence
    var_type_sequence.clear();

    for (unsigned i=0; i<upstrokeTimeMaps.size(); i++)
    {
        XSD_CREATE_WITH_FIXED_ATTR1(cp::upstrokes_map_type, temp,
                                    upstrokeTimeMaps[i],
                                    "mV");
        var_type_sequence.push_back(temp);
    }
}

void HeartConfig::SetMaxUpstrokeVelocityMaps (std::vector<double>& maxUpstrokeVelocityMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap)& max_upstroke_velocity_maps_sequence
        = mpUserParameters->PostProcessing()->MaxUpstrokeVelocityMap();

    //Erase or create a sequence
    max_upstroke_velocity_maps_sequence.clear();

    for (unsigned i=0; i<maxUpstrokeVelocityMaps.size(); i++)
    {
        XSD_CREATE_WITH_FIXED_ATTR1(cp::max_upstrokes_velocity_map_type, temp,
                                    maxUpstrokeVelocityMaps[i],
                                    "mV");


        max_upstroke_velocity_maps_sequence.push_back(temp);
    }

}

void HeartConfig::SetConductionVelocityMaps (std::vector<unsigned>& conductionVelocityMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ConductionVelocityMap)& conduction_velocity_maps_sequence
        = mpUserParameters->PostProcessing()->ConductionVelocityMap();

    //Erase or create a sequence
    conduction_velocity_maps_sequence.clear();

    for (unsigned i=0; i<conductionVelocityMaps.size(); i++)
    {
        cp::conduction_velocity_map_type temp(conductionVelocityMaps[i]);
        conduction_velocity_maps_sequence.push_back(temp);
    }
}

void HeartConfig::SetRequestedNodalTimeTraces (std::vector<unsigned>& requestedNodes)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::TimeTraceAtNode)& requested_nodes_sequence
        = mpUserParameters->PostProcessing()->TimeTraceAtNode();

    //Erase or create a sequence
    requested_nodes_sequence.clear();

    for (unsigned i=0; i<requestedNodes.size(); i++)
    {
        cp::node_number_type temp(requestedNodes[i]);
        requested_nodes_sequence.push_back(temp);
    }
}

template<unsigned SPACE_DIM>
void HeartConfig::SetPseudoEcgElectrodePositions(const std::vector<ChastePoint<SPACE_DIM> >& rPseudoEcgElectrodePositions)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::PseudoEcgElectrodePosition)& electrodes_sequence
        = mpUserParameters->PostProcessing()->PseudoEcgElectrodePosition();

    //Erase or create a sequence
    electrodes_sequence.clear();

    for (unsigned i=0; i<rPseudoEcgElectrodePositions.size(); i++)
    {
        cp::point_type temp(rPseudoEcgElectrodePositions[i].GetWithDefault(0),
                            rPseudoEcgElectrodePositions[i].GetWithDefault(1),
                            rPseudoEcgElectrodePositions[i].GetWithDefault(2));
        electrodes_sequence.push_back(temp);
    }
}


/*
 * Output visualizer
 */

void HeartConfig::EnsureOutputVisualizerExists()
{
    ENSURE_SECTION_PRESENT(mpUserParameters->Simulation()->OutputVisualizer(), cp::output_visualizer_type);
}

void HeartConfig::SetVisualizeWithMeshalyzer(bool useMeshalyzer)
{
    EnsureOutputVisualizerExists();

    mpUserParameters->Simulation()->OutputVisualizer()->meshalyzer(
        useMeshalyzer ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizeWithCmgui(bool useCmgui)
{
    EnsureOutputVisualizerExists();

    mpUserParameters->Simulation()->OutputVisualizer()->cmgui(
        useCmgui ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizeWithVtk(bool useVtk)
{
    EnsureOutputVisualizerExists();

    mpUserParameters->Simulation()->OutputVisualizer()->vtk(
        useVtk ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizeWithParallelVtk(bool useParallelVtk)
{
    EnsureOutputVisualizerExists();

    mpUserParameters->Simulation()->OutputVisualizer()->parallel_vtk(
        useParallelVtk ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizerOutputPrecision(unsigned numberOfDigits)
{
    EnsureOutputVisualizerExists();

    mpUserParameters->Simulation()->OutputVisualizer()->precision(numberOfDigits);
}


void HeartConfig::SetElectrodeParameters(bool groundSecondElectrode,
                                         unsigned index, double magnitude,
                                         double startTime, double duration )
{
    assert(index < 3);

    cp::axis_type axis = cp::axis_type::x;
    if (index==1)
    {
        axis = cp::axis_type::y;
    }
    else if (index==2)
    {
        axis = cp::axis_type::z;
    }

    XSD_CREATE_WITH_FIXED_ATTR1(cp::surface_stimulus_strength_type, strength, magnitude, "uA/cm^2");
    XSD_CREATE_WITH_FIXED_ATTR1(cp::time_type, start_time, startTime, "ms");
    XSD_CREATE_WITH_FIXED_ATTR1(cp::time_type, duration_time, duration, "ms");

    if (!IsElectrodesPresent())
    {
        cp::electrodes_type element( groundSecondElectrode ? cp::yesno_type::yes : cp::yesno_type::no,
                                     axis,
                                     strength,
                                     start_time,
                                     duration_time );
        mpUserParameters->Simulation()->Electrodes().set(element);
    }
    else
    {
        mpUserParameters->Simulation()->Electrodes()->GroundSecondElectrode(groundSecondElectrode ? cp::yesno_type::yes : cp::yesno_type::no);
        mpUserParameters->Simulation()->Electrodes()->PerpendicularToAxis(axis);
        mpUserParameters->Simulation()->Electrodes()->Strength(strength);
        mpUserParameters->Simulation()->Electrodes()->StartTime(start_time);
        mpUserParameters->Simulation()->Electrodes()->Duration(duration_time);
    }
}

void HeartConfig::GetElectrodeParameters(bool& rGroundSecondElectrode,
                                         unsigned& rIndex, double& rMagnitude,
                                         double& rStartTime, double& rDuration )
{
    if (!IsElectrodesPresent())
    {
        EXCEPTION("Attempted to get electrodes that have not been defined.");
    }
    else
    {
        rGroundSecondElectrode = (mpUserParameters->Simulation()->Electrodes()->GroundSecondElectrode()==cp::yesno_type::yes);

        cp::axis_type axis = mpUserParameters->Simulation()->Electrodes()->PerpendicularToAxis();
        if (axis==cp::axis_type::x)
        {
            rIndex = 0;
        }
        else if (axis==cp::axis_type::y)
        {
            rIndex = 1;
        }
        else
        {
            rIndex = 2;
        }

        rMagnitude = mpUserParameters->Simulation()->Electrodes()->Strength();
        rStartTime = mpUserParameters->Simulation()->Electrodes()->StartTime();
        rDuration = mpUserParameters->Simulation()->Electrodes()->Duration();
    }

}

bool HeartConfig::GetUseStateVariableInterpolation() const
{
    try
    {
        return DecideLocation( & mpUserParameters->Numerical().UseStateVariableInterpolation(),
                               & mpDefaultParameters->Numerical().UseStateVariableInterpolation(),
                               "UseStateVariableInterpolation")->get() == cp::yesno_type::yes;
    }
    catch (const Exception& e)
    {
        // It's an older version parameters & defaults (we're loading a checkpoint)
        return false;
    }
}

void HeartConfig::SetUseStateVariableInterpolation(bool useStateVariableInterpolation)
{
    if (useStateVariableInterpolation)
    {
        mpUserParameters->Numerical().UseStateVariableInterpolation().set(cp::yesno_type::yes);
    }
    else
    {
        mpUserParameters->Numerical().UseStateVariableInterpolation().set(cp::yesno_type::no);
    }
}



bool HeartConfig::HasDrugDose() const
{
    try
    {
        DecideLocation( & mpUserParameters->Physiological().ApplyDrug(),
                        & mpDefaultParameters->Physiological().ApplyDrug(),
                        "ApplyDrug")->present();
        // If there's an element
        return true;
    }
    catch (Exception &e)
    {
        // No element
        return false;
    }
}

double HeartConfig::GetDrugDose() const
{
    return DecideLocation( & mpUserParameters->Physiological().ApplyDrug(),
                           & mpDefaultParameters->Physiological().ApplyDrug(),
                           "ApplyDrug")->get().concentration();
}

void HeartConfig::SetDrugDose(double drugDose)
{
    if (!mpUserParameters->Physiological().ApplyDrug().present())
    {
        cp::apply_drug_type drug(drugDose);
        mpUserParameters->Physiological().ApplyDrug().set(drug);
    }
    else
    {
        mpUserParameters->Physiological().ApplyDrug()->concentration(drugDose);
    }
}

std::map<std::string, std::pair<double, double> > HeartConfig::GetIc50Values()
{
    std::map<std::string, std::pair<double, double> > ic50s;

    XSD_SEQUENCE_TYPE(cp::apply_drug_type::IC50)&
        ic50_seq = DecideLocation( & mpUserParameters->Physiological().ApplyDrug(),
                                   & mpDefaultParameters->Physiological().ApplyDrug(),
                                   "ApplyDrug")->get().IC50();

    for (XSD_ITERATOR_TYPE(cp::apply_drug_type::IC50) i = ic50_seq.begin();
         i != ic50_seq.end();
         ++i)
    {
        std::pair<double, double> ic50_hill(*i, i->hill());
        std::string current = i->current();
        ic50s[current] = ic50_hill;
    }

    return ic50s;
}

void HeartConfig::SetIc50Value(const std::string& rCurrentName, double ic50, double hill)
{
    if (!mpUserParameters->Physiological().ApplyDrug().present())
    {
        SetDrugDose(0.0);
    }
    XSD_SEQUENCE_TYPE(cp::apply_drug_type::IC50)& ic50_seq = mpUserParameters->Physiological().ApplyDrug()->IC50();
    if (ic50_seq.empty())
    {
        // Erase or create a sequence
        ic50_seq.clear();
    }
    bool entry_exists = false;
    cp::ic50_type ic50_elt(ic50, rCurrentName);
    ic50_elt.hill(hill);
    for (XSD_ITERATOR_TYPE(cp::apply_drug_type::IC50) i = ic50_seq.begin();
         i != ic50_seq.end();
         ++i)
    {
        if (i->current() == rCurrentName)
        {
            entry_exists = true;
            *i = ic50_elt;
            break;
        }
    }
    if (!entry_exists)
    {
        ic50_seq.push_back(ic50_elt);
    }
}



void HeartConfig::SetUseMassLumping(bool useMassLumping)
{
    mUseMassLumping = useMassLumping;
}

bool HeartConfig::GetUseMassLumping()
{
    return mUseMassLumping;
}

void HeartConfig::SetUseMassLumpingForPrecond(bool useMassLumping)
{
    mUseMassLumpingForPrecond = useMassLumping;
}

bool HeartConfig::GetUseMassLumpingForPrecond()
{
    return mUseMassLumpingForPrecond;
}

void HeartConfig::SetUseReactionDiffusionOperatorSplitting(bool useOperatorSplitting)
{
    mUseReactionDiffusionOperatorSplitting = useOperatorSplitting;
}

bool HeartConfig::GetUseReactionDiffusionOperatorSplitting()
{
    return mUseReactionDiffusionOperatorSplitting;
}

void HeartConfig::SetUseFixedNumberIterationsLinearSolver(bool useFixedNumberIterations, unsigned evaluateNumItsEveryNSolves)
{
    mUseFixedNumberIterations = useFixedNumberIterations;
    mEvaluateNumItsEveryNSolves = evaluateNumItsEveryNSolves;
}

bool HeartConfig::GetUseFixedNumberIterationsLinearSolver()
{
    return mUseFixedNumberIterations;
}

unsigned HeartConfig::GetEvaluateNumItsEveryNSolves()
{
    return mEvaluateNumItsEveryNSolves;
}

/**********************************************************************
 *                                                                    *
 *                                                                    *
 *            Utility methods for reading/transforming XML            *
 *                                                                    *
 *                                                                    *
 **********************************************************************/


void XmlTransforms::TransformArchiveDirectory(xercesc::DOMDocument* pDocument,
                                              xercesc::DOMElement* pRootElement)
{
    using namespace xercesc;
    std::vector<xercesc::DOMElement*> elts = XmlTools::FindElements(
        pRootElement,
        "ResumeSimulation/ArchiveDirectory");
    if (elts.size() > 0)
    {
        // We have an ArchiveDirectory element, so add the relative_to='chaste_test_output' attribute
        DOMElement* p_dir_elt = elts[0];
        p_dir_elt->setAttribute(X("relative_to"), X("chaste_test_output"));
    }
}

void XmlTransforms::TransformIonicModelDefinitions(xercesc::DOMDocument* pDocument,
                                                   xercesc::DOMElement* pRootElement)
{
    // Default ionic model
    std::vector<xercesc::DOMElement*> p_elt_list = XmlTools::FindElements(
        pRootElement,
        "Simulation/IonicModels/Default");
    if (p_elt_list.size() > 0)
    {
        assert(p_elt_list.size() == 1); // Asserted by schema
        XmlTools::WrapContentInElement(pDocument, p_elt_list[0], X("Hardcoded"));
        // Now do any region-specific definitions
        p_elt_list = XmlTools::FindElements(pRootElement, "Simulation/IonicModels/Region/IonicModel");
        for (unsigned i=0; i<p_elt_list.size(); i++)
        {
            XmlTools::WrapContentInElement(pDocument, p_elt_list[i], X("Hardcoded"));
        }
    }
}

void XmlTransforms::CheckForIluPreconditioner(xercesc::DOMDocument* pDocument,
                                              xercesc::DOMElement* pRootElement)
{
    std::vector<xercesc::DOMElement*> p_elt_list = XmlTools::FindElements(
        pRootElement,
        "Numerical/KSPPreconditioner");
    if (p_elt_list.size() > 0)
    {
        assert(p_elt_list.size() == 1); // Asserted by schema
        std::string text_value = X2C(p_elt_list[0]->getTextContent());
        if (text_value == "ilu")
        {
            EXCEPTION("PETSc does not have a parallel implementation of ilu, so we no longer allow it as an option.  Use bjacobi instead.");
        }
    }
}

void XmlTransforms::MoveConductivityHeterogeneities(xercesc::DOMDocument* pDocument,
                                                    xercesc::DOMElement* pRootElement)
{
    std::vector<xercesc::DOMElement*> p_elt_list = XmlTools::FindElements(
            pRootElement,
            "Simulation/ConductivityHeterogeneities");
    if (p_elt_list.size() > 0)
    {
        assert(p_elt_list.size() == 1); // Asserted by schema
        xercesc::DOMNode* p_parent = p_elt_list[0]->getParentNode();
        xercesc::DOMNode* p_child = p_parent->removeChild(p_elt_list[0]);
        std::vector<xercesc::DOMElement*> p_phys_list = XmlTools::FindElements(pRootElement, "Physiological");
        assert(p_phys_list.size() == 1); // Asserted by schema
        p_phys_list[0]->appendChild(p_child);
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation of the templated functions
/////////////////////////////////////////////////////////////////////
#define COVERAGE_IGNORE //These methods are covered above with DIM=1,2,3 but the instantiations may fail spuriously
/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template void HeartConfig::GetIonicModelRegions<3u>(std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >& , std::vector<cp::ionic_model_selection_type>&) const;
template void HeartConfig::GetStimuli<3u>(std::vector<boost::shared_ptr<AbstractStimulusFunction> >& , std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >& ) const;
template void HeartConfig::GetCellHeterogeneities<3u>(std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >& ,std::vector<double>& ,std::vector<double>& ,std::vector<double>& ,std::vector<std::map<std::string, double> >*) ;
template void HeartConfig::GetConductivityHeterogeneities<3u>(std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >& ,std::vector< c_vector<double,3> >& ,std::vector< c_vector<double,3> >& ) const;

template void HeartConfig::GetIonicModelRegions<2u>(std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >& , std::vector<cp::ionic_model_selection_type>&) const;
template void HeartConfig::GetStimuli<2u>(std::vector<boost::shared_ptr<AbstractStimulusFunction> >& , std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >& ) const;
template void HeartConfig::GetCellHeterogeneities<2u>(std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >& ,std::vector<double>& ,std::vector<double>& ,std::vector<double>& ,std::vector<std::map<std::string, double> >*) ;
template void HeartConfig::GetConductivityHeterogeneities<2u>(std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >& ,std::vector< c_vector<double,3> >& ,std::vector< c_vector<double,3> >& ) const;

template void HeartConfig::GetIonicModelRegions<1u>(std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >& , std::vector<cp::ionic_model_selection_type>&) const;
template void HeartConfig::GetStimuli<1u>(std::vector<boost::shared_ptr<AbstractStimulusFunction> >& , std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >& ) const;
template void HeartConfig::GetCellHeterogeneities<1u>(std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >& ,std::vector<double>& ,std::vector<double>& ,std::vector<double>& ,std::vector<std::map<std::string, double> >*);
template void HeartConfig::GetConductivityHeterogeneities<1u>(std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >& ,std::vector< c_vector<double,3> >& ,std::vector< c_vector<double,3> >& ) const;

template void HeartConfig::GetPseudoEcgElectrodePositions(std::vector<ChastePoint<1u> >& rPseudoEcgElectrodePositions) const;
template void HeartConfig::GetPseudoEcgElectrodePositions(std::vector<ChastePoint<2u> >& rPseudoEcgElectrodePositions) const;
template void HeartConfig::GetPseudoEcgElectrodePositions(std::vector<ChastePoint<3u> >& rPseudoEcgElectrodePositions) const;

template void HeartConfig::SetPseudoEcgElectrodePositions(const std::vector<ChastePoint<1u> >& rPseudoEcgElectrodePositions);
template void HeartConfig::SetPseudoEcgElectrodePositions(const std::vector<ChastePoint<2u> >& rPseudoEcgElectrodePositions);
template void HeartConfig::SetPseudoEcgElectrodePositions(const std::vector<ChastePoint<3u> >& rPseudoEcgElectrodePositions);
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
#undef COVERAGE_IGNORE //These methods are covered above with DIM=1,2,3 but the instantiations may fail spuriously


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HeartConfig)
