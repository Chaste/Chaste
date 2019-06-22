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

#include "CheckpointArchiveTypes.hpp"

#include "UblasCustomFunctions.hpp"

#include "AbstractChasteRegion.hpp"
#include "ArchiveLocationInfo.hpp"
#include "ChastePoint.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "HeartFileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "Version.hpp"
#include "Warnings.hpp"

#include "HeartRegionCodes.hpp"

#include "RegularStimulus.hpp"
#include "SimpleStimulus.hpp"

#include <cassert>
#include <fstream>
#include <istream>
#include <map>
#include <string>

#include <xsd/cxx/tree/exceptions.hxx>
#include "XmlTools.hpp"
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

#include <boost/current_function.hpp>
/**
 * This macro gives a friendly-ish exception method if you call a Get method that is
 * missing structure up the tree.
 *
 * @param test  the existence test code
 * @param path  the XML path you're looking for
 */
#define CHECK_EXISTS(test, path)                                                         \
    do                                                                                   \
    {                                                                                    \
        if (!test)                                                                       \
        {                                                                                \
            EXCEPTION("No XML element " << path << " found in parameters when calling '" \
                                        << BOOST_CURRENT_FUNCTION << "'");               \
        }                                                                                \
    } while (false)

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

    /**
     * Release 3.3 changed the default setting for meshalyzer visualization from
     * true to false.  Older parameters files need to retain the original default.
     *
     * @param pDocument  the DOM document containing the tree to be transformed
     * @param pRootElement  the root of the tree to be transformed
     */
    static void SetDefaultVisualizer(xercesc::DOMDocument* pDocument,
                                     xercesc::DOMElement* pRootElement);
};

//
// Default settings
//
#include "HeartConfigDefaults.hpp"

//
// Definition of static member variables
//
boost::shared_ptr<HeartConfig> HeartConfig::mpInstance;

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

    mpParameters = CreateDefaultParameters();
    //CheckTimeSteps(); // necessity of this line of code is not tested -- remove with caution!

    //initialise the member variable of the layers
    mEpiFraction = -1.0;
    mEndoFraction = -1.0;
    mMidFraction = -1.0;
    mUserAskedForCellularTransmuralHeterogeneities = false;
    // initialise to senseless values (these should be only 0, 1 and 2)
    // note: the 'minus 3' is for checking purposes as we need to add 0, 1 or 2 to this initial value
    // and UINT_MAX+1 seems to be 0
    mIndexMid = UINT_MAX - 3u;
    mIndexEpi = UINT_MAX - 3u;
    mIndexEndo = UINT_MAX - 3u;

    mUseReactionDiffusionOperatorSplitting = false;

    /// \todo #1703 This defaults should be set in HeartConfigDefaults.hpp
    mTissueIdentifiers.insert(0);
    mBathIdentifiers.insert(1);
}

HeartConfig::~HeartConfig()
{
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

    // Sometimes this method is called collectively, sometimes not,
    // in any case the below file operations only want to be performed by
    // the master - so exit. Caller takes responsibility for nice
    // exception handling.
    if (!PetscTools::AmMaster())
    {
        return;
    }

    out_stream p_parameters_file(new std::ofstream((output_dirname + "ChasteParameters.xml").c_str()));

    if (!p_parameters_file->is_open())
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
    map["cp31"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/3_1";
    map["cp31"].schema = "ChasteParameters_3_1.xsd";
    map["cp33"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/3_3";
    map["cp33"].schema = "ChasteParameters_3_3.xsd";
    map["cp34"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/3_4";
    map["cp34"].schema = "ChasteParameters_3_4.xsd";
    // We use 'cp' as prefix for the latest version to avoid having to change saved
    // versions for comparison at every release.
    map["cp"].name = "https://chaste.comlab.ox.ac.uk/nss/parameters/2017_1";
    map["cp"].schema = "ChasteParameters_2017_1.xsd";

    cp::ChasteParameters(*p_parameters_file, *mpParameters, map);

    // If we're archiving, try to save a copy of the latest schema too
    if (useArchiveLocationInfo)
    {
        CopySchema(output_dirname);
    }
}

void HeartConfig::LoadFromCheckpoint()
{
    /*
     *  This method implements the logic required by HeartConfig to be able to handle resuming a simulation via the executable.
     *
     *  When the control reaches the method mpParameters points to the file specified as resuming parameters.
     *  However SetParametersFile() will set this variable to point to the archived parameters.
     *
     *  We make a temporary copy of mpParameters so we don't lose its content.
     *  At the end of the method we update the new mpParameters with the resuming parameters.
     */
    assert(mpParameters.use_count() > 0);
    boost::shared_ptr<cp::chaste_parameters_type> p_new_parameters = mpParameters;

    /*
     *  When we unarchive a simulation, we load the old parameters file in order to inherit things such
     *  as default cell model, stimuli, heterogeneities, ... This has the side effect of inheriting the
     *  <CheckpointSimulation> element (if defined).
     *
     *  We disable checkpointing definition coming from the unarchived config file. We will enable it again
     *  if defined in the resume config file.
     */
    std::string parameters_filename_xml = ArchiveLocationInfo::GetArchiveDirectory() + "ChasteParameters.xml";
    mpParameters = ReadFile(parameters_filename_xml);
    mParametersFilePath.SetPath(parameters_filename_xml, RelativeTo::AbsoluteOrCwd);

    // Release 3.0 and earlier wrote a separate defaults file in the checkpoint
    std::string defaults_filename_xml = ArchiveLocationInfo::GetArchiveDirectory() + "ChasteDefaults.xml";
    if (FileFinder(defaults_filename_xml).Exists())
    {
        boost::shared_ptr<cp::chaste_parameters_type> p_defaults = ReadFile(defaults_filename_xml);
        MergeDefaults(mpParameters, p_defaults);
    }

    HeartConfig::Instance()->SetCheckpointSimulation(false);

    // If we are resuming a simulation, some parameters can be altered at this point.
    if (p_new_parameters->ResumeSimulation().present())
    {
        UpdateParametersFromResumeSimulation(p_new_parameters);
    }

    CheckTimeSteps(); // For consistency with SetParametersFile
}

void HeartConfig::CopySchema(const std::string& rToDirectory)
{
    // N.B. This method should only be called by the master process,
    // in a situation where it can handle EXCEPTION()s nicely, e.g.
    // TRY_IF_MASTER(CopySchema(...));

    std::string schema_name("ChasteParameters_2017_1.xsd");
    FileFinder schema_location("heart/src/io/" + schema_name, RelativeTo::ChasteSourceRoot);
    if (!schema_location.Exists())
    {
        // Try a relative path instead
        schema_location.SetPath(schema_name, RelativeTo::CWD);
        if (!schema_location.Exists())
        {
            // Warn the user
            std::string message("Unable to locate schema file " + schema_name + ". You will need to ensure it is available when resuming from the checkpoint.");
            WARN_ONCE_ONLY(message);
        }
    }
    if (schema_location.Exists())
    {
        FileFinder output_directory(rToDirectory, RelativeTo::Absolute);
        schema_location.CopyTo(output_directory);
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
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_3"] = root_dir + "ChasteParameters_3_3.xsd";
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_4"] = root_dir + "ChasteParameters_3_4.xsd";
    mSchemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2017_1"] = root_dir + "ChasteParameters_2017_1.xsd";
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
    // which returns a std::unique_ptr. We convert to a shared_ptr for easier semantics.
    try
    {
        // Make sure Xerces finalization happens
        XmlTools::Finalizer finalizer(false);
        // Parse XML to DOM
        auto p_doc = XmlTools::ReadXmlFile(rFileName, props);
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
        if (version < 3003) // Changes made in release 3.3
        {
            XmlTransforms::SetDefaultVisualizer(p_doc.get(), p_root_elt);
        }
        if (version < 3004) // Not the latest in release 3.4
        {
            XmlTools::SetNamespace(p_doc.get(), p_root_elt, "https://chaste.comlab.ox.ac.uk/nss/parameters/3_4");
        }
        if (version < 2017001) // Not the latest release
        {
            XmlTools::SetNamespace(p_doc.get(), p_root_elt, "https://chaste.comlab.ox.ac.uk/nss/parameters/2017_1");
        }
        // Parse DOM to object model
        boost::shared_ptr<cp::chaste_parameters_type> p_params(cp::ChasteParameters(*p_doc, ::xml_schema::flags::dont_initialize, props));
        // Get rid of the DOM stuff
        p_doc.reset();

        return boost::shared_ptr<cp::chaste_parameters_type>(p_params);
    }
    catch (const xml_schema::exception& e)
    {
        std::cerr << e << std::endl;
        // Make sure we don't store invalid parameters
        mpParameters.reset();
        EXCEPTION("XML parsing error in configuration file: " + rFileName);
    }
    catch (...)
    {
        // Make sure we don't store invalid parameters
        mpParameters.reset();
        throw;
    }
}

void HeartConfig::SetParametersFile(const std::string& rFileName)
{
    mpParameters = ReadFile(rFileName);
    MergeDefaults(mpParameters, CreateDefaultParameters());
    mParametersFilePath.SetPath(rFileName, RelativeTo::AbsoluteOrCwd);

    if (IsSimulationDefined())
    {
        CheckTimeSteps(); // Resume files might not have time steps defined
    }
}

FileFinder HeartConfig::GetParametersFilePath()
{
    return mParametersFilePath;
}

void HeartConfig::UpdateParametersFromResumeSimulation(boost::shared_ptr<cp::chaste_parameters_type> pResumeParameters)
{
    // Check for user foolishness
    if ((pResumeParameters->ResumeSimulation()->SpaceDimension() != HeartConfig::Instance()->GetSpaceDimension())
        || (pResumeParameters->ResumeSimulation()->Domain() != HeartConfig::Instance()->GetDomain()))
    {
        EXCEPTION("Problem type and space dimension should match when restarting a simulation.");
    }

    // New simulation duration
    HeartConfig::Instance()->SetSimulationDuration(pResumeParameters->ResumeSimulation()->SimulationDuration());

    // Stimulus definition.  For these we always replace any previous definitions (at least for now...)
    if (pResumeParameters->ResumeSimulation()->Stimuli().present())
    {
        mpParameters->Simulation()->Stimuli().set(pResumeParameters->ResumeSimulation()->Stimuli().get());
    }

    // Cell heterogeneities.  Note that while we copy the elements here, other code in CardiacSimulation actually updates
    // the loaded simulation to take account of the new settings.
    if (pResumeParameters->ResumeSimulation()->CellHeterogeneities().present())
    {
        if (!mpParameters->Simulation()->CellHeterogeneities().present())
        {
            // Original parameters had no heterogeneities, so just copy the whole element
            mpParameters->Simulation()->CellHeterogeneities().set(pResumeParameters->ResumeSimulation()->CellHeterogeneities().get());
        }
        else
        {
            // Need to append the new heterogeneity defitions to the original sequence
            XSD_SEQUENCE_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity)& new_seq = pResumeParameters->ResumeSimulation()->CellHeterogeneities()->CellHeterogeneity();
            XSD_SEQUENCE_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity)& orig_seq = mpParameters->Simulation()->CellHeterogeneities()->CellHeterogeneity();
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

    //Visualization parameters are no longer compulsory
    if (pResumeParameters->ResumeSimulation()->OutputVisualizer().present())
    {
        HeartConfig::Instance()->SetVisualizeWithParallelVtk(pResumeParameters->ResumeSimulation()->OutputVisualizer()->parallel_vtk() == cp::yesno_type::yes);
        HeartConfig::Instance()->SetVisualizeWithVtk(pResumeParameters->ResumeSimulation()->OutputVisualizer()->vtk() == cp::yesno_type::yes);
        HeartConfig::Instance()->SetVisualizeWithCmgui(pResumeParameters->ResumeSimulation()->OutputVisualizer()->cmgui() == cp::yesno_type::yes);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(pResumeParameters->ResumeSimulation()->OutputVisualizer()->meshalyzer() == cp::yesno_type::yes);
    }

    // Numerical parameters may be overridden
    {
        cp::numerical_type& r_resume = pResumeParameters->Numerical();
        cp::numerical_type& r_user = mpParameters->Numerical();
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
        cp::postprocessing_type& r_user = mpParameters->PostProcessing().get();
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
        if (!r_resume.PseudoEcgElectrodePosition().empty())
        {
            r_user.PseudoEcgElectrodePosition() = r_resume.PseudoEcgElectrodePosition();
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
    return mpParameters->Simulation().present();
}

bool HeartConfig::IsSimulationResumed() const
{
    return mpParameters->ResumeSimulation().present();
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
        CHECK_EXISTS(mpParameters->Simulation()->SpaceDimension().present(), "Simulation/SpaceDimension");
        return mpParameters->Simulation()->SpaceDimension().get();
    }
    else
    {
        return mpParameters->ResumeSimulation()->SpaceDimension();
    }
}

double HeartConfig::GetSimulationDuration() const
{
    if (IsSimulationDefined())
    {
        CHECK_EXISTS(mpParameters->Simulation()->SimulationDuration().present(), "Simulation/SimulationDuration");
        return mpParameters->Simulation()->SimulationDuration().get();
    }
    else // IsSimulationResumed
    {
        return mpParameters->ResumeSimulation()->SimulationDuration();
    }
}

cp::domain_type HeartConfig::GetDomain() const
{
    if (IsSimulationDefined())
    {
        CHECK_EXISTS(mpParameters->Simulation()->Domain().present(), "Simulation/Domain");
        return mpParameters->Simulation()->Domain().get();
    }
    else
    {
        return mpParameters->ResumeSimulation()->Domain();
    }
}

cp::ionic_model_selection_type HeartConfig::GetDefaultIonicModel() const
{
    CheckSimulationIsDefined("DefaultIonicModel");

    return mpParameters->Simulation()->IonicModels()->Default();
}

template <unsigned DIM>
void HeartConfig::GetIonicModelRegions(std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& definedRegions,
                                       std::vector<cp::ionic_model_selection_type>& ionicModels) const
{
    CheckSimulationIsDefined("IonicModelRegions");
    definedRegions.clear();
    ionicModels.clear();

    XSD_SEQUENCE_TYPE(cp::ionic_models_type::Region)& regions = mpParameters->Simulation()->IonicModels()->Region();

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
                        ChastePoint<DIM> chaste_point_a(point_a.x());
                        ChastePoint<DIM> chaste_point_b(point_b.x());
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b)));
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a(point_a.x(), point_a.y());
                        ChastePoint<DIM> chaste_point_b(point_b.x(), point_b.y());
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b)));
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a(point_a.x(), point_a.y(), point_a.z());
                        ChastePoint<DIM> chaste_point_b(point_b.x(), point_b.y(), point_b.z());
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b)));
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
                cp::point_type radii = ionic_model_region.Location().Ellipsoid()->Radii();
                switch (DIM)
                {
                    case 1:
                    {
                        ChastePoint<DIM> chaste_point_a(centre.x());
                        ChastePoint<DIM> chaste_point_b(radii.x());
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b)));
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a(centre.x(), centre.y());
                        ChastePoint<DIM> chaste_point_b(radii.x(), radii.y());
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b)));
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a(centre.x(), centre.y(), centre.z());
                        ChastePoint<DIM> chaste_point_b(radii.x(), radii.y(), radii.z());
                        definedRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b)));
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
        else if (ionic_model_region.Location().EpiLayer().present() || ionic_model_region.Location().MidLayer().present() || ionic_model_region.Location().EndoLayer().present())
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
    return mpParameters->Simulation()->Mesh().present();
}

bool HeartConfig::GetCreateMesh() const
{
    CheckSimulationIsDefined("Mesh");
    CHECK_EXISTS(IsMeshProvided(), "Simulation/Mesh");
    cp::mesh_type mesh = mpParameters->Simulation()->Mesh().get();
    return (mesh.Slab().present() || mesh.Sheet().present() || mesh.Fibre().present());
}

bool HeartConfig::GetCreateSlab() const
{
    CheckSimulationIsDefined("Mesh");
    CHECK_EXISTS(IsMeshProvided(), "Simulation/Mesh");
    cp::mesh_type mesh = mpParameters->Simulation()->Mesh().get();
    return (mesh.Slab().present());
}

bool HeartConfig::GetCreateSheet() const
{
    CheckSimulationIsDefined("Mesh");
    CHECK_EXISTS(IsMeshProvided(), "Simulation/Mesh");
    cp::mesh_type mesh = mpParameters->Simulation()->Mesh().get();
    return (mesh.Sheet().present());
}

bool HeartConfig::GetCreateFibre() const
{
    CheckSimulationIsDefined("Mesh");
    CHECK_EXISTS(IsMeshProvided(), "Simulation/Mesh");
    cp::mesh_type mesh = mpParameters->Simulation()->Mesh().get();
    return (mesh.Fibre().present());
}

bool HeartConfig::GetLoadMesh() const
{
    CheckSimulationIsDefined("Mesh");
    CHECK_EXISTS(IsMeshProvided(), "Simulation/Mesh");
    return (mpParameters->Simulation()->Mesh()->LoadMesh().present());
}

void HeartConfig::GetSlabDimensions(c_vector<double, 3>& slabDimensions) const
{
    CheckSimulationIsDefined("Slab");

    if (GetSpaceDimension() != 3 || !GetCreateSlab())
    {
        EXCEPTION("Tissue slabs can only be defined in 3D");
    }

    optional<cp::slab_type, false> slab_dimensions = mpParameters->Simulation()->Mesh()->Slab();

    slabDimensions[0] = slab_dimensions->x();
    slabDimensions[1] = slab_dimensions->y();
    slabDimensions[2] = slab_dimensions->z();
}

void HeartConfig::GetSheetDimensions(c_vector<double, 2>& sheetDimensions) const
{
    CheckSimulationIsDefined("Sheet");

    if (GetSpaceDimension() != 2 || !GetCreateSheet())
    {
        EXCEPTION("Tissue sheets can only be defined in 2D");
    }

    optional<cp::sheet_type, false> sheet_dimensions = mpParameters->Simulation()->Mesh()->Sheet();

    sheetDimensions[0] = sheet_dimensions->x();
    sheetDimensions[1] = sheet_dimensions->y();
}

void HeartConfig::GetFibreLength(c_vector<double, 1>& fibreLength) const
{
    CheckSimulationIsDefined("Fibre");

    if (GetSpaceDimension() != 1 || !GetCreateFibre())
    {
        EXCEPTION("Tissue fibres can only be defined in 1D");
    }

    optional<cp::fibre_type, false> fibre_length = mpParameters->Simulation()->Mesh()->Fibre();

    fibreLength[0] = fibre_length->x();
}

double HeartConfig::GetInterNodeSpace() const
{
    CheckSimulationIsDefined("InterNodeSpace");

    switch (GetSpaceDimension())
    {
        case 3:
            CHECK_EXISTS(GetCreateSlab(), "Simulation/Mesh/Slab");
            return mpParameters->Simulation()->Mesh()->Slab()->inter_node_space();
            break;
        case 2:
            CHECK_EXISTS(GetCreateSheet(), "Simulation/Mesh/Sheet");
            return mpParameters->Simulation()->Mesh()->Sheet()->inter_node_space();
            break;
        case 1:
            CHECK_EXISTS(GetCreateFibre(), "Simulation/Mesh/Fibre");
            return mpParameters->Simulation()->Mesh()->Fibre()->inter_node_space();
            break;
        default:
            NEVER_REACHED;
            // LCOV_EXCL_START
            return 0.0; //To fool the compiler
            // LCOV_EXCL_STOP
    }
}

std::string HeartConfig::GetMeshName() const
{
    CheckSimulationIsDefined("LoadMesh");
    CHECK_EXISTS(GetLoadMesh(), "Mesh/LoadMesh");

    return mpParameters->Simulation()->Mesh()->LoadMesh()->name();
}

cp::media_type HeartConfig::GetConductivityMedia() const
{
    CheckSimulationIsDefined("LoadMesh");
    CHECK_EXISTS(GetLoadMesh(), "Mesh/LoadMesh");

    return mpParameters->Simulation()->Mesh()->LoadMesh()->conductivity_media();
}

template <unsigned DIM>
void HeartConfig::GetStimuli(std::vector<boost::shared_ptr<AbstractStimulusFunction> >& rStimuliApplied,
                             std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rStimulatedAreas) const
{
    CheckSimulationIsDefined("Stimuli");

    if (!mpParameters->Simulation()->Stimuli().present())
    {
        // Finding no stimuli defined is allowed (although HeartConfigRelatedFactory does
        // throw an exception is no stimuli and no electrodes)
        return;
    }

    XSD_SEQUENCE_TYPE(cp::stimuli_type::Stimulus)
    stimuli = mpParameters->Simulation()->Stimuli()->Stimulus();

    for (XSD_ITERATOR_TYPE(cp::stimuli_type::Stimulus) i = stimuli.begin();
         i != stimuli.end();
         ++i)
    {
        cp::stimulus_type stimulus(*i);
        if (stimulus.Location().Cuboid().present() || stimulus.Location().Ellipsoid().present())
        {
            boost::shared_ptr<AbstractChasteRegion<DIM> > area_ptr;
            if (stimulus.Location().Cuboid().present())
            {
                cp::point_type point_a = stimulus.Location().Cuboid()->LowerCoordinates();
                cp::point_type point_b = stimulus.Location().Cuboid()->UpperCoordinates();
                switch (DIM)
                {
                    case 1:
                    {
                        ChastePoint<DIM> chaste_point_a(point_a.x());
                        ChastePoint<DIM> chaste_point_b(point_b.x());
                        area_ptr.reset(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b));
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a(point_a.x(), point_a.y());
                        ChastePoint<DIM> chaste_point_b(point_b.x(), point_b.y());
                        area_ptr.reset(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b));
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a(point_a.x(), point_a.y(), point_a.z());
                        ChastePoint<DIM> chaste_point_b(point_b.x(), point_b.y(), point_b.z());
                        area_ptr.reset(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b));
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
                cp::point_type radii = stimulus.Location().Ellipsoid()->Radii();
                switch (DIM)
                {
                    case 1:
                    {
                        ChastePoint<DIM> chaste_point_a(centre.x());
                        ChastePoint<DIM> chaste_point_b(radii.x());
                        area_ptr.reset(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b));
                        break;
                    }
                    case 2:
                    {
                        ChastePoint<DIM> chaste_point_a(centre.x(), centre.y());
                        ChastePoint<DIM> chaste_point_b(radii.x(), radii.y());
                        area_ptr.reset(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b));
                        break;
                    }
                    case 3:
                    {
                        ChastePoint<DIM> chaste_point_a(centre.x(), centre.y(), centre.z());
                        ChastePoint<DIM> chaste_point_b(radii.x(), radii.y(), radii.z());
                        area_ptr.reset(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b));
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
            rStimuliApplied.push_back(stim);
        }
        else if (stimulus.Location().EpiLayer().present() || stimulus.Location().MidLayer().present() || stimulus.Location().EndoLayer().present())
        {
            EXCEPTION("Definition of transmural layers is not yet supported for specifying stimulated areas, please use cuboids instead");
        }
        else
        {
            EXCEPTION("Invalid region type for stimulus definition");
        }
    }
}

template <unsigned DIM>
void HeartConfig::GetCellHeterogeneities(std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rCellHeterogeneityRegions,
                                         std::vector<double>& rScaleFactorGks,
                                         std::vector<double>& rScaleFactorIto,
                                         std::vector<double>& rScaleFactorGkr,
                                         std::vector<std::map<std::string, double> >* pParameterSettings)
{
    CheckSimulationIsDefined("CellHeterogeneities");

    if (!mpParameters->Simulation()->CellHeterogeneities().present())
    {
        // finding no heterogeneities defined is allowed
        return;
    }
    XSD_SEQUENCE_TYPE(cp::cell_heterogeneities_type::CellHeterogeneity)
    cell_heterogeneity
        = mpParameters->Simulation()->CellHeterogeneities()->CellHeterogeneity();

    bool user_supplied_negative_value = false;
    mUserAskedForCellularTransmuralHeterogeneities = false; // overwritten with true below if necessary
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

            ChastePoint<DIM> chaste_point_a(point_a.x(), point_a.y(), point_a.z());
            ChastePoint<DIM> chaste_point_b(point_b.x(), point_b.y(), point_b.z());

            rCellHeterogeneityRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b)));
        }
        else if (ht.Location().Ellipsoid().present())
        {
            user_asked_for_cuboids_or_ellipsoids = true;
            cp::point_type centre = ht.Location().Ellipsoid()->Centre();
            cp::point_type radii = ht.Location().Ellipsoid()->Radii();

            ChastePoint<DIM> chaste_point_a(centre.x(), centre.y(), centre.z());
            ChastePoint<DIM> chaste_point_b(radii.x(), radii.y(), radii.z());
            rCellHeterogeneityRegions.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b)));
        }
        else if (ht.Location().EpiLayer().present())
        {
            mEpiFraction = ht.Location().EpiLayer().get();

            mUserAskedForCellularTransmuralHeterogeneities = true;
            if (mEpiFraction < 0)
            {
                user_supplied_negative_value = true;
            }
            mIndexEpi = counter_of_heterogeneities;
        }
        else if (ht.Location().EndoLayer().present())
        {
            mEndoFraction = ht.Location().EndoLayer().get();

            mUserAskedForCellularTransmuralHeterogeneities = true;
            if (mEndoFraction < 0)
            {
                user_supplied_negative_value = true;
            }
            mIndexEndo = counter_of_heterogeneities;
        }
        else if (ht.Location().MidLayer().present())
        {
            mMidFraction = ht.Location().MidLayer().get();

            mUserAskedForCellularTransmuralHeterogeneities = true;
            if (mMidFraction < 0)
            {
                user_supplied_negative_value = true;
            }
            mIndexMid = counter_of_heterogeneities;
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

    if (mUserAskedForCellularTransmuralHeterogeneities)
    {
        // cuboids/ellipsoids and layers at the same time are not yet supported
        if (user_asked_for_cuboids_or_ellipsoids)
        {
            EXCEPTION("Specification of cellular heterogeneities by cuboids/ellipsoids and layers at the same time is not yet supported");
        }

        //check that the user supplied all three layers, the indexes should be 0, 1 and 2.
        // As they are initialised to a higher value, if their summation is higher than 3,
        // one (or more) is missing
        if ((mIndexMid + mIndexEndo + mIndexEpi) > 3)
        {
            EXCEPTION("Three specifications of layers must be supplied");
        }
        if (fabs((mEndoFraction + mMidFraction + mEpiFraction) - 1) > 1e-2)
        {
            EXCEPTION("Summation of epicardial, midmyocardial and endocardial fractions should be 1");
        }
        if (user_supplied_negative_value)
        {
            EXCEPTION("Fractions must be positive");
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
    return mpParameters->Physiological().ConductivityHeterogeneities().present();
}

template <unsigned DIM>
void HeartConfig::GetConductivityHeterogeneities(
    std::vector<boost::shared_ptr<AbstractChasteRegion<DIM> > >& rConductivitiesHeterogeneityAreas,
    std::vector<c_vector<double, 3> >& rIntraConductivities,
    std::vector<c_vector<double, 3> >& rExtraConductivities) const
{
    CheckSimulationIsDefined("ConductivityHeterogeneities");
    CHECK_EXISTS(GetConductivityHeterogeneitiesProvided(), "Physiological/ConductivityHeterogeneities");
    XSD_ANON_SEQUENCE_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity)& conductivity_heterogeneity = mpParameters->Physiological().ConductivityHeterogeneities()->ConductivityHeterogeneity();

    for (XSD_ANON_ITERATOR_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity) i = conductivity_heterogeneity.begin();
         i != conductivity_heterogeneity.end();
         ++i)
    {
        cp::conductivity_heterogeneity_type ht(*i);

        if (ht.Location().Cuboid().present())
        {
            cp::point_type point_a = ht.Location().Cuboid()->LowerCoordinates();
            cp::point_type point_b = ht.Location().Cuboid()->UpperCoordinates();
            ChastePoint<DIM> chaste_point_a(point_a.x(), point_a.y(), point_a.z());
            ChastePoint<DIM> chaste_point_b(point_b.x(), point_b.y(), point_b.z());
            rConductivitiesHeterogeneityAreas.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteCuboid<DIM>(chaste_point_a, chaste_point_b)));
        }
        else if (ht.Location().Ellipsoid().present())
        {
            cp::point_type centre = ht.Location().Ellipsoid()->Centre();
            cp::point_type radii = ht.Location().Ellipsoid()->Radii();
            ChastePoint<DIM> chaste_point_a(centre.x(), centre.y(), centre.z());
            ChastePoint<DIM> chaste_point_b(radii.x(), radii.y(), radii.z());
            rConductivitiesHeterogeneityAreas.push_back(boost::shared_ptr<AbstractChasteRegion<DIM> >(new ChasteEllipsoid<DIM>(chaste_point_a, chaste_point_b)));
        }
        else if (ht.Location().EpiLayer().present() || ht.Location().MidLayer().present() || ht.Location().EndoLayer().present())
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

            rIntraConductivities.push_back(Create_c_vector(intra_x, intra_y, intra_z));
        }
        else
        {
            c_vector<double, 3> intra_conductivities;
            GetIntracellularConductivities(intra_conductivities);
            rIntraConductivities.push_back(intra_conductivities);
        }

        if (ht.ExtracellularConductivities().present())
        {
            double extra_x = ht.ExtracellularConductivities()->longi();
            double extra_y = ht.ExtracellularConductivities()->trans();
            double extra_z = ht.ExtracellularConductivities()->normal();

            rExtraConductivities.push_back(Create_c_vector(extra_x, extra_y, extra_z));
        }
        else
        {
            c_vector<double, 3> extra_conductivities;
            GetExtracellularConductivities(extra_conductivities);
            rExtraConductivities.push_back(extra_conductivities);
        }
    }
}

std::string HeartConfig::GetOutputDirectory() const
{
    CheckSimulationIsDefined("Simulation/OutputDirectory");
    CHECK_EXISTS(mpParameters->Simulation()->OutputDirectory().present(), "Simulation/OutputDirectory");
    return mpParameters->Simulation()->OutputDirectory().get();
}

std::string HeartConfig::GetOutputFilenamePrefix() const
{
    CheckSimulationIsDefined("Simulation/OutputFilenamePrefix");
    CHECK_EXISTS(mpParameters->Simulation()->OutputFilenamePrefix().present(), "Simulation/OutputFilenamePrefix");
    return mpParameters->Simulation()->OutputFilenamePrefix().get();
}

bool HeartConfig::GetOutputVariablesProvided() const
{
    CheckSimulationIsDefined("OutputVariables");
    return mpParameters->Simulation()->OutputVariables().present();
}

void HeartConfig::GetOutputVariables(std::vector<std::string>& rOutputVariables) const
{
    CHECK_EXISTS(GetOutputVariablesProvided(), "Simulation/OutputVariables");
    XSD_SEQUENCE_TYPE(cp::output_variables_type::Var)& output_variables = mpParameters->Simulation()->OutputVariables()->Var();
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
    CheckSimulationIsDefined("OutputUsingOriginalNodeOrdering");
    bool result = false;
    if (mpParameters->Simulation()->OutputUsingOriginalNodeOrdering().present())
    {
        result = (mpParameters->Simulation()->OutputUsingOriginalNodeOrdering().get() == cp::yesno_type::yes);
    }
    return result;
}

bool HeartConfig::GetCheckpointSimulation() const
{
    return IsSimulationDefined() && mpParameters->Simulation()->CheckpointSimulation().present();
}

double HeartConfig::GetCheckpointTimestep() const
{
    CHECK_EXISTS(GetCheckpointSimulation(), "Simulation/CheckpointSimulation");
    return mpParameters->Simulation()->CheckpointSimulation()->timestep();
}

unsigned HeartConfig::GetMaxCheckpointsOnDisk() const
{
    CHECK_EXISTS(GetCheckpointSimulation(), "Simulation/CheckpointSimulation");
    return mpParameters->Simulation()->CheckpointSimulation()->max_checkpoints_on_disk();
}

HeartFileFinder HeartConfig::GetArchivedSimulationDir() const
{
    CheckResumeSimulationIsDefined("GetArchivedSimulationDir");

    return HeartFileFinder(mpParameters->ResumeSimulation()->ArchiveDirectory());
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 3>& rIntraConductivities) const
{
    CHECK_EXISTS(mpParameters->Physiological().IntracellularConductivities().present(), "Physiological/IntracellularConductivities");
    cp::conductivities_type intra_conductivities
        = mpParameters->Physiological().IntracellularConductivities().get();
    double intra_x_cond = intra_conductivities.longi();
    double intra_y_cond = intra_conductivities.trans();
    double intra_z_cond = intra_conductivities.normal();
    ;

    assert(intra_y_cond != DBL_MAX);
    assert(intra_z_cond != DBL_MAX);

    rIntraConductivities[0] = intra_x_cond;
    rIntraConductivities[1] = intra_y_cond;
    rIntraConductivities[2] = intra_z_cond;
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 2>& rIntraConductivities) const
{
    CHECK_EXISTS(mpParameters->Physiological().IntracellularConductivities().present(), "Physiological/IntracellularConductivities");
    cp::conductivities_type intra_conductivities
        = mpParameters->Physiological().IntracellularConductivities().get();
    double intra_x_cond = intra_conductivities.longi();
    double intra_y_cond = intra_conductivities.trans();

    assert(intra_y_cond != DBL_MAX);

    rIntraConductivities[0] = intra_x_cond;
    rIntraConductivities[1] = intra_y_cond;
}

void HeartConfig::GetIntracellularConductivities(c_vector<double, 1>& rIntraConductivities) const
{
    CHECK_EXISTS(mpParameters->Physiological().IntracellularConductivities().present(), "Physiological/IntracellularConductivities");
    cp::conductivities_type intra_conductivities
        = mpParameters->Physiological().IntracellularConductivities().get();
    double intra_x_cond = intra_conductivities.longi();

    rIntraConductivities[0] = intra_x_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 3>& rExtraConductivities) const
{
    CHECK_EXISTS(mpParameters->Physiological().ExtracellularConductivities().present(), "Physiological/ExtracellularConductivities");
    cp::conductivities_type extra_conductivities
        = mpParameters->Physiological().ExtracellularConductivities().get();
    double extra_x_cond = extra_conductivities.longi();
    double extra_y_cond = extra_conductivities.trans();
    double extra_z_cond = extra_conductivities.normal();
    ;

    assert(extra_y_cond != DBL_MAX);
    assert(extra_z_cond != DBL_MAX);

    rExtraConductivities[0] = extra_x_cond;
    rExtraConductivities[1] = extra_y_cond;
    rExtraConductivities[2] = extra_z_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 2>& rExtraConductivities) const
{
    CHECK_EXISTS(mpParameters->Physiological().ExtracellularConductivities().present(), "Physiological/ExtracellularConductivities");
    cp::conductivities_type extra_conductivities
        = mpParameters->Physiological().ExtracellularConductivities().get();
    double extra_x_cond = extra_conductivities.longi();
    double extra_y_cond = extra_conductivities.trans();

    assert(extra_y_cond != DBL_MAX);

    rExtraConductivities[0] = extra_x_cond;
    rExtraConductivities[1] = extra_y_cond;
}

void HeartConfig::GetExtracellularConductivities(c_vector<double, 1>& rExtraConductivities) const
{
    CHECK_EXISTS(mpParameters->Physiological().ExtracellularConductivities().present(), "Physiological/ExtracellularConductivities");
    cp::conductivities_type extra_conductivities
        = mpParameters->Physiological().ExtracellularConductivities().get();
    double extra_x_cond = extra_conductivities.longi();

    rExtraConductivities[0] = extra_x_cond;
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
        CHECK_EXISTS(mpParameters->Physiological().BathConductivity().present(), "Physiological/BathConductivity");
        return mpParameters->Physiological().BathConductivity().get();
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
            CHECK_EXISTS(mpParameters->Physiological().BathConductivity().present(), "Physiological/BathConductivity");
            return mpParameters->Physiological().BathConductivity().get();
        }
    }
}
const std::set<unsigned>& HeartConfig::rGetTissueIdentifiers()
{
    return mTissueIdentifiers;
}

const std::set<unsigned>& HeartConfig::rGetBathIdentifiers()
{
    return mBathIdentifiers;
}

double HeartConfig::GetSurfaceAreaToVolumeRatio() const
{
    CHECK_EXISTS(mpParameters->Physiological().SurfaceAreaToVolumeRatio().present(), "Physiological/SurfaceAreaToVolumeRatio");
    return mpParameters->Physiological().SurfaceAreaToVolumeRatio().get();
}

double HeartConfig::GetCapacitance() const
{
    CHECK_EXISTS(mpParameters->Physiological().Capacitance().present(), "Physiological/Capacitance");
    return mpParameters->Physiological().Capacitance().get();
}

double HeartConfig::GetOdeTimeStep() const
{
    CHECK_EXISTS(mpParameters->Numerical().TimeSteps().present(), "Numerical/TimeSteps");
    return mpParameters->Numerical().TimeSteps()->ode();
}

double HeartConfig::GetPdeTimeStep() const
{
    CHECK_EXISTS(mpParameters->Numerical().TimeSteps().present(), "Numerical/TimeSteps");
    return mpParameters->Numerical().TimeSteps()->pde();
}

double HeartConfig::GetPrintingTimeStep() const
{
    CHECK_EXISTS(mpParameters->Numerical().TimeSteps().present(), "Numerical/TimeSteps");
    return mpParameters->Numerical().TimeSteps()->printing();
}

bool HeartConfig::GetUseAbsoluteTolerance() const
{
    CHECK_EXISTS(mpParameters->Numerical().KSPTolerances().present(), "Numerical/KSPTolerances");
    return mpParameters->Numerical().KSPTolerances()->KSPAbsolute().present();
}

double HeartConfig::GetAbsoluteTolerance() const
{
    CHECK_EXISTS(mpParameters->Numerical().KSPTolerances().present(), "Numerical/KSPTolerances");
    if (!GetUseAbsoluteTolerance())
    {
        EXCEPTION("Absolute tolerance is not set in Chaste parameters");
    }
    return mpParameters->Numerical().KSPTolerances()->KSPAbsolute().get();
}

bool HeartConfig::GetUseRelativeTolerance() const
{
    CHECK_EXISTS(mpParameters->Numerical().KSPTolerances().present(), "Numerical/KSPTolerances");
    return mpParameters->Numerical().KSPTolerances()->KSPRelative().present();
}

double HeartConfig::GetRelativeTolerance() const
{
    CHECK_EXISTS(mpParameters->Numerical().KSPTolerances().present(), "Numerical/KSPTolerances");
    if (!GetUseRelativeTolerance())
    {
        EXCEPTION("Relative tolerance is not set in Chaste parameters");
    }
    return mpParameters->Numerical().KSPTolerances()->KSPRelative().get();
}

const char* HeartConfig::GetKSPSolver() const
{
    CHECK_EXISTS(mpParameters->Numerical().KSPSolver().present(), "Numerical/KSPSolver");
    switch (mpParameters->Numerical().KSPSolver().get())
    {
        case cp::ksp_solver_type::gmres:
            return "gmres";
        case cp::ksp_solver_type::cg:
            return "cg";
        case cp::ksp_solver_type::symmlq:
            return "symmlq";
        case cp::ksp_solver_type::chebychev:
            return "chebychev";
    }
    // LCOV_EXCL_START
    EXCEPTION("Unknown ksp solver");
    // LCOV_EXCL_STOP
}

const char* HeartConfig::GetKSPPreconditioner() const
{
    CHECK_EXISTS(mpParameters->Numerical().KSPPreconditioner().present(), "Numerical/KSPPreconditioner");
    switch (mpParameters->Numerical().KSPPreconditioner().get())
    {
        case cp::ksp_preconditioner_type::jacobi:
            return "jacobi";
        case cp::ksp_preconditioner_type::bjacobi:
            return "bjacobi";
        case cp::ksp_preconditioner_type::hypre:
            return "hypre";
        case cp::ksp_preconditioner_type::ml:
            return "ml";
        case cp::ksp_preconditioner_type::spai:
            return "spai";
        case cp::ksp_preconditioner_type::blockdiagonal:
            return "blockdiagonal";
        case cp::ksp_preconditioner_type::ldufactorisation:
            return "ldufactorisation";
        case cp::ksp_preconditioner_type::twolevelsblockdiagonal:
            return "twolevelsblockdiagonal";
        case cp::ksp_preconditioner_type::none:
            return "none";
    }
    // LCOV_EXCL_START
    EXCEPTION("Unknown ksp preconditioner");
    // LCOV_EXCL_STOP
}

DistributedTetrahedralMeshPartitionType::type HeartConfig::GetMeshPartitioning() const
{
    CHECK_EXISTS(mpParameters->Numerical().MeshPartitioning().present(), "Numerical/MeshPartitioning");
    switch (mpParameters->Numerical().MeshPartitioning().get())
    {
        case cp::mesh_partitioning_type::dumb:
            return DistributedTetrahedralMeshPartitionType::DUMB;
        case cp::mesh_partitioning_type::metis:
            return DistributedTetrahedralMeshPartitionType::METIS_LIBRARY;
        case cp::mesh_partitioning_type::parmetis:
            return DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY;
        case cp::mesh_partitioning_type::petsc:
            return DistributedTetrahedralMeshPartitionType::PETSC_MAT_PARTITION;
    }
    // LCOV_EXCL_START
    EXCEPTION("Unknown mesh partitioning type");
    // LCOV_EXCL_STOP
}

bool HeartConfig::IsAdaptivityParametersPresent() const
{
    bool IsAdaptivityParametersPresent = mpParameters->Numerical().AdaptivityParameters().present();
    if (IsAdaptivityParametersPresent)
    {
        WARNING("Use of the Adaptivity library is deprecated");
    }
    return IsAdaptivityParametersPresent;
}

/*
 * PostProcessing
 */

bool HeartConfig::IsPostProcessingSectionPresent() const
{
    return mpParameters->PostProcessing().present();
}

void HeartConfig::EnsurePostProcessingSectionPresent()
{
    ENSURE_SECTION_PRESENT(mpParameters->PostProcessing(), cp::postprocessing_type);
}

bool HeartConfig::IsPostProcessingRequested() const
{
    if (IsPostProcessingSectionPresent() == false)
    {
        return false;
    }
    else
    {
        return (IsApdMapsRequested() || IsUpstrokeTimeMapsRequested() || IsMaxUpstrokeVelocityMapRequested() || IsConductionVelocityMapsRequested() || IsAnyNodalTimeTraceRequested() || IsPseudoEcgCalculationRequested());
    }
}
bool HeartConfig::IsApdMapsRequested() const
{
    bool result = false;
    if (IsPostProcessingSectionPresent())
    {
        XSD_SEQUENCE_TYPE(cp::postprocessing_type::ActionPotentialDurationMap)& apd_maps = mpParameters->PostProcessing()->ActionPotentialDurationMap();
        result = (apd_maps.begin() != apd_maps.end());
    }
    return result;
}

void HeartConfig::GetApdMaps(std::vector<std::pair<double, double> >& apd_maps) const
{
    CHECK_EXISTS(IsApdMapsRequested(), "PostProcessing/ActionPotentialDurationMap");
    apd_maps.clear();

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ActionPotentialDurationMap)& apd_maps_sequence = mpParameters->PostProcessing()->ActionPotentialDurationMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::ActionPotentialDurationMap) i = apd_maps_sequence.begin();
         i != apd_maps_sequence.end();
         ++i)
    {
        std::pair<double, double> map(i->repolarisation_percentage(), i->threshold());

        apd_maps.push_back(map);
    }
}

bool HeartConfig::IsUpstrokeTimeMapsRequested() const
{
    bool result = false;
    if (IsPostProcessingSectionPresent())
    {
        XSD_SEQUENCE_TYPE(cp::postprocessing_type::UpstrokeTimeMap)& upstroke_map = mpParameters->PostProcessing()->UpstrokeTimeMap();
        result = (upstroke_map.begin() != upstroke_map.end());
    }
    return result;
}
void HeartConfig::GetUpstrokeTimeMaps(std::vector<double>& upstroke_time_maps) const
{
    CHECK_EXISTS(IsUpstrokeTimeMapsRequested(), "PostProcessing/UpstrokeTimeMap");
    assert(upstroke_time_maps.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::UpstrokeTimeMap)& upstroke_maps_sequence = mpParameters->PostProcessing()->UpstrokeTimeMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::UpstrokeTimeMap) i = upstroke_maps_sequence.begin();
         i != upstroke_maps_sequence.end();
         ++i)
    {
        upstroke_time_maps.push_back(i->threshold());
    }
}

bool HeartConfig::IsMaxUpstrokeVelocityMapRequested() const
{
    bool result = false;
    if (IsPostProcessingSectionPresent())
    {
        XSD_SEQUENCE_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap)& max_upstroke_velocity_map = mpParameters->PostProcessing()->MaxUpstrokeVelocityMap();
        result = (max_upstroke_velocity_map.begin() != max_upstroke_velocity_map.end());
    }
    return result;
}

void HeartConfig::GetMaxUpstrokeVelocityMaps(std::vector<double>& upstroke_velocity_maps) const
{
    CHECK_EXISTS(IsMaxUpstrokeVelocityMapRequested(), "PostProcessing/MaxUpstrokeVelocityMap");
    assert(upstroke_velocity_maps.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap)& max_upstroke_velocity_maps_sequence = mpParameters->PostProcessing()->MaxUpstrokeVelocityMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap) i = max_upstroke_velocity_maps_sequence.begin();
         i != max_upstroke_velocity_maps_sequence.end();
         ++i)
    {
        upstroke_velocity_maps.push_back(i->threshold());
    }
}

bool HeartConfig::IsConductionVelocityMapsRequested() const
{
    bool result = false;
    if (IsPostProcessingSectionPresent())
    {
        XSD_SEQUENCE_TYPE(cp::postprocessing_type::ConductionVelocityMap)& cond_vel_maps = mpParameters->PostProcessing()->ConductionVelocityMap();
        result = (cond_vel_maps.begin() != cond_vel_maps.end());
    }
    return result;
}

void HeartConfig::GetConductionVelocityMaps(std::vector<unsigned>& conduction_velocity_maps) const
{
    CHECK_EXISTS(IsConductionVelocityMapsRequested(), "PostProcessing/ConductionVelocityMap");
    assert(conduction_velocity_maps.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ConductionVelocityMap)& cond_vel_maps_sequence = mpParameters->PostProcessing()->ConductionVelocityMap();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::ConductionVelocityMap) i = cond_vel_maps_sequence.begin();
         i != cond_vel_maps_sequence.end();
         ++i)
    {
        conduction_velocity_maps.push_back(i->origin_node());
    }
}

bool HeartConfig::IsAnyNodalTimeTraceRequested() const
{
    bool result = false;
    if (IsPostProcessingSectionPresent())
    {
        XSD_SEQUENCE_TYPE(cp::postprocessing_type::TimeTraceAtNode)& requested_nodes = mpParameters->PostProcessing()->TimeTraceAtNode();
        result = (requested_nodes.begin() != requested_nodes.end());
    }
    return result;
}

void HeartConfig::GetNodalTimeTraceRequested(std::vector<unsigned>& rRequestedNodes) const
{
    CHECK_EXISTS(IsAnyNodalTimeTraceRequested(), "PostProcessing/TimeTraceAtNode");
    assert(rRequestedNodes.size() == 0);

    XSD_SEQUENCE_TYPE(cp::postprocessing_type::TimeTraceAtNode)& req_nodes = mpParameters->PostProcessing()->TimeTraceAtNode();

    for (XSD_ITERATOR_TYPE(cp::postprocessing_type::TimeTraceAtNode) i = req_nodes.begin();
         i != req_nodes.end();
         ++i)
    {
        rRequestedNodes.push_back(i->node_number());
    }
}

bool HeartConfig::IsPseudoEcgCalculationRequested() const
{
    bool result = false;
    if (IsPostProcessingSectionPresent())
    {
        XSD_SEQUENCE_TYPE(cp::postprocessing_type::PseudoEcgElectrodePosition)& electrodes = mpParameters->PostProcessing()->PseudoEcgElectrodePosition();
        result = (electrodes.begin() != electrodes.end());
    }
    return result;
}

template <unsigned SPACE_DIM>
void HeartConfig::GetPseudoEcgElectrodePositions(std::vector<ChastePoint<SPACE_DIM> >& rPseudoEcgElectrodePositions) const
{
    rPseudoEcgElectrodePositions.clear();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::PseudoEcgElectrodePosition)& electrodes = mpParameters->PostProcessing()->PseudoEcgElectrodePosition();
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

    return mpParameters->Simulation()->OutputVisualizer().present();
}

bool HeartConfig::GetVisualizeWithMeshalyzer() const
{
    if (!IsOutputVisualizerPresent())
    {
        return false;
    }
    else
    {
        return mpParameters->Simulation()->OutputVisualizer()->meshalyzer() == cp::yesno_type::yes;
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
        return mpParameters->Simulation()->OutputVisualizer()->cmgui() == cp::yesno_type::yes;
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
        return mpParameters->Simulation()->OutputVisualizer()->parallel_vtk() == cp::yesno_type::yes;
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
        return mpParameters->Simulation()->OutputVisualizer()->vtk() == cp::yesno_type::yes;
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
        return mpParameters->Simulation()->OutputVisualizer()->precision();
    }
}

bool HeartConfig::IsElectrodesPresent() const
{
    return mpParameters->Simulation()->Electrodes().present();
}

/*
 *  Set methods
 */
void HeartConfig::SetSpaceDimension(unsigned spaceDimension)
{
    mpParameters->Simulation()->SpaceDimension().set(spaceDimension);
}

void HeartConfig::SetSimulationDuration(double simulationDuration)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::time_type, time, simulationDuration, "ms");
    mpParameters->Simulation()->SimulationDuration().set(time);
}

void HeartConfig::SetDomain(const cp::domain_type& rDomain)
{
    mpParameters->Simulation()->Domain().set(rDomain);
}

void HeartConfig::SetDefaultIonicModel(const cp::ionic_models_available_type& rIonicModel)
{
    cp::ionic_model_selection_type ionic_model;
    ionic_model.Hardcoded(rIonicModel);
    cp::ionic_models_type container(ionic_model);
    mpParameters->Simulation()->IonicModels().set(container);
}

void HeartConfig::SetSlabDimensions(double x, double y, double z, double inter_node_space)
{
    if (!mpParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    cp::slab_type slab_definition(x, y, z, inter_node_space);
    mpParameters->Simulation()->Mesh()->Slab().set(slab_definition);
}

void HeartConfig::SetSheetDimensions(double x, double y, double inter_node_space)
{
    if (!mpParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    cp::sheet_type sheet_definition(x, y, inter_node_space);
    mpParameters->Simulation()->Mesh()->Sheet().set(sheet_definition);
}

void HeartConfig::SetFibreLength(double x, double inter_node_space)
{
    if (!mpParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    cp::fibre_type fibre_definition(x, inter_node_space);
    mpParameters->Simulation()->Mesh()->Fibre().set(fibre_definition);
}

void HeartConfig::SetMeshFileName(std::string meshPrefix, cp::media_type fibreDefinition)
{
    if (!mpParameters->Simulation()->Mesh().present())
    {
        XSD_CREATE_WITH_FIXED_ATTR(cp::mesh_type, mesh_to_load, "cm");
        mpParameters->Simulation()->Mesh().set(mesh_to_load);
    }

    XSD_NESTED_TYPE(cp::mesh_type::LoadMesh)
    mesh_prefix(meshPrefix, fibreDefinition);
    mpParameters->Simulation()->Mesh()->LoadMesh().set(mesh_prefix);
}

void HeartConfig::SetIonicModelRegions(std::vector<ChasteCuboid<3> >& rDefinedRegions,
                                       std::vector<cp::ionic_model_selection_type>& rIonicModels) const
{
    assert(rDefinedRegions.size() == rIonicModels.size());
    // You need to have defined a default model first...
    assert(mpParameters->Simulation()->IonicModels().present());
    XSD_SEQUENCE_TYPE(cp::ionic_models_type::Region)& regions = mpParameters->Simulation()->IonicModels()->Region();
    regions.clear();
    for (unsigned region_index = 0; region_index < rDefinedRegions.size(); region_index++)
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

void HeartConfig::SetConductivityHeterogeneities(std::vector<ChasteCuboid<3> >& rConductivityAreas,
                                                 std::vector<c_vector<double, 3> >& rIntraConductivities,
                                                 std::vector<c_vector<double, 3> >& rExtraConductivities)
{
    assert(rConductivityAreas.size() == rIntraConductivities.size());
    assert(rIntraConductivities.size() == rExtraConductivities.size());

    XSD_ANON_SEQUENCE_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity)
    heterogeneities_container;

    for (unsigned region_index = 0; region_index < rConductivityAreas.size(); region_index++)
    {
        cp::point_type point_a(rConductivityAreas[region_index].rGetLowerCorner()[0],
                               rConductivityAreas[region_index].rGetLowerCorner()[1],
                               rConductivityAreas[region_index].rGetLowerCorner()[2]);

        cp::point_type point_b(rConductivityAreas[region_index].rGetUpperCorner()[0],
                               rConductivityAreas[region_index].rGetUpperCorner()[1],
                               rConductivityAreas[region_index].rGetUpperCorner()[2]);

        XSD_CREATE_WITH_FIXED_ATTR(cp::location_type, locn, "cm");
        locn.Cuboid().set(cp::box_type(point_a, point_b));
        cp::conductivity_heterogeneity_type ht(locn);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                    rIntraConductivities[region_index][0],
                                    rIntraConductivities[region_index][1],
                                    rIntraConductivities[region_index][2],
                                    "mS/cm");

        ht.IntracellularConductivities(intra);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                    rExtraConductivities[region_index][0],
                                    rExtraConductivities[region_index][1],
                                    rExtraConductivities[region_index][2],
                                    "mS/cm");

        ht.ExtracellularConductivities(extra);

        heterogeneities_container.push_back(ht);
    }

    XSD_ANON_TYPE(cp::physiological_type, ConductivityHeterogeneities)
    heterogeneities_object;
    heterogeneities_object.ConductivityHeterogeneity(heterogeneities_container);

    mpParameters->Physiological().ConductivityHeterogeneities().set(heterogeneities_object);
}

void HeartConfig::SetConductivityHeterogeneitiesEllipsoid(std::vector<ChasteEllipsoid<3> >& rConductivityAreas,
                                                          std::vector<c_vector<double, 3> >& rIntraConductivities,
                                                          std::vector<c_vector<double, 3> >& rExtraConductivities)
{
    assert(rConductivityAreas.size() == rIntraConductivities.size());
    assert(rIntraConductivities.size() == rExtraConductivities.size());

    XSD_ANON_SEQUENCE_TYPE(cp::physiological_type, ConductivityHeterogeneities, ConductivityHeterogeneity)
    heterogeneities_container;

    for (unsigned region_index = 0; region_index < rConductivityAreas.size(); region_index++)
    {
        cp::point_type centre(rConductivityAreas[region_index].rGetCentre()[0],
                              rConductivityAreas[region_index].rGetCentre()[1],
                              rConductivityAreas[region_index].rGetCentre()[2]);

        cp::point_type radii(rConductivityAreas[region_index].rGetRadii()[0],
                             rConductivityAreas[region_index].rGetRadii()[1],
                             rConductivityAreas[region_index].rGetRadii()[2]);

        XSD_CREATE_WITH_FIXED_ATTR(cp::location_type, locn, "cm");
        locn.Ellipsoid().set(cp::ellipsoid_type(centre, radii));
        cp::conductivity_heterogeneity_type ht(locn);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                    rIntraConductivities[region_index][0],
                                    rIntraConductivities[region_index][1],
                                    rIntraConductivities[region_index][2],
                                    "mS/cm");

        ht.IntracellularConductivities(intra);

        XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                    rExtraConductivities[region_index][0],
                                    rExtraConductivities[region_index][1],
                                    rExtraConductivities[region_index][2],
                                    "mS/cm");

        ht.ExtracellularConductivities(extra);

        heterogeneities_container.push_back(ht);
    }

    XSD_ANON_TYPE(cp::physiological_type, ConductivityHeterogeneities)
    heterogeneities_object;
    heterogeneities_object.ConductivityHeterogeneity(heterogeneities_container);

    mpParameters->Physiological().ConductivityHeterogeneities().set(heterogeneities_object);
}

void HeartConfig::SetOutputDirectory(const std::string& rOutputDirectory)
{
    mpParameters->Simulation()->OutputDirectory().set(rOutputDirectory);
}

void HeartConfig::SetOutputFilenamePrefix(const std::string& rOutputFilenamePrefix)
{
    mpParameters->Simulation()->OutputFilenamePrefix().set(rOutputFilenamePrefix);
}

void HeartConfig::SetOutputVariables(const std::vector<std::string>& rOutputVariables)
{
    if (!mpParameters->Simulation()->OutputVariables().present())
    {
        cp::output_variables_type variables_requested;
        mpParameters->Simulation()->OutputVariables().set(variables_requested);
    }

    XSD_SEQUENCE_TYPE(cp::output_variables_type::Var)& var_type_sequence = mpParameters->Simulation()->OutputVariables()->Var();
    // Erase or create a sequence
    var_type_sequence.clear();

    for (unsigned i = 0; i < rOutputVariables.size(); i++)
    {
        cp::var_type temp(rOutputVariables[i]);
        var_type_sequence.push_back(temp);
    }
}

void HeartConfig::SetOutputUsingOriginalNodeOrdering(bool useOriginal)
{
    //What if it doesn't exist?
    mpParameters->Simulation()->OutputUsingOriginalNodeOrdering().set(useOriginal ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetCheckpointSimulation(bool saveSimulation, double checkpointTimestep, unsigned maxCheckpointsOnDisk)
{
    if (saveSimulation)
    {
        // Make sure values for the optional parameters have been provided
        assert(checkpointTimestep != -1.0 && maxCheckpointsOnDisk != UINT_MAX);

        XSD_CREATE_WITH_FIXED_ATTR2(cp::simulation_type::XSD_NESTED_TYPE(CheckpointSimulation),
                                    cs,
                                    checkpointTimestep,
                                    maxCheckpointsOnDisk,
                                    "ms");
        mpParameters->Simulation()->CheckpointSimulation().set(cs);
    }
    else
    {
        mpParameters->Simulation()->CheckpointSimulation().reset();
    }

    CheckTimeSteps();
}

// Physiological

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 3>& rIntraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                rIntraConductivities[0],
                                rIntraConductivities[1],
                                rIntraConductivities[2],
                                "mS/cm");

    mpParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 2>& rIntraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                rIntraConductivities[0],
                                rIntraConductivities[1],
                                0.0, "mS/cm");

    mpParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetIntracellularConductivities(const c_vector<double, 1>& rIntraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, intra,
                                rIntraConductivities[0],
                                0.0, 0.0, "mS/cm");

    mpParameters->Physiological().IntracellularConductivities().set(intra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 3>& rExtraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                rExtraConductivities[0],
                                rExtraConductivities[1],
                                rExtraConductivities[2],
                                "mS/cm");

    mpParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 2>& rExtraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                rExtraConductivities[0],
                                rExtraConductivities[1],
                                0.0, "mS/cm");

    mpParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetExtracellularConductivities(const c_vector<double, 1>& rExtraConductivities)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::conductivities_type, extra,
                                rExtraConductivities[0],
                                0.0, 0.0, "mS/cm");

    mpParameters->Physiological().ExtracellularConductivities().set(extra);
}

void HeartConfig::SetBathConductivity(double bathConductivity)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::conductivity_type, cond, bathConductivity, "mS/cm");
    mpParameters->Physiological().BathConductivity().set(cond);
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
    if (tissueIds.empty() || bathIds.empty())
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
    mTissueIdentifiers = tissueIds;
    mBathIdentifiers = bathIds;
}

void HeartConfig::SetSurfaceAreaToVolumeRatio(double ratio)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::inverse_length_type, ratio_object, ratio, "1/cm");
    mpParameters->Physiological().SurfaceAreaToVolumeRatio().set(ratio_object);
}

void HeartConfig::SetCapacitance(double capacitance)
{
    XSD_CREATE_WITH_FIXED_ATTR1(cp::capacitance_type, capacitance_object, capacitance, "uF/cm^2");
    mpParameters->Physiological().Capacitance().set(capacitance_object);
}

// Numerical
void HeartConfig::SetOdePdeAndPrintingTimeSteps(double odeTimeStep, double pdeTimeStep, double printingTimeStep)
{
    XSD_CREATE_WITH_FIXED_ATTR3(cp::time_steps_type, time_steps,
                                odeTimeStep, pdeTimeStep, printingTimeStep, "ms");
    mpParameters->Numerical().TimeSteps().set(time_steps);
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

    if (GetPdeTimeStep() > GetPrintingTimeStep())
    {
        EXCEPTION("Printing time-step should not be smaller than PDE time-step");
    }

    if (!Divides(GetPdeTimeStep(), GetPrintingTimeStep()))
    {
        EXCEPTION("Printing time-step should be a multiple of PDE time step");
    }

    if (GetOdeTimeStep() > GetPdeTimeStep())
    {
        EXCEPTION("Ode time-step should not be greater than PDE time-step");
    }

    if (GetCheckpointSimulation())
    {
        if (GetCheckpointTimestep() <= 0.0)
        {
            EXCEPTION("Checkpoint time-step should be positive");
        }

        if (!Divides(GetPrintingTimeStep(), GetCheckpointTimestep()))
        {
            EXCEPTION("Checkpoint time-step should be a multiple of printing time-step");
        }
    }
}

void HeartConfig::SetUseRelativeTolerance(double relativeTolerance)
{
    ENSURE_SECTION_PRESENT(mpParameters->Numerical().KSPTolerances(), cp::ksp_tolerances_type);
    //Remove any reference to tolerances is user parameters
    mpParameters->Numerical().KSPTolerances()->KSPAbsolute().reset();
    mpParameters->Numerical().KSPTolerances()->KSPRelative().set(relativeTolerance);
}

void HeartConfig::SetUseAbsoluteTolerance(double absoluteTolerance)
{
    ENSURE_SECTION_PRESENT(mpParameters->Numerical().KSPTolerances(), cp::ksp_tolerances_type);
    //Remove any reference to tolerances is user parameters
    mpParameters->Numerical().KSPTolerances()->KSPRelative().reset();
    mpParameters->Numerical().KSPTolerances()->KSPAbsolute().set(absoluteTolerance);
}

void HeartConfig::SetKSPSolver(const char* kspSolver, bool warnOfChange)
{
    if (warnOfChange && strcmp(GetKSPSolver(), kspSolver) != 0)
    {
        //Warn
        WARNING("Code has changed the KSP solver type from " << GetKSPSolver() << " to " << kspSolver);
    }

    /* Note that changes in these conditions need to be reflected in the Doxygen*/
    if (strcmp(kspSolver, "gmres") == 0)
    {
        mpParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::gmres);
        return;
    }
    if (strcmp(kspSolver, "cg") == 0)
    {
        mpParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::cg);
        return;
    }
    if (strcmp(kspSolver, "symmlq") == 0)
    {
        mpParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::symmlq);
        return;
    }
    if (strcmp(kspSolver, "chebychev") == 0)
    {
        mpParameters->Numerical().KSPSolver().set(cp::ksp_solver_type::chebychev);
        return;
    }

    EXCEPTION("Unknown solver type provided");
}

void HeartConfig::SetKSPPreconditioner(const char* kspPreconditioner)
{
    /* Note that changes in these conditions need to be reflected in the Doxygen*/
    if (strcmp(kspPreconditioner, "jacobi") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::jacobi);
        return;
    }
    if (strcmp(kspPreconditioner, "bjacobi") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::bjacobi);
        return;
    }
    if (strcmp(kspPreconditioner, "hypre") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::hypre);
        return;
    }
    if (strcmp(kspPreconditioner, "ml") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::ml);
        return;
    }
    if (strcmp(kspPreconditioner, "spai") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::spai);
        return;
    }
    if (strcmp(kspPreconditioner, "twolevelsblockdiagonal") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::twolevelsblockdiagonal);
        return;
    }
    if (strcmp(kspPreconditioner, "blockdiagonal") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::blockdiagonal);
        return;
    }
    if (strcmp(kspPreconditioner, "ldufactorisation") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::ldufactorisation);
        return;
    }
    if (strcmp(kspPreconditioner, "none") == 0)
    {
        mpParameters->Numerical().KSPPreconditioner().set(cp::ksp_preconditioner_type::none);
        return;
    }

    EXCEPTION("Unknown preconditioner type provided");
}

void HeartConfig::SetMeshPartitioning(const char* meshPartioningMethod)
{
    /* Note that changes in these conditions need to be reflected in the Doxygen*/
    if (strcmp(meshPartioningMethod, "dumb") == 0)
    {
        mpParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::dumb);
        return;
    }
    if (strcmp(meshPartioningMethod, "metis") == 0)
    {
        mpParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::metis);
        return;
    }
    if (strcmp(meshPartioningMethod, "parmetis") == 0)
    {
        mpParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::parmetis);
        return;
    }
    if (strcmp(meshPartioningMethod, "petsc") == 0)
    {
        mpParameters->Numerical().MeshPartitioning().set(cp::mesh_partitioning_type::petsc);
        return;
    }

    EXCEPTION("Unknown mesh partitioning method provided");
}

void HeartConfig::SetApdMaps(const std::vector<std::pair<double, double> >& apdMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ActionPotentialDurationMap)& apd_maps_sequence
        = mpParameters->PostProcessing()->ActionPotentialDurationMap();
    //Erase or create a sequence
    apd_maps_sequence.clear();

    for (unsigned i = 0; i < apdMaps.size(); i++)
    {
        XSD_CREATE_WITH_FIXED_ATTR2(cp::apd_map_type, temp,
                                    apdMaps[i].first, apdMaps[i].second,
                                    "mV");
        apd_maps_sequence.push_back(temp);
    }
}

void HeartConfig::SetUpstrokeTimeMaps(std::vector<double>& upstrokeTimeMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::UpstrokeTimeMap)& var_type_sequence
        = mpParameters->PostProcessing()->UpstrokeTimeMap();

    //Erase or create a sequence
    var_type_sequence.clear();

    for (unsigned i = 0; i < upstrokeTimeMaps.size(); i++)
    {
        XSD_CREATE_WITH_FIXED_ATTR1(cp::upstrokes_map_type, temp,
                                    upstrokeTimeMaps[i],
                                    "mV");
        var_type_sequence.push_back(temp);
    }
}

void HeartConfig::SetMaxUpstrokeVelocityMaps(std::vector<double>& maxUpstrokeVelocityMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::MaxUpstrokeVelocityMap)& max_upstroke_velocity_maps_sequence
        = mpParameters->PostProcessing()->MaxUpstrokeVelocityMap();

    //Erase or create a sequence
    max_upstroke_velocity_maps_sequence.clear();

    for (unsigned i = 0; i < maxUpstrokeVelocityMaps.size(); i++)
    {
        XSD_CREATE_WITH_FIXED_ATTR1(cp::max_upstrokes_velocity_map_type, temp,
                                    maxUpstrokeVelocityMaps[i],
                                    "mV");

        max_upstroke_velocity_maps_sequence.push_back(temp);
    }
}

void HeartConfig::SetConductionVelocityMaps(std::vector<unsigned>& conductionVelocityMaps)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::ConductionVelocityMap)& conduction_velocity_maps_sequence
        = mpParameters->PostProcessing()->ConductionVelocityMap();

    //Erase or create a sequence
    conduction_velocity_maps_sequence.clear();

    for (unsigned i = 0; i < conductionVelocityMaps.size(); i++)
    {
        cp::conduction_velocity_map_type temp(conductionVelocityMaps[i]);
        conduction_velocity_maps_sequence.push_back(temp);
    }
}

void HeartConfig::SetRequestedNodalTimeTraces(std::vector<unsigned>& requestedNodes)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::TimeTraceAtNode)& requested_nodes_sequence
        = mpParameters->PostProcessing()->TimeTraceAtNode();

    //Erase or create a sequence
    requested_nodes_sequence.clear();

    for (unsigned i = 0; i < requestedNodes.size(); i++)
    {
        cp::node_number_type temp(requestedNodes[i]);
        requested_nodes_sequence.push_back(temp);
    }
}

template <unsigned SPACE_DIM>
void HeartConfig::SetPseudoEcgElectrodePositions(const std::vector<ChastePoint<SPACE_DIM> >& rPseudoEcgElectrodePositions)
{
    EnsurePostProcessingSectionPresent();
    XSD_SEQUENCE_TYPE(cp::postprocessing_type::PseudoEcgElectrodePosition)& electrodes_sequence
        = mpParameters->PostProcessing()->PseudoEcgElectrodePosition();

    //Erase or create a sequence
    electrodes_sequence.clear();

    for (unsigned i = 0; i < rPseudoEcgElectrodePositions.size(); i++)
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
    ENSURE_SECTION_PRESENT(mpParameters->Simulation()->OutputVisualizer(), cp::output_visualizer_type);
}

void HeartConfig::SetVisualizeWithMeshalyzer(bool useMeshalyzer)
{
    EnsureOutputVisualizerExists();

    mpParameters->Simulation()->OutputVisualizer()->meshalyzer(
        useMeshalyzer ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizeWithCmgui(bool useCmgui)
{
    EnsureOutputVisualizerExists();

    mpParameters->Simulation()->OutputVisualizer()->cmgui(
        useCmgui ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizeWithVtk(bool useVtk)
{
    EnsureOutputVisualizerExists();

    mpParameters->Simulation()->OutputVisualizer()->vtk(
        useVtk ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizeWithParallelVtk(bool useParallelVtk)
{
    EnsureOutputVisualizerExists();

    mpParameters->Simulation()->OutputVisualizer()->parallel_vtk(
        useParallelVtk ? cp::yesno_type::yes : cp::yesno_type::no);
}

void HeartConfig::SetVisualizerOutputPrecision(unsigned numberOfDigits)
{
    EnsureOutputVisualizerExists();

    mpParameters->Simulation()->OutputVisualizer()->precision(numberOfDigits);
}

void HeartConfig::SetElectrodeParameters(bool groundSecondElectrode,
                                         unsigned index, double magnitude,
                                         double startTime, double duration)
{
    assert(index < 3);

    cp::axis_type axis = cp::axis_type::x;
    if (index == 1)
    {
        axis = cp::axis_type::y;
    }
    else if (index == 2)
    {
        axis = cp::axis_type::z;
    }

    XSD_CREATE_WITH_FIXED_ATTR1(cp::surface_stimulus_strength_type, strength, magnitude, "uA/cm^2");
    XSD_CREATE_WITH_FIXED_ATTR1(cp::time_type, start_time, startTime, "ms");
    XSD_CREATE_WITH_FIXED_ATTR1(cp::time_type, duration_time, duration, "ms");

    if (!IsElectrodesPresent())
    {
        cp::electrodes_type element(groundSecondElectrode ? cp::yesno_type::yes : cp::yesno_type::no,
                                    axis,
                                    strength,
                                    start_time,
                                    duration_time);
        mpParameters->Simulation()->Electrodes().set(element);
    }
    else
    {
        mpParameters->Simulation()->Electrodes()->GroundSecondElectrode(groundSecondElectrode ? cp::yesno_type::yes : cp::yesno_type::no);
        mpParameters->Simulation()->Electrodes()->PerpendicularToAxis(axis);
        mpParameters->Simulation()->Electrodes()->Strength(strength);
        mpParameters->Simulation()->Electrodes()->StartTime(start_time);
        mpParameters->Simulation()->Electrodes()->Duration(duration_time);
    }
}

void HeartConfig::GetElectrodeParameters(bool& rGroundSecondElectrode,
                                         unsigned& rIndex, double& rMagnitude,
                                         double& rStartTime, double& rDuration)
{
    if (!IsElectrodesPresent())
    {
        EXCEPTION("Attempted to get electrodes that have not been defined.");
    }
    else
    {
        rGroundSecondElectrode = (mpParameters->Simulation()->Electrodes()->GroundSecondElectrode() == cp::yesno_type::yes);

        cp::axis_type axis = mpParameters->Simulation()->Electrodes()->PerpendicularToAxis();
        if (axis == cp::axis_type::x)
        {
            rIndex = 0;
        }
        else if (axis == cp::axis_type::y)
        {
            rIndex = 1;
        }
        else
        {
            rIndex = 2;
        }

        rMagnitude = mpParameters->Simulation()->Electrodes()->Strength();
        rStartTime = mpParameters->Simulation()->Electrodes()->StartTime();
        rDuration = mpParameters->Simulation()->Electrodes()->Duration();
    }
}

bool HeartConfig::GetUseStateVariableInterpolation() const
{
    // If it's an older version parameters & defaults (we're loading a checkpoint) say 'no'
    bool result = false;
    if (mpParameters->Numerical().UseStateVariableInterpolation().present())
    {
        result = mpParameters->Numerical().UseStateVariableInterpolation().get() == cp::yesno_type::yes;
    }
    return result;
}

void HeartConfig::SetUseStateVariableInterpolation(bool useStateVariableInterpolation)
{
    if (useStateVariableInterpolation)
    {
        mpParameters->Numerical().UseStateVariableInterpolation().set(cp::yesno_type::yes);
    }
    else
    {
        mpParameters->Numerical().UseStateVariableInterpolation().set(cp::yesno_type::no);
    }
}

bool HeartConfig::HasDrugDose() const
{
    return mpParameters->Physiological().ApplyDrug().present();
}

double HeartConfig::GetDrugDose() const
{
    CHECK_EXISTS(HasDrugDose(), "Physiological/ApplyDrug");
    return mpParameters->Physiological().ApplyDrug()->concentration();
}

void HeartConfig::SetDrugDose(double drugDose)
{
    if (!mpParameters->Physiological().ApplyDrug().present())
    {
        cp::apply_drug_type drug(drugDose);
        mpParameters->Physiological().ApplyDrug().set(drug);
    }
    else
    {
        mpParameters->Physiological().ApplyDrug()->concentration(drugDose);
    }
}

std::map<std::string, std::pair<double, double> > HeartConfig::GetIc50Values()
{
    CHECK_EXISTS(HasDrugDose(), "Physiological/ApplyDrug");
    std::map<std::string, std::pair<double, double> > ic50s;

    XSD_SEQUENCE_TYPE(cp::apply_drug_type::IC50)& ic50_seq = mpParameters->Physiological().ApplyDrug()->IC50();

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
    if (!mpParameters->Physiological().ApplyDrug().present())
    {
        SetDrugDose(0.0);
    }
    XSD_SEQUENCE_TYPE(cp::apply_drug_type::IC50)& ic50_seq = mpParameters->Physiological().ApplyDrug()->IC50();
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

//
// Purkinje methods
//

bool HeartConfig::HasPurkinje()
{
    CheckSimulationIsDefined("Purkinje");
    return mpParameters->Simulation()->Purkinje().present();
}

double HeartConfig::GetPurkinjeCapacitance()
{
    CHECK_EXISTS(mpParameters->Physiological().Purkinje().present(), "Physiological/Purkinje");
    CHECK_EXISTS(mpParameters->Physiological().Purkinje()->Capacitance().present(),
                 "Physiological/Purkinje/Capacitance");
    return mpParameters->Physiological().Purkinje()->Capacitance().get();
}

void HeartConfig::SetPurkinjeCapacitance(double capacitance)
{
    ENSURE_SECTION_PRESENT(mpParameters->Physiological().Purkinje(), cp::purkinje_physiological_type);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::capacitance_type, purk_Cm, capacitance, "uF/cm^2");
    mpParameters->Physiological().Purkinje()->Capacitance().set(purk_Cm);
}

double HeartConfig::GetPurkinjeSurfaceAreaToVolumeRatio()
{
    CHECK_EXISTS(mpParameters->Physiological().Purkinje().present(), "Physiological/Purkinje");
    CHECK_EXISTS(mpParameters->Physiological().Purkinje()->SurfaceAreaToVolumeRatio().present(),
                 "Physiological/Purkinje/SurfaceAreaToVolumeRatio");
    return mpParameters->Physiological().Purkinje()->SurfaceAreaToVolumeRatio().get();
}

void HeartConfig::SetPurkinjeSurfaceAreaToVolumeRatio(double ratio)
{
    ENSURE_SECTION_PRESENT(mpParameters->Physiological().Purkinje(), cp::purkinje_physiological_type);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::inverse_length_type, purk_Am, ratio, "1/cm");
    mpParameters->Physiological().Purkinje()->SurfaceAreaToVolumeRatio().set(purk_Am);
}

double HeartConfig::GetPurkinjeConductivity()
{
    CHECK_EXISTS(mpParameters->Physiological().Purkinje().present(), "Physiological/Purkinje");
    CHECK_EXISTS(mpParameters->Physiological().Purkinje()->Conductivity().present(),
                 "Physiological/Purkinje/Conductivity");
    return mpParameters->Physiological().Purkinje()->Conductivity().get();
}

void HeartConfig::SetPurkinjeConductivity(double conductivity)
{
    ENSURE_SECTION_PRESENT(mpParameters->Physiological().Purkinje(), cp::purkinje_physiological_type);
    XSD_CREATE_WITH_FIXED_ATTR1(cp::conductivity_type, purkinje_conductivity, conductivity, "mS/cm");
    mpParameters->Physiological().Purkinje()->Conductivity().set(purkinje_conductivity);
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
        for (unsigned i = 0; i < p_elt_list.size(); i++)
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

void XmlTransforms::SetDefaultVisualizer(xercesc::DOMDocument* pDocument,
                                         xercesc::DOMElement* pRootElement)
{
    std::vector<xercesc::DOMElement*> p_sim_list = XmlTools::FindElements(pRootElement, "Simulation");
    if (p_sim_list.size() > 0)
    {
        std::vector<xercesc::DOMElement*> p_viz_list = XmlTools::FindElements(p_sim_list[0], "OutputVisualizer");
        if (p_viz_list.empty())
        {
            // Create the element and set meshalyzer (only) to on
            xercesc::DOMElement* p_viz_elt = pDocument->createElementNS(X("https://chaste.comlab.ox.ac.uk/nss/parameters/3_3"), X("OutputVisualizer"));
            p_sim_list[0]->appendChild(p_viz_elt);
            p_viz_elt->setAttribute(X("meshalyzer"), X("yes"));
        }
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation of the templated functions
/////////////////////////////////////////////////////////////////////
// LCOV_EXCL_START //These methods are covered above with DIM=1,2,3 but the instantiations may fail spuriously
/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template void HeartConfig::GetIonicModelRegions<3u>(std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >&, std::vector<cp::ionic_model_selection_type>&) const;
template void HeartConfig::GetStimuli<3u>(std::vector<boost::shared_ptr<AbstractStimulusFunction> >&, std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >&) const;
template void HeartConfig::GetCellHeterogeneities<3u>(std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<std::map<std::string, double> >*);
template void HeartConfig::GetConductivityHeterogeneities<3u>(std::vector<boost::shared_ptr<AbstractChasteRegion<3u> > >&, std::vector<c_vector<double, 3> >&, std::vector<c_vector<double, 3> >&) const;

template void HeartConfig::GetIonicModelRegions<2u>(std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >&, std::vector<cp::ionic_model_selection_type>&) const;
template void HeartConfig::GetStimuli<2u>(std::vector<boost::shared_ptr<AbstractStimulusFunction> >&, std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >&) const;
template void HeartConfig::GetCellHeterogeneities<2u>(std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<std::map<std::string, double> >*);
template void HeartConfig::GetConductivityHeterogeneities<2u>(std::vector<boost::shared_ptr<AbstractChasteRegion<2u> > >&, std::vector<c_vector<double, 3> >&, std::vector<c_vector<double, 3> >&) const;

template void HeartConfig::GetIonicModelRegions<1u>(std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >&, std::vector<cp::ionic_model_selection_type>&) const;
template void HeartConfig::GetStimuli<1u>(std::vector<boost::shared_ptr<AbstractStimulusFunction> >&, std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >&) const;
template void HeartConfig::GetCellHeterogeneities<1u>(std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<std::map<std::string, double> >*);
template void HeartConfig::GetConductivityHeterogeneities<1u>(std::vector<boost::shared_ptr<AbstractChasteRegion<1u> > >&, std::vector<c_vector<double, 3> >&, std::vector<c_vector<double, 3> >&) const;

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
// LCOV_EXCL_STOP //These methods are covered above with DIM=1,2,3 but the instantiations may fail spuriously

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(HeartConfig)
