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

#include "CardiacSimulationArchiver.hpp"  // Must go first
#include "CardiacSimulation.hpp"


boost::shared_ptr<AbstractUntemplatedCardiacProblem> CardiacSimulation::GetSavedProblem()
{
    return mSavedProblem;
}

std::string CardiacSimulation::BoolToString(bool yesNo)
{
    std::string result;
    if (yesNo)
    {
        result = "yes";
    }
    else
    {
        result = "no";
    }
    return result;
}

void CardiacSimulation::CreateResumeXmlFile(const std::string& rOutputDirectory, const std::string& rArchiveDirectory)
{
    OutputFileHandler handler(rOutputDirectory, false);
    if (PetscTools::AmMaster())
    {
        out_stream p_file = handler.OpenOutputFile("ResumeParameters.xml");
        (*p_file) << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
        (*p_file) << "<ChasteParameters xmlns='https://chaste.comlab.ox.ac.uk/nss/parameters/3_0' "
                  << "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' "
                  << "xsi:schemaLocation='https://chaste.comlab.ox.ac.uk/nss/parameters/3_0 ChasteParameters_3_0.xsd'>" << std::endl;
        (*p_file) << std::endl;
        (*p_file) << "    <ResumeSimulation>" << std::endl;
        (*p_file) << "        <ArchiveDirectory relative_to='this_file'>" << rArchiveDirectory << "</ArchiveDirectory>" << std::endl;
        (*p_file) << "        <SpaceDimension>" << HeartConfig::Instance()->GetSpaceDimension() << "</SpaceDimension>" << std::endl;
        (*p_file) << "        <SimulationDuration unit='ms'>0.0</SimulationDuration> <!-- Edit with new simulation duration. Please "
                  << "note that the simulation does not restart at t=0 but at the time where the checkpoint was created.-->" << std::endl;
        (*p_file) << "        <Domain>" << HeartConfig::Instance()->GetDomain() << "</Domain>" << std::endl;
        (*p_file) << "        <CheckpointSimulation timestep='" << HeartConfig::Instance()->GetCheckpointTimestep()
                  << "' unit='ms' max_checkpoints_on_disk='" << HeartConfig::Instance()->GetMaxCheckpointsOnDisk()
                  << "'/> <!-- This is optional; if not given, the loaded simulation will NOT itself be checkpointed -->" << std::endl;
        (*p_file) << "        <OutputVisualizer meshalyzer='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithMeshalyzer())
                  << "' vtk='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithVtk())
                  << "' parallel_vtk='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithParallelVtk())
                  << "' cmgui='" << BoolToString(HeartConfig::Instance()->GetVisualizeWithCmgui()) << "'/>" << std::endl;
        (*p_file) << "    </ResumeSimulation>" << std::endl;
        (*p_file) << std::endl;
        (*p_file) << "    <!-- These elements must exist, but their contents are ignored -->" << std::endl;
        (*p_file) << "    <Physiological/>" << std::endl;
        (*p_file) << "    <Numerical/>" << std::endl;
        (*p_file) << "</ChasteParameters>" << std::endl;
        p_file->close();
    }

    TRY_IF_MASTER(HeartConfig::Instance()->CopySchema(handler.GetOutputDirectoryFullPath()));
}

CardiacSimulation::CardiacSimulation(std::string parameterFileName,
                                     bool writeProvenanceInfo,
                                     bool saveProblemInstance)
    : mSaveProblemInstance(saveProblemInstance)
{
    // If we have been passed an XML file then parse the XML file, otherwise throw
    if (parameterFileName == "")
    {
        EXCEPTION("No XML file name given");
    }
    ReadParametersFromFile(parameterFileName);
    Run();
    HeartEventHandler::Headings();
    HeartEventHandler::Report();
    if (writeProvenanceInfo)
    {
        ExecutableSupport::SetOutputDirectory(HeartConfig::Instance()->GetOutputDirectory());
        ExecutableSupport::WriteProvenanceInfoFile();
        ExecutableSupport::WriteMachineInfoFile("machine_info");
    }
}

void CardiacSimulation::ReadParametersFromFile(std::string parameterFileName)
{
    // Ensure the singleton is in a clean state
    HeartConfig::Reset();
    try
    {
        // Try the hardcoded schema location first
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile(parameterFileName);
    }
    catch (Exception& e)
    {
        if (e.CheckShortMessageContains("Missing file parsing configuration") == "")
        {
            // Try using the schema location given in the XML
            HeartConfig::Reset();
            HeartConfig::Instance()->SetUseFixedSchemaLocation(false);
            HeartConfig::Instance()->SetParametersFile(parameterFileName);
        }
        else
        {
            throw e;
        }
    }
}


#define DOMAIN_CASE(VALUE, CLASS, DIM)    \
    case VALUE:                           \
    {                                     \
        CreateAndRun<CLASS<DIM>, DIM>();  \
        break;                            \
    }

#define DOMAIN_SWITCH(DIM)                                                     \
    switch (HeartConfig::Instance()->GetDomain())                              \
    {                                                                          \
        DOMAIN_CASE(cp::domain_type::Mono, MonodomainProblem, DIM)             \
        DOMAIN_CASE(cp::domain_type::Bi, BidomainProblem, DIM)                 \
        DOMAIN_CASE(cp::domain_type::BiWithBath, BidomainWithBathProblem, DIM) \
        default:                                                               \
            NEVER_REACHED;                                                     \
    }                                                                          \
    break
// Note that if the domain is not set correctly then the XML parser will have picked it up before now!
// Missing semi-colon after break so we can put it after the macro call.

void CardiacSimulation::Run()
{
    switch (HeartConfig::Instance()->GetSpaceDimension())
    {
        case 3:
        {
            DOMAIN_SWITCH(3);
        }
        case 2:
        {
            DOMAIN_SWITCH(2);
        }
        case 1:
        {
            DOMAIN_SWITCH(1);
        }
        default:
            // We could check for this too in the XML Schema...
            EXCEPTION("Space dimension not supported: should be 1, 2 or 3");
    }
}

// These aren't needed externally
#undef DOMAIN_SWITCH
#undef DOMAIN_CASE

