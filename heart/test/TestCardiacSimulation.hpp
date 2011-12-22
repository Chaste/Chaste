/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTCARDIACSIMULATION_HPP_
#define TESTCARDIACSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
#include <sstream>

#include "AbstractCardiacProblem.hpp"
#include "MonodomainProblem.hpp"

#include "CardiacSimulation.hpp"

#include "OutputFileHandler.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "FileFinder.hpp"
#include "PetscTools.hpp"
#include "HeartEventHandler.hpp"
#include "AbstractCardiacCell.hpp"
#include "DistributedVectorFactory.hpp"
#include "Warnings.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCardiacSimulation : public CxxTest::TestSuite
{
    void setUp()
    {
        HeartEventHandler::Reset();
    }

    void CreateOptionsFile(const OutputFileHandler& rHandler,
                           const std::string& rModelName,
                           const std::vector<std::string>& rArgs,
                           const std::string& rExtraXml="")
    {
        if (PetscTools::AmMaster())
        {
            out_stream p_optfile = rHandler.OpenOutputFile(rModelName + "-conf.xml");
            (*p_optfile) << "<?xml version='1.0'?>" << std::endl
                         << "<pycml_config>" << std::endl
                         << "<command_line_args>" << std::endl;
            for (unsigned i=0; i<rArgs.size(); i++)
            {
                (*p_optfile) << "<arg>" << rArgs[i] << "</arg>" << std::endl;
            }
            (*p_optfile) << "</command_line_args>" << std::endl
                         << rExtraXml
                         << "</pycml_config>" << std::endl;
            p_optfile->close();
        }
        PetscTools::Barrier("CreateOptionsFile");
    }
public:

    void TestMono1dSmall() throw(Exception)
    {
        // Fox2002BackwardEuler cell model
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_small.xml", true);
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "mono_1d_small", false,
                                                 "SaveMono1D", "SimulationResults", true, 1e-6));
        CardiacSimulation simulation2("heart/test/data/xml/monodomain1d_resume.xml", true);
    }

    void TestMono2dSmall() throw(Exception)
    {
        {
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
            CardiacSimulation simulation("heart/test/data/xml/monodomain2d_small.xml", false, true);
            boost::shared_ptr<AbstractUntemplatedCardiacProblem> p_problem = simulation.GetSavedProblem();
            TS_ASSERT(p_problem);
            MonodomainProblem<2,2>* p_mono_problem = dynamic_cast<MonodomainProblem<2,2>*>(p_problem.get());
            TS_ASSERT(p_mono_problem != NULL);
            DistributedVectorFactory* p_vector_factory = p_mono_problem->rGetMesh().GetDistributedVectorFactory();
            for (unsigned node_global_index = p_vector_factory->GetLow();
                 node_global_index < p_vector_factory->GetHigh();
                 node_global_index++)
            {
                AbstractCardiacCell* p_cell = p_mono_problem->GetTissue()->GetCardiacCell(node_global_index);
                TS_ASSERT_DELTA(p_cell->GetParameter("membrane_fast_sodium_current_conductance"), 23 * 0.99937539038101175, 1e-6);
                TS_ASSERT_DELTA(p_cell->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance"), 0.282/3.0, 1e-6);
            }
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), p_vector_factory->GetLocalOwnership());
            if (p_vector_factory->GetLocalOwnership() > 0)
            {
                std::stringstream msg;
                msg << "Cannot apply drug to cell at node " << p_vector_factory->GetLow() << " as it has no parameter named 'not_a_current_conductance'.";
                TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), msg.str());
            }
        }

        CardiacSimulation simulation2("heart/test/data/xml/monodomain2d_resume.xml");
        Warnings::QuietDestroy();
    }

    void TestMono3dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/monodomain3d_resume.xml");
    }

    void TestMono1dSodiumBlockBySettingNamedParameter() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_sodium_block.xml");
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "mono_1d_sodium_block", false,
                                                 "Mono1dSodiumBlock", "SimulationResults", true,
                                                 1e-3));

        // Test exception
        TS_ASSERT_THROWS_THIS(CardiacSimulation bad_param("heart/test/data/xml/bad_cell_parameter.xml"),
                              "No parameter named 'missing-parameter'.");
    }

    void TestMonoStimUsingEllipsoids() throw(Exception)
    {
        // Fox2002BackwardEuler cell model
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_stim_using_ellipsoid.xml", true);
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "monodomain1d_stim_using_ellipsoid", false,
                                                 "Mono1DStimUsingEllipsoid", "SimulationResults", true, 1e-6));
    }

    void TestBi1dSmall() throw(Exception)
    {
        { CardiacSimulation simulation("heart/test/data/xml/bidomain1d_small.xml"); }
        { CardiacSimulation simulation2("heart/test/data/xml/bidomain1d_resume.xml"); }
        {
            // The default resume file specifies a simulation duration of zero.
            // In reality the user should edit the file to specify something sensible...
            FileFinder resume_xml("SaveBi1D_checkpoints/0.1ms/ResumeParameters.xml", RelativeTo::ChasteTestOutput);
            TS_ASSERT_THROWS_THIS(CardiacSimulation simulation(resume_xml.GetAbsolutePath()),
                                  "The simulation duration must be positive, not -0.1");
        }
    }
    void TestBi2dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain2d_resume.xml");
    }
    void TestBi3dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain3d_resume.xml");
    }

    void TestBiWithBath1dSmall() throw(Exception)
    {
        { CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath1d_small.xml"); }
        { CardiacSimulation simulation2("heart/test/data/xml/bidomain_with_bath1d_resume.xml"); }
        {
            // The default resume file specifies a simulation duration of zero.
            // In reality the user should edit the file to specify something sensible...
            FileFinder resume_xml("SaveBiWithBath1D_checkpoints/0.1ms/ResumeParameters.xml", RelativeTo::ChasteTestOutput);
            TS_ASSERT_THROWS_THIS(CardiacSimulation simulation(resume_xml.GetAbsolutePath()),
                                  "The simulation duration must be positive, not -0.1");
        }
    }

    void TestBiWithBath2dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath2d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain_with_bath2d_resume.xml");
    }

    void TestBiWithBath3dSmall() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain_with_bath3d_resume.xml");
    }

    void TestBiWith2dHeterogeneousConductivities() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_heterogeneous.xml", true);
    }

    void TestCardiacSimulationBasicBidomainShort() throw(Exception)
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain_with_bath simulation
        CardiacSimulation simulation("heart/test/data/xml/base_bidomain_short.xml");

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_bidomain_short_results", false,
                                                 "BaseBidomainShort", "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationBasicMonodomainShort() throw(Exception)
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_short.xml");
        std::string foldername = "BaseMonodomainShort";

       // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_monodomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationPostprocessMonodomain() throw(Exception)
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/postprocess_monodomain_short.xml");
        std::string foldername = "PostprocessMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "postprocess_monodomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationArchiveBidomain() throw(Exception)
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_bidomain_short.xml");
        std::string foldername = "SaveBidomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT(CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_bidomain_short_results", false,
                                                foldername, "SimulationResults", true, 1e-6));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0);
    }

    // requires TestCardiacSimulationArchiveBidomain() to have been run
    void TestCardiacSimulationResumeBidomain() throw(Exception)
    {
        // run a bidomain simulation
        HeartConfig::Instance()->SetSpaceDimension(1);
        FileFinder resume_xml("heart/test/data/xml/resume_bidomain_short.xml", RelativeTo::ChasteSourceRoot);
        OutputFileHandler checkpoint_dir("SaveBidomainShort_checkpoints/0.2ms", false);
        FileFinder copied_xml = checkpoint_dir.CopyFileTo(resume_xml);
        CardiacSimulation simulation(copied_xml.GetAbsolutePath());
        std::string foldername = "SaveBidomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_bidomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationArchiveMonodomain() throw(Exception)
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_monodomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0);
    }

    // requires TestCardiacSimulationArchiveMonodomain() to have been run
    void TestCardiacSimulationResumeMonodomain() throw(Exception)
    {
        // run a monodomain simulation
        HeartConfig::Instance()->SetSpaceDimension(1);
        CardiacSimulation simulation("heart/test/data/xml/resume_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationArchiveDynamic() throw(Exception)
    {
#ifdef CHASTE_CAN_CHECKPOINT_DLLS
        // run a monodomain simulation
        {
            CardiacSimulation simulation("heart/test/data/xml/save_monodomain_dynamic.xml");
        }
        std::string foldername = "SaveMonodomainDynamic";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_monodomain_dynamic", false,
                   foldername, "SimulationResults", true));

        std::string command = "test -e " +  OutputFileHandler::GetChasteTestOutputDirectory() + foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0";
        int return_value = system(command.c_str());
        TS_ASSERT_EQUALS(return_value,0);

        //resume the simulation
        {
            CardiacSimulation simulation("heart/test/data/xml/resume_monodomain_dynamic.xml");
        }

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_dynamic", false,
                   foldername, "SimulationResults", true));
#endif // CHASTE_CAN_CHECKPOINT_DLLS
    }

    /**
     * Run TestCardiacSimulationArchiveBidomain on 4 processors to create the archive for this test,
     * and copy it to the repository using:
     *
       scons build=GccOpt_hostconfig,boost=1-33-1_4 test_suite=heart/test/TestCardiacSimulation.hpp
       cp -r /tmp/$USER/testoutput/SaveBidomainShort_checkpoints/0.2ms heart/test/data/checkpoint_migration_via_xml/
       rm heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort/progress_status.txt
       rm -rf heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort/output
     */
    void TestCardiacSimulationResumeMigration() throw(Exception)
    {
        // We can only load simulations from CHASTE_TEST_OUTPUT, so copy the archives there
        std::string source_directory = "heart/test/data/checkpoint_migration_via_xml/0.2ms/";
        // Clear the target directories
        OutputFileHandler h1("SaveBidomainShort");
        OutputFileHandler h2("SaveBidomainShort_0.2ms");
        if (PetscTools::AmMaster())
        {
            ABORT_IF_NON0(system, "cp " + source_directory + "SaveBidomainShort/* " + h1.GetOutputDirectoryFullPath());
            ABORT_IF_NON0(system, "cp " + source_directory + "SaveBidomainShort_0.2ms/* " + h2.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestCardiacSimulationResumeMigration");

        // Resume
        CardiacSimulation simulation("heart/test/data/xml/resume_migration.xml");
        // Compare results
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_bidomain_short_results", false,
                                                 "SaveBidomainShort", "SimulationResults", true, 1e-6));
    }

    void runSimulation(const std::string& rParametersFileName)
    {
        CardiacSimulation simulation(rParametersFileName);
    }

    void checkParameter(AbstractCardiacCell* pCell, unsigned globalIndex)
    {
        // Check parameter has been set in the central region
        TS_ASSERT_EQUALS(pCell->GetNumberOfParameters(), 2u);
        double expected_value;
        if (globalIndex <= 4 || globalIndex >= 16)
        {
            expected_value = 23.0;
        }
        else
        {
            expected_value = 0.0;
        }
        double actual_value = pCell->GetParameter(0);
        TS_ASSERT_EQUALS(actual_value, expected_value);
        // Check stimulus has been replaced.  It started as 0-1ms at x<=0.02, and should now be 600-601ms at x<=0.02
        if (globalIndex < 3)
        {
            TS_ASSERT_EQUALS(pCell->GetStimulus(0.0), 0.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(0.5), 0.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(1.0), 0.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(599.9), 0.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(600.0), -200000.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(600.5), -200000.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(601.0), -200000.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(601.1), 0.0);
        }
        else
        {
            // Always zero...
            TS_ASSERT_EQUALS(pCell->GetStimulus(0.0), 0.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(0.5), 0.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(600.5), 0.0);
            TS_ASSERT_EQUALS(pCell->GetStimulus(10.0), 0.0);
        }
    }

    void doTestResumeChangingSettings(const std::string& rParametersFileName)
    {
        std::string foldername = "SaveMonodomainWithParameter";

        // Save
        runSimulation(rParametersFileName);
        // Just check that the checkpoint exists
        FileFinder archive(foldername + "_checkpoints/1ms/" + foldername + "_1ms/archive.arch.0", RelativeTo::ChasteTestOutput);
        TS_ASSERT(archive.Exists());

        { // Load
            CardiacSimulation simulation("heart/test/data/xml/resume_monodomain_changing_parameter.xml", false, true);

            boost::shared_ptr<AbstractUntemplatedCardiacProblem> p_problem = simulation.GetSavedProblem();
            TS_ASSERT(p_problem);
            MonodomainProblem<1,1>* p_mono_problem = dynamic_cast<MonodomainProblem<1,1>*>(p_problem.get());
            TS_ASSERT(p_mono_problem != NULL);
            DistributedVectorFactory* p_vector_factory = p_mono_problem->rGetMesh().GetDistributedVectorFactory();
            for (unsigned node_global_index = p_vector_factory->GetLow();
                 node_global_index < p_vector_factory->GetHigh();
                 node_global_index++)
            {
                AbstractCardiacCell* p_cell = p_mono_problem->GetTissue()->GetCardiacCell(node_global_index);
                checkParameter(p_cell, node_global_index);
            }
            // compare the files, using the CompareFilesViaHdf5DataReader() method
            TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_changing_parameter_results", false,
                                                     foldername, "SimulationResults", true));
        }
    }

    void TestResumeChangingSettings() throw(Exception)
    {
        doTestResumeChangingSettings("heart/test/data/xml/save_monodomain_with_parameter.xml");
        doTestResumeChangingSettings("heart/test/data/xml/save_monodomain_with_parameter_append.xml");
    }

    void TestCardiacSimulationPatchwork() throw(Exception)
    {
        OutputFileHandler handler("DynamicallyLoadedModel");
        if (PetscTools::AmMaster())
        {
            // Copy CellML file into output dir
            FileFinder cellml_file("heart/dynamic/luo_rudy_1991_dyn.cellml", RelativeTo::ChasteSourceRoot);
            ABORT_IF_NON0(system, "cp " + cellml_file.GetAbsolutePath() + " " + handler.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestCardiacSimulationPatchwork");

        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_patchwork.xml");
        std::string foldername = "Patchwork";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT(CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "patchwork_results", false,
                                                foldername, "SimulationResults", true, 1e-5));
    }

    void TestCardiacSimulationKirsten() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_kirsten.xml");
        std::string foldername = "Kirsten";
        TS_ASSERT(CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "Kirsten", false,
                                                          foldername, "SimulationResults", true));
    }

    void TestTransmuralCellularheterogeneities() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/ChasteParametersCellHeterogeneities.xml");
        std::string foldername = "ChasteResults_heterogeneities";

        TS_ASSERT( CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "transmural_heterogeneities_results", false,
                   foldername, "SimulationResults", true));


    }

    void TestElectrodes() throw(Exception)
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath2d_electrodes.xml");
        std::string foldername = "ChasteResults_electrodes";

        TS_ASSERT( CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "electrodes_results", false,
                   foldername, "SimulationResults", true, 5e-5));
    }

    void TestExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/monodomain8d_small.xml"),
                              "Space dimension not supported: should be 1, 2 or 3");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/bidomain8d_small.xml"),
                              "Space dimension not supported: should be 1, 2 or 3");

        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/base_monodomain_frankenstein.xml"),
                              "XML parsing error in configuration file: heart/test/data/xml/base_monodomain_frankenstein.xml");

        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("no file"),
                              "Missing file parsing configuration file: no file");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation(""), "No XML file name given");

        FileFinder model("file_does_not_exist.so", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/missing_dynamic_model.xml"),
                              "Dynamically loadable cell model '" + model.GetAbsolutePath() + "' does not exist.");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath2d_noelectrodes.xml"),
                              "Simulation needs a stimulus (either <Stimuli> or <Electrodes>).");

#ifndef CHASTE_CAN_CHECKPOINT_DLLS
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/dynamic_checkpoint.xml"),
                              "Checkpointing is not compatible with dynamically loaded cell models on Boost<1.37.");
#endif
        // Coverage - using CVODE should throw
#ifdef CHASTE_CVODE
        OutputFileHandler handler_cvode("DynamicallyLoadedModelCvode");
        if (PetscTools::AmMaster())
        {
            // Copy CellML file into output dir
            FileFinder cellml_file("heart/dynamic/luo_rudy_1991_dyn.cellml", RelativeTo::ChasteSourceRoot);
            ABORT_IF_NON0(system, "cp " + cellml_file.GetAbsolutePath() + " " + handler_cvode.GetOutputDirectoryFullPath());
        }
        PetscTools::Barrier("TestExceptions");
        std::vector<std::string> args;
        args.push_back("--cvode");
        CreateOptionsFile(handler_cvode, "luo_rudy_1991_dyn", args);
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/dynamic_cvode_model.xml"),
                              "CVODE cannot be used as a cell model solver in tissue simulations: do not use the --cvode flag.");
#endif
    }
};

#endif /*TESTCARDIACSIMULATION_HPP_*/

