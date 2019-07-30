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
#include "AbstractCardiacCellInterface.hpp"
#include "DistributedVectorFactory.hpp"
#include "Warnings.hpp"
#include "CellMLToSharedLibraryConverter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCardiacSimulation : public CxxTest::TestSuite
{
    void setUp()
    {
        HeartEventHandler::Reset();
    }
//
//    void CreateOptionsFile(const OutputFileHandler& rHandler,
//                           const std::string& rModelName,
//                           const std::vector<std::string>& rArgs,
//                           const std::string& rExtraXml="")
//    {
//        if (PetscTools::AmMaster())
//        {
//            out_stream p_optfile = rHandler.OpenOutputFile(rModelName + "-conf.xml");
//            (*p_optfile) << "<?xml version='1.0'?>" << std::endl
//                         << "<pycml_config>" << std::endl
//                         << "<command_line_args>" << std::endl;
//            for (unsigned i=0; i<rArgs.size(); i++)
//            {
//                (*p_optfile) << "<arg>" << rArgs[i] << "</arg>" << std::endl;
//            }
//            (*p_optfile) << "</command_line_args>" << std::endl
//                         << rExtraXml
//                         << "</pycml_config>" << std::endl;
//            p_optfile->close();
//        }
//        PetscTools::Barrier("CreateOptionsFile");
//    }
public:

    void TestMono1dSmall()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        // Fox2002BackwardEuler cell model
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_small.xml", true);
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "mono_1d_small", false,
                                                 "SaveMono1D", "SimulationResults", true, 1e-6));

        /* If the above fails, and you are happy the new results are correct, uncomment the following line,
         * run the test, and then do
         cp /tmp/$USER/testoutput/SaveMono1D/SimulationResults.h5 heart/test/data/cardiac_simulations/mono_1d_small.h5
         */
        //assert(0);

        CardiacSimulation simulation2("heart/test/data/xml/monodomain1d_resume.xml", true);
    }

    void TestMono2dSmall()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        //Clear any warnings from previous tests
        Warnings::QuietDestroy();
        {
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
                AbstractCardiacCellInterface* p_cell = p_mono_problem->GetTissue()->GetCardiacCell(node_global_index);
                TS_ASSERT_DELTA(p_cell->GetParameter("membrane_fast_sodium_current_conductance"), 23 * 0.99937539038101175, 1e-6);
                TS_ASSERT_DELTA(p_cell->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance"), 0.282/3.0, 1e-6);
            }

            if (p_vector_factory->GetLocalOwnership() > 0)
            {
                TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), p_vector_factory->GetLocalOwnership());
                std::stringstream msg;
                msg << "Cannot apply drug to cell at node " << p_vector_factory->GetLow() << " as it has no parameter named 'not_a_current_conductance'.";
                TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), msg.str());
            }
            else
            {
                //There is now a warning on the top processes about the fact that no cells were assigned
                TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
                std::stringstream msg;
                msg << "No cells were assigned to process "<< PetscTools::GetMyRank();
                msg << " in AbstractCardiacTissue constructor. Advice: Make total number of processors no greater than number of nodes in the mesh";
                TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), msg.str());
            }
        }

        //Check that the post-processed file is there and remove it
        FileFinder pseudoecg("SaveMono2D/output/PseudoEcgFromElectrodeAt_0.05_0.05_0.dat", RelativeTo::ChasteTestOutput);

        if (PetscTools::AmMaster())
        {
            TS_ASSERT(pseudoecg.Exists()); // Only master tests. This prevents master from removing file before other processes have seen it
            pseudoecg.Remove();
            TS_ASSERT(pseudoecg.Exists() == false);
        }

        //Check that archive which has just been produced can be read
        CardiacSimulation simulation2("heart/test/data/xml/monodomain2d_resume.xml");

        //Check that the post-processed file is back after the simulation has been restarted
        TS_ASSERT(pseudoecg.Exists());  // (Should check that it's bigger than the one we deleted)

        Warnings::QuietDestroy();
    }

    /* Do the same as before but ask for post-processing after the simulation has been run and checkpointed. */
    void TestMono2dSmallAddPostprocessingOnResume()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        //Clear any warnings from previous tests
        Warnings::QuietDestroy();
        {
            CardiacSimulation simulation("heart/test/data/xml/monodomain2d_small2.xml", false, true);
        }

        //Check that the post-processed file is not produced in the original simulation
        FileFinder pseudoecg("SaveMono2D2/output/PseudoEcgFromElectrodeAt_0.05_0.05_0.dat", RelativeTo::ChasteTestOutput);
        TS_ASSERT(pseudoecg.Exists() == false);

        //Check that archive which has just been produced can be read
        CardiacSimulation simulation2("heart/test/data/xml/monodomain2d_resume2.xml");

        //Check that the post-processed file is present after the simulation has been restarted
        TS_ASSERT(pseudoecg.Exists());  // (Should check that it's bigger than the one we deleted)

        Warnings::QuietDestroy();
    }

    void TestMono3dSmall()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        CardiacSimulation simulation("heart/test/data/xml/monodomain3d_small.xml");
        //Check that archive which has just been produced can be read
        CardiacSimulation simulation2("heart/test/data/xml/monodomain3d_resume.xml");
    }

    void TestMono1dSodiumBlockBySettingNamedParameter()
    {
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_sodium_block.xml");
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "mono_1d_sodium_block", false,
                                                 "Mono1dSodiumBlock", "SimulationResults", true,
                                                 2.5e-3));

        // Test exception
        TS_ASSERT_THROWS_THIS(CardiacSimulation bad_param("heart/test/data/xml/bad_cell_parameter.xml"),
                              "No parameter named 'missing-parameter'.");
    }

    void TestMonoStimUsingEllipsoids()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        // Fox2002BackwardEuler cell model
        CardiacSimulation simulation("heart/test/data/xml/monodomain1d_stim_using_ellipsoid.xml", true);
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "monodomain1d_stim_using_ellipsoid", false,
                                                 "Mono1DStimUsingEllipsoid", "SimulationResults", true, 1e-6));
    }

    void TestBi1dSmall()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
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
    void TestBi2dSmall()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_small.xml");
        //Check that archive which has just been produced can be read
        CardiacSimulation simulation2("heart/test/data/xml/bidomain2d_resume.xml");
    }
    void TestBi3dSmall()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
       CardiacSimulation simulation("heart/test/data/xml/bidomain3d_small.xml");
        //Check that archive which has just been produced can be read
        CardiacSimulation simulation2("heart/test/data/xml/bidomain3d_resume.xml");
    }

    void TestBiWithBath1dSmall()
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

    void TestBiWithBath2dSmall()
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath2d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain_with_bath2d_resume.xml");
    }

    void TestBiWithBath3dSmall()
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath3d_small.xml");
        CardiacSimulation simulation2("heart/test/data/xml/bidomain_with_bath3d_resume.xml");
    }

    void TestBiWith2dHeterogeneousConductivities()
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain2d_heterogeneous.xml", true);
    }

    void TestCardiacSimulationBasicBidomainShort()
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain_with_bath simulation
        CardiacSimulation simulation("heart/test/data/xml/base_bidomain_short.xml");

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_bidomain_short_results", false,
                                                 "BaseBidomainShort", "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationBasicMonodomainShort()
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_short.xml");

       // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_monodomain_short_results", false,
                                                 "BaseMonodomainShort", "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationPostprocessMonodomain()
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/postprocess_monodomain_short.xml");
        std::string foldername = "PostprocessMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "postprocess_monodomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));
        {
            // look for the existence of post-processing files
            TS_ASSERT(FileFinder(foldername + "/output/Apd_90_minus_30_Map.dat", RelativeTo::ChasteTestOutput).Exists());
            TS_ASSERT(FileFinder(foldername + "/output/ConductionVelocityFromNode10.dat", RelativeTo::ChasteTestOutput).Exists());
            TS_ASSERT(FileFinder(foldername + "/output/ConductionVelocityFromNode20.dat", RelativeTo::ChasteTestOutput).Exists());
            TS_ASSERT(FileFinder(foldername + "/output/MaxUpstrokeVelocityMap_minus_30.dat", RelativeTo::ChasteTestOutput).Exists());
            TS_ASSERT(FileFinder(foldername + "/output/UpstrokeTimeMap_minus_30.dat", RelativeTo::ChasteTestOutput).Exists());
            TS_ASSERT(FileFinder(foldername + "/output/PseudoEcgFromElectrodeAt_0.05_0.05_0.dat", RelativeTo::ChasteTestOutput).Exists());
        }
    }

    void TestCardiacSimulationArchiveBidomain()
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_bidomain_short.xml");
        std::string foldername = "SaveBidomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT(CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_bidomain_short_results", false,
                                                foldername, "SimulationResults", true, 1e-5));

        FileFinder file(foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0",
                        RelativeTo::ChasteTestOutput);
        TS_ASSERT(file.Exists());
        /* If you want to update the .h5 results for this test for any reason, you need to stop the following test adding to them.
         * So uncomment the assert(), run the test, and then do:
         cp /tmp/chaste/testoutput/SaveBidomainShort/SimulationResults.h5 heart/test/data/cardiac_simulations/save_bidomain_short_results.h5
         */
        //assert(0);
    }

    // requires TestCardiacSimulationArchiveBidomain() to have been run
    void TestCardiacSimulationResumeBidomain()
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
                                                 foldername, "SimulationResults", true, 1e-5));
        //assert(0);
    }

    void TestCardiacSimulationArchiveMonodomain()
    {
        // Fox2002BackwardEuler cell model
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_monodomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));

        FileFinder file(foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0",
                        RelativeTo::ChasteTestOutput);
        TS_ASSERT(file.Exists());
        /* If you want to update the .h5 results for this test for any reason, you need to stop the following test adding to them.
         * So uncomment the assert(), run the test, and then do:
         cp /tmp/chaste/testoutput/SaveMonodomainShort/SimulationResults.h5 heart/test/data/cardiac_simulations/save_monodomain_short_results.h5
         */
        //assert(0);
    }

    // requires TestCardiacSimulationArchiveMonodomain() to have been run
    void TestCardiacSimulationResumeMonodomain()
    {
        // run a monodomain simulation
        HeartConfig::Instance()->SetSpaceDimension(1);
        CardiacSimulation simulation("heart/test/data/xml/resume_monodomain_short.xml");
        std::string foldername = "SaveMonodomainShort";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_short_results", false,
                                                 foldername, "SimulationResults", true, 1e-6));
    }

    void TestCardiacSimulationArchiveDynamic()
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

        FileFinder file(foldername + "_checkpoints/0.2ms/" + foldername + "_0.2ms/archive.arch.0",
                        RelativeTo::ChasteTestOutput);
        TS_ASSERT(file.Exists());

        /* If you want to update the .h5 results for this test for any reason, you need to stop the following lines adding to them.
         * So uncomment the assert(), run the test, and then do:
         cp /tmp/chaste/testoutput/SaveMonodomainDynamic/SimulationResults.h5 heart/test/data/cardiac_simulations/save_monodomain_dynamic.h5
         */
        //assert(0);

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
     * Note: from Chaste release 3.1 onward we no longer support Boost 1.33.
     * The earliest version of Boost supported is 1.34
     * Run TestCardiacSimulationArchiveBidomain on 4 processors to create the archive for this test,
     * and copy it to the repository using:
     *
       scons build=GccOpt_hostconfig,boost=1-34_4 test_suite=heart/test/TestCardiacSimulation.hpp
       cp -r /tmp/$USER/testoutput/SaveBidomainShort_checkpoints/0.2ms heart/test/data/checkpoint_migration_via_xml/
       rm -f heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort/progress_status.txt
       rm -f heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort_0.2ms/mesh.ncl
       rm -f heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort_0.2ms/ChasteParameters_?_?xsd
       rm -rf heart/test/data/checkpoint_migration_via_xml/0.2ms/SaveBidomainShort/output
     */
    void TestCardiacSimulationResumeMigration()
    {
        // We can only load simulations from CHASTE_TEST_OUTPUT, so copy the archives there
        std::string source_directory = "heart/test/data/checkpoint_migration_via_xml/0.2ms/";
        std::string folder_1 = "SaveBidomainShort";
        std::string folder_2 = "SaveBidomainShort_0.2ms";

        FileFinder to_directory1(OutputFileHandler::GetChasteTestOutputDirectory() + folder_1, RelativeTo::Absolute);
        FileFinder to_directory2(OutputFileHandler::GetChasteTestOutputDirectory() + folder_2, RelativeTo::Absolute);

        FileFinder from_directory1(source_directory + folder_1, RelativeTo::ChasteSourceRoot);
        FileFinder from_directory2(source_directory + folder_2, RelativeTo::ChasteSourceRoot);

        TRY_IF_MASTER(
            to_directory1.Remove();
            to_directory2.Remove();
            TS_ASSERT_EQUALS(to_directory1.Exists(), false);
            TS_ASSERT_EQUALS(to_directory2.Exists(), false);
            from_directory1.CopyTo(to_directory1);
            from_directory2.CopyTo(to_directory2);
        );

        // Resume
        CardiacSimulation simulation("heart/test/data/xml/resume_migration.xml");

        // Compare results
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_bidomain_short_results", false,
                                                 "SaveBidomainShort", "SimulationResults", true, 1e-5));
    }

    void runSimulation(const std::string& rParametersFileName)
    {
        CardiacSimulation simulation(rParametersFileName);
    }

    void checkParameter(AbstractCardiacCellInterface* pCell, unsigned globalIndex)
    {
        // Check one of the parameters has been set in the central region
        TS_ASSERT_EQUALS(pCell->GetNumberOfParameters(), 3u);
        double expected_value;
        if (globalIndex <= 4 || globalIndex >= 16)
        {
            expected_value = 23.0;
        }
        else
        {
            expected_value = 0.0;
        }
        TS_ASSERT_DELTA(pCell->GetParameter("membrane_fast_sodium_current_conductance"), expected_value, 1e-12);
        TS_ASSERT_DELTA(pCell->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance"), 0.282, 1e-12);
        TS_ASSERT_DELTA(pCell->GetParameter("membrane_L_type_calcium_current_conductance"), 0.09, 1e-12);

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
                AbstractCardiacCellInterface* p_cell = p_mono_problem->GetTissue()->GetCardiacCell(node_global_index);
                checkParameter(p_cell, node_global_index);
            }
            // compare the files, using the CompareFilesViaHdf5DataReader() method
            TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "resume_monodomain_changing_parameter_results", false,
                                                     foldername, "SimulationResults", true));
        }
    }

    void TestResumeChangingSettings()
    {
        doTestResumeChangingSettings("heart/test/data/xml/save_monodomain_with_parameter.xml");
        doTestResumeChangingSettings("heart/test/data/xml/save_monodomain_with_parameter_append.xml");
    }

    void TestCardiacSimulationPatchwork()
    {
        OutputFileHandler handler("DynamicallyLoadedModel");
        FileFinder cellml_file("heart/dynamic/luo_rudy_1991_dyn.cellml", RelativeTo::ChasteSourceRoot);
        handler.CopyFileTo(cellml_file);

        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_patchwork.xml");
        std::string foldername = "Patchwork";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT(CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "patchwork_results", false,
                                                foldername, "SimulationResults", true, 1e-5));
    }

    void TestCardiacSimulationKirsten()
    {
        if (PetscTools::GetNumProcs() > 2u)
        {
            // There are only 2 layers of nodes in this simulation -- z length is equal to space step.
            TS_TRACE("This test is not suitable for more than 2 processes.");
            return;
        }

        CardiacSimulation simulation("heart/test/data/xml/base_monodomain_tt06_region.xml");
        std::string foldername = "Kirsten";
        TS_ASSERT(CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "Kirsten", false,
                                                          foldername, "SimulationResults", true, 5e-4)); // lower tolerance as comparing with non-backward-euler results.
    }

    void TestTransmuralCellularheterogeneities()
    {
        CardiacSimulation simulation("heart/test/data/xml/ChasteParametersCellHeterogeneities.xml");
        std::string foldername = "ChasteResults_heterogeneities";

        TS_ASSERT( CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "transmural_heterogeneities_results", false,
                   foldername, "SimulationResults", true));
    }

    void TestElectrodes()
    {
        CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath2d_electrodes.xml");
        std::string foldername = "ChasteResults_electrodes";

        TS_ASSERT( CompareFilesViaHdf5DataReaderGlobalNorm("heart/test/data/cardiac_simulations", "electrodes_results", false,
                   foldername, "SimulationResults", true, 1e-4));
    }

    void TestExceptions()
    {
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/monodomain8d_small.xml"),
                              "Space dimension not supported: should be 1, 2 or 3");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/bidomain8d_small.xml"),
                              "Space dimension not supported: should be 1, 2 or 3");

#ifndef __APPLE__
        ///\todo Passing error is fatal on Mac OSX
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/base_monodomain_frankenstein.xml"),
                              "XML parsing error in configuration file: heart/test/data/xml/base_monodomain_frankenstein.xml");
#endif

        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("no file"),
                              "Missing file parsing configuration file: no file");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation(""), "No XML file name given");

#ifdef __APPLE__
        FileFinder model("file_does_not_exist.dylib", RelativeTo::ChasteSourceRoot);
#else
        FileFinder model("file_does_not_exist.so", RelativeTo::ChasteSourceRoot);
#endif
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/missing_dynamic_model.xml"),
                              "Dynamically loadable cell model '" + model.GetAbsolutePath() + "' does not exist.");
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/bidomain_with_bath2d_noelectrodes.xml"),
                              "Simulation needs a stimulus (either <Stimuli> or <Electrodes>).");

#ifndef CHASTE_CAN_CHECKPOINT_DLLS
        TS_ASSERT_THROWS_THIS(CardiacSimulation simulation("heart/test/data/xml/dynamic_checkpoint.xml"),
                              "Checkpointing is not compatible with dynamically loaded cell models on Mac OS X.");
#endif
    }

    void TestDynamicallyLoadingCvodeCell()
    {
        // Coverage - using native CVODE cells should no longer throw
#ifdef CHASTE_CVODE
        OutputFileHandler handler_cvode("DynamicallyLoadedModelCvode");
        FileFinder cellml_file("heart/dynamic/luo_rudy_1991_dyn.cellml", RelativeTo::ChasteSourceRoot);
        handler_cvode.CopyFileTo(cellml_file);

        std::vector<std::string> args;
        args.push_back("--cvode");

        CellMLToSharedLibraryConverter::CreateOptionsFile(handler_cvode, "luo_rudy_1991_dyn", args);
        CardiacSimulation simulation("heart/test/data/xml/dynamic_cvode_model.xml");
#else
        std::cout << "CVODE is not enabled.\n";
#endif
    }
};

#endif /*TESTCARDIACSIMULATION_HPP_*/

