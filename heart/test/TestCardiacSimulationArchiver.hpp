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

#ifndef TESTCARDIACSIMULATIONARCHIVER_HPP_
#define TESTCARDIACSIMULATIONARCHIVER_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp" // Needs to be before other Chaste code
#include "CardiacSimulationArchiver.hpp"

#include "Exception.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "ArchiveOpener.hpp"
#include "ChasteSyscalls.hpp"

#include "AbstractCardiacCellInterface.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "FoxModel2002BackwardEuler.hpp"
#include "FaberRudy2000.hpp"

#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "Electrodes.hpp"
#include "SimpleBathProblemSetup.hpp"

#include "FileComparison.hpp"
#include "NumericFileComparison.hpp"

/// For checkpoint migration tests
#define ABS_TOL 1e-6

/*
 * NB There are some tests in here that are only run when the boost version is 1.34
 * (i.e. on chaste-bob), so don't be too surprised if it fails there for just some
 * of the builds and doesn't elsewhere.
 *
 * If you intended to change the archive format, then you should have incremented the
 * BOOST_CLASS_VERSION(SomeClass, 1) to maintain backwards compatibility if possible.
 * You can then scroll down to lines ~426 and follow the instructions to generate some
 * new reference archives for the current version of the classes.
 */
class TestCardiacSimulationArchiver : public CxxTest::TestSuite
{
private:

    std::vector<double> mSolutionReplicated1d2ms;///<Used to test differences between tests

    void setUp()
    {
        HeartConfig::Instance()->Reset();
    }

public:

    /*
     *  Simple bidomain simulation to test against in TestArchivingWithHelperClass below
     */
    void TestSimpleBidomain1D()
    {
        OutputFileHandler handler("BidomainSimple1d",true); // This test was accruing output - true to wipe directory.

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // Check some voltages
        ReplicatableVector solution_replicated(bidomain_problem.GetSolution());
        double atol=5e-3;
        //Evens are the Vm voltages
        TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol); // V at 1
        TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol); // V at 5
        TS_ASSERT_DELTA(solution_replicated[20], -83.3582, atol); // V at 10
        //Odds are Phi_e
        TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol); // Phi_e at 0
        TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
        TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
        TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
        TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);

        // Save for checking in the next test
        for (unsigned index=0; index<solution_replicated.GetSize(); index++)
        {
            mSolutionReplicated1d2ms.push_back(solution_replicated[index]);
        }
    }

    void TestArchivingWithHelperClass()
    {
        std::string archive_dir("bidomain_problem_archive_helper");

        // Save
        {
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("BiProblemArchiveHelper");
            HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

            PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
            BidomainProblem<1> bidomain_problem( &cell_factory );

            bidomain_problem.Initialise();
            HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
            bidomain_problem.Solve();

            CardiacSimulationArchiver<BidomainProblem<1> >::Save(bidomain_problem, archive_dir, false);
        }

        // Load and run, outputting to a different directory.
        {
            OutputFileHandler handler("BidomainSimple1d_moved",true); // This test was accruing output - true to wipe directory.

            BidomainProblem<1> *p_bidomain_problem;
            p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load(archive_dir);

            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1d_moved");
            p_bidomain_problem->Solve();

            // check some voltages
            ReplicatableVector solution_replicated(p_bidomain_problem->GetSolution());
            TS_ASSERT_EQUALS(solution_replicated.GetSize(), mSolutionReplicated1d2ms.size()); //This in to make sure that the first test in the suite has been run!
            for (unsigned index=0; index<solution_replicated.GetSize(); index++)
            {
                // Shouldn't differ from the original run at all
                TS_ASSERT_DELTA(solution_replicated[index], mSolutionReplicated1d2ms[index],  5e-11);
            }

            Hdf5DataReader reader("BidomainSimple1d_moved", "BidomainLR91_1d", true);
            // This new file should go from 1.0 to 2.0, in printing time steps of 0.1.
            std::vector<double> times = reader.GetUnlimitedDimensionValues();
            TS_ASSERT_EQUALS(times.size(),11u);
            for (unsigned i=0; i<times.size(); i++)
            {
                TS_ASSERT_DELTA(times[i], 1.0 + 0.1*(double)(i), 1e-9);
            }

            // The last entry of the .h5 file over time should be the solution at present time.
            // There are 11 nodes in this mesh.
            for (unsigned node_idx = 0; node_idx < solution_replicated.GetSize()/2; node_idx++) // solution is V and phi_e
            {
                std::vector<double> v_over_time = reader.GetVariableOverTime("V", node_idx);
                std::vector<double> phi_over_time = reader.GetVariableOverTime("Phi_e", node_idx);
                TS_ASSERT_DELTA(v_over_time[v_over_time.size()-1u], solution_replicated[2*node_idx], 1e-4);
                TS_ASSERT_DELTA(phi_over_time[phi_over_time.size()-1u], solution_replicated[2*node_idx+1], 1e-4);
            }

            // Free memory
            delete p_bidomain_problem;
        }

        // Load and run, outputting to the same directory
        {
            BidomainProblem<1> *p_bidomain_problem;
            p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load(archive_dir);

            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            p_bidomain_problem->Solve();

            // check some voltages
            ReplicatableVector solution_replicated(p_bidomain_problem->GetSolution());
            TS_ASSERT_EQUALS(solution_replicated.GetSize(), mSolutionReplicated1d2ms.size()); //This in to make sure that the first test in the suite has been run!
            for (unsigned index=0; index<solution_replicated.GetSize(); index++)
            {
                //Shouldn't differ from the original run at all
                TS_ASSERT_DELTA(solution_replicated[index], mSolutionReplicated1d2ms[index],  5e-11);
            }
            // check output file contains results for the whole simulation
            TS_ASSERT(CompareFilesViaHdf5DataReader("BiProblemArchiveHelper", "BidomainLR91_1d", true,
                                                    "BidomainSimple1d", "BidomainLR91_1d", true));

            // Free memory
            delete p_bidomain_problem;
        }

        // Check we get an error if we try and extend the same h5 file again with results that aren't
        // consistent in time - original file went 0->1, block above extended it 0->2,
        // now we are asking for it to go 0->2,1->2 which doesn't make sense.
        {
            BidomainProblem<1> *p_bidomain_problem;
            p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load(archive_dir);

            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            TS_ASSERT_THROWS_CONTAINS(p_bidomain_problem->Solve(),
                                  "with results from time = 1, but it already contains results up to time = 2.");

            // Free memory
            delete p_bidomain_problem;
        }

        {
            // Coverage
            TS_ASSERT_THROWS_CONTAINS(CardiacSimulationArchiver<BidomainProblem<1> >::Load("absent_directory"),
                                      "Checkpoint directory does not exist: ");
            std::string empty_dir(archive_dir + "_empty");
            OutputFileHandler handler(empty_dir); // Ensure folder is empty
            TS_ASSERT_THROWS_CONTAINS(CardiacSimulationArchiver<BidomainProblem<1> >::Load(empty_dir),
                                      "Unable to open archive information file: ");
        }
        // Load into the wrong dimension
        {
            //When you try to load into the wrong dimension, the mesh reader (unarchiver)
            //will throw an exception.  Boost will eat this exception and set the mesh pointer
            //to NULL.  Your mileage may vary depending on Boost version.
            //We detect a NULL pointer and turn it back into an exception.
            TS_ASSERT_THROWS_THIS(CardiacSimulationArchiver<BidomainProblem<2> >::Load(archive_dir),
                "Failed to load from checkpoint because the dimensions of the archive do not match the object it's being read into.");

        }
    }

    /**
     *  Test used to generate data for the acceptance test resume_bidomain. We run the same simulation as in save_bidomain
     *  and archive it. resume_bidomain will load it and resume the simulation.
     *
     *  If the archiving format changes, both acceptance tests (save_bidomain and resume_bidomain) and the second part of this
     *  test will fail. Do the following to fix them.
     *
     *  Updated save_bidomain/ChasteResults_10ms_arch_0.chaste
     *    This can be easily done with texttest GUI. Run save_bidomain test to find that ChasteResults_10ms_arch_0.chaste
     *    contains differences (but NO other file). Highlight it and use the Save button on the window top left part.
     *
     * Note: from Chaste release 3.1 onward we no longer support Boost 1.33.
     * The earliest version of Boost supported in 1.34
     *
scons build=GccOpt_hostconfig,boost=1-34  test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp -r /tmp/$USER/testoutput/save_bidomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_bidomain/
cp -r /tmp/$USER/testoutput/SaveBidomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_bidomain/
     */
    void TestGenerateResultsForResumeBidomain()
    {
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/save_bidomain/ChasteParameters.xml");
        // We reset the mesh filename to include the relative path
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_1626_elements");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetMeshName(),
                         "mesh/test/data/cube_1626_elements");
        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(),
                         cp::ionic_models_available_type::Fox2002BackwardEuler);

        HeartConfig::Instance()->SetOutputDirectory("SaveBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("Bidomain_3d");

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<CellFoxModel2002FromCellMLBackwardEuler, 3> cell_factory(-80000.0, 1.0);
        BidomainProblem<3> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, "save_bidomain", false);

        // For coverage, pretend we're resuming
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/resume_bidomain/ChasteParameters.xml");

        if (PetscTools::IsSequential())
        {
            FileFinder archive_dir("apps/texttest/chaste/resume_bidomain/save_bidomain", RelativeTo::ChasteSourceRoot);
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, "archive.arch");
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            BidomainProblem<3> *p_bidomain_problem;
            TS_ASSERT_THROWS_NOTHING((*p_arch) >> p_bidomain_problem);
            //Double check that the acceptance test archive really is valid
            TS_ASSERT_EQUALS(p_bidomain_problem->mMeshFilename, "");
            TS_ASSERT_EQUALS(p_bidomain_problem->mPrintOutput, true);
            TS_ASSERT_EQUALS(p_bidomain_problem->mNodesToOutput.size(), 0u);
            TS_ASSERT_EQUALS(p_bidomain_problem->mCurrentTime, 10.0);
            delete p_bidomain_problem;
        }
    }


    /**
     *  Test used to generate data for the acceptance test resume_monodomain. We run the same simulation as in save_monodomain
     *  and archive it. resume_monodomain will load it and resume the simulation.
     *
     *  If the archiving format changes, both acceptance tests (save_monodomain and resume_monodomain) and the second part of this
     *  test will fail. Do the following to fix them.
     *
     *  Updated save_monodomain/ChasteResults_10ms_arch_0.chaste
     *    This can be easily done with texttest GUI. Run save_monodomain test to find that ChasteResults_10ms_arch_0.chaste
     *    contains differences (but NO other file). Highlight it and use the Save button on the window top left part.
     *
     * Note: from Chaste release 3.1 onward we no longer support Boost 1.33.
     * The earliest version of Boost supported in 1.34
     *
scons build=GccOpt_hostconfig,boost=1-34  test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp -r /tmp/$USER/testoutput/save_monodomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_monodomain/
cp -r /tmp/$USER/testoutput/SaveMonodomain/ ~/eclipse/workspace/Chaste/apps/texttest/chaste/resume_monodomain/
 */
    void TestGenerateResultsForResumeMonodomain()
    {
        HeartConfig::Instance()->SetUseFixedSchemaLocation(true);
        HeartConfig::Instance()->SetParametersFile("apps/texttest/chaste/save_monodomain/ChasteParameters.xml");

        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 10.0);
        TS_ASSERT(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().present());
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetDefaultIonicModel().Hardcoded().get(),
                         cp::ionic_models_available_type::Fox2002BackwardEuler);

        HeartConfig::Instance()->SetOutputDirectory("SaveMonodomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("Monodomain_2d");

        // This cell factory should apply the same stimulus described in the xml config file.
        PlaneStimulusCellFactory<CellFoxModel2002FromCellMLBackwardEuler, 2> cell_factory(-600000.0, 1.0);
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        CardiacSimulationArchiver<MonodomainProblem<2> >::Save(monodomain_problem, "save_monodomain", false);

        if (PetscTools::IsSequential())
        {
            FileFinder archive_dir("apps/texttest/chaste/resume_monodomain/save_monodomain", RelativeTo::ChasteSourceRoot);
            // We can't use CardiacSimulationArchiver::Load for this, as it can't read files in the repository
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, "archive.arch");
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            MonodomainProblem<2> *p_monodomain_problem;
            TS_ASSERT_THROWS_NOTHING((*p_arch) >> p_monodomain_problem);
            //Double check that the acceptance test archive really is valid
            TS_ASSERT_EQUALS(p_monodomain_problem->mMeshFilename, "");
            TS_ASSERT_EQUALS(p_monodomain_problem->mPrintOutput, true);
            TS_ASSERT_EQUALS(p_monodomain_problem->mNodesToOutput.size(), 0u);
            TS_ASSERT_EQUALS(p_monodomain_problem->mCurrentTime, 10.0);
            delete p_monodomain_problem;
        }
    }

    /***********************************************************************
     *                                                                     *
     *         Below this point are the checkpoint migration tests         *
     *                                                                     *
     * Note that archives are created for these using a variety of numbers *
     * of processes: 1, 2, and 3, so we get a good range of migration test *
     * cases.  The continuous build runs this testsuite with 1,2 processes *
     * and the nightly Parallel_Continuous tests use 3.                    *
     ***********************************************************************/

private:
    // Helper functions for the migration tests defined below.
    template<class Problem, unsigned DIM>
    Problem* DoMigrateAndBasicTests(const FileFinder& rArchiveDirectory,
                                    const std::string& rRefArchiveDir, // relative to CHASTE_TEST_OUTPUT
                                    const std::string& rNewArchiveDir, // relative to CHASTE_TEST_OUTPUT
                                    const unsigned totalNumCells,
                                    bool isDistributedMesh=false,
                                    double currentTime=0.0)
    {
        // Do the migration to sequential
        Problem* p_problem = CardiacSimulationArchiver<Problem>::Load(rArchiveDirectory);

        // Some basic tests that we have the right data
        TS_ASSERT_EQUALS(p_problem->mMeshFilename, "");
        TS_ASSERT_EQUALS(p_problem->mPrintOutput, true);
        TS_ASSERT_EQUALS(p_problem->mNodesToOutput.size(), 0u);
        TS_ASSERT_EQUALS(p_problem->mCurrentTime, currentTime);
        DistributedVectorFactory* p_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
        TS_ASSERT(p_factory->GetOriginalFactory());
        TS_ASSERT_EQUALS(p_factory->GetOriginalFactory()->GetProblemSize(), totalNumCells);
        unsigned local_num_cells = p_factory->GetLocalOwnership();
        TS_ASSERT_EQUALS(p_problem->GetTissue()->rGetCellsDistributed().size(), local_num_cells);
        TS_ASSERT_EQUALS(p_problem->rGetMesh().GetNumAllNodes(), totalNumCells);
        TS_ASSERT_EQUALS(p_problem->rGetMesh().GetNumNodes(), totalNumCells);
        TS_ASSERT_EQUALS(&(p_problem->rGetMesh()), p_problem->GetTissue()->pGetMesh());
        // Check the mesh is/isn't the parallel variety
        const DistributedTetrahedralMesh<DIM,DIM>* p_dist_mesh = dynamic_cast<const DistributedTetrahedralMesh<DIM,DIM>*>(p_problem->GetTissue()->pGetMesh());
        if (isDistributedMesh)
        {
            TS_ASSERT(p_dist_mesh != NULL);
        }
        else
        {
            TS_ASSERT(p_dist_mesh == NULL);
        }

        // All real cells should be at initial conditions if we're starting at t=0.
        if (currentTime == 0.0)
        {
            std::vector<double> inits;
            for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
            {
                AbstractCardiacCell* p_cell = static_cast<AbstractCardiacCell*>(p_problem->GetTissue()->GetCardiacCell(i));
                FakeBathCell* p_fake_cell = dynamic_cast<FakeBathCell*>(p_cell);
                if (p_fake_cell == NULL)
                {
                    if (inits.empty())
                    {
                        inits = p_cell->GetInitialConditions();
                    }
                    std::vector<double> state = p_cell->GetStateVariables();
                    TS_ASSERT_EQUALS(state.size(), inits.size());
                    for (unsigned j=0; j<state.size(); j++)
                    {
                        TS_ASSERT_DELTA(state[j], inits[j], 1e-10);
                    }
                }
            }
        }

        // Save it to an archive for the current number of processes
        CardiacSimulationArchiver<Problem>::Save(*p_problem, rNewArchiveDir);

        // Compare with the global archive from the previous test
        OutputFileHandler handler(rNewArchiveDir, false);
        std::string ref_archive = handler.GetChasteTestOutputDirectory() + rRefArchiveDir + "/archive.arch";
        std::string my_archive = handler.GetOutputDirectoryFullPath() + "archive.arch";
        TS_ASSERT(FileComparison(ref_archive, my_archive).CompareFiles());
        //THIS WON'T WORK WITH DIFFERENT VERSIONS OF BOOST:
#ifndef BOOST_VERSION
        TS_FAIL("This test needs to know the version of Boost with which it was compiled.");
#endif
#if BOOST_VERSION >= 103400 && BOOST_VERSION < 103500
        if (PetscTools::IsSequential())
        {
            /*
             * This would differ because we get extra copies of the ODE solver and intracellular stimulus objects
             * (one per original process which had them).
             * If this fails you probably just need to copy a new reference_0_archive file from "my_archive.0",
             * (running with Boost 1.34 - this code isn't executed otherwise!)
             * but do check that's the case!
             *
             * (Assuming 1.34)
scons test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp /tmp/$USER/testoutput/TestLoadAsSequential/archive.arch.0 ./heart/test/data/checkpoint_migration/reference_0_archive
cp /tmp/$USER/testoutput/TestLoadAsSequentialWithBath/archive.arch.0 ./heart/test/data/checkpoint_migration_with_bath/reference_0_archive
cp /tmp/$USER/testoutput/TestBcsOnNonMasterOnly/archive.arch.0 ./heart/test/data/checkpoint_migration_bcs_on_non_master_only/reference_0_archive
cp /tmp/$USER/testoutput/TestMigrateAfterSolve/archive.arch.0 ./heart/test/data/checkpoint_migration_after_solve/reference_0_archive
scons test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp /tmp/$USER/testoutput/TestLoadAsSequentialWithBathAndDistributedMesh/archive.arch.0 ./heart/test/data/checkpoint_migration_with_bath_and_distributed_mesh/reference_0_archive
             */
            //Not a collective call in general, but we only do this one in serial
            FileComparison comp(rArchiveDirectory.GetAbsolutePath() + "reference_0_archive", my_archive + ".0");
            TS_ASSERT(comp.CompareFiles());
        }
#endif

        // Return the problem for further tests
        return p_problem;
    }

    /**
     * Write results to a subfolder of the original output dir, to avoid over-writing
     * original results & archives.
     *
     * Also change the end time if we're resuming a simulation, and copy the previous
     * results h5 file into the new subfolder in that case, so that it can be extended.
     */
    void SetOutputDirAndEndTime(const FileFinder& rDirToCopyFrom,
                                const std::string& rRefArchiveDir, // relative to CHASTE_TEST_OUTPUT
                                const std::string& rSubDir,
                                double endTime)
    {
        // Change output directory to avoid over-writing original results & archives
        std::string output_dir = rRefArchiveDir + "/" + rSubDir;
        HeartConfig::Instance()->SetOutputDirectory(output_dir);
        if (endTime > 0.0)
        {
            HeartConfig::Instance()->SetSimulationDuration(endTime);
            // Also need to copy the old results into the new output folder
            OutputFileHandler handler2(output_dir, false); // must be collective

            FileFinder original_file(rDirToCopyFrom.GetAbsolutePath() + "simulation.h5", RelativeTo::Absolute);
            handler2.CopyFileTo(original_file); // must be collective.
        }
    }

    template<class Problem>
    void DoSimulationsAfterMigrationAndCompareResults(Problem* pProblem, // the loaded problem
                                                      const FileFinder& rArchiveDirectory, // where the archive was stored
                                                      const std::string& rRefArchiveDir, // relative to CHASTE_TEST_OUTPUT
                                                      const std::string& rNewArchiveDir, // relative to CHASTE_TEST_OUTPUT
                                                      unsigned numVars, // mono or bidomain?
                                                      double endTime=0.0) // if >0, we're resuming a simulation
    {
        // Simulate this problem
        SetOutputDirAndEndTime(rArchiveDirectory, rRefArchiveDir, "mig1", endTime);
        pProblem->Solve();
        // Replicate the results vector
        ReplicatableVector migrated_soln_1(pProblem->GetSolution());
        // and destroy the problem, so we don't get confusion from 2 problems at the same time
        delete pProblem;

        // Compare the results with simulating the archive from the previous test
        Problem* p_orig_problem = CardiacSimulationArchiver<Problem>::Load(rRefArchiveDir);
        FileFinder ref_archive_dir(rRefArchiveDir, RelativeTo::ChasteTestOutput);
        SetOutputDirAndEndTime(ref_archive_dir, rRefArchiveDir, "orig", endTime);
        p_orig_problem->Solve();
        ReplicatableVector orig_soln(p_orig_problem->GetSolution());
        TS_ASSERT_EQUALS(migrated_soln_1.GetSize(), orig_soln.GetSize());
        for (unsigned i=0; i<migrated_soln_1.GetSize(); i++)
        {
            double tol;
            if (PetscTools::GetNumProcs() < 4)
            {
                tol = ABS_TOL;
            }
            else
            {
                // This is horrible, but it seems that the change in partitioning when using a constructed
                // mesh leads to slight differences in PETSc's linear solve (we think)...
                tol = 300*ABS_TOL;
            }
            TS_ASSERT_DELTA(orig_soln[i], migrated_soln_1[i], tol);
        }
        delete p_orig_problem;

        // Now try loading the migrated simulation that we saved above
        Problem* p_problem = CardiacSimulationArchiver<Problem>::Load(rNewArchiveDir);
        SetOutputDirAndEndTime(rArchiveDirectory, rRefArchiveDir, "mig2", endTime);
        // Change end time?
        if (endTime > 0.0)
        {
            HeartConfig::Instance()->SetSimulationDuration(endTime);
        }
        p_problem->Solve();
        // and again compare the results
        ReplicatableVector migrated_soln_2(p_problem->GetSolution());
        TS_ASSERT_EQUALS(migrated_soln_1.GetSize(), migrated_soln_2.GetSize());
        for (unsigned i=0; i<migrated_soln_1.GetSize(); i++)
        {
            TS_ASSERT_DELTA(migrated_soln_2[i], migrated_soln_1[i], ABS_TOL);
        }
        delete p_problem;
    }

public:
    void TestMigrationExceptions()
    {
        FileFinder archive_dir("heart/test/data/checkpoint_migration_exception/", RelativeTo::ChasteSourceRoot);
        TS_ASSERT_THROWS_CONTAINS(CardiacSimulationArchiver<BidomainProblem<3> >::Load(archive_dir),
                                  "Cannot load main archive file: ");
        TS_ASSERT_THROWS_CONTAINS(CardiacSimulationArchiver<BidomainProblem<3> >::Load("absent_directory"),
                                  "Checkpoint directory does not exist: ");
    }

    /**
     * Run this in parallel (build=_3) to create the archive for TestLoadAsSequential.
     *
scons build=GccOpt_hostconfig,boost=1-34_3 test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp /tmp/$USER/testoutput/TestCreateArchiveForLoadAsSequential/?* ./heart/test/data/checkpoint_migration/
     *
     * Sets up a simulation and archives it without solving at all.
     *
     * When running sequentially, this creates an archive we can compare with
     * that produced by the next test.
     *
     * Generates a 3d cube mesh with 125 nodes, corners at (0,0,0) and (1,1,1)
     * with nodal spacing of 0.2cm.
     */
    void TestCreateArchiveForLoadAsSequential()
    {
        std::string directory = "TestCreateArchiveForLoadAsSequential";
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSlabDimensions(1, 1, 1, 0.2);
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("simulation");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(ABS_TOL);

        PlaneStimulusCellFactory<CellFoxModel2002FromCellMLBackwardEuler, 3> cell_factory(-80000.0, 1.0);
        BidomainProblem<3> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();

        CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, directory);

        // This is quite a nasty exception to cover, as you need to be able to write to all files
        // in the directory except archive.info.  It also has to be done in parallel too, to
        // cover the exception replication.
        OutputFileHandler handler(directory, false);
        std::string info_path = handler.GetOutputDirectoryFullPath() + "archive.info";
        if (PetscTools::AmMaster())
        {
            chmod(info_path.c_str(), CHASTE_READONLY);
        }
        TS_ASSERT_THROWS_CONTAINS(CardiacSimulationArchiver<BidomainProblem<3> >::Save(bidomain_problem, directory, false),
                                  "Unable to open archive information file");
        if (PetscTools::AmMaster())
        {
            chmod(info_path.c_str(), CHASTE_READ_WRITE);
        }
    }

    /**
     * #1159 - the first part of migrating a checkpoint to a different number of processes.
     *
     * The original LoadAsSequential and LoadFromSequential methods are now gone, since
     * we can do everything with a single Load method.  But the tests are still named after
     * the original methods.
     */
    void TestLoadAsSequential()
    {
        FileFinder source_directory("heart/test/data/checkpoint_migration/", RelativeTo::ChasteSourceRoot);
        std::string new_directory = "TestLoadAsSequential";
        std::string ref_archive_dir = "TestCreateArchiveForLoadAsSequential";

        // Do the migration
        const unsigned num_cells = 216u;
        BidomainProblem<3>* p_problem = DoMigrateAndBasicTests<BidomainProblem<3>,3>(source_directory, ref_archive_dir, new_directory, num_cells, true);

        // All cells at x=0 should have a SimpleStimulus(-80000, 1).
        DistributedVectorFactory* p_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
        for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
        {
            AbstractCardiacCellInterface* p_cell = p_problem->GetTissue()->GetCardiacCell(i);
            double x = p_problem->rGetMesh().GetNode(i)->GetPoint()[0];

            if (x*x < 1e-10)
            {
                // Stim exists
                TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), -80000.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(1.0), -80000.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(1.001), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(-1e-10), 0.0, 1e-10);
            }
            else
            {
                // No stim
                TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(1.0), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(1.001), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(-1e-10), 0.0, 1e-10);
            }
        }

        // Test bccs - none defined in this problem
        TS_ASSERT(! p_problem->mpDefaultBoundaryConditionsContainer);
        TS_ASSERT(! p_problem->mpBoundaryConditionsContainer);

        DoSimulationsAfterMigrationAndCompareResults(p_problem, source_directory, ref_archive_dir, new_directory, 2);
    }


    /**
     * Run this in parallel (build=_3) to create the archive for TestLoadAsSequentialWithBath.
     *
scons build=GccOpt_hostconfig,boost=1-34_3 test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp  /tmp/$USER/testoutput/TestCreateArchiveForLoadAsSequentialWithBath/?* ./heart/test/data/checkpoint_migration_with_bath/
cp  /tmp/$USER/testoutput/TestCreateArchiveForLoadAsSequentialWithBathAndDistributedMesh/?* ./heart/test/data/checkpoint_migration_with_bath_and_distributed_mesh/
     */
    void TestCreateArchiveForLoadAsSequentialWithBath()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        HeartConfig::Instance()->SetOutputFilenamePrefix("simulation");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(ABS_TOL);

        // boundary flux for Phi_e. -10e3 is under threshold, -14e3 crashes the cell model
        HeartConfig::Instance()->SetOdeTimeStep(0.001);  // ms
        double boundary_flux = -11.0e3;
        double duration = 1.9; // of the stimulus, in ms

        HeartConfig::Instance()->SetElectrodeParameters(false/*don't ground*/, 0/*x*/,
                                                        boundary_flux, 0.0, duration);

        {
            std::string directory = "TestCreateArchiveForLoadAsSequentialWithBath";
            HeartConfig::Instance()->SetOutputDirectory(directory);

            TetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<TetrahedralMesh<2,2> >(
                "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);
            ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;

            BidomainProblem<2> bidomain_problem( &cell_factory, true );

            bidomain_problem.SetMesh(p_mesh);
            bidomain_problem.Initialise();

            CardiacSimulationArchiver<BidomainProblem<2> >::Save(bidomain_problem, directory, false);

            delete p_mesh;
        }

        // And now with a distributed mesh, for coverage
        {
            std::string directory = "TestCreateArchiveForLoadAsSequentialWithBathAndDistributedMesh";
            HeartConfig::Instance()->SetOutputDirectory(directory);

            DistributedTetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<DistributedTetrahedralMesh<2,2> >(
                "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);
            ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;

            BidomainProblem<2> bidomain_problem( &cell_factory, true );

            bidomain_problem.SetMesh(p_mesh);
            bidomain_problem.Initialise();

            CardiacSimulationArchiver<BidomainProblem<2> >::Save(bidomain_problem, directory, false);

            delete p_mesh;
        }
    }

    void TestLoadAsSequentialWithBath()
    {
        FileFinder source_directory("heart/test/data/checkpoint_migration_with_bath/", RelativeTo::ChasteSourceRoot);
        std::string new_archive_dir = "TestLoadAsSequentialWithBath";
        std::string ref_archive_dir = "TestCreateArchiveForLoadAsSequentialWithBath";

        BidomainProblem<2>* p_problem;
        // Do the migration to sequential
        const unsigned num_cells = 221u;
        p_problem = DoMigrateAndBasicTests<BidomainProblem<2>,2>(source_directory, ref_archive_dir, new_archive_dir, num_cells, false);

        // All cells should have no stimulus.
        DistributedVectorFactory* p_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
        for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
        {
            AbstractCardiacCellInterface* p_cell = p_problem->GetTissue()->GetCardiacCell(i);
            AbstractStimulusFunction* p_stim = p_cell->GetStimulusFunction().get();
            ZeroStimulus* p_zero_stim = dynamic_cast<ZeroStimulus*>(p_stim);
            TS_ASSERT(p_zero_stim != NULL);
            TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), 0.0, 1e-10);
            TS_ASSERT_DELTA(p_cell->GetStimulus(1.0), 0.0, 1e-10);
        }

        // Test bccs
        TS_ASSERT( ! p_problem->mpDefaultBoundaryConditionsContainer);
        TS_ASSERT( ! p_problem->mpBoundaryConditionsContainer);
        // Problem's BCC will be created first time BidomainProblem::Solve() is called,
        // so check the Electrodes' BCC directly
        boost::shared_ptr<BoundaryConditionsContainer<2,2,2> > p_bcc = p_problem->mpElectrodes->GetBoundaryConditionsContainer();
        // We have neumann boundary conditions from the electrodes
        TS_ASSERT(p_bcc->AnyNonZeroNeumannConditions());
        for (BoundaryConditionsContainer<2,2,2>::NeumannMapIterator it = p_bcc->BeginNeumann();
             it != p_bcc->EndNeumann();
             ++it)
        {
            ChastePoint<2> centroid(it->first->CalculateCentroid());
            // Negative flux at x=0
            if (fabs(centroid[0] - 0.0) < 1e-6)
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(it->first, centroid, 1), -11e3, 1e-10);
            }
            // Positive flux at x=0.1
            else if (fabs(centroid[0] - 0.1) < 1e-6)
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(it->first, centroid, 1), +11e3, 1e-10);
            }
            else
            {
                TS_FAIL("Unexpected Neumann BC found.");
            }
        }
        // There are no dirichlet BCs.
        TS_ASSERT( ! p_bcc->HasDirichletBoundaryConditions());
        // Now check that all relevant boundary elements have neumann conditions
        for (AbstractTetrahedralMesh<2,2>::BoundaryElementIterator iter = p_problem->rGetMesh().GetBoundaryElementIteratorBegin();
             iter != p_problem->rGetMesh().GetBoundaryElementIteratorEnd();
             iter++)
        {
            ChastePoint<2> centroid((*iter)->CalculateCentroid());
            if (fabs(centroid[0] - 0.0) < 1e-6)
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter,centroid,1), -11e3, 1e-10);
            }
            else if (fabs(centroid[0] - 0.1) < 1e-6)
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter,centroid,1), +11e3, 1e-10);
            }
            else
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter,centroid,1), 0, 1e-10);
            }
            // No neumann stimulus applied to V
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter,centroid,0), 0, 1e-10);
        }

        DoSimulationsAfterMigrationAndCompareResults(p_problem, source_directory, ref_archive_dir, new_archive_dir, 2);

        // And now a shorter test with a distributed mesh, for coverage
        source_directory = FileFinder("heart/test/data/checkpoint_migration_with_bath_and_distributed_mesh/", RelativeTo::ChasteSourceRoot);
        new_archive_dir = "TestLoadAsSequentialWithBathAndDistributedMesh";
        ref_archive_dir = "TestCreateArchiveForLoadAsSequentialWithBathAndDistributedMesh";

        BidomainProblem<2>* p_problem2;
        // Do the migration to sequential
        p_problem2 = DoMigrateAndBasicTests<BidomainProblem<2>,2>(source_directory, ref_archive_dir, new_archive_dir, num_cells, true);
        delete p_problem2;
    }

private:
    template<class Problem, unsigned DIM>
    Problem* DoMigrateFromSequentialAndBasicTests(
            const FileFinder& rArchiveDirectory,
            const std::string& rRefArchiveDir,
            const std::string& rNewArchiveDir,
            const unsigned totalNumCells,
            bool isDistributedMesh,
            bool isConstructedMesh,
            double currentTime=0.0)
    {
        // Do the migration
        Problem* p_problem = CardiacSimulationArchiver<Problem>::Load(rArchiveDirectory);

        // Some basic tests that we have the right data
        TS_ASSERT_EQUALS(p_problem->mMeshFilename, "");
        TS_ASSERT_EQUALS(p_problem->mPrintOutput, true);
        TS_ASSERT_EQUALS(p_problem->mNodesToOutput.size(), 0u);
        TS_ASSERT_EQUALS(p_problem->mCurrentTime, currentTime);
        DistributedVectorFactory* p_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
        TS_ASSERT(p_factory->GetOriginalFactory());
        TS_ASSERT_EQUALS(p_factory->GetOriginalFactory()->GetProblemSize(), totalNumCells);
        TS_ASSERT_EQUALS(p_factory->GetOriginalFactory()->GetNumProcs(), 1u);
        TS_ASSERT_EQUALS(p_problem->GetTissue()->rGetCellsDistributed().size(), p_factory->GetLocalOwnership());
        TS_ASSERT_EQUALS(p_problem->rGetMesh().GetNumAllNodes(), totalNumCells);
        TS_ASSERT_EQUALS(p_problem->rGetMesh().GetNumNodes(), totalNumCells);
        TS_ASSERT_EQUALS(&(p_problem->rGetMesh()), p_problem->GetTissue()->pGetMesh());
        // Check the mesh is/isn't the parallel variety
        const DistributedTetrahedralMesh<DIM,DIM>* p_dist_mesh = dynamic_cast<const DistributedTetrahedralMesh<DIM,DIM>*>(p_problem->GetTissue()->pGetMesh());
        if (isDistributedMesh)
        {
            TS_ASSERT(p_dist_mesh != NULL);
        }
        else
        {
            TS_ASSERT(p_dist_mesh == NULL);
        }

        if (currentTime == 0.0)
        {
            // All cells should be at initial conditions.
            if (p_factory->GetHigh() > p_factory->GetLow())
            {
                std::vector<double> inits;
                for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
                {
                    AbstractCardiacCellInterface* p_cell = p_problem->GetTissue()->GetCardiacCell(i);
                    FakeBathCell* p_fake_cell = dynamic_cast<FakeBathCell*>(p_cell);
                    if (p_fake_cell == NULL)
                    {
                        if (inits.empty())
                        {
                            inits = static_cast<AbstractCardiacCell*>(p_cell)->GetInitialConditions();
                        }
                        std::vector<double> state = p_cell->GetStdVecStateVariables();
                        TS_ASSERT_EQUALS(state.size(), inits.size());
                        for (unsigned j=0; j<state.size(); j++)
                        {
                            TS_ASSERT_DELTA(state[j], inits[j], 1e-10);
                        }
                    }
                }
            }
        }

        // Save it to a normal archive
        CardiacSimulationArchiver<Problem>::Save(*p_problem, rNewArchiveDir);

        // Compare with the archive from the previous test
        OutputFileHandler handler(rNewArchiveDir, false);
        std::string ref_archive = handler.GetChasteTestOutputDirectory() + rRefArchiveDir + "/archive.arch";
        std::string my_archive = handler.GetOutputDirectoryFullPath() + "archive.arch";
        ///\todo When #1199 is done, re-instate this if
//        if (PetscTools::IsParallel())
//        {
            // If we're loading on 1 process, then the archives will differ in one digit:
            // DistributedTetrahedralMesh::mMetisPartitioning !
            FileComparison comparer(ref_archive, my_archive, false);
            TS_ASSERT(comparer.CompareFiles());
//        }
        if (!isConstructedMesh || PetscTools::GetNumProcs() < 4)
        {
            // This works because the original archive was created by a single process,
            // but only if the mesh was read from disk: a constructed mesh doesn't use
            // a DUMB partition, so the archives differ if the number of processes doesn't
            // divide the number of rows in the mesh.
            std::stringstream proc_id;
            proc_id << PetscTools::GetMyRank();
            std::string suffix = "." + proc_id.str();
            // We can't do a straight FileComparison:
            // because we may have done a little simulation already, so there may be differences
            // between serial and parallel results, below the numerical solver tolerance.

            std::string file1 = ref_archive + suffix;
            std::string file2 = my_archive + suffix;
            NumericFileComparison same_data(file1, file2);
            TS_ASSERT(same_data.CompareFiles(ABS_TOL, 1)); //Tolerance of comparison.  Ignore first line (first 1 lines).
        }

        return p_problem;
    }

public:
    /**
     * Run this in sequential to create the archive for TestLoadFromSequential.
scons build=GccOpt_hostconfig,boost=1-34 test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp /tmp/$USER/testoutput/TestCreateArchiveForLoadFromSequential/?* ./heart/test/data/checkpoint_migration_from_seq/
     *
     * Sets up a simulation and archives it without solving at all.
     *
     * Generates a 3d cube mesh with 125 nodes, corners at (0,0,0) and (1,1,1)
     * with nodal spacing of 0.2cm.
     */
    void TestCreateArchiveForLoadFromSequential()
    {
        std::string directory = "TestCreateArchiveForLoadFromSequential";
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSlabDimensions(1, 1, 1, 0.2);
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("simulation");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(ABS_TOL);

        PlaneStimulusCellFactory<CellFaberRudy2000FromCellML, 3> cell_factory(-25500.0, 2.0);
        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        CardiacSimulationArchiver<MonodomainProblem<3> >::Save(monodomain_problem, directory);
    }

    /**
     * #1159 - the second part of migrating a checkpoint to a different number of processes.
     */
    void TestLoadFromSequential()
    {
        FileFinder source_directory("heart/test/data/checkpoint_migration_from_seq/", RelativeTo::ChasteSourceRoot);
        std::string new_archive_dir = "TestLoadFromSequential";
        std::string ref_archive_dir = "TestCreateArchiveForLoadFromSequential";

        // Loading from a sequential archive should work just as well running sequentially as in parallel -
        // if running sequentially it's just the same as a normal load.
        const unsigned num_cells = 216u;
        MonodomainProblem<3>* p_problem = DoMigrateFromSequentialAndBasicTests<MonodomainProblem<3>,3>(source_directory, ref_archive_dir, new_archive_dir, num_cells, true, true);

        // All cells at x=0 should have a SimpleStimulus(-25500, 2).
        DistributedVectorFactory* p_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
        for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
        {
            AbstractCardiacCellInterface* p_cell = p_problem->GetTissue()->GetCardiacCell(i);
            double x = p_problem->rGetMesh().GetNode(i)->GetPoint()[0];

            if (x*x < 1e-10)
            {
                // Stim exists
                TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), -25500.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(2.0), -25500.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(2.001), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(-1e-10), 0.0, 1e-10);
            }
            else
            {
                // No stim
                TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(2.0), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(2.001), 0.0, 1e-10);
                TS_ASSERT_DELTA(p_cell->GetStimulus(-1e-10), 0.0, 1e-10);
            }
        }

        // Test bccs - none defined in this problem
        TS_ASSERT(! p_problem->mpDefaultBoundaryConditionsContainer);
        TS_ASSERT(! p_problem->mpBoundaryConditionsContainer);

        DoSimulationsAfterMigrationAndCompareResults(p_problem, source_directory, ref_archive_dir, new_archive_dir, 1);
    }

    /**
     * Run this in sequential to create the archive for TestLoadFromSequentialWithBath.
scons build=GccOpt_hostconfig,boost=1-34  test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp /tmp/$USER/testoutput/TestCreateArchiveForLoadFromSequentialWithBath/?* ./heart/test/data/checkpoint_migration_from_seq_with_bath/
     *
     * Sets up a simulation and archives it, solving for one PDE step first to set up default BCs.
     */
    void TestCreateArchiveForLoadFromSequentialWithBath()
    {
        std::string directory = "TestCreateArchiveForLoadFromSequentialWithBath";
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSimulationDuration(0.2);
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("simulation");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(ABS_TOL);

        DistributedTetrahedralMesh<2,2>* p_mesh = Load2dMeshAndSetCircularTissue<DistributedTetrahedralMesh<2,2> >(
            "mesh/test/data/2D_0_to_1mm_400_elements", 0.05, 0.05, 0.02);
        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;

        // boundary flux for Phi_e. -10e3 is under threshold, -14e3 crashes the cell model
        HeartConfig::Instance()->SetOdeTimeStep(0.001);  // ms
        double boundary_flux = -11.0e3;
        double duration = 1.9; // of the stimulus, in ms
        HeartConfig::Instance()->SetElectrodeParameters(true/*do ground*/, 0/*x*/,
                                                        boundary_flux, 0.0, duration);

        BidomainProblem<2> bidomain_problem( &cell_factory, true );
        bidomain_problem.SetMesh(p_mesh);

        // We solve for a small period of time so BCCs are created (next test wants to check it)
        HeartConfig::Instance()->SetSimulationDuration(0.1);
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // We increase the simulation time so next test finds something left to simulate in the archive
        HeartConfig::Instance()->SetSimulationDuration(0.2);

        CardiacSimulationArchiver<BidomainProblem<2> >::Save(bidomain_problem, directory, false);

        delete p_mesh;
    }

    void TestLoadFromSequentialWithBath()
    {
        FileFinder source_directory("heart/test/data/checkpoint_migration_from_seq_with_bath/", RelativeTo::ChasteSourceRoot);
        std::string new_archive_dir = "TestLoadFromSequentialWithBath";
        std::string ref_archive_dir = "TestCreateArchiveForLoadFromSequentialWithBath";

        // Loading from a sequential archive should work just as well running sequentially as in parallel -
        // if running sequentially it's essentially just the same as a normal load.
        const unsigned num_cells = 221u;
        BidomainProblem<2>* p_problem = DoMigrateFromSequentialAndBasicTests<BidomainProblem<2>,2>(source_directory, ref_archive_dir, new_archive_dir, num_cells, true, false, 0.1);

        // All cells should have a ZeroStimulus.
        DistributedVectorFactory* p_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
        for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
        {
            AbstractCardiacCellInterface* p_cell = p_problem->GetTissue()->GetCardiacCell(i);
            AbstractStimulusFunction* p_stim = p_cell->GetStimulusFunction().get();
            ZeroStimulus* p_zero_stim = dynamic_cast<ZeroStimulus*>(p_stim);
            TS_ASSERT(p_zero_stim != NULL);
            TS_ASSERT_DELTA(p_cell->GetStimulus(0.0), 0.0, 1e-10);
            TS_ASSERT_DELTA(p_cell->GetStimulus(1.0), 0.0, 1e-10);
        }

        // Test bccs
        TS_ASSERT(p_problem->mpDefaultBoundaryConditionsContainer); /// \todo: see todo in BidomainProblem.cpp:344
        TS_ASSERT(p_problem->mpBoundaryConditionsContainer);
        boost::shared_ptr<BoundaryConditionsContainer<2,2,2> > p_bcc = p_problem->mpBoundaryConditionsContainer;
        // We have neumann and dirichlet boundary conditions from the electrodes (at least on some processes)
        if (PetscTools::IsSequential())
        {
            TS_ASSERT(p_bcc->AnyNonZeroNeumannConditions());
        }
        TS_ASSERT(p_bcc->HasDirichletBoundaryConditions());
        for (BoundaryConditionsContainer<2,2,2>::NeumannMapIterator it = p_bcc->BeginNeumann();
             it != p_bcc->EndNeumann();
             ++it)
        {
            ChastePoint<2> centroid(it->first->CalculateCentroid());
            // Negative flux at x=0
            if (fabs(centroid[0] - 0.0) < 1e-6)
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(it->first, centroid, 1), -11e3, 1e-10);
            }
            else
            {
                TS_FAIL("Unexpected Neumann BC found.");
            }
        }
        // Now check that all relevant boundary elements have neumann conditions
        for (DistributedTetrahedralMesh<2,2>::BoundaryElementIterator iter = p_problem->rGetMesh().GetBoundaryElementIteratorBegin();
             iter != p_problem->rGetMesh().GetBoundaryElementIteratorEnd();
             iter++)
        {
            ChastePoint<2> centroid((*iter)->CalculateCentroid());
            if (fabs(centroid[0] - 0.0) < 1e-6)
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter,centroid,1), -11e3, 1e-10);
            }
            else
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter,centroid,1), 0, 1e-10);
            }
            // No neumann stimulus applied to V
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter,centroid,0), 0, 1e-10);
        }
        // Check dirichlet conditions exist only at x=0.1
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter=p_problem->rGetMesh().GetNodeIteratorBegin();
             iter != p_problem->rGetMesh().GetNodeIteratorEnd();
             ++iter)
        {
            Node<2>* p_node = &(*iter);
            double x = p_node->rGetLocation()[0];
            if (fabs(x - 0.1) < 1e-6)
            {
                TS_ASSERT(p_bcc->HasDirichletBoundaryCondition(p_node, 1));
                TS_ASSERT_DELTA(p_bcc->GetDirichletBCValue(p_node, 1), 0, 1e-10); // grounded
            }
            else
            {
                TS_ASSERT(!p_bcc->HasDirichletBoundaryCondition(p_node, 1));
            }
            // No dirichlet condition on V
            TS_ASSERT(!p_bcc->HasDirichletBoundaryCondition(p_node, 0));
        }

        DoSimulationsAfterMigrationAndCompareResults(p_problem, source_directory, ref_archive_dir, new_archive_dir, 2, 0.2);
    }

    /*
     * The following tests that there are no problems when running and resuming a simulation with postprocessing.
     */
    void TestSolveAndResumeWithPostprocessing()
    {
        setUp();
        HeartConfig::Instance()->SetOutputDirectory("SolveAndResumeWithPostprocessing");
        const std::string archive_dir = "SolveAndResumeWithPostprocessing/saved_simulation";
        // Solve with some postprocessing switched on
        {
            HeartConfig::Instance()->SetSimulationDuration(400.0);
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 1.0);
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");

            std::vector<double> upstroke_maps;
            upstroke_maps.push_back(0.0);
            HeartConfig::Instance()->SetUpstrokeTimeMaps(upstroke_maps);

            std::vector<std::pair<double,double> > apd_maps;
            apd_maps.push_back(std::pair<double, double>(90.0, 0.0)); // APD90, 0 mV threshold
            HeartConfig::Instance()->SetApdMaps(apd_maps);

            PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory(-5.0e5, 1.0);
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();
            monodomain_problem.Solve();

            CardiacSimulationArchiver<MonodomainProblem<1> >::Save(monodomain_problem, archive_dir, true); // Clear just in case
        }
        {
            MonodomainProblem<1>* p_problem = CardiacSimulationArchiver<MonodomainProblem<1> >::Load(archive_dir);
            HeartConfig::Instance()->SetSimulationDuration(500.0);
            // All other settings inherited, including output folder and maps, so the writer will try and overwrite.
            p_problem->Solve();
            delete p_problem;
        }
        {
            // Do it all again but this time in one go, for comparison.
            HeartConfig::Instance()->SetOutputDirectory("SolveWithPostprocessing");
            HeartConfig::Instance()->SetSimulationDuration(500.0);

            PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory(-5.0e5, 1.0);
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();
            monodomain_problem.Solve();
        }
        // Compare
        TS_ASSERT(CompareFilesViaHdf5DataReader("SolveAndResumeWithPostprocessing", "results", true,
                                                "SolveWithPostprocessing", "results", true));
    }

    /**
     * The next pair of tests cover the case where the master process (that creates the archive)
     * has no boundary conditions, but at least one other process does.
     * We set a zero dirichlet boundary condition on the right end of a parallel 1d mesh.
     *
scons build=GccOpt_hostconfig,boost=1-34_2  test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp /tmp/$USER/testoutput/TestCreateArchiveForBcsOnNonMasterOnly/?* ./heart/test/data/checkpoint_migration_bcs_on_non_master_only/
     */
    void TestCreateArchiveForBcsOnNonMasterOnly()
    {
        std::string directory = "TestCreateArchiveForBcsOnNonMasterOnly";
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSimulationDuration(0.1); //ms
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("simulation");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(ABS_TOL);

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_10_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);
        // Set the whole mesh to be bath
        DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
        for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
        {
            mesh.GetNode(i)->SetRegion(HeartRegionCode::GetValidBathId());
        }
        for (DistributedTetrahedralMesh<1,1>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            it->SetAttribute(HeartRegionCode::GetValidBathId());
        }
        mesh.SetMeshHasChangedSinceLoading();

        ZeroStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;

        BidomainProblem<1> bidomain_problem( &cell_factory, true );
        bidomain_problem.SetMesh(&mesh);

        // Boundary condition - fix phi=0 on RHS
        unsigned rhs_index = 10u;
        if (p_factory->IsGlobalIndexLocal(rhs_index))
        {
            boost::shared_ptr<BoundaryConditionsContainer<1,1,2> > p_bcc(new BoundaryConditionsContainer<1,1,2>);
            ConstBoundaryCondition<1>* p_zero_stim = new ConstBoundaryCondition<1>(0.0);
            Node<1>* p_node = mesh.GetNode(rhs_index);
            p_bcc->AddDirichletBoundaryCondition(p_node, p_zero_stim, 1u);
            bidomain_problem.SetBoundaryConditionsContainer(p_bcc);
        }

        bidomain_problem.Initialise();

        // And Save
        CardiacSimulationArchiver<BidomainProblem<1> >::Save(bidomain_problem, directory);
    }

    void TestBcsOnNonMasterOnly()
    {
        FileFinder source_directory("heart/test/data/checkpoint_migration_bcs_on_non_master_only/", RelativeTo::ChasteSourceRoot);
        std::string new_archive_dir = "TestBcsOnNonMasterOnly";
        std::string ref_archive_dir = "TestCreateArchiveForBcsOnNonMasterOnly";

        const unsigned num_nodes = 11u;
        BidomainProblem<1>* p_problem = DoMigrateAndBasicTests<BidomainProblem<1>,1>(source_directory, ref_archive_dir, new_archive_dir, num_nodes, true);

        // Test the bccs - zero dirichlet on RHS only
        TS_ASSERT(! p_problem->mpDefaultBoundaryConditionsContainer);
        boost::shared_ptr<BoundaryConditionsContainer<1,1,2> > p_bcc = p_problem->mpBoundaryConditionsContainer;
        // If running on 2 procs we now don't migrate, so process 0 doesn't have a BCC.
        if (PetscTools::GetNumProcs() != 2)
        {
            TS_ASSERT(p_bcc);
            TS_ASSERT(! p_bcc->AnyNonZeroNeumannConditions());
            TS_ASSERT( p_bcc->HasDirichletBoundaryConditions());
            DistributedVectorFactory* p_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
            for (unsigned i=p_factory->GetLow(); i<p_factory->GetHigh(); i++)
            {
                Node<1>* p_node = p_problem->rGetMesh().GetNode(i);
                TS_ASSERT(! p_bcc->HasDirichletBoundaryCondition(p_node, 0u));
                if (i == num_nodes-1)
                {
                    TS_ASSERT(p_bcc->HasDirichletBoundaryCondition(p_node, 1u));
                    TS_ASSERT_EQUALS(p_bcc->GetDirichletBCValue(p_node, 1u), 0u);
                }
                else
                {
                    TS_ASSERT(! p_bcc->HasDirichletBoundaryCondition(p_node, 1u));
                }
            }
        }

        DoSimulationsAfterMigrationAndCompareResults(p_problem, source_directory, ref_archive_dir, new_archive_dir, 2);
    }

    /**
     * The next pair of tests covers the case where we migrate after a solve.
     * This means we will have default boundary conditions included in the archive.
     * Saving the simulation in parallel, and using a parallel mesh, means we also cover the
     * direct call of mpBoundaryConditionsContainer->MergeFromArchive in LoadExtraArchive.
     *
     * Run this in parallel (build=_2) to create the archive for TestMigrateAfterSolve.
scons build=GccOpt_hostconfig,boost=1-34_2 test_suite=heart/test/TestCardiacSimulationArchiver.hpp
cp /tmp/$USER/testoutput/TestCreateArchiveForMigrateAfterSolve/archive/?* ./heart/test/data/checkpoint_migration_after_solve/
     */
    void TestCreateArchiveForMigrateAfterSolve()
    {
        std::string directory = "TestCreateArchiveForMigrateAfterSolve";
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSimulationDuration(0.1); //ms
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("simulation");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(ABS_TOL);

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Stimulate at LHS
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;

        MonodomainProblem<1> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();

        // Solve...
        monodomain_problem.Solve();

        // ...and Save
        CardiacSimulationArchiver<MonodomainProblem<1> >::Save(monodomain_problem, directory + "/archive");

        // Copy the simulation.h5 file into the archive folder, for convenience
        OutputFileHandler handler(directory, false);

        TRY_IF_MASTER(
            FileFinder source_file = handler.FindFile("simulation.h5");
            FileFinder destination(handler.GetOutputDirectoryFullPath() + "archive/", RelativeTo::Absolute);
            source_file.CopyTo(destination);
        );
    }

    void TestMigrateAfterSolve()
    {
        FileFinder source_directory("heart/test/data/checkpoint_migration_after_solve/", RelativeTo::ChasteSourceRoot);
        std::string new_archive_dir = "TestMigrateAfterSolve";
        std::string ref_archive_dir = "TestCreateArchiveForMigrateAfterSolve/archive";

        const unsigned num_nodes = 101u;
        MonodomainProblem<1>* p_problem = DoMigrateAndBasicTests<MonodomainProblem<1>,1>(source_directory, ref_archive_dir, new_archive_dir, num_nodes, true, 0.1);

        // Test the bccs - zero neumann on the boundary
        TS_ASSERT( p_problem->mpDefaultBoundaryConditionsContainer);
        boost::shared_ptr<BoundaryConditionsContainer<1,1,1> > p_bcc = p_problem->mpBoundaryConditionsContainer;
        TS_ASSERT(p_bcc);
        TS_ASSERT_EQUALS(p_bcc, p_problem->mpDefaultBoundaryConditionsContainer);
        TS_ASSERT(! p_bcc->AnyNonZeroNeumannConditions());
        TS_ASSERT(! p_bcc->HasDirichletBoundaryConditions());
        // Check all conditions are on the boundary
        for (BoundaryConditionsContainer<1,1,1>::NeumannMapIterator it = p_bcc->BeginNeumann();
             it != p_bcc->EndNeumann();
             ++it)
        {
            ChastePoint<1> centroid(it->first->CalculateCentroid());
            if (fabs(centroid[0] - 0.0) < 1e-6 || fabs(centroid[0] - 1.0) < 1e-6)
            {
                TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(it->first, centroid, 0), 0.0, 1e-10);
            }
            else
            {
                TS_FAIL("Unexpected Neumann BC found.");
            }
        }
        // Now check that all relevant boundary elements have neumann conditions
        for (DistributedTetrahedralMesh<1,1>::BoundaryElementIterator it = p_problem->rGetMesh().GetBoundaryElementIteratorBegin();
             it != p_problem->rGetMesh().GetBoundaryElementIteratorEnd();
             it++)
        {
            ChastePoint<1> centroid((*it)->CalculateCentroid());
            TS_ASSERT(p_bcc->HasNeumannBoundaryCondition(*it, 0));
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*it, centroid, 0), 0, 1e-10);
        }

        DoSimulationsAfterMigrationAndCompareResults(p_problem, source_directory, ref_archive_dir, new_archive_dir, 1, 0.2);
    }

    /*
     *  Check that we can read for a permuted mesh (or permuted archive) and then correctly record that it was permuted
     */
    void TestPermutedBidomain1D()
    {
        std::string archive_dir("ArchiveBidomainPermuted");

        std::vector<unsigned> permutation(11);
        permutation[0]=0;
        permutation[1]=10; // 1->10
        permutation[2]=2;
        permutation[3]=3;
        permutation[4]=4;
        permutation[5]=1;  // 5->1
        permutation[6]=6;
        permutation[7]=7;
        permutation[8]=8;
        permutation[9]=9;
        permutation[10]=5; // 10->5
        {
            // Run permuted simulation and archive it
            OutputFileHandler handler("PermutedBidomain1d",true);
            //In order to apply a permutation, we need to have a binary mesh
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1mm_10_elements_binary");
            mesh_reader.SetNodePermutation(permutation);

            DistributedTetrahedralMesh<1,1> permuted_mesh;
            permuted_mesh.ConstructFromMeshReader(mesh_reader);

            std::vector<unsigned> stored_perm = permuted_mesh.rGetNodePermutation();
            TS_ASSERT_EQUALS(stored_perm.size(), 11u);

            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            HeartConfig::Instance()->SetOutputDirectory("PermutedBidomain1d");
            HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

            PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
            BidomainProblem<1> bidomain_problem( &cell_factory );
            bidomain_problem.SetMesh(&permuted_mesh);
            bidomain_problem.Initialise();
            bidomain_problem.Solve();

            // Check some voltages
            ReplicatableVector solution_replicated(bidomain_problem.GetSolution());
            double atol=5e-3;
            TS_ASSERT_DELTA(solution_replicated[2], 25.3148 , atol); // V at 1 (is the original 5)
            TS_ASSERT_DELTA(solution_replicated[10],-83.3582, atol); // V at 5 (is the original 10)
            TS_ASSERT_DELTA(solution_replicated[20], 22.8117, atol); // V at 10 (is the original 1)

            CardiacSimulationArchiver<BidomainProblem<1> >::Save(bidomain_problem, archive_dir, false);
        }
        {
            // Unarchive and check that the archiver has recorded a permutation
            BidomainProblem<1> *p_bidomain_problem;
            p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<1> >::Load(archive_dir);

            // This part will pass (regardless of whether the permutation has been saved)
            ReplicatableVector solution_replicated(p_bidomain_problem->GetSolution());
            double atol=5e-3;
            TS_ASSERT_DELTA(solution_replicated[2], 25.3148 , atol); // V at 1 (is the original 5)
            TS_ASSERT_DELTA(solution_replicated[10],-83.3582, atol); // V at 5 (is the original 10)
            TS_ASSERT_DELTA(solution_replicated[20], 22.8117, atol); // V at 10 (is the original 1)

            std::vector<unsigned> archived_permutation = p_bidomain_problem->rGetMesh().rGetNodePermutation();
            TS_ASSERT_EQUALS(archived_permutation.size(), 11u);
            for (unsigned i=0; i<archived_permutation.size(); i++)  //We use "size" rather than 11 in case the previous test failed
            {
                TS_ASSERT_EQUALS(archived_permutation[i], permutation[i]);
            }
            delete p_bidomain_problem;
        }
    }
};

#endif /*TESTCARDIACSIMULATIONARCHIVER_HPP_*/
