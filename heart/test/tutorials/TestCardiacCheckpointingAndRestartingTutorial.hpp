

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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTCARDIACCHECKPOINTINGANDRESTARTINGTUTORIAL_HPP_
#define TESTCARDIACCHECKPOINTINGANDRESTARTINGTUTORIAL_HPP_

/* HOW_TO_TAG Cardiac/Problem definition
 * Save ('checkpoint') and reload simulations
 */

/*
 * = Checkpointing and restarting cardiac simulations =
 *
 * In this tutorial we show how to save and reload cardiac simulations
 *
 * EMPTYLINE
 *
 * `CardiacSimulationArchiver` is the main class that takes care of checkpointing
 * cardiac simulations.
 */
#include <cxxtest/TestSuite.h>
#include "CardiacSimulationArchiver.hpp"
#include "BidomainProblem.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"


class TestCardiacCheckpointingAndRestartingTutorial : public CxxTest::TestSuite
{
public:
    /* First, the checkpointing test. */
    void TestCheckpointing()
    {
        /* We set up exactly the same simulation as in UserTutorials/AnotherBidomainSimulation */
        HeartConfig::Instance()->Reset();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML,2> cell_factory(-2000000);
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainCheckpointingTutorial");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements", cp::media_type::Orthotropic);

        double scale = 2;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75*scale, 0.19*scale));
        BidomainProblem<2> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        /* To save the entire simulation, use the `CardiacSimulationArchiver` class, as shown in the following.
         * Note the `BidomainProblem<2>` as the template parameter. The output directory is relative to
         * CHASTE_TEST_OUTPUT. */
        CardiacSimulationArchiver<BidomainProblem<2> >::Save(bidomain_problem, "BidomainCheckpointingTutorial/saved_simulation");
    }

    /* This is how to restart the test. */
    void TestRestarting()
    {
        /* To restart from the saved simulation directory we  use the `CardiacSimulationArchiver` class, as shown in the following.
         * Note the `BidomainProblem<2>` as the template parameter again.  The dimension (2) must match the one given in the
         * saved archive directory.
         * The output directory is again relative to CHASTE_TEST_OUTPUT. */
        BidomainProblem<2>* p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<2> >::Load("BidomainCheckpointingTutorial/saved_simulation");

        /* The simulation duration has to be amended.
         * Note that the duration is always given with respect to the origin of the first solve.
         * This means that we are running from {{{t=5 ms}}} (the end of the previous simulation) to {{{t=10 ms}}}.
         * The output files are concatenated so that they appear to be made by a single simulation running from
         * {{{t=0 ms}}} to {{{t=10 ms}}}.
         * Note: loading an archive also loads `HeartConfig` options, so `HeartConfig` calls such as this one must appear
         * ''after'' CardiacSimulationArchiver::Load().
         */
        HeartConfig::Instance()->SetSimulationDuration(10); //ms

        /* One point of checkpointing and restarting is that there may be something which we want to change
         * during the course of experiment.  Here we change the conductivity. */
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(3.0, 0.3));

        p_bidomain_problem->Solve();

        /* Note that the pointer p_bidomain_problem exists in the scope of this test and that the object
         * which was unarchived was created on the CardiacSimulationArchiver::Load() line above.  We are therefore
         * responsible for deleting the memory.
         */
        delete p_bidomain_problem;
    }
};

/* == Notes ==
 *
 * * Making a checkpoint does add a significant overhead at present, in particular because the mesh is
 *   written out to disk at each checkpoint. This is to ensure that each checkpoint directory contains everything
 *   needed to resume the simulation. The mesh written out will be in permuted form if it was
 *   partitioned for a parallel simulation.
 * * To make this process slightly more efficient, Chaste will copy the original mesh files if the mesh was loaded
 *   from disk and hasn't been modified (e.g. by permuting).  Because of this, if you modify the mesh in memory,
 *   e.g. by setting element attributes in [UserTutorials/BidomainWithBath bidomain-with-bath simulations], then
 *   you need to inform Chaste by calling `mesh.SetMeshHasChangedSinceLoading()`, so your modifications aren't lost.
 * * Meshes written in checkpoints use a binary form of the !Triangle/Tetgen mesh format. This makes checkpoints
 *   significantly smaller but will cause portability problems if checkpoints are moved between little-endian systems
 *   (e.g. x86) and big-endian systems (e.g. PowerPC).
 * * Checkpoints may be resumed on any number of processes — you are not restricted to the number on which it was
 *   saved. However, the mesh will not be re-partitioned if loaded on a different number of processes, so the parallel
 *   efficiency of the simulation may be significantly reduced in this case.
 * * Resuming a checkpoint will attempt to extend the original results HDF5 file if it exists (specified by
 *   `HeartConfig::SetOutputDirectory` and `HeartConfig::SetOutputFilenamePrefix`), so that the file contains the complete
 *   simulation results. If this file does not exist a new file will be created containing just the results from the
 *   resume time.
 *
 * * When checkpointing, the `progress_status.txt` file only reports the percentage of time to
 *   go until the ''next'' checkpoint, not until the end of the simulation.  This makes it slightly less
 *   useful; however, the presence of the checkpoint directories (`1ms`, `2ms`, etc.) provides overall
 *   progress information instead.
 * * On older versions of Boost (before 1.37) there are some restrictions on what may be successfully checkpointed.
 *   Firstly, checkpointing a simulation that uses a dynamically loaded cell model is impossible.
 *   Secondly, if you are using a "non-standard" cell model, i.e. one that is hard-coded but not one of those available
 *   via the XML configuration file, then you may get an error "terminate called after throwing an instance of
 *   'boost::archive::archive_exception' what():  unregistered class".  To resolve this you unfortunately have to
 *   include the header for the relevant cell model in `CardiacSimulationArchiver.cpp`; after the inclusion of
 *   `HeartConfigRelatedCellFactory.hpp` is a sensible location.
 * * Related to the above point, on all Boost versions if you are saving and loading a simulation in different
 *   source or test files, ensure that you have the same list of includes in both cases, or the serialization code
 *   may not know about all classes when loading, and give the "unregistered class" error.
 *
 */

#endif /*TESTCARDIACCHECKPOINTINGANDRESTARTINGTUTORIAL_HPP_*/
