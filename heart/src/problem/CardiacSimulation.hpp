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

#ifndef CARDIACSIMULATION_HPP_
#define CARDIACSIMULATION_HPP_

#include <vector>
#include <memory>

#include "UblasIncludes.hpp"

#include "AbstractCardiacProblem.hpp"
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "BidomainWithBathProblem.hpp"
#include "CardiacSimulationArchiver.hpp"
#include "PetscTools.hpp"
#include "TimeStepper.hpp"
#include "Exception.hpp"

#include "HeartConfig.hpp"
#include "HeartConfigRelatedCellFactory.hpp"
#include "HeartFileFinder.hpp"

#include "TetrahedralMesh.hpp"
#include "NonCachedTetrahedralMesh.hpp"
#include "ChastePoint.hpp"
#include "ChasteCuboid.hpp"
#include "MeshalyzerMeshWriter.hpp"
#include "TrianglesMeshWriter.hpp"

#include "OrthotropicConductivityTensors.hpp"
#include "PostProcessingWriter.hpp"

#include "OutputDirectoryFifoQueue.hpp"
#include "ExecutableSupport.hpp"

/**
 * A class which encapsulates the executable functionality.
 *
 * Takes in a chaste parameters XML file and runs the relevant simulation.
 *
 * The XML Schema, which describes what is allowed in the XML configuration file,
 * can be found at heart/src/io/ChasteParameters_2_0.xsd (for Chaste release 2.0).
 * It contains documentation describing what settings are available.  The
 * documentation of the HeartConfig class may also be of use.
 */
class CardiacSimulation
{
private:
    /**
     * Read parameters from the HeartConfig XML file.
     *
     * @param parameterFileName a string containing the chaste simulation parameters XML file name.
     */
    void ReadParametersFromFile(std::string parameterFileName);

    /**
     * Templated method which creates and runs a cardiac simulation, based on the
     * XML file passed to our constructor.
     */
    template<class Problem, unsigned SPACE_DIM>
    void CreateAndRun()
    {
        boost::shared_ptr<Problem> p_problem;

        if (HeartConfig::Instance()->IsSimulationDefined())
        {
            HeartConfigRelatedCellFactory<SPACE_DIM> cell_factory;
            p_problem.reset(new Problem(&cell_factory));

            p_problem->Initialise();
        }
        else // (HeartConfig::Instance()->IsSimulationResumed())
        {
            p_problem.reset(CardiacSimulationArchiver<Problem>::Load(HeartConfig::Instance()->GetArchivedSimulationDir()));
            // Any changes to parameters that normally only take effect at problem creation time...
            HeartConfigRelatedCellFactory<SPACE_DIM> cell_factory;
            cell_factory.SetMesh(&(p_problem->rGetMesh()));
            AbstractCardiacTissue<SPACE_DIM, SPACE_DIM>* p_tissue = p_problem->GetTissue();
            DistributedVectorFactory* p_vector_factory = p_problem->rGetMesh().GetDistributedVectorFactory();
            for (unsigned node_global_index = p_vector_factory->GetLow();
                 node_global_index < p_vector_factory->GetHigh();
                 node_global_index++)
            {
                // Overwrite any previous stimuli if new ones are defined
                cell_factory.SetCellIntracellularStimulus(p_tissue->GetCardiacCell(node_global_index), node_global_index);
                // Modify cell model parameters
                cell_factory.SetCellParameters(p_tissue->GetCardiacCell(node_global_index), node_global_index);
            }
        }

        if (HeartConfig::Instance()->GetCheckpointSimulation())
        {
            // Create the checkpoints directory and set up a fifo queue of subdirectory names
            OutputDirectoryFifoQueue directory_queue(HeartConfig::Instance()->GetOutputDirectory() + "_checkpoints/",
                                                     HeartConfig::Instance()->GetMaxCheckpointsOnDisk());

            TimeStepper checkpoint_stepper(p_problem->GetCurrentTime(), HeartConfig::Instance()->GetSimulationDuration(), HeartConfig::Instance()->GetCheckpointTimestep());
            while ( !checkpoint_stepper.IsTimeAtEnd() )
            {
                // Solve checkpoint timestep
                HeartConfig::Instance()->SetSimulationDuration(checkpoint_stepper.GetNextTime());
                p_problem->Solve();

                // Create directory that will contain archive and partial results for this checkpoint timestep.
                std::stringstream checkpoint_id;
                checkpoint_id << HeartConfig::Instance()->GetSimulationDuration() << "ms/";
                std::string checkpoint_dir_basename = directory_queue.CreateNextDir(checkpoint_id.str());

                // Archive simulation (in a subdirectory of checkpoint_dir_basename).
                char time_stamp[60]; // Write out without scientific notation in time - ticket 2861

                // Here we cope with times of up to 16 significant figures, which should be OK without floating-point crazy decimal places like 1.000000000000000000000000000000000001
                std::sprintf(time_stamp,"%0.16g",HeartConfig::Instance()->GetSimulationDuration());

                std::stringstream archive_foldername;
                archive_foldername << HeartConfig::Instance()->GetOutputDirectory() << "_" << time_stamp << "ms";

                CardiacSimulationArchiver<Problem>::Save(*(p_problem.get()), checkpoint_dir_basename + archive_foldername.str(), false);

                // Put a copy of the partial results aside (in a subdirectory of checkpoint_dir_basename).
                OutputFileHandler checkpoint_dir_basename_handler(checkpoint_dir_basename, false);
                OutputFileHandler partial_output_dir_handler(HeartConfig::Instance()->GetOutputDirectory(), false);

                TRY_IF_MASTER(
                        partial_output_dir_handler.FindFile("").CopyTo(checkpoint_dir_basename_handler.FindFile(""));
                );

                // Create an XML file to help in resuming
                CreateResumeXmlFile(checkpoint_dir_basename, archive_foldername.str());

                // Advance time stepper
                checkpoint_stepper.AdvanceOneTimeStep();
            }
        }
        else
        {
            p_problem->Solve();
        }
        if (mSaveProblemInstance)
        {
            mSavedProblem = p_problem;
        }
    }

    /**
     * Run the simulation.
     * This method basically contains switches on the problem type and space dimension,
     * and calls CreateAndRun() to do the work.
     */
    void Run();

    /**
     * Write a ResumeParameters.xml file to the checkpoint directory, to help users in resuming
     * a checkpointed simulation.  If the contents of rOutputDirectory are copied to CHASTE_TEST_OUTPUT,
     * and ResumeParameters.xml edited to specify a sensible simulation duration, then it can be used
     * as the input parameters file to resume from the given checkpoint.
     *
     * @param rOutputDirectory  the directory to put the XML file in
     * @param rArchiveDirectory  the relative path from this directory to the archive directory
     */
    void CreateResumeXmlFile(const std::string& rOutputDirectory, const std::string& rArchiveDirectory);

    /**
     * Convert a boolean to a "yes" or "no" string.
     * @param yesNo
     * @return a "yes" or "no" string.
     */
    std::string BoolToString(bool yesNo);
public:
    /**
     * Constructor
     *
     * This also runs the simulation immediately.
     *
     * @param parameterFileName  The name of the chaste parameters xml file to use to run a simulation.
     * @param writeProvenanceInfo  Whether to write provanence and machine information files.
     * @param saveProblemInstance  Whether to save a copy of the problem instance for examination by tests.
     */
    CardiacSimulation(std::string parameterFileName,
                      bool writeProvenanceInfo=false,
                      bool saveProblemInstance=false);

    /**
     * @return the saved problem instance, if any.  Will return an empty pointer if the
     * instance was not saved.
     */
    boost::shared_ptr<AbstractUntemplatedCardiacProblem> GetSavedProblem();
private:
    /** Whether to save a copy of the problem instance for examination by tests. */
    bool mSaveProblemInstance;

    /** The saved problem instance, if any. */
    boost::shared_ptr<AbstractUntemplatedCardiacProblem> mSavedProblem;
};

#endif /*CARDIACSIMULATION_HPP_*/
