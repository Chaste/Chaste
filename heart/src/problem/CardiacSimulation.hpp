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

#ifndef CARDIACSIMULATION_HPP_
#define CARDIACSIMULATION_HPP_

#include <vector>
#include <ctime>
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

#include "Hdf5ToMeshalyzerConverter.hpp"
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
        std::auto_ptr<Problem> p_problem;

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
                std::stringstream archive_foldername;
                archive_foldername << HeartConfig::Instance()->GetOutputDirectory() << "_" << HeartConfig::Instance()->GetSimulationDuration() << "ms";
                CardiacSimulationArchiver<Problem>::Save(*(p_problem.get()), checkpoint_dir_basename + archive_foldername.str(), false);

                // Put a copy of the partial results aside (in a subdirectory of checkpoint_dir_basename).
                OutputFileHandler checkpoint_dir_basename_handler(checkpoint_dir_basename, false);
                OutputFileHandler partial_output_dir_handler(HeartConfig::Instance()->GetOutputDirectory(), false);
                if (PetscTools::AmMaster())
                {
                    ABORT_IF_NON0(system, "cp -r " + partial_output_dir_handler.GetOutputDirectoryFullPath() + " " + checkpoint_dir_basename_handler.GetOutputDirectoryFullPath());
                }

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
     * Get the saved problem instance, if any.  Will return an empty pointer if the
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
