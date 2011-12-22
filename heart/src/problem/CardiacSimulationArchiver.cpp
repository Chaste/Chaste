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

#include <fstream>

// Must be included before any other serialization headers
#include "CheckpointArchiveTypes.hpp"
#include "CardiacSimulationArchiver.hpp"

#include "Exception.hpp"
#include "ArchiveOpener.hpp"
#include "OutputFileHandler.hpp"
#include "ArchiveLocationInfo.hpp"
#include "DistributedVectorFactory.hpp"
#include "PetscTools.hpp"
#include "FileFinder.hpp"

/*
 * Archiving extravaganza:
 *
 * Various objects are archived through a pointer to an abstract class.  All the potential
 * concrete classes need to be included here, so they are registered with Boost.  If not,
 * Boost (<1.37) won't be able to find the archiving methods of the concrete class and will
 * throw the following exception:
 *
 *       terminate called after throwing an instance of 'boost::archive::archive_exception
 *       what():  unregistered class
 *
 * No member variable is defined to be of any of these clases, so removing them won't
 * produce any compiler error.  The exception above will occur at runtime.
 *
 * This might not be even necessary in certain cases, if the file is included implicitly by
 * another header file or by the test itself.  It's safer though.
 */
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include "BidomainWithBathProblem.hpp"
#include "HeartConfigRelatedCellFactory.hpp"


template<class PROBLEM_CLASS>
void CardiacSimulationArchiver<PROBLEM_CLASS>::Save(PROBLEM_CLASS& rSimulationToArchive,
                                                    const std::string& rDirectory,
                                                    bool clearDirectory)
{
    // Clear directory if requested (and make sure it exists)
    OutputFileHandler handler(rDirectory, clearDirectory);

    // Nest the archive writing, so the ArchiveOpener goes out of scope before
    // the method ends.
    {
        // Open the archive files
        FileFinder dir(rDirectory, RelativeTo::ChasteTestOutput);
        ArchiveOpener<boost::archive::text_oarchive, std::ofstream> archive_opener(dir, "archive.arch");
        boost::archive::text_oarchive* p_main_archive = archive_opener.GetCommonArchive();

        // And save
        PROBLEM_CLASS* const p_simulation_to_archive = &rSimulationToArchive;
        (*p_main_archive) & p_simulation_to_archive;
    }

    // Write the info file
    if (PetscTools::AmMaster())
    {
        std::string info_path = handler.GetOutputDirectoryFullPath() + "archive.info";
        std::ofstream info_file(info_path.c_str());
        if (!info_file.is_open())
        {
            // Avoid deadlock...
            PetscTools::ReplicateBool(true);
            EXCEPTION("Unable to open archive information file: " + info_path);
        }
        PetscTools::ReplicateBool(false);
        unsigned archive_version = 0; /// \todo #1026 get a real version number!
        info_file << PetscTools::GetNumProcs() << " " << archive_version;
    }
    else
    {
        bool master_threw = PetscTools::ReplicateBool(false);
        if (master_threw)
        {
            EXCEPTION("Unable to open archive information file");
        }
    }
    // Make sure everything is written before any process continues.
    PetscTools::Barrier("CardiacSimulationArchiver::Save");
}

template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::Load(const std::string& rDirectory)
{
    FileFinder directory(rDirectory, RelativeTo::ChasteTestOutput);
    return CardiacSimulationArchiver<PROBLEM_CLASS>::Migrate(directory);
}

template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::Load(const FileFinder& rDirectory)
{
    return CardiacSimulationArchiver<PROBLEM_CLASS>::Migrate(rDirectory);
}


template<class PROBLEM_CLASS>
PROBLEM_CLASS* CardiacSimulationArchiver<PROBLEM_CLASS>::Migrate(const FileFinder& rDirectory)
{
    // Check the directory exists
    std::string dir_path = rDirectory.GetAbsolutePath();
    if (!rDirectory.IsDir() || !rDirectory.Exists())
    {
        EXCEPTION("Checkpoint directory does not exist: " + dir_path);
    }
    assert(*(dir_path.end()-1) == '/'); // Paranoia
    // Load the info file
    std::string info_path = dir_path + "archive.info";
    std::ifstream info_file(info_path.c_str());
    if (!info_file.is_open())
    {
        EXCEPTION("Unable to open archive information file: " + info_path);
    }
    unsigned num_procs, archive_version;
    info_file >> num_procs >> archive_version;

    PROBLEM_CLASS *p_unarchived_simulation;

    // Avoid the DistributedVectorFactory throwing a 'wrong number of processes' exception when loading,
    // and make it get the original DistributedVectorFactory from the archive so we can compare against
    // num_procs.
    DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(false);
    // Put what follows in a try-catch to make sure we reset this
    try
    {
        // Figure out which process-specific archive to load first.  If we're loading on the same number of
        // processes, we must load our own one, or the mesh gets confused.  Otherwise, start with 0 to make
        // sure it exists.
        unsigned initial_archive = num_procs == PetscTools::GetNumProcs() ? PetscTools::GetMyRank() : 0u;

        // Load the master and initial process-specific archive files.
        // This will also set up ArchiveLocationInfo for us.
        ArchiveOpener<boost::archive::text_iarchive, std::ifstream> archive_opener(rDirectory, "archive.arch", initial_archive);
        boost::archive::text_iarchive* p_main_archive = archive_opener.GetCommonArchive();
        (*p_main_archive) >> p_unarchived_simulation;

        // Work out how many more process-specific files to load
        DistributedVectorFactory* p_factory = p_unarchived_simulation->rGetMesh().GetDistributedVectorFactory();
        assert(p_factory != NULL);
        unsigned original_num_procs = p_factory->GetOriginalFactory()->GetNumProcs();
        assert(original_num_procs == num_procs); // Paranoia

        // Merge in the extra data
        for (unsigned archive_num=0; archive_num<original_num_procs; archive_num++)
        {
            if (archive_num != initial_archive)
            {
                std::string archive_path = ArchiveLocationInfo::GetProcessUniqueFilePath("archive.arch", archive_num);
                std::ifstream ifs(archive_path.c_str());
                boost::archive::text_iarchive archive(ifs);
                p_unarchived_simulation->LoadExtraArchive(archive, archive_version);
            }
        }
    }
    catch (Exception &e)
    {
        DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(true);
        throw e;
    }

    // Done.
    DistributedVectorFactory::SetCheckNumberOfProcessesOnLoad(true);
    return p_unarchived_simulation;
}

//
// Explicit instantiation
//

template class CardiacSimulationArchiver<MonodomainProblem<1> >;
template class CardiacSimulationArchiver<MonodomainProblem<2> >;
template class CardiacSimulationArchiver<MonodomainProblem<3> >;

template class CardiacSimulationArchiver<BidomainProblem<1> >;
template class CardiacSimulationArchiver<BidomainProblem<2> >;
template class CardiacSimulationArchiver<BidomainProblem<3> >;

template class CardiacSimulationArchiver<BidomainWithBathProblem<1> >;
template class CardiacSimulationArchiver<BidomainWithBathProblem<2> >;
template class CardiacSimulationArchiver<BidomainWithBathProblem<3> >;
