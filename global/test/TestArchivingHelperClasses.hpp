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

#ifndef TESTARCHIVINGHELPERCLASSES_HPP_
#define TESTARCHIVINGHELPERCLASSES_HPP_

#include <climits>

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/foreach.hpp>
#include <boost/version.hpp>

#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"
#include "ChasteSyscalls.hpp"
#include "FileFinder.hpp"
#include "GetCurrentWorkingDirectory.hpp"
#include "OutputFileHandler.hpp"
#include "PosixPathFixer.hpp"
#include "ProcessSpecificArchive.hpp"
#include "PetscSetupAndFinalize.hpp"

// Save typing, and allow the use of these in cxxtest macros
typedef ArchiveOpener<boost::archive::text_iarchive, std::ifstream> InputArchiveOpener;
typedef ArchiveOpener<boost::archive::text_oarchive, std::ofstream> OutputArchiveOpener;

class TestArchivingHelperClasses : public CxxTest::TestSuite
{
public:
    void TestArchiveLocationInfoMethods()
    {
        // These throw because we are getting things before they are set.
        TS_ASSERT_THROWS_THIS(ArchiveLocationInfo::GetArchiveDirectory(),
                              "ArchiveLocationInfo::mDirAbsPath has not been set");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshFilename(), "mesh"); //default value

        // To test exceptions (default value is now "mesh".)
        ArchiveLocationInfo::SetMeshFilename("");
        TS_ASSERT_THROWS_THIS(ArchiveLocationInfo::GetMeshFilename(),
                              "ArchiveLocationInfo::mMeshFilename has not been set");

        ArchiveLocationInfo::SetMeshPathname("archive_dir", "mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveDirectory(),
                         OutputFileHandler::GetChasteTestOutputDirectory() + "archive_dir/");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshFilename(), "mesh_name");

        // With absolute path...
        ArchiveLocationInfo::SetMeshPathname(ArchiveLocationInfo::GetArchiveDirectory(), "mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveDirectory(),
                         OutputFileHandler::GetChasteTestOutputDirectory() + "archive_dir/");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshFilename(), "mesh_name");

        FileFinder dir("new_archive_dir", RelativeTo::CWD);
        ArchiveLocationInfo::SetArchiveDirectory(dir);
        ArchiveLocationInfo::SetMeshFilename("new_mesh_name");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveDirectory(),
                         GetCurrentWorkingDirectory() + "/new_archive_dir/");
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetMeshFilename(), "new_mesh_name");

        // This gives the absolute path, since it isn't relative to CHASTE_TEST_OUTPUT
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveRelativePath(),
                         GetCurrentWorkingDirectory() + "/new_archive_dir/");
        TS_ASSERT(!ArchiveLocationInfo::GetIsDirRelativeToChasteTestOutput());
        FileFinder dir2("relative_archive_dir", RelativeTo::ChasteTestOutput);
        ArchiveLocationInfo::SetArchiveDirectory(dir2);
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetArchiveRelativePath(), "relative_archive_dir/");
        TS_ASSERT(ArchiveLocationInfo::GetIsDirRelativeToChasteTestOutput());
    }

    void TestArchiveLocationInfoProcessUniqueNaming()
    {
        FileFinder dir("new_archive_dir", RelativeTo::CWD);
        ArchiveLocationInfo::SetArchiveDirectory(dir);

        std::stringstream expected_filepath;
        expected_filepath << ArchiveLocationInfo::GetArchiveDirectory() << "fred";
        expected_filepath << "." << PetscTools::GetMyRank();

        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetProcessUniqueFilePath("fred"), expected_filepath.str());

        std::string expected2 = ArchiveLocationInfo::GetArchiveDirectory() + "fred.12";
        TS_ASSERT_EQUALS(ArchiveLocationInfo::GetProcessUniqueFilePath("fred", 12), expected2);
    }

    void TestProcessSpecificArchive()
    {
        TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_oarchive>::Get(),
                              "A ProcessSpecificArchive has not been set up.");
        TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_iarchive>::Get(),
                              "A ProcessSpecificArchive has not been set up.");

        // Set up an output archive pointer
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string arch_path = ArchiveLocationInfo::GetProcessUniqueFilePath("test.arch");
        std::ofstream ofs(arch_path.c_str());
        boost::archive::text_oarchive* p_arch = new boost::archive::text_oarchive(ofs);

        // Test the ProcessSpecificArchive Get and Set methods with this
        ProcessSpecificArchive<boost::archive::text_oarchive>::Set(p_arch);
        TS_ASSERT(ProcessSpecificArchive<boost::archive::text_oarchive>::Get() == p_arch);
        delete p_arch;

        // Set up an input archive pointer
        std::ifstream ifs(arch_path.c_str());
        boost::archive::text_iarchive* p_arch2 = new boost::archive::text_iarchive(ifs);

        // Test the ProcessSpecificArchive Get and Set methods with this
        ProcessSpecificArchive<boost::archive::text_iarchive>::Set(p_arch2);
        TS_ASSERT(ProcessSpecificArchive<boost::archive::text_iarchive>::Get() == p_arch2);
        delete p_arch2;

        // Clean up
        ProcessSpecificArchive<boost::archive::text_iarchive>::Set(NULL);
        ProcessSpecificArchive<boost::archive::text_oarchive>::Set(NULL);
    }

    std::string mArchiveDir;

    void TestArchiveOpenerReadAndWrite()
    {
        // Should this test fail with an exception involving
        // apps/texttest/chaste/resume_bidomain/save_bidomain
        // then look at TestCardiacSimulationArchiver
        mArchiveDir = "archiving_helpers";
        FileFinder archive_dir(mArchiveDir, RelativeTo::ChasteTestOutput);
        std::string archive_file = "archive_opener.arch";
        const unsigned test_int = 123;

        // Write
        {
            OutputArchiveOpener archive_opener_out(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = archive_opener_out.GetCommonArchive();
            boost::archive::text_oarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_oarchive>::Get();

            (*p_arch) & test_int; // All can write to the common archive - non-masters will write to /dev/null.
            (*p_process_arch) & test_int;

            // archive_opener_out will do a PetscTools::Barrier when it is destructed
        }

        // Read
        {
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_oarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_iarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");

            InputArchiveOpener archive_opener_in(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = archive_opener_in.GetCommonArchive();
            boost::archive::text_iarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_iarchive>::Get();

            unsigned test_int1, test_int2;
            (*p_arch) & test_int1;
            (*p_process_arch) & test_int2;

            TS_ASSERT_EQUALS(test_int1, test_int);
            TS_ASSERT_EQUALS(test_int2, test_int);
        }

        // Cover the case of an archive in the chaste folder (i.e. a path relative to the working directory)
        if (PetscTools::IsSequential())
        {
            // Read
            FileFinder save_bidomain_dir("apps/texttest/chaste/resume_bidomain/save_bidomain", RelativeTo::ChasteSourceRoot);
            InputArchiveOpener archive_opener_relative(save_bidomain_dir, "archive.arch");
        }

        PetscTools::Barrier(); // Make sure all processes have finished this test before proceeding
    }

    // This test relies on TestArchiveOpenerReadAndWrite succeeding
    void TestArchiveOpenerExceptions()
    {
        OutputFileHandler handler(mArchiveDir, false);
        handler.SetArchiveDirectory();
        FileFinder archive_dir_finder(mArchiveDir, RelativeTo::ChasteTestOutput);
        std::string archive_base_name = "archive_opener.arch";

        // Remove the process-specific archive for this process
        FileFinder(ArchiveLocationInfo::GetProcessUniqueFilePath(archive_base_name)).Remove();
        TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener archive_opener_in(archive_dir_finder, archive_base_name),
                                  "Cannot load secondary archive file: ");
        PetscTools::Barrier("TestArchiveOpenerExceptions-1");

        // Remove the main archive
        if (PetscTools::AmMaster())
        {
            ABORT_IF_THROWS(handler.FindFile(archive_base_name).Remove());
        }
        PetscTools::Barrier("TestArchiveOpenerExceptions-2");
        TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener archive_opener_in(archive_dir_finder, archive_base_name),
                                  "Cannot load main archive file: ");

// Remove write permissions on the archive dir
//Note: changing *directory* permissions and other attributes does not work on Windows
//See http://support.microsoft.com/kb/326549
#ifndef _MSC_VER
        if (PetscTools::AmMaster())
        {
            chmod(handler.GetOutputDirectoryFullPath().c_str(), CHASTE_READONLY);
        }
        PetscTools::Barrier("TestArchiveOpenerExceptions-3");

        /*
         * Now neither the master nor the slaves can write to their output files.
         * This avoids hitting a PetscBarrier() in the ~ArchiveOpener() because they
         * all throw an error first.
         *
         * If this test starts hanging it is because these TS_ASSERT_THROWS_CONTAINS
         * are not being thrown (rather than a real parallel calling problem).
         */
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_THROWS_CONTAINS(OutputArchiveOpener archive_opener_out(archive_dir_finder, archive_base_name),
                                      "Failed to open main archive file for writing: ");
        }
        else
        {
            TS_ASSERT_THROWS_CONTAINS(OutputArchiveOpener archive_opener_out(archive_dir_finder, archive_base_name),
                                      "Failed to open secondary archive file for writing: ");
        }
        PetscTools::Barrier("TestArchiveOpenerExceptions-4");
        if (PetscTools::AmMaster())
        {
            // Restore permissions on the folder before allowing processes to continue.
            chmod(handler.GetOutputDirectoryFullPath().c_str(), CHASTE_READ_WRITE_EXECUTE);
        }
#endif // _MSC_VER
        PetscTools::Barrier("TestArchiveOpenerExceptions-5");
    }

    void TestSpecifyingSecondaryArchive()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "specific_secondary.arch";
        const unsigned test_int = 321;
        const unsigned proc_id = PetscTools::GetMyRank();

        // Writing when specifying the secondary archive doesn't make sense
        {
            TS_ASSERT_THROWS_THIS(OutputArchiveOpener archive_opener_out(archive_dir, archive_file, UINT_MAX),
                                  "Specifying the secondary archive file ID doesn't make sense when writing.");
        }

        // Write normally so we can test reading
        {
            OutputArchiveOpener archive_opener_out(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = archive_opener_out.GetCommonArchive();
            boost::archive::text_oarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_oarchive>::Get();

            (*p_arch) & test_int; // All can write to the common archive - non-masters will write to /dev/null.
            (*p_process_arch) & proc_id;

            // archive_opener_out will do a PetscTools::Barrier when it is destructed
        }

        // Read
        {
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_oarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");
            TS_ASSERT_THROWS_THIS(ProcessSpecificArchive<boost::archive::text_iarchive>::Get(),
                                  "A ProcessSpecificArchive has not been set up.");

            InputArchiveOpener archive_opener_in(archive_dir, archive_file, 0);
            boost::archive::text_iarchive* p_arch = archive_opener_in.GetCommonArchive();
            boost::archive::text_iarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_iarchive>::Get();

            unsigned test_int1, test_int2;
            (*p_arch) & test_int1;
            (*p_process_arch) & test_int2;

            TS_ASSERT_EQUALS(test_int1, test_int);
            TS_ASSERT_EQUALS(test_int2, 0u);
        }
    }

    void TestOpenFutureBoostArchive()
    {
        //Check testout/archive/specific_secondary.arch
        FileFinder archive_dir("global/test/data", RelativeTo::ChasteSourceRoot);
        std::string archive_file = "future_boost.arch";
        // future_boost has got archive version 18 in it
        // 1.33 => 3
        // 1.34 => 4
        // 1.36 => 5
        // 1.37 => 5
        // 1.40 => 5
        // 1.42 => 7
        // 1.46 => 9
        // 1.48 => 9
        // 1.49 => 9
        // 1.51 => 9
        // 1.52 => ??
        // 1.53 => 10
        // 1.54 => 10
        // 1.55 => 10
        // 1.56 => 11
        // 1.57 => 11
        // 1.58 => 12
        // 1.59 => 13
        // 1.60 => 14
        // 1.61 => 14
        // 1.62 => 14
        // 1.63 => 14
        // 1.64 => 15
        // 1.65 => 15
        // 1.66 => 16
        // 1.67 => 16
        // 1.68 => 17
        // 1.69 => 17

#ifndef BOOST_VERSION
        TS_FAIL("This test needs to know the version of Boost with which it was compiled.");
        return;
#endif
        //#if BOOST_VERSION >= 999999
        //        InputArchiveOpener archive_opener_in(archive_dir, archive_file, 0);
        //        boost::archive::text_iarchive* p_arch = archive_opener_in.GetCommonArchive();
        //        boost::archive::text_iarchive* p_process_arch = ProcessSpecificArchive<boost::archive::text_iarchive>::Get();
        //
        //        const unsigned test_int = 321;
        //        unsigned test_int1, test_int2;
        //        (*p_arch) & test_int1;
        //        (*p_process_arch) & test_int2;
        //
        //        TS_ASSERT_EQUALS(test_int1, test_int);
        //        TS_ASSERT_EQUALS(test_int2, 0u);
        //#else
        //Current Boost can't read this archive...
        TS_ASSERT_THROWS_CONTAINS(InputArchiveOpener archive_opener_in(archive_dir, archive_file, 0),
                                  "Could not open Boost archive '");
        //#endif
    }
};

#endif /*TESTARCHIVINGHELPERCLASSES_HPP_*/
