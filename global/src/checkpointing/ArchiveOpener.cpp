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

// Must be included before any other serialization headers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <sstream>
#include <fstream>

#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "ProcessSpecificArchive.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"

/**
 * Specialization for input archives.
 * @param rDirectory
 * @param rFileNameBase
 * @param procId
 */
template<>
ArchiveOpener<boost::archive::text_iarchive, std::ifstream>::ArchiveOpener(
        const FileFinder& rDirectory,
        const std::string& rFileNameBase,
        unsigned procId)
    : mpCommonStream(nullptr),
      mpPrivateStream(nullptr),
      mpCommonArchive(nullptr),
      mpPrivateArchive(nullptr)
{
    // Figure out where things live
    ArchiveLocationInfo::SetArchiveDirectory(rDirectory);
    std::string private_path = ArchiveLocationInfo::GetProcessUniqueFilePath(rFileNameBase, procId);
    std::stringstream common_path;
    common_path << ArchiveLocationInfo::GetArchiveDirectory() << rFileNameBase;

    // Try to open the main archive for replicated data
    mpCommonStream = new std::ifstream(common_path.str().c_str(), std::ios::binary);
    if (!mpCommonStream->is_open())
    {
        delete mpCommonStream;
        EXCEPTION("Cannot load main archive file: " + common_path.str());
    }

    try
    {
        mpCommonArchive = new boost::archive::text_iarchive(*mpCommonStream);
    }
    catch (boost::archive::archive_exception& boost_exception)
    {
        if (boost_exception.code == boost::archive::archive_exception::unsupported_version)
        {
            // This is forward compatibility issue.  We can't open the archive because it's been written by a more recent Boost.
            delete mpCommonArchive;
            delete mpCommonStream;
            EXCEPTION("Could not open Boost archive '" + common_path.str() + "' because it was written by a more recent Boost.  Check process-specific archives too");
        }
        else
        {
            // We don't understand the exception, so we shouldn't continue
            throw boost_exception; // LCOV_EXCL_LINE
        }
    }

    // Try to open the secondary archive for distributed data
    mpPrivateStream = new std::ifstream(private_path.c_str(), std::ios::binary);
    if (!mpPrivateStream->is_open())
    {
        delete mpPrivateStream;
        delete mpCommonArchive;
        delete mpCommonStream;
        EXCEPTION("Cannot load secondary archive file: " + private_path);
    }
    mpPrivateArchive = new boost::archive::text_iarchive(*mpPrivateStream);
    ProcessSpecificArchive<boost::archive::text_iarchive>::Set(mpPrivateArchive);
}

template<>
ArchiveOpener<boost::archive::text_iarchive, std::ifstream>::~ArchiveOpener()
{
    ProcessSpecificArchive<boost::archive::text_iarchive>::Set(nullptr);
    delete mpPrivateArchive;
    delete mpPrivateStream;
    delete mpCommonArchive;
    delete mpCommonStream;
}

/**
 * Specialization for output archives.
 * @param rDirectory
 * @param rFileNameBase
 * @param procId
 */
template<>
ArchiveOpener<boost::archive::text_oarchive, std::ofstream>::ArchiveOpener(
        const FileFinder& rDirectory,
        const std::string& rFileNameBase,
        unsigned procId)
    : mpCommonStream(nullptr),
      mpPrivateStream(nullptr),
      mpCommonArchive(nullptr),
      mpPrivateArchive(nullptr)
{
    // Check for user error
    if (procId != PetscTools::GetMyRank())
    {
        EXCEPTION("Specifying the secondary archive file ID doesn't make sense when writing.");
    }

    // Figure out where things live
    ArchiveLocationInfo::SetArchiveDirectory(rDirectory);
    if (ArchiveLocationInfo::GetIsDirRelativeToChasteTestOutput())
    {
        // Ensure the directory exists
        OutputFileHandler handler(ArchiveLocationInfo::GetArchiveRelativePath(), false);
    }
    std::string private_path = ArchiveLocationInfo::GetProcessUniqueFilePath(rFileNameBase);
    std::stringstream common_path;
    common_path << ArchiveLocationInfo::GetArchiveDirectory() << rFileNameBase;

    // Create master archive for replicated data
    if (PetscTools::AmMaster())
    {
        mpCommonStream = new std::ofstream(common_path.str().c_str(), std::ios::binary | std::ios::trunc);
        if (!mpCommonStream->is_open())
        {
            delete mpCommonStream;
            EXCEPTION("Failed to open main archive file for writing: " + common_path.str());
        }
    }
    else
    {
        // Non-master processes need to go through the serialization methods, but not write any data
#ifdef _MSC_VER
        mpCommonStream = new std::ofstream("NUL", std::ios::binary | std::ios::trunc);
#else
        mpCommonStream = new std::ofstream("/dev/null", std::ios::binary | std::ios::trunc);
#endif
        // LCOV_EXCL_START
        if (!mpCommonStream->is_open())
        {
            delete mpCommonStream;
            EXCEPTION("Failed to open dummy archive file '/dev/null' for writing");
        }
        // LCOV_EXCL_STOP
    }
    mpCommonArchive = new boost::archive::text_oarchive(*mpCommonStream);

    // Create secondary archive for distributed data
    mpPrivateStream = new std::ofstream(private_path.c_str(), std::ios::binary | std::ios::trunc);
    if (!mpPrivateStream->is_open())
    {
        delete mpPrivateStream;
        delete mpCommonArchive;
        delete mpCommonStream;
        EXCEPTION("Failed to open secondary archive file for writing: " + private_path);
    }
    mpPrivateArchive = new boost::archive::text_oarchive(*mpPrivateStream);
    ProcessSpecificArchive<boost::archive::text_oarchive>::Set(mpPrivateArchive);
}

template<>
ArchiveOpener<boost::archive::text_oarchive, std::ofstream>::~ArchiveOpener()
{
    ProcessSpecificArchive<boost::archive::text_oarchive>::Set(nullptr);
    delete mpPrivateArchive;
    delete mpPrivateStream;
    delete mpCommonArchive;
    delete mpCommonStream;

    /* In a parallel setting, make sure all processes have finished writing before
     * continuing, to avoid nasty race conditions.
     * For example, many tests will write an archive then immediately read it back
     * in, which could easily break without this.
     */
    PetscTools::Barrier("~ArchiveOpener");
}
