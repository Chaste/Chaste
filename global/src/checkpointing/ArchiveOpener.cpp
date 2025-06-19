/*

Copyright (c) 2005-2025, University of Oxford.
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
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>
#include <sstream>

#include "ArchiveLocationInfo.hpp"
#include "ArchiveOpener.hpp"
#include "Exception.hpp"
#include "OutputFileHandler.hpp"
#include "ProcessSpecificArchive.hpp"

/**
 * A helper class to add the same members and accessor functionality to all template
 * specializations of the ArchiveOpener class.
 *
 * Was added, because class member specialization is not allowed and class
 * specialization requires the members to be defined for all specializations.
 * Added to reduce code duplication
 */
template <class Archive, class Stream>
class ArchiveOpenerBase
{
private:
    friend class TestArchivingHelperClasses;

public:
    ArchiveOpenerBase() : mpCommonStream(nullptr),
                          mpPrivateStream(nullptr),
                          mpCommonArchive(nullptr),
                          mpPrivateArchive(nullptr) {}
    /**
     * @return the main archive for replicated data.
     */
    Archive* GetCommonArchive()
    {
        assert(mpCommonArchive != NULL);
        return mpCommonArchive;
    }

    ~ArchiveOpenerBase()
    {
        delete mpPrivateArchive;
        delete mpPrivateStream;
        delete mpCommonArchive;
        delete mpCommonStream;
    }

protected:
    /** The file stream for the main archive. */
    Stream* mpCommonStream;

    /** The file stream for the secondary archive. */
    Stream* mpPrivateStream;

    /** The main archive. */
    Archive* mpCommonArchive;

    /** The secondary archive. */
    Archive* mpPrivateArchive;
};

/**
 * @brief Partial class specialization to specialize class members for input archives
 *
 * @tparam InputArchive Type of the input archive type, which can vary between text and binary input archives from boost::archive
 */
template <class InputArchive>
class ArchiveOpener<InputArchive, std::ifstream> : public ArchiveOpenerBase<InputArchive, std::ifstream>
{
private:
    friend class TestArchivingHelperClasses;

public:
    /**
     * Specialization for input archives.
     * @param rDirectory
     * @param rFileNameBase
     * @param procId
     */
    ArchiveOpener(
        const FileFinder& rDirectory,
        const std::string& rFileNameBase,
        unsigned procId)
            : ArchiveOpenerBase<InputArchive, std::ifstream>()
    {
        // Figure out where things live
        ArchiveLocationInfo::SetArchiveDirectory(rDirectory);
        std::string private_path = ArchiveLocationInfo::GetProcessUniqueFilePath(rFileNameBase, procId);
        std::stringstream common_path;
        common_path << ArchiveLocationInfo::GetArchiveDirectory() << rFileNameBase;

        // Try to open the main archive for replicated data
        this->mpCommonStream = new std::ifstream(common_path.str().c_str(), std::ios::binary);
        if (!this->mpCommonStream->is_open())
        {
            delete this->mpCommonStream;
            EXCEPTION("Cannot load main archive file: " + common_path.str());
        }

        try
        {
            this->mpCommonArchive = new InputArchive(*this->mpCommonStream);
        }
        catch (boost::archive::archive_exception& boost_exception)
        {
            if (boost_exception.code == boost::archive::archive_exception::unsupported_version)
            {
                // This is forward compatibility issue.  We can't open the archive because it's been written by a more recent Boost.
                delete this->mpCommonArchive;
                delete this->mpCommonStream;
                EXCEPTION("Could not open Boost archive '" + common_path.str() + "' because it was written by a more recent Boost. Check process-specific archives too");
            }
            else
            {
                // We don't understand the exception, so we shouldn't continue
                throw boost_exception; // LCOV_EXCL_LINE
            }
        }

        // Try to open the secondary archive for distributed data
        this->mpPrivateStream = new std::ifstream(private_path.c_str(), std::ios::binary);
        if (!this->mpPrivateStream->is_open())
        {
            delete this->mpPrivateStream;
            delete this->mpCommonArchive;
            delete this->mpCommonStream;
            EXCEPTION("Cannot load secondary archive file: " + private_path);
        }
        this->mpPrivateArchive = new InputArchive(*this->mpPrivateStream);
        ProcessSpecificArchive<InputArchive>::Set(this->mpPrivateArchive);
    }

    ~ArchiveOpener()
    {
        ProcessSpecificArchive<InputArchive>::Set(nullptr);
    }
};

/**
 * @brief Partial class specialization to specialize class members for output archives
 *
 * @tparam OutputArchive Type of the output archive type, which can vary between text and binary output archives from boost::archive
 */
template <class OutputArchive>
class ArchiveOpener<OutputArchive, std::ofstream> : public ArchiveOpenerBase<OutputArchive, std::ofstream>
{
private:
    friend class TestArchivingHelperClasses;

public:
    /**
     * Specialization for output archives.
     * @param rDirectory
     * @param rFileNameBase
     * @param procId
     */
    ArchiveOpener(
        const FileFinder& rDirectory,
        const std::string& rFileNameBase,
        unsigned procId)
            : ArchiveOpenerBase<OutputArchive, std::ofstream>()
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
            this->mpCommonStream = new std::ofstream(common_path.str().c_str(), std::ios::binary | std::ios::trunc);
            if (!this->mpCommonStream->is_open())
            {
                delete this->mpCommonStream;
                EXCEPTION("Failed to open main archive file for writing: " + common_path.str());
            }
        }
        else
        {
            // Non-master processes need to go through the serialization methods, but not write any data
#ifdef _MSC_VER
            this->mpCommonStream = new std::ofstream("NUL", std::ios::binary | std::ios::trunc);
#else
            this->mpCommonStream = new std::ofstream("/dev/null", std::ios::binary | std::ios::trunc);
#endif
            // LCOV_EXCL_START
            if (!this->mpCommonStream->is_open())
            {
                delete this->mpCommonStream;
                EXCEPTION("Failed to open dummy archive file '/dev/null' for writing");
            }
            // LCOV_EXCL_STOP
        }
        this->mpCommonArchive = new OutputArchive(*this->mpCommonStream);

        // Create secondary archive for distributed data
        this->mpPrivateStream = new std::ofstream(private_path.c_str(), std::ios::binary | std::ios::trunc);
        if (!this->mpPrivateStream->is_open())
        {
            delete this->mpPrivateStream;
            delete this->mpCommonArchive;
            delete this->mpCommonStream;
            EXCEPTION("Failed to open secondary archive file for writing: " + private_path);
        }
        this->mpPrivateArchive = new OutputArchive(*this->mpPrivateStream);
        ProcessSpecificArchive<OutputArchive>::Set(this->mpPrivateArchive);
    }

    ~ArchiveOpener()
    {
        ProcessSpecificArchive<OutputArchive>::Set(nullptr);

        /* In a parallel setting, make sure all processes have finished writing before
         * continuing, to avoid nasty race conditions.
         * For example, many tests will write an archive then immediately read it back
         * in, which could easily break without this.
         */
        PetscTools::Barrier("~ArchiveOpener");
    }
};

// Explicit instantiation
template class ArchiveOpener<boost::archive::text_iarchive, std::ifstream>;
template class ArchiveOpener<boost::archive::text_oarchive, std::ofstream>;
template class ArchiveOpener<boost::archive::binary_iarchive, std::ifstream>;
template class ArchiveOpener<boost::archive::binary_oarchive, std::ofstream>;