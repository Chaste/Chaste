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
#ifndef ARCHIVEOPENER_HPP_
#define ARCHIVEOPENER_HPP_

#include <string>
#include <cassert>

#include "PetscTools.hpp"
#include "FileFinder.hpp"

/**
 * A convenience class to assist with managing archives for parallel checkpointing.
 *
 * When checkpointing a parallel simulation, there are two kinds of data that need to be saved:
 * replicated (same for every process) and distributed (different on each process).  We wish to
 * write these to separate archive files.  This class hides the complexity of doing so, such
 * that all a user needs to do is create an instance of this class, call GetCommonArchive, and
 * read from/write to the returned archive.  When done, just destroy the instance (e.g. by
 * closing the scope).
 *
 * Internally the class uses ProcessSpecificArchive<Archive> to store the secondary archive.
 *
 * Note also that implementations of this templated class only exist for text archives, i.e.
 * Archive = boost::archive::text_iarchive (with Stream = std::ifstream), or
 * Archive = boost::archive::text_oarchive (with Stream = std::ofstream).
 */
template <class Archive, class Stream>
class ArchiveOpener
{
private:
    friend class TestArchivingHelperClasses;
public:

    /**
     * Open the archives for this process, either for reading or writing depending on the
     * template parameter Archive.
     *
     * Note that when writing, only the master process writes to the main archive.  For other
     * processes the main archive is a dummy, writing to /dev/null.
     *
     * @note Must be called collectively, i.e. by all processes!
     *
     * If writing, and rDirectory is relative to CHASTE_TEST_OUTPUT, it will be created if it
     * doesn't exist.
     *
     * @param rDirectory  folder containing archive files.
     * @param rFileNameBase  base name of archive files.  This will be used for the main archive
     *     (for replicated data) with ".n" (where n is the process index) being appended for
     *     the secondary archive.
     * @param procId  this can be specified to read a specific secondary archive, rather than
     *     this process' default.  Should not be used for writing!
     */
    ArchiveOpener(const FileFinder& rDirectory,
                  const std::string& rFileNameBase,
                  unsigned procId=PetscTools::GetMyRank());

    /**
     * Close the opened archives.
     */
    ~ArchiveOpener();

    /**
     * @return the main archive for replicated data.
     */
    Archive* GetCommonArchive()
    {
        assert(mpCommonArchive != NULL);
        return mpCommonArchive;
    }

private:

    /** The file stream for the main archive. */
    Stream* mpCommonStream;

    /** The file stream for the secondary archive. */
    Stream* mpPrivateStream;

    /** The main archive. */
    Archive* mpCommonArchive;

    /** The secondary archive. */
    Archive* mpPrivateArchive;
};

#endif /*ARCHIVEOPENER_HPP_*/
