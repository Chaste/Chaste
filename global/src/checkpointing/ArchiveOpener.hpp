/*

Copyright (C) University of Oxford, 2005-2012

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
     * Get the main archive for replicated data.
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
