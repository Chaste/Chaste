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
#ifndef PROCESSSPECIFICARCHIVE_HPP_
#define PROCESSSPECIFICARCHIVE_HPP_

#include <string>
#include <sstream>
#include <cassert>
#include <iostream>

#include "Exception.hpp"
#include "PetscTools.hpp"

/**
 * When checkpointing a parallel simulation, there are two kinds of data that need to be saved:
 * replicated (same for every process) and distributed (different on each process).  We wish to
 * write these to separate locations, so that the replicated data is only written to disk once,
 * and to make it easier to re-load on a different number of processes (in which case the
 * distributed data will need to be re-distributed).  However, the Boost Serialization library
 * expects to be writing to just one archive.
 *
 * This class provides access to a secondary archive in which to store the distributed data.
 * When opening an archive in a (potentially) parallel setting, using either the ArchiveOpener
 * or CardiacSimulationArchiver, the Set method will be called to specify the archive.  Classes
 * which need to save distributed data can then use the Get method to access and write to/read
 * from this archive.
 *
 * Note that because this class stores just a pointer to the archive, whatever object owns the
 * archive must ensure it exists for the duration of the serialization process, and call
 * Set(NULL) prior to closing the archive for safety.
 *
 * Note also that implementations of this templated class only exist for text archives, i.e.
 * Archive = boost::archive::text_iarchive or Archive = boost::archive::text_oarchive.
 */
template <class Archive>
class ProcessSpecificArchive
{
private:

    /** The secondary archive for this process. */
    static Archive* mpArchive;

public:

    /** Retrieve the stored secondary archive for this process. */
    static Archive* Get(void)
    {
        if (mpArchive == NULL)
        {
            EXCEPTION("A ProcessSpecificArchive has not been set up.");
        }
        return mpArchive;
    }

    /**
     * Set the secondary archive for this process.
     *
     * @param pArchive  the archive to use.
     */
    static void Set(Archive* pArchive)
    {
        mpArchive = pArchive;
    }
};

#endif /*PROCESSSPECIFICARCHIVE_HPP_*/
