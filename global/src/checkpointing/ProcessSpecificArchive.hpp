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
#ifndef PROCESSSPECIFICARCHIVE_HPP_
#define PROCESSSPECIFICARCHIVE_HPP_

#include "Exception.hpp"

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

    /** @return the stored secondary archive for this process. */
    static Archive* Get(void)
    {
        if (mpArchive == nullptr)
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
