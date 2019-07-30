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

#ifndef DISTRIBUTEDVECTORFACTORY_HPP_
#define DISTRIBUTEDVECTORFACTORY_HPP_

#include "ChasteSerialization.hpp"
#include <petscvec.h>

#include "DistributedVector.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

/**
 * Factory for creating PETSc vectors distributed across processes.
 *
 * Replacement for the vector creation portions of DistributedVector (which
 * was implemented using static methods and data), the factory class allows
 * several patterns of PETSc vector length (and distributions among
 * processes) to co-exist.
 *
 * All vectors created by a factory instance will have the same base size
 * and parallelisation pattern.
 */
 class DistributedVectorFactory
{
private:

    // Data global to all vectors created by this factory

    /** The first entry owned by the current processor. */
    unsigned mLo;

    /** One above the last entry owned by the current processor. */
    unsigned mHi;

    /** The problem size, i.e. the number of nodes in the mesh (the number of unknowns may be larger in a Stripe). */
    unsigned mProblemSize;

    /** How many processes this factory is expecting. */
    unsigned mNumProcs;

    /** Whether we've checked that PETSc is initialised. */
    bool mPetscStatusKnown;

    /** A cached vector of #mLo values from each process. */
    std::vector<unsigned> mGlobalLows;

    /**
     * Whether, when loading an instance from an archive, to check that the
     * current number of processes matches that used in creating the archive.
     */
    static bool msCheckNumberOfProcessesOnLoad;

    /**
     * If this instance was loaded from an archive, this points to a factory with
     * the settings from the archive, which may not be the same as this instance.
     */
    DistributedVectorFactory* mpOriginalFactory;

    /**
     * Double check (in debug code) that PETSc has been initialised properly.
     */
    void CheckForPetsc();

    /**
     * Helper method for the constructors.
     *
     * @param vec the sample PETSc vector from which to calculate ownerships
     */
    void CalculateOwnership(Vec vec);

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Nothing to do - all done in load_construct_data
    }

public:

    /**
     * Set the problem with an existing PETSc vector -- must have stride=1.
     *
     * @param vec is a PETSc vector which we want to use as the pattern for future vectors produced by this factory
     */
    DistributedVectorFactory(Vec vec);

    /**
     * Set the problem size specifying distribution over local processor.
     *
     * @param size  the problem size
     * @param local defaults to PETSc's default
     */
    DistributedVectorFactory(unsigned size, PetscInt local=PETSC_DECIDE);

    /**
     * Constructor for use in archiving.
     * Note that this constructor is only called when the number of processes is different from the original.
     * Therefore, the orignal local node ownership cannot be used, and a new even partition will be applied.
     *
     * @param pOriginalFactory  see #mpOriginalFactory
     */
    DistributedVectorFactory(DistributedVectorFactory* pOriginalFactory);

    /**
     * Constructor intended for use in archiving.  Allows complete manual
     * specification of the factory.
     *
     * @param lo  first index owned by this process
     * @param hi  one beyond last index owned by this process
     * @param size  total size of vectors
     * @param numProcs  the number of processes expected (defaults to the current number)
     */
    DistributedVectorFactory(unsigned lo, unsigned hi, unsigned size,
                             unsigned numProcs=PetscTools::GetNumProcs());

    /** Destructor deletes #mpOriginalFactory if it exists. */
    ~DistributedVectorFactory();

    /**
     * Create a PETSc vector of the problem size.
     *
     * @return new PETSc vector
     */
    Vec CreateVec();

    /**
     * Create a striped PETSc vector of size: stride * problem size.
     *
     * @param stride
     * @return new PETSc vector
     */
    Vec CreateVec(unsigned stride);

    /**
     * Create a distributed vector which wraps a given petsc vector.
     *
     * @param vec is the vector
     * @param readOnly flag to pass on to request that the distributed vector will only be used read-only
     * @return the distributed vector
     */
    DistributedVector CreateDistributedVector(Vec vec, bool readOnly=false);

    /**
     * Test if the given global index is owned by the current process, i.e. is local to it.
     *
     * @param globalIndex a global index
     * @return true if the global index can be accessed by this process
     */
    bool IsGlobalIndexLocal(unsigned globalIndex);

    /**
     * @return The number of elements in the vector owned by the local process.
     */
    unsigned GetLocalOwnership() const
    {
        return mHi - mLo;
    }

    /**
     * @return #mHi - The next index above the top one owned by the process.
     */
    unsigned GetHigh() const
    {
        return mHi;
    }

    /**
     * @return #mLo - The lowest index owned by the process.
     */
    unsigned GetLow() const
    {
        return mLo;
    }

    /**
     * @return The number of elements in (normal) vectors created by CreateVec().
     */
    unsigned GetProblemSize() const
    {
        return mProblemSize;
    }

    /**
     * @return #mNumProcs - how many processes this factory is expecting.
     */
    unsigned GetNumProcs() const
    {
        return mNumProcs;
    }

    /**
     * Set whether, when loading an instance from an archive, to check that the
     * current number of processes matches that used in creating the archive.
     *
     * @param checkNumberOfProcessesOnLoad
     */
    static void SetCheckNumberOfProcessesOnLoad(bool checkNumberOfProcessesOnLoad=true)
    {
        msCheckNumberOfProcessesOnLoad = checkNumberOfProcessesOnLoad;
    }

    /**
     * Determine whether, when loading an instance from an archive, to check that the
     * current number of processes matches that used in creating the archive.
     * @return true when check is required
     */
    static bool CheckNumberOfProcessesOnLoad()
    {
        return msCheckNumberOfProcessesOnLoad;
    }

    /**
     * If #msCheckNumberOfProcessesOnLoad is not set, and this factory was loaded from
     * an archive, then return a factory with the settings from the archive, which may
     * not be the same as ours - if running on a different number of processes from
     * the original, we will have used PETSC_DECIDE to set the local ownership on load.
     * @return the original factory
     */
    DistributedVectorFactory* GetOriginalFactory()
    {
        return mpOriginalFactory;
    }

    /**
     * Set method for #mpOriginalFactory, used by archiving (load_construct_data).
     * @param pOriginalFactory  see #mpOriginalFactory
     */
    void SetOriginalFactory(DistributedVectorFactory* pOriginalFactory)
    {
        mpOriginalFactory = pOriginalFactory;
    }

    /**
     * Set #mLo and #mHi from another vector factory. Used by archiving.
     * @param pFactory  the factory to set from.
     */
    void SetFromFactory(DistributedVectorFactory* pFactory);

    /**
     * @return the mLo value from each process in a vector. This is calculated on the first call
     * and cached for later use.
     */
    std::vector<unsigned> &rGetGlobalLows();

//    /**
//     * For debugging.
//     */
//    void Dump(std::string msg=std::string())
//    {
//        std::cout << "DVF(" << this << "): " << msg;
//        std::cout << " lo=" << mLo << " hi=" << mHi << " size=" << mProblemSize << " np=" << mNumProcs;
//        std::cout << " [running np=" << PetscTools::GetNumProcs() << " rank=" << PetscTools::GetMyRank() << "]\n" << std::flush;
//    }
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(DistributedVectorFactory)

namespace boost
{
namespace serialization
{

template<class Archive>
inline void save_construct_data(
    Archive & ar, const DistributedVectorFactory * t, const unsigned int file_version)
{
    unsigned num_procs, lo, hi, size;
    hi = t->GetHigh();
    ar << hi;
    lo = t->GetLow();
    ar << lo;
    size = t->GetProblemSize();
    ar << size;
    num_procs = PetscTools::GetNumProcs();
    ar << num_procs;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate an instance (using existing constructor).
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, DistributedVectorFactory * t, const unsigned int file_version)
{
    unsigned num_procs, lo, hi, size;
    ar >> hi;
    ar >> lo;
    ar >> size;
    ar >> num_procs;
    DistributedVectorFactory* p_original_factory = new DistributedVectorFactory(lo, hi, size, num_procs);

    if (!DistributedVectorFactory::CheckNumberOfProcessesOnLoad())
    {
        ::new(t)DistributedVectorFactory(p_original_factory);
    }
    else
    {
        if (num_procs != PetscTools::GetNumProcs())
        {
            // We need to have a ::new here, or Boost will try to free non-allocated memory
            // when the exception is thrown.  However, if we use the ::new line after this if,
            // then PETSc complains about wrong sizes.
            ::new(t)DistributedVectorFactory(size);
            EXCEPTION("This archive was written for a different number of processors");
        }
        ::new(t)DistributedVectorFactory(size, hi-lo);
        t->SetOriginalFactory(p_original_factory);
    }
}

}
} // namespace ...

#endif /*DISTRIBUTEDVECTORFACTORY_HPP_*/
