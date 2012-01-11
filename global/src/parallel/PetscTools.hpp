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

#ifndef PETSCTOOLS_HPP_
#define PETSCTOOLS_HPP_

/**
 * @file
 * Contains the PetscTools class.
 */

#include <string>
#include <vector>
#include <cstdlib> // For EXIT_FAILURE

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

/** For use in tests that do not work when run in parallel. */
#define EXIT_IF_PARALLEL if(PetscTools::IsParallel()){TS_TRACE("This test does not pass in parallel yet.");return;}
/** For use in tests that should ONLY be run in parallel. */
#define EXIT_IF_SEQUENTIAL if(PetscTools::IsSequential()){TS_TRACE("This test is not meant to be executed in sequential.");return;}

#ifndef NDEBUG
// Uncomment this to trace calls to PetscTools::Barrier
//#define DEBUG_BARRIERS
#endif

/**
 * A helper class of static methods.
 */
class PetscTools
{
private:

    /** Whether PETSc has been initialised. */
    static bool mPetscIsInitialised;

    /** The total number of processors. */
    static unsigned mNumProcessors;

#ifdef DEBUG_BARRIERS
    /** Used to debug number of barriers. */
    static unsigned mNumBarriers;
#endif

    /** Which processors we are. */
    static unsigned mRank;

    /** Whether to pretend that we're just running many master processes independently. */
    static bool mIsolateProcesses;

    /** Private method makes sure that (if this is the first use within a test) then PETSc has been probed. */
    static inline void CheckCache()
    {
        if (mNumProcessors == 0)
        {
            ResetCache();
        }
    }

public:

    /**
     * As a convention, we consider processor 0 the master process.
     */
    static const unsigned MASTER_RANK=0;

    /**
     * Reset our cached values: whether PETSc is initialised,
     * how many processors there are, and which one we are.
     */
    static void ResetCache();

    /**
     * Just returns whether there is one process or not.
     */
    static bool IsSequential();

    /**
     * Just returns whether there is more than one process.
     */
    static bool IsParallel();

    /**
     * Returns total number of processors.
     */
    static unsigned GetNumProcs();

    /**
     * Return our rank.
     *
     * If PETSc has not been initialized, returns 0.
     */
    static unsigned GetMyRank();

    /**
     * Just returns whether it is the master process or not.
     *
     * If not running in parallel, always returns true.
     */
    static bool AmMaster();

    /**
     * Just returns whether it is the right-most process or not.
     *
     * If not running in parallel, always returns true.
     */
    static bool AmTopMost();

    /**
     * If MPI is set up, perform a barrier synchronisation.
     * If not, it's a noop.
     *
     * @param callerId  only used in debug mode; printed before & after the barrier call
     */
    static void Barrier(const std::string callerId="");

    /**
     * Call at the start of a block of code that should be executed by each process in turn.
     */
    static void BeginRoundRobin();

    /**
     * Call at the end of a block of code that should be executed by each process in turn.
     */
    static void EndRoundRobin();

    /**
     * Where work can be split between isolated processes, it would be nice to be able to do
     * so easily without worrying about collective calls made inside classes such as OutputFileHandler
     * leading to deadlock.
     * This method attempts to enable this behaviour.  If the flag is set then AmMaster always
     * returns true, Barrier becomes a no-op, and ReplicateBool doesn't replicate.
     *
     * @param isolate  whether to consider processes as isolated
     */
    static void IsolateProcesses(bool isolate=true);

    /**
     * Create a vector of the specified size. SetFromOptions is called.
     *
     * @param size  the size of the vector
     * @param localSize  the local number of items owned by this process
     * @param ignoreOffProcEntries whether to ignore entries destined to be stored on a separate processor when assembling (eliminates global reductions).
     */
    static Vec CreateVec(int size, int localSize=PETSC_DECIDE, bool ignoreOffProcEntries = true);

    /**
     * Create a Vec from the given data.
     *
     * @param data  some data
     */
    static Vec CreateVec(std::vector<double> data);

    /**
     * Create a vector of the specified size with all values set to be the given
     * constant. SetFromOptions is called.
     *
     * @param size  the size of the vector
     * @param value  the value to set each entry
     */
    static Vec CreateAndSetVec(int size, double value);

    /**
     * Set up a matrix - set the size using the given parameters. The number of local rows
     * and columns is by default PETSC_DECIDE. SetFromOptions is called.
     *
     * @param rMat the matrix
     * @param numRows the number of rows in the matrix
     * @param numColumns the number of columns in the matrix
     * @param rowPreallocation the max number of nonzero entries expected on a row
     *   A value of 0 is allowed: no preallocation is then done and the user must
     *   preallocate the memory for the matrix themselves.
     * @param numLocalRows the number of local rows (defaults to PETSC_DECIDE)
     * @param numLocalColumns the number of local columns (defaults to PETSC_DECIDE)
     * @param ignoreOffProcEntries tells PETSc to drop off-processor entries
     */
    static void SetupMat(Mat& rMat, int numRows, int numColumns,
                         unsigned rowPreallocation,
                         int numLocalRows=PETSC_DECIDE,
                         int numLocalColumns=PETSC_DECIDE,
                         bool ignoreOffProcEntries=true);

    /**
     * Boolean AND of a flags between processes
     *
     * @param flag is set to true on this process.
     * @return true if any process' flag is set to true.
     */
    static bool ReplicateBool(bool flag);

    /**
     * Ensure exceptions are handled cleanly in parallel code, by causing all processes to
     * throw if any one does.
     *
     * @param flag is set to true if this process has thrown.
     */
    static void ReplicateException(bool flag);

    /**
     * Dumps a given PETSc object to disk.
     *
     * @param rMat a matrix
     * @param rOutputFileFullPath where to dump the matrix to disk
     */
    static void DumpPetscObject(const Mat& rMat, const std::string& rOutputFileFullPath);

    /**
     * Dumps a given PETSc object to disk.
     *
     * @param rVec a vector
     * @param rOutputFileFullPath where to dump the vector to disk
     */
    static void DumpPetscObject(const Vec& rVec, const std::string& rOutputFileFullPath);

    /**
     * Read a previously dumped PETSc object from disk.
     *
     * @param rMat a matrix
     * @param rOutputFileFullPath where to read the matrix from
     * @param rParallelLayout If provided, rMat will have the same parallel layout. Its content is irrelevant.
     */
    static void ReadPetscObject(Mat& rMat, const std::string& rOutputFileFullPath, Vec rParallelLayout=NULL);

    /**
     * Read a previously dumped PETSc object from disk.
     *
     * @param rVec a vector
     * @param rOutputFileFullPath where to read the matrix from
     * @param rParallelLayout If provided, rMat will have the same parallel layout. Its content is irrelevant.
     */
    static void ReadPetscObject(Vec& rVec, const std::string& rOutputFileFullPath, Vec rParallelLayout=NULL);

    /**
     * Destroy method
     * Note that PETSc 3.1 and previous destroy based on a PETSc object
     * but PETSc 3.2 and later destroy based on a pointer to a PETSc object
     * 
     * @param rVec a reference to the PETSc object
     */
     static inline void Destroy(Vec& rVec)
     {
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        VecDestroy(&rVec);
#else
        VecDestroy(rVec);
#endif   
     }
     
    /**
     * Destroy method
     * Note that PETSc 3.1 and previous destroy based on a PETSc object
     * but PETSc 3.2 and later destroy based on a pointer to a PETSc object
     * 
     * @param rMat a reference to the PETSc object
     */
     static inline void Destroy(Mat& rMat)
     {
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        MatDestroy(&rMat);
#else
        MatDestroy(rMat);
#endif   
     }        
};

#endif /*PETSCTOOLS_HPP_*/
