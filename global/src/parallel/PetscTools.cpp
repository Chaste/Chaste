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

#include "PetscTools.hpp"
#include "Exception.hpp"
#include "Warnings.hpp"
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstring> // For strcmp etc. Needed in gcc-4.3

bool PetscTools::mPetscIsInitialised = false;
unsigned PetscTools::mNumProcessors = 0;
unsigned PetscTools::mRank = 0;
bool PetscTools::mIsolateProcesses = false;

#ifndef NDEBUG
// Uncomment this to trace calls to PetscTools::Barrier
//#define DEBUG_BARRIERS
#endif

#ifdef DEBUG_BARRIERS
static unsigned mNumBarriers = 0u;
#endif

void PetscTools::ResetCache()
{
    PetscBool is_there;
    PetscInitialized(&is_there);
    if (is_there)
    {
        mPetscIsInitialised = true;

        PetscInt num_procs;
        MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
        mNumProcessors = (unsigned) num_procs;

        PetscInt my_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
        mRank = (unsigned) my_rank;
    }
    else
    {
        // No PETSc
        mPetscIsInitialised = false;
        mNumProcessors = 1;
        mRank = 0;
    }
}

// Information methods

bool PetscTools::IsInitialised()
{
    CheckCache();
    return mPetscIsInitialised;
}

bool PetscTools::IsSequential()
{
    CheckCache();
    return (mIsolateProcesses || mNumProcessors == 1);
}

bool PetscTools::IsParallel()
{
    CheckCache();
    return (mNumProcessors > 1);
}

bool PetscTools::IsIsolated()
{
    return mIsolateProcesses;
}

unsigned PetscTools::GetNumProcs()
{
    CheckCache();
    return mNumProcessors;
}

unsigned PetscTools::GetMyRank()
{
    CheckCache();
    return mRank;
}

bool PetscTools::AmMaster()
{
    CheckCache();
    return (mRank == MASTER_RANK || mIsolateProcesses);
}

bool PetscTools::AmTopMost()
{
    CheckCache();
    return (mRank == mNumProcessors - 1 || mIsolateProcesses);
}

// Little utility methods

void PetscTools::Barrier(const std::string callerId)
{
    CheckCache();
    if (mPetscIsInitialised && !mIsolateProcesses)
    {
#ifdef DEBUG_BARRIERS
        // "Before" is alphabetically before "Post" so that one can sort the output on process/event/barrier
        std::cout << "DEBUG: proc " << PetscTools::GetMyRank() << ": Before " << "Barrier " << mNumBarriers << " \""<< callerId <<  "\"." << std::endl << std::flush;
#endif
        PetscBarrier(PETSC_NULL);
#ifdef DEBUG_BARRIERS
        std::cout << "DEBUG: proc " << PetscTools::GetMyRank() << ": Post " << "Barrier " << mNumBarriers++ << " \""<< callerId <<  "\"." << std::endl << std::flush;
#endif
    }
}

void PetscTools::BeginRoundRobin()
{
    Barrier("PetscTools::RoundRobin"); // We want barriers both before all and after all, just in case
    const unsigned my_rank = GetMyRank();
    for (unsigned turn=0; turn<my_rank; turn++)
    {
        Barrier("PetscTools::RoundRobin");
    }
}

void PetscTools::EndRoundRobin()
{
    const unsigned num_procs = GetNumProcs();
    for (unsigned turn=GetMyRank(); turn<num_procs; turn++)
    {
        Barrier("PetscTools::RoundRobin");
    }
}

void PetscTools::IsolateProcesses(bool isolate)
{
    mIsolateProcesses = isolate;
}

MPI_Comm PetscTools::GetWorld()
{
    if (mIsolateProcesses)
    {
        return PETSC_COMM_SELF;
    }
    else
    {
        return PETSC_COMM_WORLD;
    }
}

bool PetscTools::ReplicateBool(bool flag)
{
    CheckCache();
    unsigned my_flag = (unsigned) flag;
    unsigned anyones_flag_is_true = my_flag;
    if (mPetscIsInitialised && !mIsolateProcesses)
    {
        MPI_Allreduce(&my_flag, &anyones_flag_is_true, 1, MPI_UNSIGNED, MPI_MAX, PETSC_COMM_WORLD);
    }
    return (anyones_flag_is_true == 1);
}

void PetscTools::ReplicateException(bool flag)
{
    bool anyones_error=ReplicateBool(flag);
    if (flag)
    {
        // Return control to exception thrower
        return;
    }
    if (anyones_error)
    {
        EXCEPTION("Another process threw an exception; bailing out.");
    }
}

// Vector & matrix creation routines

Vec PetscTools::CreateVec(int size, int localSize, bool ignoreOffProcEntries)
{
    assert(size >= 0); // There is one test where we create a zero-sized vector
    Vec ret;
    VecCreate(PETSC_COMM_WORLD, &ret);
    VecSetSizes(ret, localSize, size); // localSize usually defaults to PETSC_DECIDE
    VecSetFromOptions(ret);

    if (ignoreOffProcEntries)
    {
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
        VecSetOption(ret, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
        VecSetOption(ret, VEC_IGNORE_OFF_PROC_ENTRIES);
#endif
    }

    return ret;
}

Vec PetscTools::CreateVec(std::vector<double> data)
{
    assert(data.size() > 0);
    Vec ret = CreateVec(data.size());

    double* p_ret;
    VecGetArray(ret, &p_ret);
    int lo, hi;
    VecGetOwnershipRange(ret, &lo, &hi);

    for (int global_index=lo; global_index<hi; global_index++)
    {
        int local_index = global_index - lo;
        p_ret[local_index] = data[global_index];
    }
    VecRestoreArray(ret, &p_ret);

    return ret;
}

Vec PetscTools::CreateAndSetVec(int size, double value)
{
    assert(size > 0);
    Vec ret = CreateVec(size);

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    VecSet(&value, ret);
#else
    VecSet(ret, value);
#endif

    return ret;
}

void PetscTools::SetupMat(Mat& rMat, int numRows, int numColumns,
                          unsigned rowPreallocation,
                          int numLocalRows,
                          int numLocalColumns,
                          bool ignoreOffProcEntries,
                          bool newAllocationError)
{
    assert(numRows > 0);
    assert(numColumns > 0);
    if ((int) rowPreallocation>numColumns)
    {
        WARNING("Preallocation failure: requested number of nonzeros per row greater than number of columns");//+rowPreallocation+">"+numColumns);
        rowPreallocation = numColumns;
    }

#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    MatCreate(PETSC_COMM_WORLD,numLocalRows,numLocalColumns,numRows,numColumns,&rMat);
#else //New API
    MatCreate(PETSC_COMM_WORLD,&rMat);
    MatSetSizes(rMat,numLocalRows,numLocalColumns,numRows,numColumns);
#endif

    if (PetscTools::IsSequential())
    {
        MatSetType(rMat, MATSEQAIJ);
        if (rowPreallocation > 0)
        {
            MatSeqAIJSetPreallocation(rMat, rowPreallocation, PETSC_NULL);
        }
    }
    else
    {
        MatSetType(rMat, MATMPIAIJ);
        if (rowPreallocation > 0)
        {
            MatMPIAIJSetPreallocation(rMat, rowPreallocation, PETSC_NULL, rowPreallocation, PETSC_NULL);
        }
    }

    MatSetFromOptions(rMat);

    if (ignoreOffProcEntries)//&& IsParallel())
    {
        if (rowPreallocation == 0)
        {
            // We aren't allowed to do non-zero allocation after setting MAT_IGNORE_OFF_PROC_ENTRIES
            WARNING("Ignoring MAT_IGNORE_OFF_PROC_ENTRIES flag because we might set non-zeroes later");
        }
        else
        {
#if (PETSC_VERSION_MAJOR == 3) //PETSc 3.x.x
            MatSetOption(rMat, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
            MatSetOption(rMat, MAT_IGNORE_OFF_PROC_ENTRIES);
#endif
        }
    }
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 3) //PETSc 3.3 or later
    if (newAllocationError == false)
    {
        MatSetOption(rMat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    }
#endif
}

void PetscTools::DumpPetscObject(const Mat& rMat, const std::string& rOutputFileFullPath)
{
    PetscViewer view;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PetscViewerFileType type = PETSC_FILE_CREATE;
#else
    PetscFileMode type = FILE_MODE_WRITE;
#endif

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, rOutputFileFullPath.c_str(), type, &view);
    MatView(rMat, view);
    PetscViewerDestroy(PETSC_DESTROY_PARAM(view));
}

void PetscTools::DumpPetscObject(const Vec& rVec, const std::string& rOutputFileFullPath)
{
    PetscViewer view;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PetscViewerFileType type = PETSC_FILE_CREATE;
#else
    PetscFileMode type = FILE_MODE_WRITE;
#endif

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, rOutputFileFullPath.c_str(), type, &view);
    VecView(rVec, view);
    PetscViewerDestroy(PETSC_DESTROY_PARAM(view));
}

void PetscTools::ReadPetscObject(Mat& rMat, const std::string& rOutputFileFullPath, Vec rParallelLayout)
{
    /*
     * PETSc (as of 3.1) doesn't provide any method for loading a Mat object with a user-defined parallel
     * layout, i.e. there's no equivalent to VecLoadIntoVector for Mat's.
     *
     * It seems to be in their future plans though: http://lists.mcs.anl.gov/pipermail/petsc-users/2010-March/006062.html
     */

    PetscViewer view;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PetscViewerFileType type = PETSC_FILE_RDONLY;
#else
    PetscFileMode type = FILE_MODE_READ;
#endif

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, rOutputFileFullPath.c_str(),
                          type, &view);

#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
    MatCreate(PETSC_COMM_WORLD,&rMat);
    MatSetType(rMat,MATMPIAIJ);
    MatLoad(rMat,view);
#else
    MatLoad(view, MATMPIAIJ, &rMat);
#endif

    PetscViewerDestroy(PETSC_DESTROY_PARAM(view));

    if (rParallelLayout != nullptr)
    {
        /*
         * The idea is to copy rMat into a matrix that has the appropriate
         * parallel layout. Inefficient...
         */
        PetscInt num_rows, num_local_rows;

        VecGetSize(rParallelLayout, &num_rows);
        VecGetLocalSize(rParallelLayout, &num_local_rows);

        Mat temp_mat;
        /// \todo: #1082 work out appropriate nz allocation.
        PetscTools::SetupMat(temp_mat, num_rows, num_rows, 100, num_local_rows, num_local_rows, false);

        MatCopy(rMat, temp_mat, DIFFERENT_NONZERO_PATTERN);

        PetscTools::Destroy(rMat);
        rMat = temp_mat;
    }
}

void PetscTools::ReadPetscObject(Vec& rVec, const std::string& rOutputFileFullPath, Vec rParallelLayout)
{
    PetscViewer view;
#if (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 2) //PETSc 2.2
    PetscViewerFileType type = PETSC_FILE_RDONLY;
#else
    PetscFileMode type = FILE_MODE_READ;
#endif

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, rOutputFileFullPath.c_str(),
                          type, &view);
    if (rParallelLayout == nullptr)
    {
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        VecCreate(PETSC_COMM_WORLD,&rVec);
        VecSetType(rVec,VECMPI);
        VecLoad(rVec,view);
#else
        VecLoad(view, VECMPI, &rVec);
#endif
    }
    else
    {
        VecDuplicate(rParallelLayout, &rVec);
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2) //PETSc 3.2 or later
        VecLoad(rVec,view);
#else
        VecLoadIntoVector(view, rVec);
#endif
    }
    PetscViewerDestroy(PETSC_DESTROY_PARAM(view));
}

bool PetscTools::HasParMetis()
{
#ifdef __INTEL_COMPILER
    //Old versions of the intel compiler can result in a PETSC ERROR and the program aborting if parmetis is checked for before
    //some other calls to Petsc are made. This nasty hack ensures that the HasParMetis method always works on Intel
    Mat mat;
    PetscTools::SetupMat(mat, 1, 1, 1);
    PetscTools::Destroy(mat);
#endif

    MatPartitioning part;
    MatPartitioningCreate(PETSC_COMM_WORLD, &part);


    // We are expecting an error from PETSC on systems that don't have the interface, so suppress it
    // in case it aborts
    PetscPushErrorHandler(PetscIgnoreErrorHandler, nullptr);

#if (PETSC_VERSION_MAJOR == 2 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR < 2))
    PetscErrorCode parmetis_installed_error = MatPartitioningSetType(part,MAT_PARTITIONING_PARMETIS);
#else
    PetscErrorCode parmetis_installed_error = MatPartitioningSetType(part,MATPARTITIONINGPARMETIS);
#endif

    // Stop supressing error
    PetscPopErrorHandler();

    // Note that this method probably leaks memory inside PETSc because if MatPartitioningCreate fails
    // then there isn't a proper handle to destroy.
    MatPartitioningDestroy(PETSC_DESTROY_PARAM(part));

    // Get out of jail free card for Windows where the latest configuration of the integration machine shows that our implementation doesn't work as expected.
#ifdef _MSC_VER
//\todo #2016 (or similar ticket).  The method NodePartitioner::PetscMatrixPartitioning is not working in parallel
    if (parmetis_installed_error == 0)
    {
        WARN_ONCE_ONLY("The PETSc/parMETIS interface is correctly installed but does not yet work in Windows so matrix-based partitioning will be turned off.");
    }
    return false;
#endif

    return (parmetis_installed_error == 0);
}
