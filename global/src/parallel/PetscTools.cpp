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
//unsigned PetscTools::mMaxNumNonzerosIfMatMpiAij = 54;

#ifdef DEBUG_BARRIERS
unsigned PetscTools::mNumBarriers = 0u;
#endif

void PetscTools::ResetCache()
{
    PetscTruth is_there;
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

bool PetscTools::IsSequential()
{
    CheckCache();
    return (mNumProcessors == 1);
}

bool PetscTools::IsParallel()
{
    CheckCache();
    return (mNumProcessors > 1);
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
    return (mRank == mNumProcessors - 1);
}

// Little utility methods

void PetscTools::Barrier(const std::string callerId)
{
    CheckCache();
    if (mPetscIsInitialised && !mIsolateProcesses)
    {
#ifdef DEBUG_BARRIERS
        std::cout << "DEBUG: proc " << PetscTools::GetMyRank() << ": Pre " << callerId << " Barrier " << mNumBarriers << "." << std::endl << std::flush;
#endif
        PetscBarrier(PETSC_NULL);
#ifdef DEBUG_BARRIERS
        std::cout << "DEBUG: proc " << PetscTools::GetMyRank() << ": Post " << callerId << " Barrier " << mNumBarriers++ << "." << std::endl << std::flush;
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
#if (PETSC_VERSION_MAJOR == 3)
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
                          bool ignoreOffProcEntries)
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
            ///\todo #1216 Fix the 0.7 magic number
            MatMPIAIJSetPreallocation(rMat, rowPreallocation, PETSC_NULL, (PetscInt) (rowPreallocation*0.7), PETSC_NULL);
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
#if (PETSC_VERSION_MAJOR == 3)
            MatSetOption(rMat, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
#else
            MatSetOption(rMat, MAT_IGNORE_OFF_PROC_ENTRIES);
#endif
        }
    }
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
    PetscViewerDestroy(view);
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
    PetscViewerDestroy(view);
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
    MatLoad(view, MATMPIAIJ, &rMat);
    PetscViewerDestroy(view);

    if (rParallelLayout != NULL)
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

        MatDestroy(rMat);
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
    if (rParallelLayout == NULL)
    {
        VecLoad(view, VECMPI, &rVec);
    }
    else
    {
        VecDuplicate(rParallelLayout, &rVec);
        VecLoadIntoVector(view, rVec);
    }
    PetscViewerDestroy(view);
}
