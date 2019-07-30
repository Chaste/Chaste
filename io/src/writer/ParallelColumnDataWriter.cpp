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

#include "ParallelColumnDataWriter.hpp"
#include "Exception.hpp"
#include "DistributedVectorFactory.hpp"

ParallelColumnDataWriter::ParallelColumnDataWriter(const std::string& rDirectory,
                                                   const std::string& rBaseName,
                                                   bool cleanDirectory)
    : ColumnDataWriter(rDirectory, rBaseName, cleanDirectory),
      mConcentrated(nullptr)
{
    int num_procs;
    MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
    if (num_procs==1)
    {
        mIsParallel = false;
    }
    else
    {
        mIsParallel = true;
    }
}

void ParallelColumnDataWriter::PutVector(int variableID, Vec petscVector)
{
    int size;
    VecGetSize(petscVector,&size);

    if (size != mFixedDimensionSize)
    {
        EXCEPTION("Size of vector does not match FixedDimensionSize.");
    }

    // Construct the appropriate "scatter" object to concentrate the vector on the master
    if (mConcentrated==nullptr)
    {
        VecScatterCreateToZero(petscVector, &mToMaster, &mConcentrated);
    }

//    int size2;
//    VecGetSize(mConcentrated, &size2);
//    std::cout << "Vector size=" << size << "," << size2 << std::endl << std::flush;

//PETSc-3.x.x or PETSc-2.3.3
#if ((PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
    VecScatterBegin(mToMaster, petscVector, mConcentrated, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(mToMaster, petscVector, mConcentrated, INSERT_VALUES, SCATTER_FORWARD);
#else
    VecScatterBegin(petscVector, mConcentrated, INSERT_VALUES, SCATTER_FORWARD, mToMaster);
    VecScatterEnd(petscVector, mConcentrated, INSERT_VALUES, SCATTER_FORWARD, mToMaster);
#endif

//    std::cout << "Done scatter" << std::endl << std::flush;

    if (PetscTools::AmMaster())
    {
        double *concentrated_vector;
        VecGetArray(mConcentrated, &concentrated_vector);
        for (int i=0; i<size; i++)
        {
            ColumnDataWriter::PutVariable(variableID, concentrated_vector[i], i);
        }
        VecRestoreArray(mConcentrated, &concentrated_vector);
    }
}

void ParallelColumnDataWriter::PutVectorStripe(int variableId, DistributedVector::Stripe& rStripe)
{
    // Put the stripe into its own 'unstriped' vector
    DistributedVectorFactory* p_factory = rStripe.GetFactory();
    Vec unstriped_petsc = p_factory->CreateVec();
    DistributedVector unstriped = p_factory->CreateDistributedVector(unstriped_petsc);
    for (DistributedVector::Iterator index = unstriped.Begin();
         index!= unstriped.End();
         ++index)
    {
        unstriped[index] = rStripe[index];
    }

    // Put the unstriped vector
    ParallelColumnDataWriter::PutVector(variableId, unstriped_petsc);
    PetscTools::Destroy(unstriped_petsc);
}

void ParallelColumnDataWriter::EndDefineMode()
{
    if (PetscTools::AmMaster())
    {
        ColumnDataWriter::EndDefineMode();
    }
    else
    {
        mIsInDefineMode = false;
    }
}

/**
 * There are two ways of calling PutVariable:
 * 1) All processes call it as a collective operation from the user's code.
 *    This only makes sense if they are writing the unlimited dimension (time) variable.
 *    It is actually a no-op if any non-master process attempts to write anything at all.
 * 2) The master calls the equivalent method in the parent class after concentrating
 *      the data into a single Vec (ie. from the method PutVector() above).
 */
void ParallelColumnDataWriter::PutVariable(int variableID, double variableValue, long dimensionPosition)
{
    if (PetscTools::AmMaster())
    {
        // Master process is allowed to write
        ColumnDataWriter::PutVariable(variableID, variableValue, dimensionPosition);
    }
}

ParallelColumnDataWriter::~ParallelColumnDataWriter()
{
    if (mConcentrated != nullptr)
    {
        VecScatterDestroy(PETSC_DESTROY_PARAM(mToMaster));
        PetscTools::Destroy(mConcentrated);
    }
    Close();
}

void ParallelColumnDataWriter::AdvanceAlongUnlimitedDimension()
{
    // Paranoia
    PetscTools::Barrier("ParallelColumnDataWriter::AdvanceAlongUnlimitedDimension");

    if (PetscTools::AmMaster())
    {
        ColumnDataWriter::DoAdvanceAlongUnlimitedDimension();
    }
}

void ParallelColumnDataWriter::Close()
{
    // Paranoia
    PetscTools::Barrier("ParallelColumnDataWriter::Close");

    if (PetscTools::AmMaster())
    {
        ColumnDataWriter::Close();
    }
}
