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

#include "ParallelColumnDataWriter.hpp"
#include "Exception.hpp"
#include "DistributedVectorFactory.hpp"

ParallelColumnDataWriter::ParallelColumnDataWriter(const std::string& rDirectory,
                                                   const std::string& rBaseName,
                                                   bool cleanDirectory)
    : ColumnDataWriter::ColumnDataWriter(rDirectory, rBaseName, cleanDirectory),
      mConcentrated(NULL)
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
    if (mConcentrated==NULL)
    {
        VecScatterCreateToZero(petscVector, &mToMaster, &mConcentrated);
    }

//    int size2;
//    VecGetSize(mConcentrated, &size2);
//    std::cout << "Vector size=" << size << "," << size2 << std::endl << std::flush;

//PETSc-3.x.x or PETSc-2.3.3
#if ( (PETSC_VERSION_MAJOR == 3) || (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR == 3 && PETSC_VERSION_SUBMINOR == 3)) //2.3.3 or 3.x.x
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
    VecDestroy(unstriped_petsc);
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
    if (mConcentrated != NULL)
    {
        VecScatterDestroy(mToMaster);
        VecDestroy(mConcentrated);
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
