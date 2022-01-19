/*

Copyright (c) 2005-2022, University of Oxford.
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
#include "DistributedBoxCollection.hpp"
#include "Exception.hpp"
#include "MathsCustomFunctions.hpp"
#include "Warnings.hpp"

// Static member for "fudge factor" is instantiated here
template<unsigned DIM>
const double DistributedBoxCollection<DIM>::msFudge = 5e-14;

template<unsigned DIM>
DistributedBoxCollection<DIM>::DistributedBoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize, bool isPeriodicInX, bool isPeriodicInY, bool isPeriodicInZ, int localRows)
    : mBoxWidth(boxWidth),
      mIsPeriodicInX(isPeriodicInX),
      mIsPeriodicInY(isPeriodicInY),
      mIsPeriodicInZ(isPeriodicInZ),
      mAreLocalBoxesSet(false),
      mCalculateNodeNeighbours(true)
{
    // If the domain size is not 'divisible' (i.e. fmod(width, box_size) > 0.0) we swell the domain to enforce this.
    for (unsigned i=0; i<DIM; i++)
    {
        double r = fmod((domainSize[2*i+1]-domainSize[2*i]), boxWidth);
        if (r > 0.0)
        {
            domainSize[2*i+1] += boxWidth - r;
        }
    }

    mDomainSize = domainSize;

    // Set the periodicity across procs flag
    mIsPeriodicAcrossProcs = (DIM==1 && mIsPeriodicInX) || (DIM==2 && mIsPeriodicInY) || (DIM==3 && mIsPeriodicInZ);

    // Calculate the number of boxes in each direction.
    mNumBoxesEachDirection = scalar_vector<unsigned>(DIM, 0u);

    for (unsigned i=0; i<DIM; i++)
    {
        double counter = mDomainSize(2*i);
        while (counter + msFudge < mDomainSize(2*i+1))
        {
            mNumBoxesEachDirection(i)++;
            counter += mBoxWidth;
        }
    }

    // Make sure there are enough boxes for the number of processes.
    if (mNumBoxesEachDirection(DIM-1) < PetscTools::GetNumProcs())
    {
        WARNING("There are more processes than convenient for the domain/mesh/box size.  The domain size has been swollen.")
        mDomainSize[2*DIM - 1] += (PetscTools::GetNumProcs() - mNumBoxesEachDirection(DIM-1))*mBoxWidth;
        mNumBoxesEachDirection(DIM-1) = PetscTools::GetNumProcs();
    }

    // Make a distributed vector factory to split the rows of boxes between processes.
    mpDistributedBoxStackFactory = new DistributedVectorFactory(mNumBoxesEachDirection(DIM-1), localRows);

    // Calculate how many boxes in a row / face. A useful piece of data in the class.
    mNumBoxes = 1u;
    for (unsigned dim=0; dim<DIM; dim++)
    {
        mNumBoxes *= mNumBoxesEachDirection(dim);
    }

    mNumBoxesInAFace = mNumBoxes / mNumBoxesEachDirection(DIM-1);

    unsigned num_local_boxes = mNumBoxesInAFace * GetNumLocalRows();

    mMinBoxIndex = mpDistributedBoxStackFactory->GetLow() * mNumBoxesInAFace;
    mMaxBoxIndex = mpDistributedBoxStackFactory->GetHigh() * mNumBoxesInAFace - 1;

    // Create the correct number of boxes and set up halos
    mBoxes.resize(num_local_boxes);
    SetupHaloBoxes();
}

template<unsigned DIM>
DistributedBoxCollection<DIM>::~DistributedBoxCollection()
{
    delete mpDistributedBoxStackFactory;
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::EmptyBoxes()
{
    for (unsigned i=0; i<mBoxes.size(); i++)
    {
        mBoxes[i].ClearNodes();
    }
    for (unsigned i=0; i<mHaloBoxes.size(); i++)
    {
        mHaloBoxes[i].ClearNodes();
    }
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::SetupHaloBoxes()
{
    // We don't need to do this if not parallel
    if ( PetscTools::GetNumProcs() == 1 )
    {
        return;
    }

    // Get top-most and bottom-most value of Distributed Box Stack.
    unsigned hi = mpDistributedBoxStackFactory->GetHigh();
    unsigned lo = mpDistributedBoxStackFactory->GetLow();

    // If I am not the top-most process, add halo structures above.
    if (!PetscTools::AmTopMost())
    {
        for (unsigned i=0; i < mNumBoxesInAFace; i++)
        {
            Box<DIM> new_box;
            mHaloBoxes.push_back(new_box);

            unsigned global_index = hi * mNumBoxesInAFace + i;
            mHaloBoxesMapping[global_index] = mHaloBoxes.size()-1;
            mHalosRight.push_back(global_index - mNumBoxesInAFace);
        }
    }
    // Otherwise if I am the top most and periodic in y (2d) or z (3d) add halo boxes for
    // the base process
    else if ( mIsPeriodicAcrossProcs )
    {
        for (unsigned i=0; i < mNumBoxesInAFace; i++)
        {
            Box<DIM> new_box;
            mHaloBoxes.push_back(new_box);

            mHaloBoxesMapping[i] = mHaloBoxes.size()-1;

            mHalosRight.push_back( (hi-1)*mNumBoxesInAFace + i );
        }
    }

    // If I am not the bottom-most process, add halo structures below.
    if (!PetscTools::AmMaster())
    {
        for (unsigned i=0; i< mNumBoxesInAFace; i++)
        {
            Box<DIM> new_box;
            mHaloBoxes.push_back(new_box);

            unsigned global_index = (lo - 1) * mNumBoxesInAFace + i;
            mHaloBoxesMapping[global_index] = mHaloBoxes.size() - 1;

            mHalosLeft.push_back(global_index  + mNumBoxesInAFace);
        }
    }
    // Otherwise if I am the bottom most and periodic in y (2d) or z (3d) add halo boxes for
    // the top process
    else if ( mIsPeriodicAcrossProcs )
    {
        for (unsigned i=0; i < mNumBoxesInAFace; i++)
        {
            Box<DIM> new_box;
            mHaloBoxes.push_back(new_box);

            unsigned global_index = (mNumBoxesEachDirection(DIM-1) - 1) * mNumBoxesInAFace + i;
            mHaloBoxesMapping[global_index] = mHaloBoxes.size()-1;

            mHalosLeft.push_back( i );
        }

    }
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::UpdateHaloBoxes()
{
    mHaloNodesLeft.clear();
    for (unsigned i=0; i<mHalosLeft.size(); i++)
    {
        for (typename std::set<Node<DIM>* >::iterator iter=this->rGetBox(mHalosLeft[i]).rGetNodesContained().begin();
                iter!=this->rGetBox(mHalosLeft[i]).rGetNodesContained().end();
                iter++)
        {
            mHaloNodesLeft.push_back((*iter)->GetIndex());
        }
    }

    // Send right
    mHaloNodesRight.clear();
    for (unsigned i=0; i<mHalosRight.size(); i++)
    {
        for (typename std::set<Node<DIM>* >::iterator iter=this->rGetBox(mHalosRight[i]).rGetNodesContained().begin();
                iter!=this->rGetBox(mHalosRight[i]).rGetNodesContained().end();
                iter++)
        {
            mHaloNodesRight.push_back((*iter)->GetIndex());
        }
    }
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::GetNumLocalRows() const
{
    return mpDistributedBoxStackFactory->GetHigh() - mpDistributedBoxStackFactory->GetLow();
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::IsBoxOwned(unsigned globalIndex)
{
    return (!(globalIndex<mMinBoxIndex) && !(mMaxBoxIndex<globalIndex));
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::IsHaloBox(unsigned globalIndex)
{
    bool is_halo_right = ((globalIndex > mMaxBoxIndex) && !(globalIndex > mMaxBoxIndex + mNumBoxesInAFace));
    bool is_halo_left = ((globalIndex < mMinBoxIndex) && !(globalIndex < mMinBoxIndex - mNumBoxesInAFace));

    // Also need to check for periodic boxes
    if ( mIsPeriodicAcrossProcs )
    {
        if ( PetscTools::AmTopMost() )
        {
            is_halo_right = (globalIndex < mNumBoxesInAFace);
        }
        if ( PetscTools::AmMaster() )
        {
            is_halo_left = (!(globalIndex < (mNumBoxesEachDirection[DIM-1]-1)*mNumBoxesInAFace) && (globalIndex < mNumBoxesEachDirection[DIM-1]*mNumBoxesInAFace ));
        }
    }

    return (PetscTools::IsParallel() && (is_halo_right || is_halo_left));
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::IsInteriorBox(unsigned globalIndex)
{
    bool is_on_boundary = !(globalIndex < mMaxBoxIndex - mNumBoxesInAFace) || (globalIndex < mMinBoxIndex + mNumBoxesInAFace);

    return (PetscTools::IsSequential() || !(is_on_boundary));
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::CalculateGlobalIndex(c_vector<unsigned, DIM> gridIndices)
{
    ///\todo #2308 etc. We need to make allowance for periodicity here...
    unsigned global_index;

    switch (DIM)
    {
        case 1:
        {
            global_index = gridIndices(0);
            break;
        }
        case 2:
        {
            global_index = gridIndices(0) +
                           gridIndices(1) * mNumBoxesEachDirection(0);
            break;
        }
        case 3:
        {
            global_index = gridIndices(0) +
                           gridIndices(1) * mNumBoxesEachDirection(0) +
                           gridIndices(2) * mNumBoxesEachDirection(0) * mNumBoxesEachDirection(1);
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }
    return global_index;
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::CalculateContainingBox(Node<DIM>* pNode)
{
    // Get the location of the node
    c_vector<double, DIM> location = pNode->rGetLocation();
    return CalculateContainingBox(location);
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::CalculateContainingBox(c_vector<double, DIM>& rLocation)
{
    // The node must lie inside the boundary of the box collection
    for (unsigned i=0; i<DIM; i++)
    {
        if ((rLocation[i] < mDomainSize(2*i)) || !(rLocation[i] < mDomainSize(2*i+1)))
        {
            EXCEPTION("The point provided is outside all of the boxes");
        }
    }

    // Compute the containing box index in each dimension
    c_vector<unsigned, DIM> containing_box_indices = scalar_vector<unsigned>(DIM, 0u);
    for (unsigned i=0; i<DIM; i++)
    {
        double box_counter = mDomainSize(2*i);
        while (!((box_counter + mBoxWidth) > rLocation[i] + msFudge))
        {
            containing_box_indices[i]++;
            box_counter += mBoxWidth;
        }
    }

    // Use these to compute the index of the containing box
    unsigned containing_box_index = 0;
    for (unsigned i=0; i<DIM; i++)
    {
        unsigned temp = 1;
        for (unsigned j=0; j<i; j++)
        {
            temp *= mNumBoxesEachDirection(j);
        }
        containing_box_index += temp*containing_box_indices[i];
    }

    // This index must be less than the total number of boxes
    assert(containing_box_index < mNumBoxes);

    return containing_box_index;
}

template<unsigned DIM>
c_vector<unsigned, DIM> DistributedBoxCollection<DIM>::CalculateGridIndices(unsigned globalIndex)
{
    c_vector<unsigned, DIM> grid_indices;

    switch (DIM)
    {
        case 1:
        {
            grid_indices(0) = globalIndex;
            break;
        }
        case 2:
        {
            unsigned num_x = mNumBoxesEachDirection(0);
            grid_indices(0) = globalIndex % num_x;
            grid_indices(1) = (globalIndex - grid_indices(0)) / num_x;
            break;
        }
        case 3:
        {
            unsigned num_x = mNumBoxesEachDirection(0);
            unsigned num_xy = mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1);
            grid_indices(0) = globalIndex % num_x;
            grid_indices(1) = (globalIndex % num_xy - grid_indices(0)) / num_x;
            grid_indices(2) = globalIndex / num_xy;
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    return grid_indices;
}

template<unsigned DIM>
Box<DIM>& DistributedBoxCollection<DIM>::rGetBox(unsigned boxIndex)
{
    // Check first for local ownership
    if (!(boxIndex<mMinBoxIndex) && !(mMaxBoxIndex<boxIndex))
    {
        return mBoxes[boxIndex-mMinBoxIndex];
    }

    // If normal execution reaches this point then the box does not belong to the process so we will check for a halo box
    return rGetHaloBox(boxIndex);
}

template<unsigned DIM>
Box<DIM>& DistributedBoxCollection<DIM>::rGetHaloBox(unsigned boxIndex)
{
    assert(IsHaloBox(boxIndex));

    unsigned local_index = mHaloBoxesMapping.find(boxIndex)->second;

    return mHaloBoxes[local_index];
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::GetNumBoxes()
{
    return mNumBoxes;
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::GetNumLocalBoxes()
{
    return mBoxes.size();
}

template<unsigned DIM>
c_vector<double, 2*DIM> DistributedBoxCollection<DIM>::rGetDomainSize() const
{
    return mDomainSize;
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::GetAreLocalBoxesSet() const
{
    return mAreLocalBoxesSet;
}

template<unsigned DIM>
double DistributedBoxCollection<DIM>::GetBoxWidth() const
{
    return mBoxWidth;
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::GetIsPeriodicInX() const
{
    return mIsPeriodicInX;
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::GetIsPeriodicInY() const
{
    return mIsPeriodicInY;
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::GetIsPeriodicInZ() const
{
    return mIsPeriodicInZ;
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::GetIsPeriodicAcrossProcs() const
{
    return mIsPeriodicAcrossProcs;
}

template<unsigned DIM>
c_vector<bool,DIM> DistributedBoxCollection<DIM>::GetIsPeriodicAllDims() const
{
    c_vector<bool, DIM> periodic_dims;
    periodic_dims(0) = mIsPeriodicInX;
    if (DIM > 1)
    {
        periodic_dims(1) = mIsPeriodicInY;
    }
    if (DIM>2)
    {
        periodic_dims(2) = mIsPeriodicInZ;
    }
    return periodic_dims;
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::GetNumRowsOfBoxes() const
{
    return mpDistributedBoxStackFactory->GetHigh() - mpDistributedBoxStackFactory->GetLow();
}

template<unsigned DIM>
int DistributedBoxCollection<DIM>::LoadBalance(std::vector<int> localDistribution)
{
    MPI_Status status;

    int proc_right = (PetscTools::AmTopMost()) ? MPI_PROC_NULL : (int)PetscTools::GetMyRank() + 1;
    int proc_left = (PetscTools::AmMaster()) ? MPI_PROC_NULL : (int)PetscTools::GetMyRank() - 1;

    // A variable that will return the new number of rows.
    int new_rows = localDistribution.size();

    /**
     * Shift information on distribution of nodes to the right, so processes can manage their left/bottom/back boundary (1d/2d/3d)
     */
    unsigned rows_on_left_process = 0;
    std::vector<int> node_distr_on_left_process;

    unsigned num_local_rows = localDistribution.size();

    MPI_Send(&num_local_rows, 1, MPI_UNSIGNED, proc_right, 123, PETSC_COMM_WORLD);
    MPI_Recv(&rows_on_left_process, 1, MPI_UNSIGNED, proc_left, 123, PETSC_COMM_WORLD, &status);

    node_distr_on_left_process.resize(rows_on_left_process > 0 ? rows_on_left_process : 1);

    MPI_Send(&localDistribution[0], num_local_rows, MPI_INT, proc_right, 123, PETSC_COMM_WORLD);
    MPI_Recv(&node_distr_on_left_process[0], rows_on_left_process, MPI_INT, proc_left, 123, PETSC_COMM_WORLD, &status);

    /**
     * Calculate change in balance of loads by shifting the left/bottom boundary in either direction
     */
    int local_load = 0;
    for (unsigned i=0; i<localDistribution.size(); i++)
    {
        local_load += localDistribution[i];
    }
    int load_on_left_proc = 0;
    for (unsigned i=0; i<node_distr_on_left_process.size(); i++)
    {
        load_on_left_proc += node_distr_on_left_process[i];
    }

    if (!PetscTools::AmMaster())
    {
        // Calculate (Difference in load with a shift) - (Difference in current loads) for a move left and right of the boundary
        // This code uses integer arithmetic in order to avoid the rounding errors associated with doubles
        int local_to_left_sq = (local_load - load_on_left_proc) * (local_load - load_on_left_proc);
        int delta_left =  ( (local_load + node_distr_on_left_process[node_distr_on_left_process.size() - 1]) - (load_on_left_proc - node_distr_on_left_process[node_distr_on_left_process.size() - 1]) );
        delta_left = delta_left*delta_left - local_to_left_sq;

        int delta_right = ( (local_load - localDistribution[0]) - (load_on_left_proc + localDistribution[0]));
        delta_right = delta_right*delta_right - local_to_left_sq;

        // If a delta is negative we should accept that change. If both are negative choose the largest change.
        int local_change = 0;
        bool move_left = (!(delta_left > 0) && (node_distr_on_left_process.size() > 1));
        if (move_left)
        {
            local_change = 1;
        }

        bool move_right = !(delta_right > 0) && (localDistribution.size() > 2);
        if (move_right)
        {
            local_change = -1;
        }

        if (move_left && move_right)
        {
            local_change = (fabs((double)delta_right) > fabs((double)delta_left)) ? -1 : 1;
        }

        // Update the number of local rows.
        new_rows += local_change;

        // Send the result of the calculation back to the left processes.
        MPI_Send(&local_change, 1, MPI_INT, proc_left, 123, PETSC_COMM_WORLD);
    }

    // Receive changes from right hand process.
    int remote_change = 0;
    MPI_Recv(&remote_change, 1, MPI_INT, proc_right, 123, PETSC_COMM_WORLD, &status);

    // Update based on change or right/top boundary
    new_rows -= remote_change;

    return new_rows;
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::SetupLocalBoxesHalfOnly()
{
    if (mAreLocalBoxesSet)
    {
        EXCEPTION("Local Boxes Are Already Set");
    }
    else
    {
        switch (DIM)
        {
            case 1:
            {
                // We only need to look for neighbours in the current and successive boxes plus some others for halos
                mLocalBoxes.clear();

                // Iterate over the global box indices
                for (unsigned global_index = mMinBoxIndex; global_index<mMaxBoxIndex+1; global_index++)
                {
                    std::set<unsigned> local_boxes;

                    // Insert the current box
                    local_boxes.insert(global_index);

                    // Set some bools to find out where we are
                    bool right = (global_index==mNumBoxesEachDirection(0)-1);
                    bool left = (global_index == 0);
                    bool proc_left = (global_index == mpDistributedBoxStackFactory->GetLow());

                    // If we're not at the right-most box, then insert the box to the right
                    if (!right)
                    {
                        local_boxes.insert(global_index+1);
                    }
                    // If we're on a left process boundary and not on process 0, insert the (halo) box to the left
                    if (proc_left && !left)
                    {
                        local_boxes.insert(global_index-1);
                    }

                    mLocalBoxes.push_back(local_boxes);
                }
                break;
            }
            case 2:
            {
                // Clear the local boxes
                mLocalBoxes.clear();

                // Now we need to work out the min and max y indices
                unsigned j_start = CalculateGridIndices(mMinBoxIndex)(1);
                unsigned j_end = CalculateGridIndices(mMaxBoxIndex)(1);
                // Determine the number of boxes in each direction
                unsigned nI = mNumBoxesEachDirection(0);
                unsigned nJ = mNumBoxesEachDirection(1);
                // Locally store the bottom processor row
                unsigned bottom_proc = mpDistributedBoxStackFactory->GetLow();
                unsigned top_proc = mpDistributedBoxStackFactory->GetHigh();

                // Outer loop: over y
                for ( unsigned j = j_start; j < (j_end+1); j++ )
                {
                    // Inner loop: over x
                    for ( unsigned i = 0; i < nI; i++ )
                    {
                        std::set<unsigned> local_boxes;
                        /* We want to add (for non-boundary)
                        (i-1,j+1) (i,j+1) (i+1,j+1)
                                    (i,j)   (i+1,j)   */

                        int j_mod = j; // This is used to account for y periodic boundaries
                        // If we are on the bottom of the processor, we may need to add the row below
                        int dj = -1 * (int)(j == bottom_proc && (j > 0 || (mIsPeriodicInY && top_proc < nJ)) );
                        int j_mod_2 = 0;
                        // The min ensures we don't go above the top boundary
                        for (; dj < std::min((int)nJ-j_mod,(int)2); dj++ )
                        {
                            // We need to change to the top row if we are in the condition where dj == -1 and j == 0 which is only
                            // when we are periodic in y and the top row is on a different processor to the bottom row
                            if ( mIsPeriodicInY && j == 0 && dj < 1 && top_proc < nJ )
                            {
                                j_mod_2 = ( dj < 0 ) ? nJ : 0;
                            }

                            // The -1*dj ensures we get the upper left, the max ensures we don't hit the left boundary,
                            // the min ensures we don't hit the right boundary
                            int boxi = std::max((int) i-1*std::abs(dj),(int) 0);
                            for ( ; boxi < std::min((int)i+2,(int)nI); boxi++ )
                            {
                                local_boxes.insert( (j_mod+dj+j_mod_2)*nI + boxi );
                            }
                            // Add in the x periodicity at the right boundary, we then want: (0,j) and (0,j+1)
                            if ( i==(nI-1) && mIsPeriodicInX )
                            {
                                local_boxes.insert( (j_mod+dj+j_mod_2)*nI );
                            }
                            // If the y boundary is periodic, the new j level is the bottom row; we use jmod to adjust
                            // for this to adjust for dj = 1
                            if ( mIsPeriodicInY && j == (nJ-1) && dj == 0 )
                            {
                                j_mod = -1;
                            }
                        }
                        // Now add the left-upper box if x periodic and on the boundary
                        if ( mIsPeriodicInX && i == 0 )
                        {
                            if ( j < (nJ-1) )
                            {
                                local_boxes.insert( (j+2)*nI - 1 );
                                if ( j==0 && mIsPeriodicInY && top_proc < nJ )
                                {
                                    local_boxes.insert( nI*nJ - 1 );
                                }
                            }
                            else if ( mIsPeriodicInY )
                            {
                                local_boxes.insert( nI-1 );
                            }
                            // If wee are on the bottom of the process need to add
                            // the box on the right in the row below
                            if ( j == bottom_proc && j > 0 )
                            {
                                local_boxes.insert( j*nI - 1 );
                            }
                        }

                        // Add to the local boxes
                        mLocalBoxes.push_back(local_boxes);
                    }
                }

                break;
            }
            case 3:
            {
                // Clear the local boxes
                mLocalBoxes.clear();

                // Now we need to work out the min and max z indices
                unsigned k_start = CalculateGridIndices(mMinBoxIndex)(2);
                unsigned k_end = CalculateGridIndices(mMaxBoxIndex)(2);

                // Determine the number of boxes in each direction
                unsigned nI = mNumBoxesEachDirection(0);
                unsigned nJ = mNumBoxesEachDirection(1);
                unsigned nK = mNumBoxesEachDirection(2);

                // Work out the bottom/top processor
                unsigned bottom_proc = mpDistributedBoxStackFactory->GetLow();
                unsigned top_proc = mpDistributedBoxStackFactory->GetHigh();

                // Outer loop: over z
                for ( unsigned k = k_start; k <= k_end; k++ )
                {
                    // Middle loop: over y
                    for ( unsigned j = 0; j < nJ; j++ )
                    {
                        // Inner loop: over x
                        for ( unsigned i = 0; i < nI; i++ )
                        {
                            std::set<unsigned> local_boxes;

                            // Same z level
                            unsigned z_offset = k*nI*nJ;
                            // (See case dim=2 for commented code of X-Y implementation)
                            int j_mod = (int)j;
                            for ( int dj = 0; dj < std::min((int)nJ-j_mod,2); dj++ )
                            {
                                for ( int boxi = std::max((int)i-1*dj,0);
                                            boxi < std::min((int)i+2,(int)nI); boxi++ )
                                {
                                    local_boxes.insert( z_offset + (j_mod+dj)*nI + boxi );
                                }
                                if ( i==(nI-1) && mIsPeriodicInX )
                                {
                                    local_boxes.insert( z_offset + (j_mod+dj)*nI );
                                }
                                if ( mIsPeriodicInY && j == (nJ-1) )
                                {
                                    j_mod = -1;
                                }
                            }
                            // Now add the left-upper box if x periodic and on the boundary
                            if ( mIsPeriodicInX && i == 0 )
                            {
                                if ( j < (nJ-1) )
                                {
                                    local_boxes.insert( z_offset + (j+2)*nI - 1 );
                                }
                                else if ( mIsPeriodicInY )
                                {
                                    local_boxes.insert( z_offset + nI-1 );
                                }
                            }

                            // Add all the surrounding boxes from the z level above if required
                            std::vector<unsigned> k_offset;
                            if ( k < (nK-1) )
                            {
                                k_offset.push_back(k+1);
                            }
                            if ( k == bottom_proc && k > 0 )
                            {
                                // Need the level below if we are on the lowest processor
                                k_offset.push_back(k-1);
                            }
                            if ( mIsPeriodicInZ && k == (nK-1) )
                            {
                                k_offset.push_back(0);
                            }
                            else if ( mIsPeriodicInZ && k==0 && (top_proc < nK) )
                            {
                                k_offset.push_back(nK-1);
                            }

                            for ( std::vector<unsigned>::iterator k_offset_it = k_offset.begin(); k_offset_it != k_offset.end(); ++k_offset_it )
                            {
                                z_offset = (*k_offset_it)*nI*nJ;
                                // Periodicity adjustments
                                int pX = (int) mIsPeriodicInX;
                                int pY = (int) mIsPeriodicInY;
                                for ( int boxi = std::max((int)i-1,-1*pX); boxi < std::min((int)i+2,(int)nI+pX); boxi++ )
                                {
                                    for ( int boxj = std::max((int)j-1,-1*pY); boxj < std::min((int)j+2,(int)nJ+pY); boxj++ )
                                    {
                                        int box_to_add = z_offset + (boxj*(int)nI) + boxi;
                                        // Check for x periodicity
                                        // i==(nI-1) only when we are periodic and on the right,
                                        // i==-1 only when we are periodic and on the left
                                        box_to_add += ( (unsigned)(boxi==-1) - (unsigned)(boxi==(int)nI) )*nI;
                                        // Check for y periodicity
                                        // j==nJ only when we are periodic and on the right,
                                        // j==-1 only when we are periodic and on the left
                                        box_to_add += ( (unsigned)(boxj==-1) - (unsigned)(boxj==(int)nJ) )*nI*nJ;
                                        local_boxes.insert( box_to_add );
                                    }
                                }
                            }

                            // Add to the local boxes
                            mLocalBoxes.push_back(local_boxes);
                        }
                    }
                }
                break;
            }
            default:
                NEVER_REACHED;
        }
        mAreLocalBoxesSet=true;
    }
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::SetupAllLocalBoxes()
{
    mAreLocalBoxesSet = true;
    switch (DIM)
    {
        case 1:
        {
            for (unsigned i=mMinBoxIndex; i<mMaxBoxIndex+1; i++)
            {
                std::set<unsigned> local_boxes;

                local_boxes.insert(i);

                // add the two neighbours
                if (i!=0)
                {
                    local_boxes.insert(i-1);
                }
                if (i+1 != mNumBoxesEachDirection(0))
                {
                    local_boxes.insert(i+1);
                }

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        case 2:
        {
            mLocalBoxes.clear();

            unsigned M = mNumBoxesEachDirection(0);
            unsigned N = mNumBoxesEachDirection(1);

            std::vector<bool> is_xmin(N*M); // far left
            std::vector<bool> is_xmax(N*M); // far right
            std::vector<bool> is_ymin(N*M); // bottom
            std::vector<bool> is_ymax(N*M); // top

            for (unsigned i=0; i<M*N; i++)
            {
                is_xmin[i] = (i%M==0);
                is_xmax[i] = ((i+1)%M==0);
                is_ymin[i] = (i%(M*N)<M);
                is_ymax[i] = (i%(M*N)>=(N-1)*M);
            }

            for (unsigned i=mMinBoxIndex; i<mMaxBoxIndex+1; i++)
            {
                std::set<unsigned> local_boxes;

                local_boxes.insert(i);

                // add the box to the left
                if (!is_xmin[i])
                {
                    local_boxes.insert(i-1);
                }
                else // Add Periodic Box if needed
                {
                    if (mIsPeriodicInX)
                    {
                        local_boxes.insert(i+M-1);
                    }
                }

                // add the box to the right
                if (!is_xmax[i])
                {
                    local_boxes.insert(i+1);
                }
                else // Add Periodic Box if needed
                {
                    if (mIsPeriodicInX)
                    {
                        local_boxes.insert(i-M+1);
                    }
                }

                // add the one below
                if (!is_ymin[i])
                {
                    local_boxes.insert(i-M);
                }
                else // Add periodic box if needed
                {
                    if (mIsPeriodicInY)
                    {
                        local_boxes.insert(i+(N-1)*M);
                    }
                }

                // add the one above
                if (!is_ymax[i])
                {
                    local_boxes.insert(i+M);
                }
                else // Add periodic box if needed
                {
                    if (mIsPeriodicInY)
                    {
                        local_boxes.insert(i-(N-1)*M);
                    }
                }

                // add the four corner boxes

                if ( (!is_xmin[i]) && (!is_ymin[i]) )
                {
                    local_boxes.insert(i-1-M);
                }
                if ( (!is_xmin[i]) && (!is_ymax[i]) )
                {
                    local_boxes.insert(i-1+M);
                }
                if ( (!is_xmax[i]) && (!is_ymin[i]) )
                {
                    local_boxes.insert(i+1-M);
                }
                if ( (!is_xmax[i]) && (!is_ymax[i]) )
                {
                    local_boxes.insert(i+1+M);
                }

                // Add periodic corner boxes if needed
                if (mIsPeriodicInX)
                {
                    if ( (is_xmin[i]) && (!is_ymin[i]) )
                    {
                        local_boxes.insert(i-1);
                    }
                    if ( (is_xmin[i]) && (!is_ymax[i]) )
                    {
                        local_boxes.insert(i-1+2*M);
                    }
                    if ( (is_xmax[i]) && (!is_ymin[i]) )
                    {
                        local_boxes.insert(i+1-2*M);
                    }
                    if ( (is_xmax[i]) && (!is_ymax[i]) )
                    {
                        local_boxes.insert(i+1);
                    }
                }
                if(mIsPeriodicInY)
                {
                    if( (is_ymin[i]) && !(is_xmin[i]) )
                    {
                        local_boxes.insert(i+(N-1)*M-1);
                    }
                    if( (is_ymin[i]) && !(is_xmax[i]) )
                    {
                        local_boxes.insert(i+(N-1)*M+1);
                    }
                    if( (is_ymax[i]) && !(is_xmin[i]) )
                    {
                        local_boxes.insert(i-(N-1)*M-1);
                    }
                    if( (is_ymax[i]) && !(is_xmax[i]) )
                    {
                        local_boxes.insert(i-(N-1)*M+1);
                    }
                }
                if(mIsPeriodicInX && mIsPeriodicInY)
                {
                    if( i==0 ) // Lower left corner
                    {
                        local_boxes.insert(M*N-1); // Add upper right corner
                    }
                    else if( i==(M-1) ) // Lower right corner
                    {
                        local_boxes.insert(M*(N-1)); // Add upper left corner
                    }
                    else if( i==(M*(N-1)) ) // Upper left corner
                    {
                        local_boxes.insert(M-1); // Add lower right corner
                    }
                    else if( i==(M*N-1) ) // Upper right corner
                    {
                        local_boxes.insert(0); // Lower left corner
                    }
                }

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        case 3:
        {
            mLocalBoxes.clear();

            unsigned M = mNumBoxesEachDirection(0);
            unsigned N = mNumBoxesEachDirection(1);
            unsigned P = mNumBoxesEachDirection(2);

            std::vector<bool> is_xmin(N*M*P); // far left
            std::vector<bool> is_xmax(N*M*P); // far right
            std::vector<bool> is_ymin(N*M*P); // nearest
            std::vector<bool> is_ymax(N*M*P); // furthest
            std::vector<bool> is_zmin(N*M*P); // bottom layer
            std::vector<bool> is_zmax(N*M*P); // top layer

            for (unsigned i=0; i<M*N*P; i++)
            {
                is_xmin[i] = (i%M==0);
                is_xmax[i] = ((i+1)%M==0);
                is_ymin[i] = (i%(M*N)<M);
                is_ymax[i] = (i%(M*N)>=(N-1)*M);
                is_zmin[i] = (i<M*N);
                is_zmax[i] = (i>=M*N*(P-1));
            }

            for (unsigned i=mMinBoxIndex; i<mMaxBoxIndex+1; i++)
            {
                std::set<unsigned> local_boxes;

                // add itself as a local box
                local_boxes.insert(i);

                // now add all 26 other neighbours.....

                // add the box left
                if (!is_xmin[i])
                {
                    local_boxes.insert(i-1);

                    // plus some others towards the left
                    if (!is_ymin[i])
                    {
                        local_boxes.insert(i-1-M);
                    }

                    if (!is_ymax[i])
                    {
                        local_boxes.insert(i-1+M);
                    }

                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i-1-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i-1+M*N);
                    }
                }

                // add the box to the right
                if (!is_xmax[i])
                {
                    local_boxes.insert(i+1);

                    // plus some others towards the right
                    if (!is_ymin[i])
                    {
                        local_boxes.insert(i+1-M);
                    }

                    if (!is_ymax[i])
                    {
                        local_boxes.insert(i+1+M);
                    }

                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i+1-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i+1+M*N);
                    }
                }

                // add the boxes next along the y axis
                if (!is_ymin[i])
                {
                    local_boxes.insert(i-M);

                    // and more in this plane
                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i-M-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i-M+M*N);
                    }
                }

                // add the boxes next along the y axis
                if (!is_ymax[i])
                {
                    local_boxes.insert(i+M);

                    // and more in this plane
                    if (!is_zmin[i])
                    {
                        local_boxes.insert(i+M-M*N);
                    }

                    if (!is_zmax[i])
                    {
                        local_boxes.insert(i+M+M*N);
                    }
                }

                // add the box directly above
                if (!is_zmin[i])
                {
                    local_boxes.insert(i-N*M);
                }

                // add the box directly below
                if (!is_zmax[i])
                {
                    local_boxes.insert(i+N*M);
                }

                // finally, the 8 corners are left

                if ( (!is_xmin[i]) && (!is_ymin[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i-1-M-M*N);
                }

                if ( (!is_xmin[i]) && (!is_ymin[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i-1-M+M*N);
                }

                if ( (!is_xmin[i]) && (!is_ymax[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i-1+M-M*N);
                }

                if ( (!is_xmin[i]) && (!is_ymax[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i-1+M+M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymin[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i+1-M-M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymin[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i+1-M+M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymax[i]) && (!is_zmin[i]) )
                {
                    local_boxes.insert(i+1+M-M*N);
                }

                if ( (!is_xmax[i]) && (!is_ymax[i]) && (!is_zmax[i]) )
                {
                    local_boxes.insert(i+1+M+M*N);
                }

                // Now add the periodic boxes if any periodicity
                if (mIsPeriodicInX && ( is_xmin[i] || is_xmax[i]) )
                {
                    // We are repeating the same steps on each z level, so easiest is
                    // to make a vector of the z levels we want
                    std::vector< int > z_i_offsets(1,0); // Add the current z box level
                    if ( !is_zmin[i] )
                    {
                        z_i_offsets.push_back(-M*N); // Add the z box level below
                    }
                    if ( !is_zmax[i] )
                    {
                        z_i_offsets.push_back(M*N); // Add the z box level above
                    }

                    // If we are on the left, add the nine on the right
                    if ( is_xmin[i] )
                    {
                        // Loop over the z levels
                        for ( std::vector<int>::iterator it = z_i_offsets.begin(); it != z_i_offsets.end(); it++ )
                        {
                            local_boxes.insert( i + (*it) + (M-1) ); // The right-most box on the same row
                            // We also need to check for y boundaries
                            if ( !is_ymin[i] )
                            {
                                local_boxes.insert( i + (*it) - 1 ); // The right-most box one row below
                            }
                            if ( !is_ymax[i] )
                            {
                                local_boxes.insert( i + (*it) + (2*M-1) ); // The right-most box one row above
                            }
                        }
                    }

                    // If we are on the right, add the nine on the left
                    else if ( is_xmax[i] )
                    {
                        // Loop over the z levels
                        for ( std::vector<int>::iterator it = z_i_offsets.begin(); it != z_i_offsets.end(); it++ )
                        {
                            local_boxes.insert( i + (*it) - (M-1) ); // The left-most box on the same row
                            // We also need to check for y boundaries
                            if ( !is_ymin[i] )
                            {
                                local_boxes.insert( i + (*it) - (2*M-1) ); // The left-most box one row below
                            }
                            if ( !is_ymax[i] )
                            {
                                local_boxes.insert( i + (*it) + 1 ); // The left-most box one row below
                            }
                        }
                    }
                }

                if ( mIsPeriodicInY && (is_ymax[i] || is_ymin[i]) )
                {
                    // We consider the current and upper z level and create a vector of the
                    // opposite box indices that need to be added
                    std::vector<unsigned> opp_box_i(0);
                    if ( is_ymin[i] )
                    {
                        opp_box_i.push_back(i + (N-1)*M); // Current z level
                        if ( !is_zmin[i] )
                        {
                            opp_box_i.push_back( i - M ); // z level below
                        }
                        if ( !is_zmax[i] )
                        {
                            opp_box_i.push_back(i + 2*M*N - M); // z level above
                        }
                    }
                    else if ( is_ymax[i] )
                    {
                        opp_box_i.push_back( i - (N-1)*M ); // Current z level
                        if ( !is_zmin[i] )
                        {
                            opp_box_i.push_back( i - 2*M*N + M ); // z level below
                        }
                        if ( !is_zmax[i] )
                        {
                            opp_box_i.push_back( i + M ); // z level above
                        }
                    }

                    // Now we add the different boxes, checking for left and right
                    for ( std::vector<unsigned>::iterator it_opp_box = opp_box_i.begin(); it_opp_box != opp_box_i.end(); it_opp_box++ )
                    {
                        local_boxes.insert( *it_opp_box );
                        if ( !is_xmin[i] )
                        {
                            local_boxes.insert( *it_opp_box - 1 );
                        }
                        if ( !is_xmax[i] )
                        {
                            local_boxes.insert( *it_opp_box + 1 );
                        }
                    }
                }

                if ( mIsPeriodicInX && mIsPeriodicInY )
                {
                    // Need to add the corners
                    if ( is_xmin[i] && is_ymin[i] )
                    {
                        // Current z level
                        local_boxes.insert(i+M*N-1);
                        if ( !is_zmax[i] )
                        {
                            // Upper z level
                            local_boxes.insert(i+2*M*N-1);
                        }
                        if ( !is_zmin[i] )
                        {
                            // Lower z level
                            local_boxes.insert(i-1);
                        }
                    }
                    if ( is_xmax[i] && is_ymin[i] )
                    {
                        local_boxes.insert(i + (N-2)*M + 1);
                        if ( !is_zmax[i] )
                        {
                            local_boxes.insert(i + M*N + (N-2)*M + 1);
                        }
                        if ( !is_zmin[i] )
                        {
                            // Lower z level
                            local_boxes.insert(i-2*M+1);
                        }
                    }
                    if ( is_xmin[i] && is_ymax[i] )
                    {
                        local_boxes.insert(i + (N-2)*M - 1);
                        if (!is_zmax[i])
                        {
                            // Upper z level
                            local_boxes.insert(i - 2*M - 1);
                        }
                        if (!is_zmin[i])
                        {
                            // Lower z level
                            local_boxes.insert(i-2*(N-1)*M-1);
                        }

                    }
                    if ( is_xmax[i] && is_ymax[i] )
                    {
                        if (!is_zmax[i])
                        {
                            // Upper z level
                            local_boxes.insert(i + 1);
                        }
                        if (!is_zmin[i])
                        {
                            // Lower z level
                            local_boxes.insert(i - 2*M*N + 1);
                        }

                    }
                }

                if (mIsPeriodicInZ && (is_zmin[i] || is_zmax[i]))
                {
                    if ( is_zmin[i] )
                    {
                        // We need to add the top level
                        unsigned above_box = i+(P-1)*M*N;
                        local_boxes.insert(above_box);
                        if (!is_xmin[i])
                        {
                            // Also add the boxes to the left and at the top
                            local_boxes.insert(above_box-1);
                            if (!is_ymax[i])
                            {
                                local_boxes.insert(above_box+M-1);
                            }

                            if (!is_ymin[i])
                            {
                                local_boxes.insert(above_box-M-1);
                            }
                        }
                        else if ( mIsPeriodicInX )
                        {
                            // Add x periodic box in top layer
                            local_boxes.insert(above_box+M-1);
                            if ( !is_ymin[i] )
                            {
                                local_boxes.insert(above_box+2*M-1);
                            }
                            else if ( mIsPeriodicInY )
                            {
                                // Add the xy periodic box in top layer if at y min or max
                                local_boxes.insert(above_box+M*N-1);
                            }
                            if ( !is_ymax[i] )
                            {
                                local_boxes.insert(above_box+2*M-1);
                            }
                            else if ( mIsPeriodicInY )
                            {
                                // Add the xy periodic box in top layer if at y min or max
                                local_boxes.insert(above_box-M*(N-2)-1);
                            }
                        }

                        if (!is_xmax[i])
                        {
                            // Also add the boxes to the left and at the top
                            local_boxes.insert(above_box+1);
                            if (!is_ymax[i])
                            {
                                local_boxes.insert(above_box+M+1);
                            }
                            if (!is_ymin[i])
                            {
                                local_boxes.insert(above_box-M+1);
                            }
                        }
                        else if ( mIsPeriodicInX )
                        {
                            // Add x periodic box in top layer
                            local_boxes.insert(above_box - M+1);
                            if ( mIsPeriodicInY )
                            {
                                if ( is_ymin[i] )
                                {
                                    // Add the xy periodic box in top layer if at y min or max
                                    local_boxes.insert(above_box+M*(N-2)+1);
                                }
                                else if ( is_ymax[i] )
                                {
                                    // Add the xy periodic box in top layer if at y min or max
                                    local_boxes.insert(above_box-M*N+1);
                                }
                            }
                        }


                        if (!is_ymax[i])
                        {
                            local_boxes.insert(above_box+M);
                        }
                        else if ( mIsPeriodicInY )
                        {
                            local_boxes.insert(above_box-M*(N-1));
                            if ( !is_xmin[i] )
                            {
                                local_boxes.insert(above_box-M*(N-1)-1);
                            }
                            if ( !is_xmax[i] )
                            {
                                local_boxes.insert(above_box-M*(N-1)+1);
                            }
                        }
                        if (!is_ymin[i])
                        {
                            local_boxes.insert(above_box-M);
                        }
                        else if ( mIsPeriodicInY )
                        {
                            local_boxes.insert(above_box+M*(N-1));
                            if ( !is_xmin[i] )
                            {
                                local_boxes.insert(above_box+M*(N-1)-1);
                            }
                            if ( !is_xmax[i] )
                            {
                                local_boxes.insert(above_box+M*(N-1)+1);
                            }
                        }
                    }
                    else if ( is_zmax[i] )
                    {
                        // We need to add the bottom level
                        unsigned below_box = i-(P-1)*M*N;
                        local_boxes.insert(below_box);
                        if (!is_xmin[i])
                        {
                            // Also add the boxes to the left and at the top
                            local_boxes.insert(below_box-1);
                            if (!is_ymax[i])
                            {
                                local_boxes.insert(below_box+M-1);
                            }

                            if (!is_ymin[i])
                            {
                                local_boxes.insert(below_box-M-1);
                            }
                        }
                        else if ( mIsPeriodicInX )
                        {
                            // Add x periodic box in top layer
                            local_boxes.insert(below_box+M-1);
                            if ( mIsPeriodicInY )
                            {
                                if ( is_ymin[i] )
                                {
                                    // Add the xy periodic box in top layer if at y min or max
                                    local_boxes.insert(below_box+M*N-1);
                                }
                                else if ( is_ymax[i] )
                                {
                                    // Add the xy periodic box in top layer if at y min or max
                                    local_boxes.insert(below_box-M*(N-2)-1);
                                }
                            }
                        }

                        if (!is_xmax[i])
                        {
                            // Also add the boxes to the left and at the top
                            local_boxes.insert(below_box+1);
                            if (!is_ymax[i])
                            {
                                local_boxes.insert(below_box+M+1);
                            }
                            if (!is_ymin[i])
                            {
                                local_boxes.insert(below_box-M+1);
                            }
                        }
                        else if ( mIsPeriodicInX )
                        {
                            // Add x periodic box in top layer
                            local_boxes.insert(below_box - M+1);
                            if ( mIsPeriodicInY )
                            {
                                if ( is_ymin[i] )
                                {
                                    // Add the xy periodic box in top layer if at y min or max
                                    local_boxes.insert(below_box+M*(N-2)+1);
                                }
                                else if ( is_ymax[i] )
                                {
                                    // Add the xy periodic box in top layer if at y min or max
                                    local_boxes.insert(below_box-M*N+1);
                                }
                            }
                        }


                        if (!is_ymax[i])
                        {
                            local_boxes.insert(below_box+M);
                        }
                        else if ( mIsPeriodicInY )
                        {
                            local_boxes.insert(below_box-M*(N-1));
                            if ( !is_xmin[i] )
                            {
                                local_boxes.insert(below_box-M*(N-1)+1);
                            }
                            if (!is_xmax[i])
                            {
                                local_boxes.insert(below_box-M*(N-1)-1);
                            }
                        }
                        if (!is_ymin[i])
                        {
                            local_boxes.insert(below_box-M);
                        }
                        else if ( mIsPeriodicInY )
                        {
                            local_boxes.insert(below_box+M*(N-1));
                            if ( !is_xmin[i] )
                            {
                                local_boxes.insert(below_box+M*(N-1)+1);
                            }
                            if (!is_xmax[i])
                            {
                                local_boxes.insert(below_box+M*(N-1)-1);
                            }
                        }
                    }
                }

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        default:
            NEVER_REACHED;
    }
}

template<unsigned DIM>
std::set<unsigned>& DistributedBoxCollection<DIM>::rGetLocalBoxes(unsigned boxIndex)
{
    // Make sure the box is locally owned
    assert(!(boxIndex < mMinBoxIndex) && !(mMaxBoxIndex<boxIndex));
    return mLocalBoxes[boxIndex-mMinBoxIndex];
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::IsOwned(Node<DIM>* pNode)
{
    unsigned index = CalculateContainingBox(pNode);

    return IsBoxOwned(index);
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::IsOwned(c_vector<double, DIM>& location)
{
    unsigned index = CalculateContainingBox(location);

    return IsBoxOwned(index);
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::GetProcessOwningNode(Node<DIM>* pNode)
{
    unsigned box_index = CalculateContainingBox(pNode);
    unsigned containing_process = PetscTools::GetMyRank();

    if (box_index > mMaxBoxIndex)
    {
        containing_process++;
    }
    else if (box_index < mMinBoxIndex)
    {
        containing_process--;
    }

    // Need a special case for periodicity
    if (  mIsPeriodicAcrossProcs )
    {
        if (PetscTools::AmMaster() && box_index > ((mNumBoxesEachDirection[DIM-1]-1)*mNumBoxesInAFace-1))
        {
            // It needs to move to the top process
            containing_process = PetscTools::GetNumProcs()-1;
        }
        else if (PetscTools::AmTopMost() && box_index < mNumBoxesInAFace)
        {
            // It needs to move to the bottom process
            containing_process = 0;
        }
    }

    return containing_process;
}

template<unsigned DIM>
std::vector<unsigned>& DistributedBoxCollection<DIM>::rGetHaloNodesRight()
{
    return mHaloNodesRight;
}

template<unsigned DIM>
std::vector<unsigned>& DistributedBoxCollection<DIM>::rGetHaloNodesLeft()
{
    return mHaloNodesLeft;
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::SetCalculateNodeNeighbours(bool calculateNodeNeighbours)
{
    mCalculateNodeNeighbours = calculateNodeNeighbours;
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    rNodePairs.clear();

    // Create an empty neighbours set for each node
    for (unsigned i=0; i<rNodes.size(); i++)
    {
        // Get the box containing this node as only nodes on this process have NodeAttributes
        // and therefore Neighbours setup.
        unsigned box_index = CalculateContainingBox(rNodes[i]);

        if (IsBoxOwned(box_index))
        {
            rNodes[i]->ClearNeighbours();
        }
    }

    for (unsigned box_index=mMinBoxIndex; box_index<=mMaxBoxIndex; box_index++)
    {
        AddPairsFromBox(box_index, rNodePairs);
    }

    if (mCalculateNodeNeighbours)
    {
        for (unsigned i = 0; i < rNodes.size(); i++)
        {
            // Get the box containing this node as only nodes on this process have NodeAttributes
            // and therefore Neighbours setup.
            unsigned box_index = CalculateContainingBox(rNodes[i]);

            if (IsBoxOwned(box_index))
            {
                rNodes[i]->RemoveDuplicateNeighbours();
            }
        }
    }
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::CalculateInteriorNodePairs(std::vector<Node<DIM>*>& rNodes, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    rNodePairs.clear();

    // Create an empty neighbours set for each node
    for (unsigned i=0; i<rNodes.size(); i++)
    {
        // Get the box containing this node as only nodes on this process have NodeAttributes
        // and therefore Neighbours setup.
        unsigned box_index = CalculateContainingBox(rNodes[i]);

        if (IsBoxOwned(box_index))
        {
            rNodes[i]->ClearNeighbours();
            rNodes[i]->SetNeighboursSetUp(false);
        }
    }

    for (unsigned box_index=mMinBoxIndex; box_index<=mMaxBoxIndex; box_index++)
    {
        if (IsInteriorBox(box_index))
        {
            AddPairsFromBox(box_index, rNodePairs);
        }
    }

    if (mCalculateNodeNeighbours)
    {
        for (unsigned i = 0; i < rNodes.size(); i++)
        {
            // Get the box containing this node as only nodes on this process have NodeAttributes
            // and therefore Neighbours setup.
            unsigned box_index = CalculateContainingBox(rNodes[i]);

            if (IsBoxOwned(box_index))
            {
                rNodes[i]->RemoveDuplicateNeighbours();
                rNodes[i]->SetNeighboursSetUp(true);
            }
        }
    }
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::CalculateBoundaryNodePairs(std::vector<Node<DIM>*>& rNodes, std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    for (unsigned box_index=mMinBoxIndex; box_index<=mMaxBoxIndex; box_index++)
    {
        if (!IsInteriorBox(box_index))
        {
            AddPairsFromBox(box_index, rNodePairs);
        }
    }

    if (mCalculateNodeNeighbours)
    {
        for (unsigned i = 0; i < rNodes.size(); i++)
        {
            // Get the box containing this node as only nodes on this process have NodeAttributes
            // and therefore Neighbours setup.
            unsigned box_index = CalculateContainingBox(rNodes[i]);

            if (IsBoxOwned(box_index))
            {
                rNodes[i]->RemoveDuplicateNeighbours();
                rNodes[i]->SetNeighboursSetUp(true);
            }
        }
    }
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::AddPairsFromBox(unsigned boxIndex,
                                                    std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs)
{
    // Get the box
    Box<DIM>& r_box = rGetBox(boxIndex);

    // Get the set of nodes in this box
    const std::set< Node<DIM>* >& r_contained_nodes = r_box.rGetNodesContained();

    // Get the local boxes to this box
    const std::set<unsigned>& local_boxes_indices = rGetLocalBoxes(boxIndex);

    // Loop over all the local boxes
    for (std::set<unsigned>::iterator box_iter = local_boxes_indices.begin();
         box_iter != local_boxes_indices.end();
         box_iter++)
    {
        Box<DIM>* p_neighbour_box;

        // Establish whether box is locally owned or halo.
        if (IsBoxOwned(*box_iter))
        {
            p_neighbour_box = &mBoxes[*box_iter - mMinBoxIndex];
        }
        else // Assume it is a halo.
        {
            p_neighbour_box = &mHaloBoxes[mHaloBoxesMapping[*box_iter]];
        }
        assert(p_neighbour_box);

        // Get the set of nodes contained in this box
        std::set< Node<DIM>* >& r_contained_neighbour_nodes = p_neighbour_box->rGetNodesContained();

        // Loop over these nodes
        for (typename std::set<Node<DIM>*>::iterator neighbour_node_iter = r_contained_neighbour_nodes.begin();
             neighbour_node_iter != r_contained_neighbour_nodes.end();
             ++neighbour_node_iter)
        {
            // Get the index of the other node
            unsigned other_node_index = (*neighbour_node_iter)->GetIndex();

            // Loop over nodes in this box
            for (typename std::set<Node<DIM>*>::iterator node_iter = r_contained_nodes.begin();
                 node_iter != r_contained_nodes.end();
                 ++node_iter)
            {
                unsigned node_index = (*node_iter)->GetIndex();

                // If we're in the same box, then take care not to store the node pair twice
                if (*box_iter != boxIndex || other_node_index > node_index)
                {
                    rNodePairs.push_back(std::pair<Node<DIM>*, Node<DIM>*>((*node_iter), (*neighbour_node_iter)));
                    if (mCalculateNodeNeighbours)
                    {
                        (*node_iter)->AddNeighbour(other_node_index);
                        (*neighbour_node_iter)->AddNeighbour(node_index);
                    }
                }

            }
        }
    }
}

template<unsigned DIM>
std::vector<int> DistributedBoxCollection<DIM>::CalculateNumberOfNodesInEachStrip()
{
    std::vector<int> cell_numbers(mpDistributedBoxStackFactory->GetHigh() - mpDistributedBoxStackFactory->GetLow(), 0);

    for (unsigned global_index=mMinBoxIndex; global_index<=mMaxBoxIndex; global_index++)
    {
        c_vector<unsigned, DIM> coords = CalculateGridIndices(global_index);
        unsigned location_in_vector = coords[DIM-1] - mpDistributedBoxStackFactory->GetLow();
        unsigned local_index = global_index - mMinBoxIndex;
        cell_numbers[location_in_vector] += mBoxes[local_index].rGetNodesContained().size();
    }

    return cell_numbers;
}

///////// Explicit instantiation///////

template class DistributedBoxCollection<1>;
template class DistributedBoxCollection<2>;
template class DistributedBoxCollection<3>;
