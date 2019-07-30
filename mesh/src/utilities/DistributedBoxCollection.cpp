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
#include "DistributedBoxCollection.hpp"
#include "Exception.hpp"
#include "MathsCustomFunctions.hpp"
#include "Warnings.hpp"

// Static member for "fudge factor" is instantiated here
template<unsigned DIM>
const double DistributedBoxCollection<DIM>::msFudge = 5e-14;

template<unsigned DIM>
DistributedBoxCollection<DIM>::DistributedBoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize, bool isPeriodicInX, int localRows)
    : mBoxWidth(boxWidth),
      mIsPeriodicInX(isPeriodicInX),
      mAreLocalBoxesSet(false),
      mCalculateNodeNeighbours(true)
{
    // Periodicity only works in 2d
    if (isPeriodicInX)
    {
        assert(DIM==2);    // LCOV_EXCL_LINE
    }

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
                // We only need to look for neighbours in the current box and half the neighbouring boxes plus some others for halos
                mLocalBoxes.clear();

                for (unsigned global_index = mMinBoxIndex; global_index<mMaxBoxIndex+1; global_index++)
                {
                    std::set<unsigned> local_boxes;

                    // Set up bools to find out where we are
                    bool left = (global_index%mNumBoxesEachDirection(0) == 0);
                    bool right = (global_index%mNumBoxesEachDirection(0) == mNumBoxesEachDirection(0)-1);
                    bool top = !(global_index < mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1) - mNumBoxesEachDirection(0));
                    bool bottom = (global_index < mNumBoxesEachDirection(0));
                    bool bottom_proc = (CalculateGridIndices(global_index)[1] == mpDistributedBoxStackFactory->GetLow());

                    // Insert the current box
                    local_boxes.insert(global_index);

                    // If we're on the bottom of the process boundary, but not the bottom of the domain add boxes below
                    if (!bottom && bottom_proc)
                    {
                        local_boxes.insert(global_index - mNumBoxesEachDirection(0));
                        if (!left)
                        {
                            local_boxes.insert(global_index - mNumBoxesEachDirection(0) - 1);
                        }
                        if (!right)
                        {
                            local_boxes.insert(global_index - mNumBoxesEachDirection(0) + 1);
                        }
                    }

                    // If we're not at the top of the domain insert boxes above
                    if (!top)
                    {
                        local_boxes.insert(global_index + mNumBoxesEachDirection(0));

                        if (!right)
                        {
                            local_boxes.insert(global_index + mNumBoxesEachDirection(0) + 1);
                        }
                        if (!left)
                        {
                            local_boxes.insert(global_index + mNumBoxesEachDirection(0) - 1);
                        }
                        else if ((global_index % mNumBoxesEachDirection(0) == 0) && (mIsPeriodicInX)) // If we're on the left edge but its periodic include the box on the far right and up one.
                        {
                            local_boxes.insert(global_index +  2 * mNumBoxesEachDirection(0) - 1);
                        }
                    }

                    // If we're not on the far right hand side inseryt box to the right
                    if (!right)
                    {
                        local_boxes.insert(global_index + 1);
                    }
                    // If we're on the right edge but it's periodic include the box on the far left of the domain
                    else if ((global_index % mNumBoxesEachDirection(0) == mNumBoxesEachDirection(0)-1) && (mIsPeriodicInX))
                    {
                        local_boxes.insert(global_index - mNumBoxesEachDirection(0) + 1);
                        // If we're also not on the top-most row, then insert the box above- on the far left of the domain
                        if (global_index < mBoxes.size() - mNumBoxesEachDirection(0))
                        {
                            local_boxes.insert(global_index + 1);
                        }
                    }

                    mLocalBoxes.push_back(local_boxes);
                }
                break;
            }
            case 3:
            {
                // We only need to look for neighbours in the current box and half the neighbouring boxes plus some others for halos
                mLocalBoxes.clear();
                unsigned num_boxes_xy = mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1);


                for (unsigned global_index = mMinBoxIndex; global_index<mMaxBoxIndex+1; global_index++)
                {
                    std::set<unsigned> local_boxes;

                    // Set some bools to find out where we are
                    bool top = !(global_index % num_boxes_xy < num_boxes_xy - mNumBoxesEachDirection(0));
                    bool bottom = (global_index % num_boxes_xy < mNumBoxesEachDirection(0));
                    bool left = (global_index % mNumBoxesEachDirection(0) == 0);
                    bool right = (global_index % mNumBoxesEachDirection(0) == mNumBoxesEachDirection(0) - 1);
                    bool front = (global_index < num_boxes_xy);
                    bool back = !(global_index < num_boxes_xy*mNumBoxesEachDirection(2) - num_boxes_xy);
                    bool proc_front = (CalculateGridIndices(global_index)[2] == mpDistributedBoxStackFactory->GetLow());
                    bool proc_back = (CalculateGridIndices(global_index)[2] == mpDistributedBoxStackFactory->GetHigh()-1);

                    // Insert the current box
                    local_boxes.insert(global_index);

                    // If we're not on the front face, add appropriate boxes on the closer face
                    if (!front)
                    {
                        // If we're not on the top of the domain
                        if (!top)
                        {
                            local_boxes.insert( global_index - num_boxes_xy + mNumBoxesEachDirection(0) );
                            if (!left)
                            {
                                local_boxes.insert( global_index - num_boxes_xy + mNumBoxesEachDirection(0) - 1);
                            }
                            if (!right)
                            {
                                local_boxes.insert( global_index - num_boxes_xy + mNumBoxesEachDirection(0) + 1);
                            }
                        }
                        if (!right)
                        {
                            local_boxes.insert( global_index - num_boxes_xy + 1);
                        }

                        // If we are on the front of the process we have to add extra boxes as they are halos.
                        if (proc_front)
                        {
                            local_boxes.insert( global_index - num_boxes_xy );

                            if (!left)
                            {
                                local_boxes.insert( global_index - num_boxes_xy - 1);
                            }
                            if (!bottom)
                            {
                                local_boxes.insert( global_index - num_boxes_xy - mNumBoxesEachDirection(0));

                                if (!left)
                                {
                                    local_boxes.insert( global_index - num_boxes_xy - mNumBoxesEachDirection(0) - 1);
                                }
                                if (!right)
                                {
                                    local_boxes.insert( global_index - num_boxes_xy - mNumBoxesEachDirection(0) + 1);
                                }
                            }

                        }
                    }
                    if (!right)
                    {
                        local_boxes.insert( global_index + 1);
                    }
                    // If we're not on the very top add boxes above
                    if (!top)
                    {
                        local_boxes.insert( global_index + mNumBoxesEachDirection(0));

                        if (!right)
                        {
                            local_boxes.insert( global_index + mNumBoxesEachDirection(0) + 1);
                        }
                        if (!left)
                        {
                            local_boxes.insert( global_index + mNumBoxesEachDirection(0) - 1);
                        }
                    }

                    // If we're not on the back add boxes behind
                    if (!back)
                    {
                        local_boxes.insert(global_index + num_boxes_xy);

                        if (!right)
                        {
                            local_boxes.insert(global_index + num_boxes_xy + 1);
                        }
                        if (!top)
                        {
                            local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0));
                            if (!right)
                            {
                                local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0) + 1);
                            }
                            if (!left)
                            {
                                local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0) - 1);
                            }
                        }
                        // If we are on the back proc we should make sure we get everything in the face further back
                        if (proc_back)
                        {
                            if (!left)
                            {
                                local_boxes.insert(global_index + num_boxes_xy - 1);
                            }
                            if (!bottom)
                            {
                                local_boxes.insert(global_index + num_boxes_xy - mNumBoxesEachDirection(0));

                                if (!left)
                                {
                                    local_boxes.insert(global_index + num_boxes_xy - mNumBoxesEachDirection(0) - 1);
                                }
                                if (!right)
                                {
                                    local_boxes.insert(global_index + num_boxes_xy - mNumBoxesEachDirection(0) + 1);
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

                // add the one above
                if (!is_ymax[i])
                {
                    local_boxes.insert(i+M);
                }

                // add the four corner boxes

                if ((!is_xmin[i]) && (!is_ymin[i]))
                {
                    local_boxes.insert(i-1-M);
                }
                if ((!is_xmin[i]) && (!is_ymax[i]))
                {
                    local_boxes.insert(i-1+M);
                }
                if ((!is_xmax[i]) && (!is_ymin[i]))
                {
                    local_boxes.insert(i+1-M);
                }
                if ((!is_xmax[i]) && (!is_ymax[i]))
                {
                    local_boxes.insert(i+1+M);
                }

                // Add periodic corner boxes if needed
                if (mIsPeriodicInX)
                {
                    if ((is_xmin[i]) && (!is_ymin[i]))
                    {
                        local_boxes.insert(i-1);
                    }
                    if ((is_xmin[i]) && (!is_ymax[i]))
                    {
                        local_boxes.insert(i-1+2*M);
                    }
                    if ((is_xmax[i]) && (!is_ymin[i]))
                    {
                        local_boxes.insert(i+1-2*M);
                    }
                    if ((is_xmax[i]) && (!is_ymax[i]))
                    {
                        local_boxes.insert(i+1);
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

                if ((!is_xmin[i]) && (!is_ymin[i]) && (!is_zmin[i]))
                {
                    local_boxes.insert(i-1-M-M*N);
                }

                if ((!is_xmin[i]) && (!is_ymin[i]) && (!is_zmax[i]))
                {
                    local_boxes.insert(i-1-M+M*N);
                }

                if ((!is_xmin[i]) && (!is_ymax[i]) && (!is_zmin[i]))
                {
                    local_boxes.insert(i-1+M-M*N);
                }

                if ((!is_xmin[i]) && (!is_ymax[i]) && (!is_zmax[i]))
                {
                    local_boxes.insert(i-1+M+M*N);
                }

                if ((!is_xmax[i]) && (!is_ymin[i]) && (!is_zmin[i]))
                {
                    local_boxes.insert(i+1-M-M*N);
                }

                if ((!is_xmax[i]) && (!is_ymin[i]) && (!is_zmax[i]))
                {
                    local_boxes.insert(i+1-M+M*N);
                }

                if ((!is_xmax[i]) && (!is_ymax[i]) && (!is_zmin[i]))
                {
                    local_boxes.insert(i+1+M-M*N);
                }

                if ((!is_xmax[i]) && (!is_ymax[i]) && (!is_zmax[i]))
                {
                    local_boxes.insert(i+1+M+M*N);
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
