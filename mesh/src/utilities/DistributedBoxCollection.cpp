/*

Copyright (c) 2005-2013, University of Oxford.
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


/////////////////////////////////////////////////////////////////////////////
// DistributedBoxCollection methods
/////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
DistributedBoxCollection<DIM>::DistributedBoxCollection(double boxWidth, c_vector<double, 2*DIM> domainSize, bool isPeriodicInX, int localRows)
    : mDomainSize(domainSize),
      mBoxWidth(boxWidth),
      mIsPeriodicInX(isPeriodicInX),
      mAreLocalBoxesSet(false)
{
    // Periodicity only works in 2d and in serial.
    if (isPeriodicInX)
    {
        assert(DIM==2);
        assert(PetscTools::IsSequential());
    }

    // We insist that the user provide a box width that divides the domainSize in each direction.
    ///\todo Consider adapting domain size to the width of the boxes.
    for (unsigned i=0; i<DIM; i++)
    {
        if (!Divides(boxWidth, (domainSize[2*i+1]-domainSize[2*i])))
        {
            EXCEPTION("The domainSize must be divisible by the boxWidth in each coordinate");
        }
    }

    /*
     * Start by calculating the number of boxes in each direction and total number of boxes.
     * Also create a helper vector of coefficients, whose first entry is 1 and whose i-th
     * entry (for i>1) is the i-th partial product of the vector mNumBoxesEachDirection. This
     * vector of coefficients will be used in the next code block to compute how many boxes
     * along in each dimension each box, given its index.
     */
    unsigned num_boxes = 1;
    std::vector<unsigned> coefficients;
    coefficients.push_back(1);

    for (unsigned i=0; i<DIM; i++)
    {
        mNumBoxesEachDirection(i) = (unsigned) floor((domainSize(2*i+1) - domainSize(2*i))/boxWidth + mFudge);
        num_boxes *= mNumBoxesEachDirection(i);
        coefficients.push_back(coefficients[i]*mNumBoxesEachDirection(i));
    }

    /*
     * Set up a PETSc distributed vector (0,1,2,....mNumBixesEachDirection(0)) to divide the
     * stacks of boxes among the processes.
     */
    std::vector<double> stacks_vector;
    for (unsigned i=0;i<mNumBoxesEachDirection(DIM-1);i++)
    {
        stacks_vector.push_back(i);
    }

    DistributedVectorFactory factory(mNumBoxesEachDirection(DIM-1), localRows);

    Vec petsc_vec = PetscTools::CreateVec(stacks_vector.size(), localRows);

    // Add data from mNumBoxesEachDirection(DIM-1)
    double* p_ret;
    VecGetArray(petsc_vec, &p_ret);
    int lo, hi;
    VecGetOwnershipRange(petsc_vec, &lo, &hi);

    for (int global_index=lo; global_index<hi; global_index++)
    {
        int local_index = global_index - lo;
        p_ret[local_index] = stacks_vector[global_index];
    }
    VecRestoreArray(petsc_vec, &p_ret);

    mpDistributedBoxStacks = new DistributedVector(petsc_vec, &factory);

    // Set up MPI information.
    mProcRight  = (PetscTools::AmTopMost()) ?   MPI_PROC_NULL   :   (int)PetscTools::GetMyRank()    +1;
    mProcLeft   = (PetscTools::AmMaster())  ?   MPI_PROC_NULL   :   (int)PetscTools::GetMyRank()    -1;

    mMinBoxIndex = UINT_MAX;
    mMaxBoxIndex = 0;

    for (DistributedVector::Iterator stack_index = mpDistributedBoxStacks->Begin();
         stack_index != mpDistributedBoxStacks->End();
         ++stack_index)
    {
        int num_boxes_in_plane = 1; // The number of boxes in each layer of boxes.
        std::vector<unsigned> counter_conditions(1, 1u); // The conditions to increment each counter.

        for (int i=0; i<(int)DIM-1; i++)
        {
            num_boxes_in_plane *= mNumBoxesEachDirection(i);
            counter_conditions.push_back(mNumBoxesEachDirection(i));
        }

        std::vector<unsigned> counters(DIM-1, 0);
        for (int i=0; i<num_boxes_in_plane; /*No increment here */)
        {
            c_vector<double, 2*DIM> box_coords;
            c_vector<unsigned, DIM> box_indices;
            for (int d=0; d<(int)DIM-1; d++)
            {
                box_coords[2*d] = (double)domainSize[2*d]+mBoxWidth*counters[d];
                box_coords[2*d+1] = (double)domainSize[2*d]+mBoxWidth*(1+counters[d]);
                box_indices[d] = counters[d];
            }
            box_coords[2*DIM-2]=domainSize[2*DIM-2]+mBoxWidth*(*mpDistributedBoxStacks)[stack_index];
            box_coords[2*DIM-1]=domainSize[2*DIM-2]+mBoxWidth*((*mpDistributedBoxStacks)[stack_index]+1);
            box_indices[DIM-1] = (unsigned)(*mpDistributedBoxStacks)[stack_index];

            Box<DIM> new_box(box_coords);
            mBoxes.push_back(new_box);

            unsigned global_index = CalculateGlobalIndex(box_indices);
            mBoxesMapping[global_index] = mBoxes.size()-1;
            mMinBoxIndex = (global_index < mMinBoxIndex) ? global_index : mMinBoxIndex;
            mMaxBoxIndex = (mMaxBoxIndex < global_index) ? global_index : mMaxBoxIndex;

            /* Increment counters */
            i++;
            for (int var = 0; var < (int)DIM-1; var++)
            {
              if (i%counter_conditions[var] == 0)
              {
                  counters[var] = ((counters[var]+1)%mNumBoxesEachDirection(var));
              }
            }
            /* Incremented counters */
        }
    }

    /**
     * Update the local domain region.
     */
    unsigned num_rows = mpDistributedBoxStacks->GetHigh() - mpDistributedBoxStacks->GetLow();
    for (int d=0; d<(int)DIM-1; d++)
    {
        mLocalDomainSize[2*d] = domainSize[2*d];
        mLocalDomainSize[2*d+1] = domainSize[2*d+1];
    }
    mLocalDomainSize[2*DIM-2] = domainSize[2*DIM-2] + mBoxWidth*(*mpDistributedBoxStacks)[mpDistributedBoxStacks->Begin()];
    mLocalDomainSize[2*DIM-1] = mLocalDomainSize[2*DIM-2] + mBoxWidth*num_rows;

    // Check we have set up the right number of total boxes set up.
    unsigned local_boxes=mBoxes.size();
    unsigned expected_boxes=1;
    for (unsigned i=0;i<DIM;i++)
    {
        expected_boxes*=mNumBoxesEachDirection[i];
    }
    unsigned total_boxes;

    MPI_Allreduce(&local_boxes,&total_boxes,1,MPI_UNSIGNED,MPI_SUM,PETSC_COMM_WORLD);

    mNumBoxes=total_boxes;
    PetscTools::Destroy(petsc_vec);
}

template<unsigned DIM>
DistributedBoxCollection<DIM>::~DistributedBoxCollection()
{
    delete mpDistributedBoxStacks;
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::SetupHaloBoxes()
{
    // Get top-most and bottom-most value of Distributed Box Stack.
    unsigned Hi=mpDistributedBoxStacks->GetHigh();
    unsigned Lo=mpDistributedBoxStacks->GetLow();

    c_vector<double, 2*DIM> box_coords;

    std::vector<unsigned> counter_conditions(1,1);
    unsigned num_boxes_in_plane = 1;
    for (int d=0; d<(int)DIM-1; d++)
    {
        num_boxes_in_plane*=mNumBoxesEachDirection(d);
        counter_conditions.push_back(mNumBoxesEachDirection(d));
    }

    // If I am not the top-most process, add halo structures to the right.
    if (!PetscTools::AmTopMost())
    {
        std::vector<unsigned> counters(DIM-1, 0);
        box_coords[2*DIM-2] = mDomainSize[2*DIM-2]+Hi*mBoxWidth;
        box_coords[2*DIM-1] = mDomainSize[2*DIM-2]+(Hi+1)*mBoxWidth;

        for (unsigned i=0; i< num_boxes_in_plane; /* Increment later*/)
        {
            for (int d=0; d<(int)DIM-1; d++)
            {
                box_coords[2*d] = mDomainSize[2*d] + counters[d]*mBoxWidth;
                box_coords[2*d+1] = mDomainSize[2*d] + (1+counters[d])*mBoxWidth;
            }

            // Add the halo box
            Box<DIM> new_box(box_coords);
            mHaloBoxes.push_back(new_box);

            // Add the corresponding halo on this process to mHalosRight
            c_vector<unsigned, DIM> coords;
            for (int d=0; d<(int)DIM-1; d++)
            {
                coords[d] = counters[d];
            }
            coords[DIM-1]= Hi;
            mHaloBoxesMapping[CalculateGlobalIndex(coords)]=mHaloBoxes.size()-1;

            coords[DIM-1] = Hi-1;
            mHalosRight.push_back(CalculateGlobalIndex(coords));



            /* Increment counters */
            i++;
            for (int var = 0; var < (int)DIM-1; var++)
            {
              if (i%counter_conditions[var] == 0)
              {
                  counters[var] = ((counters[var]+1)%mNumBoxesEachDirection(var));
              }
            }
            /* Incremented counters */
        }
    }
    // If I am not the bottom-most process, add halo structures to the left.
    if (!PetscTools::AmMaster())
    {
        std::vector<unsigned> counters(DIM-1, 0);
        box_coords[2*DIM-2] = mDomainSize[2*DIM-2]+(Lo-1)*mBoxWidth;
        box_coords[2*DIM-1] = mDomainSize[2*DIM-2]+Lo*mBoxWidth;

        for (unsigned i=0; i< num_boxes_in_plane; /* Increment later*/)
        {
            for (int d=0; d<(int)DIM-1; d++)
            {
                box_coords[2*d] = mDomainSize[2*d] + counters[d]*mBoxWidth;
                box_coords[2*d+1] = mDomainSize[2*d] + (1+counters[d])*mBoxWidth;
            }

            Box<DIM> new_box(box_coords);
            mHaloBoxes.push_back(new_box);

            // Add the corresponding halo on this process to mHalosLeft
            c_vector<unsigned, DIM> coords;
            for (int d=0; d<(int)DIM-1; d++)
            {
                coords[d] = counters[d];
            }
            coords[DIM-1]= Lo-1;
            mHaloBoxesMapping[CalculateGlobalIndex(coords)]=mHaloBoxes.size()-1;

            coords[DIM-1] = Lo;
            mHalosLeft.push_back(CalculateGlobalIndex(coords));

            /* Increment counters */
            i++;
            for (int var = 0; var < (int)DIM-1; var++)
            {
              if (i%counter_conditions[var] == 0)
              {
                  counters[var] = ((counters[var]+1)%mNumBoxesEachDirection(var));
              }
            }
            /* Incremented counters */
        }
    }
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::UpdateHaloBoxes()
{
    mHaloNodesLeft.clear();
    for (unsigned i=0;i<mHalosLeft.size();i++)
    {
        for (typename std::set<Node<DIM>* >::iterator iter=this->rGetBox(mHalosLeft[i]).rGetNodesContained().begin();
                iter!=this->rGetBox(mHalosLeft[i]).rGetNodesContained().end();
                iter++)
        {
            mHaloNodesLeft.insert((*iter)->GetIndex());
        }
    }

    // Send right
    mHaloNodesRight.clear();
    for (unsigned i=0;i<mHalosRight.size();i++)
    {
        for (typename std::set<Node<DIM>* >::iterator iter=this->rGetBox(mHalosRight[i]).rGetNodesContained().begin();
                iter!=this->rGetBox(mHalosRight[i]).rGetNodesContained().end();
                iter++)
        {
            mHaloNodesRight.insert((*iter)->GetIndex());
        }
    }
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::GetBoxOwnership(unsigned globalIndex)
{
    return (!(globalIndex<mMinBoxIndex) && !(mMaxBoxIndex<globalIndex));
}

template<unsigned DIM>
bool DistributedBoxCollection<DIM>::GetHaloBoxOwnership(unsigned globalIndex)
{
    assert(globalIndex < mNumBoxes);

    unsigned num_boxes_in_face = 1;
    for (int d=0; d<(int)DIM-1; d++)
    {
        num_boxes_in_face*=mNumBoxesEachDirection(d);
    }
    bool is_parallel = PetscTools::IsParallel();

    bool is_halo_right = ((globalIndex > mMaxBoxIndex) && !(globalIndex > mMaxBoxIndex + num_boxes_in_face));
    bool is_halo_left = ((globalIndex < mMinBoxIndex) && !(globalIndex < mMinBoxIndex - num_boxes_in_face));

    return (is_parallel && (is_halo_right || is_halo_left));
}

template<unsigned DIM>
unsigned DistributedBoxCollection<DIM>::CalculateGlobalIndex(c_vector<unsigned, DIM> coordinateIndices)
{
    unsigned containing_box_index = 0;
    for (unsigned i=0; i<DIM; i++)
    {
        unsigned temp = 1;
        for (unsigned j=0; j<i; j++)
        {
            temp *= mNumBoxesEachDirection(j);
        }
        containing_box_index += temp*coordinateIndices[i];
    }

    return containing_box_index;
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
        if ( (rLocation[i] < mDomainSize(2*i)) || !(rLocation[i] < mDomainSize(2*i+1)) )
        {
            EXCEPTION("The point provided is outside all of the boxes");
        }
    }

    // Compute the containing box index in each dimension
    c_vector<unsigned, DIM> containing_box_indices;
    for (unsigned i=0; i<DIM; i++)
    {
        containing_box_indices[i] = (unsigned) floor((rLocation[i] - mDomainSize(2*i) + mFudge)/mBoxWidth);
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
c_vector<unsigned, DIM> DistributedBoxCollection<DIM>::CalculateCoordinateIndices(unsigned globalIndex)
{
    c_vector<unsigned, DIM> indices;

    switch(DIM)
    {
        case 1:
        {
            indices[0]=globalIndex;
            break;
        }
        case 2:
        {
            unsigned remainder=globalIndex % mNumBoxesEachDirection(0);
            indices[0]=remainder;
            indices[1]=(unsigned)(globalIndex/mNumBoxesEachDirection(0));
            break;
        }

        case 3:
        {
            unsigned remainder1=globalIndex % (mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1));
            unsigned remainder2=remainder1 % mNumBoxesEachDirection(0);
            indices[0]=remainder2;
            indices[1]=((globalIndex-indices[0])/mNumBoxesEachDirection(0))%mNumBoxesEachDirection(1);
            indices[2]=((globalIndex-indices[0]-mNumBoxesEachDirection(0)*indices[1])/(mNumBoxesEachDirection(0)*mNumBoxesEachDirection(1)));
            break;
        }
        default:
        {
            NEVER_REACHED;
        }
    }

    return indices;
}

template<unsigned DIM>
Box<DIM>& DistributedBoxCollection<DIM>::rGetBox(unsigned boxIndex)
{
    assert(!(boxIndex<mMinBoxIndex) && !(mMaxBoxIndex<boxIndex));
    return mBoxes[boxIndex-mMinBoxIndex];
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
c_vector<double, 2*DIM> DistributedBoxCollection<DIM>::GetDomainSize() const
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
unsigned DistributedBoxCollection<DIM>::GetNumRowsOfBoxes() const
{
    return mpDistributedBoxStacks->GetHigh() - mpDistributedBoxStacks->GetLow();
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::SetupLocalBoxesHalfOnly()
{
    if(mAreLocalBoxesSet)
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
                    bool proc_left = (global_index == mpDistributedBoxStacks->GetLow());

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
                    bool bottom_proc = (CalculateCoordinateIndices(global_index)[1] == mpDistributedBoxStacks->GetLow());

                    // Insert the current box
                    local_boxes.insert(global_index);

                    // If we're on the bottom of the process boundary, but not the bottom of the domain add boxes below
                    if(!bottom && bottom_proc)
                    {
                        local_boxes.insert(global_index - mNumBoxesEachDirection(0));
                        if(!left)
                        {
                            local_boxes.insert(global_index - mNumBoxesEachDirection(0) - 1);
                        }
                        if(!right)
                        {
                            local_boxes.insert(global_index - mNumBoxesEachDirection(0) + 1);
                        }
                    }

                    // If we're not at the top of the domain insert boxes above
                    if(!top)
                    {
                        local_boxes.insert(global_index + mNumBoxesEachDirection(0));

                        if(!right)
                        {
                            local_boxes.insert(global_index + mNumBoxesEachDirection(0) + 1);
                        }
                        if(!left)
                        {
                            local_boxes.insert(global_index + mNumBoxesEachDirection(0) - 1);
                        }
                        else if ( (global_index % mNumBoxesEachDirection(0) == 0) && (mIsPeriodicInX) ) // If we're on the left edge but its periodic include the box on the far right and up one.
                        {
                            local_boxes.insert(global_index +  2 * mNumBoxesEachDirection(0) - 1);
                        }
                    }

                    // If we're not on the far right hand side inseryt box to the right
                    if(!right)
                    {
                        local_boxes.insert(global_index + 1);
                    }
                    // If we're on the right edge but it's periodic include the box on the far left of the domain
                    else if ( (global_index % mNumBoxesEachDirection(0) == mNumBoxesEachDirection(0)-1) && (mIsPeriodicInX) )
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
                    bool proc_front = (CalculateCoordinateIndices(global_index)[2] == mpDistributedBoxStacks->GetLow());
                    bool proc_back = (CalculateCoordinateIndices(global_index)[2] == mpDistributedBoxStacks->GetHigh()-1);

                    // Insert the current box
                    local_boxes.insert(global_index);

                    // If we're not on the front face, add appropriate boxes on the closer face
                    if(!front)
                    {
                        // If we're not on the top of the domain
                        if(!top)
                        {
                            local_boxes.insert( global_index - num_boxes_xy + mNumBoxesEachDirection(0) );
                            if(!left)
                            {
                                local_boxes.insert( global_index - num_boxes_xy + mNumBoxesEachDirection(0) - 1);
                            }
                            if(!right)
                            {
                                local_boxes.insert( global_index - num_boxes_xy + mNumBoxesEachDirection(0) + 1);
                            }
                        }
                        if(!right)
                        {
                            local_boxes.insert( global_index - num_boxes_xy + 1);
                        }

                        // If we are on the front of the process we have to add extra boxes as they are halos.
                        if(proc_front)
                        {
                            local_boxes.insert( global_index - num_boxes_xy );

                            if(!left)
                            {
                                local_boxes.insert( global_index - num_boxes_xy - 1);
                            }
                            if(!bottom)
                            {
                                local_boxes.insert( global_index - num_boxes_xy - mNumBoxesEachDirection(0));

                                if(!left)
                                {
                                    local_boxes.insert( global_index - num_boxes_xy - mNumBoxesEachDirection(0) - 1);
                                }
                                if(!right)
                                {
                                    local_boxes.insert( global_index - num_boxes_xy - mNumBoxesEachDirection(0) + 1);
                                }
                            }

                        }
                    }
                    if(!right)
                    {
                        local_boxes.insert( global_index + 1);
                    }
                    // If we're not on the very top add boxes above
                    if(!top)
                    {
                        local_boxes.insert( global_index + mNumBoxesEachDirection(0));

                        if(!right)
                        {
                            local_boxes.insert( global_index + mNumBoxesEachDirection(0) + 1);
                        }
                        if(!left)
                        {
                            local_boxes.insert( global_index + mNumBoxesEachDirection(0) - 1);
                        }
                    }

                    // If we're not on the back add boxes behind
                    if(!back)
                    {
                        local_boxes.insert(global_index + num_boxes_xy);

                        if(!right)
                        {
                            local_boxes.insert(global_index + num_boxes_xy + 1);
                        }
                        if(!top)
                        {
                            local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0));
                            if(!right)
                            {
                                local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0) + 1);
                            }
                            if(!left)
                            {
                                local_boxes.insert(global_index + num_boxes_xy + mNumBoxesEachDirection(0) - 1);
                            }
                        }
                        // If we are on the back proc we should make sure we get everything in the face further back
                        if(proc_back)
                        {
                            if(!left)
                            {
                                local_boxes.insert(global_index + num_boxes_xy - 1);
                            }
                            if(!bottom)
                            {
                                local_boxes.insert(global_index + num_boxes_xy - mNumBoxesEachDirection(0));

                                if(!left)
                                {
                                    local_boxes.insert(global_index + num_boxes_xy - mNumBoxesEachDirection(0) - 1);
                                }
                                if(!right)
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
                    if(mIsPeriodicInX)
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
                    if(mIsPeriodicInX)
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

                // Add Periodic Corner Boxes if needed
                if(mIsPeriodicInX)
                {
                    if( (is_xmin[i]) && (!is_ymin[i]) )
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

                mLocalBoxes.push_back(local_boxes);
            }
            break;
        }
        default:
            NEVER_REACHED;
    }
}

template<unsigned DIM>
std::set<unsigned> DistributedBoxCollection<DIM>::GetLocalBoxes(unsigned boxIndex)
{
    // Make sure the box is locally owned
    assert(!(boxIndex < mMinBoxIndex) && !(mMaxBoxIndex<boxIndex));
    return mLocalBoxes[boxIndex-mMinBoxIndex];
}

template<unsigned DIM>
c_vector<double, 2*DIM> DistributedBoxCollection<DIM>::GetLocalDomainSize()
{
    return mLocalDomainSize;
}

template<unsigned DIM>
void DistributedBoxCollection<DIM>::CalculateNodePairs(std::vector<Node<DIM>*>& rNodes, std::set<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs, std::map<unsigned, std::set<unsigned> >& rNodeNeighbours)
{
    rNodePairs.clear();
    rNodeNeighbours.clear();

    // Create an empty neighbours set for each node
    for (unsigned i=0; i<rNodes.size(); i++)
    {
        unsigned node_index = rNodes[i]->GetIndex();

        // Get the box containing this node
        unsigned box_index = CalculateContainingBox(rNodes[i]);

        if (GetBoxOwnership(box_index))
        {
            rNodeNeighbours[node_index] = std::set<unsigned>();
        }

    }
    for (unsigned i=0; i<rNodes.size(); i++)
    {
        unsigned node_index = rNodes[i]->GetIndex();

        // Get the box containing this node
        unsigned box_index = CalculateContainingBox(rNodes[i]);

        // Only calculate for owned nodes.
        if(GetBoxOwnership(box_index))
        {
            // Get the local boxes to this node
            std::set<unsigned> local_boxes_indices = GetLocalBoxes(box_index);

            // Loop over all the local boxes
            for (std::set<unsigned>::iterator box_iter = local_boxes_indices.begin();
                 box_iter != local_boxes_indices.end();
                 box_iter++)
            {
                Box<DIM>* p_box;

                // Establish whether box is locally owned or halo.
                if (GetBoxOwnership(*box_iter))
                {
                    p_box = &mBoxes[mBoxesMapping[*box_iter]];
                }
                else // Assume it is a halo.
                {
                    p_box = &mHaloBoxes[mHaloBoxesMapping[*box_iter]];
                }
                assert(p_box);

                // Get the set of nodes contained in this box
                std::set< Node<DIM>* >& r_contained_nodes = p_box->rGetNodesContained();

                // Loop over these nodes
                for (typename std::set<Node<DIM>*>::iterator node_iter = r_contained_nodes.begin();
                     node_iter != r_contained_nodes.end();
                     ++node_iter)
                {
                    // Get the index of the other node
                    unsigned other_i = (*node_iter)->GetIndex();

                    // If we're in the same box, then take care not to store the node pair twice
                    if (*box_iter == box_index)
                    {
                        if (other_i > node_index)
                        {
                            rNodePairs.insert(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[i], (*node_iter)));
                            rNodeNeighbours[node_index].insert(other_i);
                            rNodeNeighbours[other_i].insert(node_index);
                        }
                    }
                    else
                    {
                        rNodePairs.insert(std::pair<Node<DIM>*, Node<DIM>*>(rNodes[node_index], (*node_iter)));
                        rNodeNeighbours[node_index].insert(other_i);
                        rNodeNeighbours[other_i].insert(node_index);
                    }
                }
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DistributedBoxCollection<1>;
template class DistributedBoxCollection<2>;
template class DistributedBoxCollection<3>;
