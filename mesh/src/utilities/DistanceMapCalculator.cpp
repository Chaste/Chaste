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

#include "DistanceMapCalculator.hpp"
#include "DistributedTetrahedralMesh.hpp" // For dynamic cast

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::DistanceMapCalculator(
            AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh)
    : mrMesh(rMesh),
      mWorkOnEntireMesh(true),
      mNumHalosPerProcess(nullptr),
      mRoundCounter(0u),
      mPopCounter(0u),
      mTargetNodeIndex(UINT_MAX),
      mSingleTarget(false)
{
    mNumNodes = mrMesh.GetNumNodes();

    DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* p_distributed_mesh = dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(&mrMesh);
    if (PetscTools::IsSequential() || p_distributed_mesh == nullptr)
    {
        // It's a non-distributed mesh
        mLo = 0;
        mHi = mNumNodes;
    }
    else
    {
        // It's a parallel (distributed) mesh (p_parallel_mesh is non-null and we are running in parallel)
        mWorkOnEntireMesh = false;
        mLo = mrMesh.GetDistributedVectorFactory()->GetLow();
        mHi = mrMesh.GetDistributedVectorFactory()->GetHigh();

        // Get local halo information
        p_distributed_mesh->GetHaloNodeIndices(mHaloNodeIndices);

        // Share information on the number of halo nodes
        unsigned my_size = mHaloNodeIndices.size();
        mNumHalosPerProcess = new unsigned[PetscTools::GetNumProcs()];
        MPI_Allgather(&my_size, 1, MPI_UNSIGNED, mNumHalosPerProcess, 1, MPI_UNSIGNED, PETSC_COMM_WORLD);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::ComputeDistanceMap(
        const std::vector<unsigned>& rSourceNodeIndices,
        std::vector<double>& rNodeDistances)
{
    rNodeDistances.resize(mNumNodes);
    for (unsigned index=0; index<mNumNodes; index++)
    {
        rNodeDistances[index] = DBL_MAX;
    }
    assert(mActivePriorityNodeIndexQueue.empty());

    if (mSingleTarget)
    {
        assert(rSourceNodeIndices.size() == 1);
        unsigned source_node_index = rSourceNodeIndices[0];

        // We need to make sure this is local, so that we can use the geometry
        if (mLo<=source_node_index && source_node_index<mHi)
        {
            double heuristic_correction = norm_2(mrMesh.GetNode(source_node_index)->rGetLocation()-mTargetNodePoint);
            PushLocal(heuristic_correction, source_node_index);
            rNodeDistances[source_node_index] = heuristic_correction;
        }
     }
    else
    {
        for (unsigned source_index=0; source_index<rSourceNodeIndices.size(); source_index++)
        {
            unsigned source_node_index = rSourceNodeIndices[source_index];
            PushLocal(0.0, source_node_index);
            rNodeDistances[source_node_index] = 0.0;
        }
    }

    bool non_empty_queue = true;
    mRoundCounter = 0;
    mPopCounter = 0;
    while (non_empty_queue)
    {
        bool termination = WorkOnLocalQueue(rNodeDistances);

        // Sanity - check that we aren't doing this very many times
        if (mRoundCounter++ > 10 * PetscTools::GetNumProcs())
        {
            // This line will be hit if there's a parallel distributed mesh with a really bad partition
            NEVER_REACHED;
        }
        if (mSingleTarget && PetscTools::ReplicateBool(termination))
        {
            // A single process found the target already
            break;
        }
        non_empty_queue = UpdateQueueFromRemote(rNodeDistances);
    }

    if (mWorkOnEntireMesh == false)
    {
        // Update all processes with the best values from everywhere
        // Take a local copy
        std::vector<double> local_distances = rNodeDistances;

        // Share it back into the vector
        MPI_Allreduce( &local_distances[0], &rNodeDistances[0], mNumNodes, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::UpdateQueueFromRemote(std::vector<double>& rNodeDistances)
{
    if (mWorkOnEntireMesh)
    {
        // This update does nowt
        return !mActivePriorityNodeIndexQueue.empty();
    }
    for (unsigned bcast_process=0; bcast_process<PetscTools::GetNumProcs(); bcast_process++)
    {
        // Process packs cart0/cart1/...euclid/index into a 1-d array
        double* dist_exchange = new double[ mNumHalosPerProcess[bcast_process] ];
        unsigned* index_exchange = new unsigned[ mNumHalosPerProcess[bcast_process] ];
        if (PetscTools::GetMyRank() == bcast_process)
        {
            // Broadcaster fills the array
            for (unsigned index=0; index<mHaloNodeIndices.size();index++)
            {
                dist_exchange[index] = rNodeDistances[mHaloNodeIndices[index]];
                index_exchange[index] = mHaloNodeIndices[index];
            }
        }

        /*
         * Broadcast - this is can be done by casting indices to double and
         * packing everything into a single array. That would be better for
         * latency, but this is probably more readable.
         */
        MPI_Bcast(dist_exchange, mNumHalosPerProcess[bcast_process], MPI_DOUBLE,
                  bcast_process, PETSC_COMM_WORLD);
        MPI_Bcast(index_exchange, mNumHalosPerProcess[bcast_process], MPI_UNSIGNED,
                  bcast_process, PETSC_COMM_WORLD);
        if (PetscTools::GetMyRank() != bcast_process)
        {
            // Receiving process take updates
            for (unsigned index=0; index<mNumHalosPerProcess[bcast_process];index++)
            {
                unsigned global_index=index_exchange[index];
                // Is it a better answer?
                if (dist_exchange[index] < rNodeDistances[global_index]*(1.0-2*DBL_EPSILON) )
                {
                    // Copy across - this may be unnecessary when PushLocal isn't going to push because it's not local
                    rNodeDistances[global_index] = dist_exchange[index];
                    PushLocal(rNodeDistances[global_index], global_index);
                }
            }
        }
        delete [] dist_exchange;
        delete [] index_exchange;
    }
    // Is any queue non-empty?
    bool non_empty_queue = PetscTools::ReplicateBool(!mActivePriorityNodeIndexQueue.empty());
    return(non_empty_queue);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::WorkOnLocalQueue(std::vector<double>& rNodeDistances)
{
    unsigned pop_stop = mNumNodes/(PetscTools::GetNumProcs()*20);
    while (!mActivePriorityNodeIndexQueue.empty())
    {
        // Get the next index in the queue
        unsigned current_node_index = mActivePriorityNodeIndexQueue.top().second;
        double distance_when_queued = -mActivePriorityNodeIndexQueue.top().first;
        mActivePriorityNodeIndexQueue.pop();

        // Only act on nodes which haven't been acted on already
        // (It's possible that a better distance has been found and already been dealt with)
        if (distance_when_queued == rNodeDistances[current_node_index])
        {
            mPopCounter++;
            Node<SPACE_DIM>* p_current_node = mrMesh.GetNode(current_node_index);
            double current_heuristic = 0.0;
            if (mSingleTarget)
            {
                 current_heuristic=norm_2(p_current_node->rGetLocation()-mTargetNodePoint);
            }

            // Loop over the elements containing the given node
            for (typename Node<SPACE_DIM>::ContainingElementIterator element_iterator = p_current_node->ContainingElementsBegin();
                element_iterator != p_current_node->ContainingElementsEnd();
                ++element_iterator)
            {
                // Get a pointer to the container element
                Element<ELEMENT_DIM, SPACE_DIM>* p_containing_element = mrMesh.GetElement(*element_iterator);

                // Loop over the nodes of the element
                for (unsigned node_local_index=0;
                   node_local_index<p_containing_element->GetNumNodes();
                   node_local_index++)
                {
                    Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
                    unsigned neighbour_node_index = p_neighbour_node->GetIndex();

                    // Avoid revisiting the active node
                    if (neighbour_node_index != current_node_index)
                    {

                        double neighbour_heuristic=0.0;
                        if (mSingleTarget)
                        {
                             neighbour_heuristic=norm_2(p_neighbour_node->rGetLocation()-mTargetNodePoint);
                        }
                        // Test if we have found a shorter path from the source to the neighbour through current node
                        double updated_distance = rNodeDistances[current_node_index] +
                                                  norm_2(p_neighbour_node->rGetLocation() - p_current_node->rGetLocation())
                                                  - current_heuristic + neighbour_heuristic;
                        if (updated_distance < rNodeDistances[neighbour_node_index] * (1.0-2*DBL_EPSILON))
                        {
                            rNodeDistances[neighbour_node_index] = updated_distance;
                            PushLocal(updated_distance, neighbour_node_index);
                        }
                    }
                }
            }
            if (mSingleTarget)
            {
                if (current_node_index == mTargetNodeIndex)
                {
                    // Premature termination if there is a single goal in mind (and we found it)
                    return true;
                }
                if (mPopCounter%pop_stop == 0)
                {
                    // Premature termination -- in case the work has been done
                    return false;
                }
            }
        }
     }
     return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DistanceMapCalculator<ELEMENT_DIM, SPACE_DIM>::SingleDistance(unsigned sourceNodeIndex, unsigned targetNodeIndex)
{
    std::vector<unsigned> source_node_index_vector;
    source_node_index_vector.push_back(sourceNodeIndex);

    // We are re-using the mCachedDistances vector...
    mTargetNodeIndex = targetNodeIndex; // For premature termination
    mSingleTarget = true;

    // Set up information on target (for A* guidance)
    c_vector<double, SPACE_DIM> target_point = zero_vector<double>(SPACE_DIM);
    if (mrMesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(mTargetNodeIndex))
    {
        // Owner of target node sets target_point (others have zero)
        target_point = mrMesh.GetNode(mTargetNodeIndex)->rGetLocation();
    }

    // Communicate for wherever to everyone
    MPI_Allreduce(&target_point[0], &mTargetNodePoint[0], SPACE_DIM, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    //mTargetNodePoint;
    std::vector<double> distances;
    ComputeDistanceMap(source_node_index_vector, distances);

    ///\todo #1414 premature termination when we find the correct one (parallel)

    // Reset target, so we don't terminate early next time.
    mSingleTarget = false;

    // Make sure that there isn't a non-empty queue from a previous calculation
    if (!mActivePriorityNodeIndexQueue.empty())
    {
        mActivePriorityNodeIndexQueue = std::priority_queue<std::pair<double, unsigned> >();
    }

    return distances[targetNodeIndex];
}

// Explicit instantiation
template class DistanceMapCalculator<1, 1>;
template class DistanceMapCalculator<1, 2>;
template class DistanceMapCalculator<2, 2>;
template class DistanceMapCalculator<1, 3>;
//template class DistanceMapCalculator<2, 3>;
template class DistanceMapCalculator<3, 3>;
