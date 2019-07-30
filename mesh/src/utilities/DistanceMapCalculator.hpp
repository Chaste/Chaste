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
#ifndef DISTANCEMAPCALCULATOR_HPP_
#define DISTANCEMAPCALCULATOR_HPP_

#include <vector>
#include <queue>

#include "UblasIncludes.hpp"
#include "AbstractTetrahedralMesh.hpp"

/**
 * This class provides functionalities to compute a distance map in a given mesh
 * from a given surface, specifying the distance from each node to the surface.
 *
 * The mesh is specified in the constructor, and the ComputeDistanceMap computes
 * (and returns by reference) the map.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class DistanceMapCalculator
{
private:
    friend class TestDistanceMapCalculator;

    /** The mesh*/
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& mrMesh;
    /** Number of nodes in the mesh*/
    unsigned mNumNodes;
    /** Local cache of the nodes owned by this process, from mesh's DistributedVectorFactory*/
    unsigned mLo;
    /** Local cache of the nodes owned by this process, from mesh's DistributedVectorFactory*/
    unsigned mHi;
    /** Whether we should work on the entire mesh.  True if sequential.  True is the mesh is a plain TetrahedralMesh.*/
    bool mWorkOnEntireMesh;
    /** (Only used when mWorkOnEntrireMesh == false).  This forms an array of with the number of halo nodes known by each process.*/
    unsigned *mNumHalosPerProcess;
    /** (Only used when mWorkOnEntrireMesh == false).  This is a local cache of halo node indices.*/
    std::vector<unsigned> mHaloNodeIndices;
    /** Used to check parallel implementation*/
    unsigned mRoundCounter;
    /** Used to check implementation for number of queue pops per calculation*/
    unsigned mPopCounter;
    /** Used in the calculation of point-to-point distances.*/
    unsigned mTargetNodeIndex;
    /** True for point-to-point distances.*/
    bool mSingleTarget;
    /** Also used in the calculation of point-to-point distances with A* heuristic -- this requires a parallel communication*/
    c_vector<double, SPACE_DIM> mTargetNodePoint;

    /**
     * Queue of nodes to be processed (initialised with the nodes defining the surface)
     * Priorities (given as the first in the pair for lexographical ordering) are
     * initialised to -best_distance_to_source so that nodes closest to the source
     * are dealt with first.
     */
    std::priority_queue<std::pair<double, unsigned> > mActivePriorityNodeIndexQueue;

    /**
     * Work on the Queue of node indices (grass-fire across the mesh)
     *
     * @param rNodeDistances distance map computed
     * @return true when a single target has been found
     * @return false when there is work remaining or the queue is flushed
     *
     */
    bool WorkOnLocalQueue(std::vector<double>& rNodeDistances);

    /**
     * Update the local Queue of node indices using data that are from the halo nodes of remote processes.
     *
     * @param rNodeDistances distance map computed
     *
     * @return true when this update was active => there are non-empty queues left to work on
     * @return false without working or side-effects if we don't have a true distributed mesh
     *
     */
    bool UpdateQueueFromRemote(std::vector<double>& rNodeDistances);

    /**
     * Push a node index onto the queue.  In the parallel case this will only push a
     * locally-owned (not halo) node.  Halo nodes will be updated, but never pushed to the local queue
     * @param priority  Current priority/distance of this node.
     * @param nodeIndex  A global node index.
     */
    void PushLocal(double priority, unsigned nodeIndex)
    {

        if (mLo<=nodeIndex && nodeIndex<mHi)
        {
            //Push a negative priority so that the lowest one (nearest the surface) is popped first
            mActivePriorityNodeIndexQueue.push(std::pair<double, unsigned>(-priority, nodeIndex));
        }
    }

public:

    /**
     * Constructor
     *
     * @param rMesh the mesh for which to compute maps
     */
    DistanceMapCalculator(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

     /**
     * Destructor - cleans up mNumHalosPerProcess (which is normally set to NULL anyway).
     */
    ~DistanceMapCalculator()
    {
        delete [] mNumHalosPerProcess;
    }

    /**
     *  Generates a distance map of all the nodes of the mesh to the given source
     *
     *  @param rSourceNodeIndices set of node indices defining the source set or surface
     *         If the vector of source nodes is empty then the results will be a vector of node distance
     *         which are all of size DBL_MAX
     *  @param rNodeDistances distance map computed. The method will resize it if it's not big enough.
     *
     */
    void ComputeDistanceMap(const std::vector<unsigned>& rSourceNodeIndices,
                            std::vector<double>& rNodeDistances);

    /**
     *  @return calculated single point-to-point distance
     *
     *  @param sourceNodeIndex node index for source of distance computation.  Calculations will be cached so
     *  that multiple point-to-point distance computations will get faster.
     *  @param destinationNodeIndex target destination node
     */
    double SingleDistance(unsigned sourceNodeIndex, unsigned destinationNodeIndex);
};

#endif /*DISTANCEMAPCALCULATOR_HPP_*/
