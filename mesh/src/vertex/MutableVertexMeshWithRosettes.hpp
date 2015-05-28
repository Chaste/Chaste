/*

Copyright (c) 2005-2015, University of Oxford.
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
#ifndef MUTABLEVERTEXMESHWITHROSETTES_HPP_
#define MUTABLEVERTEXMESHWITHROSETTES_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "RandomNumberGenerator.hpp"
#include "MutableVertexMesh.hpp"

/**
 * A mutable vertex-based mesh class, which inherits from VertexMesh and allows for local
 * remeshing. This is implemented through simple operations including node merging, neighbour
 * exchange ("T1 swap"), node/edge merging in the case of intersections ("T3 swap") and
 * removal of small triangular elements ("T2 swap").
 *
 * MutableVertexMesh is used as a member of the VertexBasedCellPopulation class to represent
 * the junctional network of cells that forms the basis of simulations of off-lattice
 * vertex-based models.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableVertexMeshWithRosettes : public MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>
{
protected:

    /** Allows test class to access methods */
    friend class TestMutableVertexMeshWithRosettes;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /** The probability that, instead of a T1 swap, the relevant nodes merge to form a protorosette */
    double mProtorosetteFormationProbability;

    /** The probability that, in a given timestep, a protorosette node resolves into two rank-3 nodes */
    double mProtorosetteResolutionProbabilityPerTimestep;

    /** The probability that, in a given timestep, a rosette node resolves into two lower-rank nodes */
    double mRosetteResolutionProbabilityPerTimestep;

    /** A globally accessible random number generator */
    RandomNumberGenerator* mpGen;

    /**
    * Overridden helper method called by IdentifySwapType().
    *
    * Handles the case where a swap involves a junction with more than three adjacent elements.
    * This is implemented in a separate method to allow child classes to override this behaviour
    * and implement junction remodelling with high-order nodes.
    *
    * @param pNodeA one of the nodes to perform the swap with
    * @param pNodeB the other node to perform the swap
    */
    void HandleHighOrderJunctions(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
    * Overridden helper method called by IdentifySwapType().
    *
    * Handles different cases where the nodes involved in a potential swap are both
    * contained in three or fewer elements.  This is implemented in a separate method
    * to allow child classes to override the standard behaviour.
    *
    * @param pNodeA one of the nodes to perform the swap with
    * @param pNodeB the other node to perform the swap
    * @param elem_indices indices of elements touching nodes potentially involved in swap
    * @param case_num the case corresponding to location within IdentifySwapType()
    */
    void HandleAdditionalRemodellingBehaviour(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, std::set<unsigned> elem_indices, unsigned case_num);

    /**
     * Overridden helper method called by ReMesh().
     *
     * Handles any additional ReMeshing operations beyond CheckForSwapsFromShortEdges() and CheckForIntersections().
     * This is implemented in a separate method to allow child classes to override the standard behaviour (see #2664).
     */
    void HandleAdditionalReMeshingBehaviour();

    /**
     * Helper method for ReMesh(), called by HandleHighOrderJunctions().
     *
     * Merge a node of a non-rosette cell with the central node
     * of a rosette, by replacing the node in the non-rosette
     * cell with that of the rosette centre, keeping the rosette
     * centre in the same position.
     *
     * @param pNodeA one of the nodes to perform the merge with
     * @param pNodeB the other node to perform the merge with
     */
    void PerformRosetteRankIncrease(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh(), called by CheckForRosettes().
     *
     * Split protorosette in random direction.
     * Create new node and redistribute nodes along line joining
     * centres of randomly selected cell and cell opposite it.
     *
     * @param pProtorosetteNode node at centre of protorosette
     */
    void PerformProtorosetteResolution(Node<SPACE_DIM>* pProtorosetteNode);

    /**
     * Helper method for ReMesh(), called by CheckForRosettes().
     *
     * Split rosette by removing one cell at random.
     * Create new node positioned along line joining
     * rosette node and centre of randomly selected cell.
     * Rosette node will remain unmoved.
     *
     * @param pRosetteNode node at centre of rosette
     */
    void PerformRosetteRankDecrease(Node<SPACE_DIM>* pRosetteNode);

    /**
     * Helper method for ReMesh().
     *
     * Check whether the mesh contains rosettes or protorosettes, and implement resolution events
     * if necessary.
     */
    void CheckForRosettes();

    /**
     * Serialize the mesh.
     *
     * Note that if you are calling this method (from subclasses) you should archive your
     * member variables FIRST. So that this method can call a ReMesh
     * (to convert from TrianglesMeshReader input format into your native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mProtorosetteFormationProbability;
        archive & mProtorosetteResolutionProbabilityPerTimestep;
        archive & mRosetteResolutionProbabilityPerTimestep;
        archive & boost::serialization::base_object<MutableVertexMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param cellRearrangementRatio ratio between the minimum threshold distance for element
     *                                rearrangement node separation after remeshing (defaults to 1.5)
     * @param protorosetteFormationProbability the probability of a protorosette formation event happening instead of
     *                                a T1 swap
     * @param protorosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a protorosette
     *                                will resolve (similar to the completion of a T1 swap)
     * @param rosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a rosette will
     *                                resolve (reduce the number of cells sharing a common vertex by 1)
     */
    MutableVertexMeshWithRosettes(std::vector<Node<SPACE_DIM>*> nodes,
                                  std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
                                  double cellRearrangementThreshold=0.01,
                                  double t2Threshold=0.001,
                                  double cellRearrangementRatio=1.5,
                                  double protorosetteFormationProbability=1.0,
                                  double protorosetteResolutionProbabilityPerTimestep=0.0,
                                  double rosetteResolutionProbabilityPerTimestep=0.0);

    /**
     * Default constructor for use by serializer.
     */
    MutableVertexMeshWithRosettes();

    /**
     * Destructor.
     */
    virtual ~MutableVertexMeshWithRosettes();

    /**
     * Set method for mProtoRosetteFormationProbability.
     *
     * @param protorosetteFormationProbability the new value of mProtoRosetteFormationProbability
     */
    void SetProtorosetteFormationProbability(double protorosetteFormationProbability);

    /**
     * Set method for mProtoRosetteResolutionProbabilityPerTimestep.
     *
     * @param protorosetteResolutionProbabilityPerTimestep the new value of mProtoRosetteResolutionProbabilityPerTimestep
     */
    void SetProtorosetteResolutionProbabilityPerTimestep(double protorosetteResolutionProbabilityPerTimestep);

    /**
     * Set method for mRosetteResolutionProbabilityPerTimestep.
     *
     * @param rosetteResolutionProbabilityPerTimestep the new value of mRosetteResolutionProbabilityPerTimestep
     */
    void SetRosetteResolutionProbabilityPerTimestep(double rosetteResolutionProbabilityPerTimestep);

    /**
     * Get method for mProtoRosetteFormationProbability.
     *
     * @return mProtoRosetteFormationProbability
     */
    double GetProtorosetteFormationProbability();

    /**
     * Get method for mProtoRosetteResolutionProbabilityPerTimestep.
     *
     * @return mProtoRosetteResolutionProbabilityPerTimestep
     */
    double GetProtorosetteResolutionProbabilityPerTimestep();

    /**
     * Get method for mRosetteResolutionProbabilityPerTimestep.
     *
     * @return mRosetteResolutionProbabilityPerTimestep
     */
    double GetRosetteResolutionProbabilityPerTimestep();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableVertexMeshWithRosettes)

#endif /*MUTABLEVERTEXMESHWITHROSETTES_HPP_*/

