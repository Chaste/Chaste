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

#ifndef EDGEREMAPINFO_HPP_
#define EDGEREMAPINFO_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>

/**
 * Storage class contains a mapping to the old local edge indices and status of the new edges.
 */
class EdgeRemapInfo {
private:
    /**
     * Contains a mapping to the old local edge indices. Negative value means a new edge
     */
    std::vector<long int> mEdgesMapping;

    /**
     * Status
     * 0 Edge has not changed
     * 1 Edge has been split between two elements
     * 2 Completely new edge was created
     * 3 Edge above or below the current edge was deleted
     * 4 Edge above has been merged into the current one
     */
    std::vector<unsigned int> mEdgeStatus;

    /**
     * Determines how close the new node on the split edges is to the previous (lower) node
     * Value of 0 means the new node is at the same position as the lower node,
     * and value of 1 means that its at the upper node of the edge to be split.
     */
    std::vector<double> mSplitProportions;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mEdgesMapping;
        archive & mEdgeStatus;
        archive & mSplitProportions;
    }
public:
    /**
     * Default constructor. Does nothing.
     */
    EdgeRemapInfo();

    /**
     * Constructor for edge remapping.
     * @param edgesMapping the map between the new edge indices and their local index in the element prior
     * to rearrangement
     * @param edgesStatus status of the edges in the element
     */
    EdgeRemapInfo(const std::vector<long int> &edgesMapping, const std::vector<unsigned int> &edgesStatus);

    /**
     * Here for boost serialization
     */
    ~EdgeRemapInfo();

    /**
     * Contains a mapping to the old local edges index. Negative value means a new edge
     * @return edge map
     */
    std::vector<long int> GetEdgesMapping() const;

    /**
     * @return vector containing the status of each edge
     */
    std::vector<unsigned int> GetEdgesStatus() const;

    /**
     * Gets split proportions. Used in VerteBasedPopulationSrn class
     */
    std::vector<double> GetSplitProportions() const;

    /**
     * Sets split proportions. Used in VertexMeshOperationRecorder class
     * @param thetas
     */
    void SetSplitProportions(const std::vector<double> thetas);
};


#endif //EDGEREMAPINFO_HPP_
