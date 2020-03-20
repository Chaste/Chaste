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

#ifndef VERTEXBASEDPOPULATIONSRN_HPP_
#define VERTEXBASEDPOPULATIONSRN_HPP_

#include <VertexBasedCellPopulation.hpp>
#include "ChasteSerialization.hpp"
#include "SrnCellModel.hpp"
#include "EdgeRemapInfo.hpp"
#include "EdgeOperation.hpp"
#include "VertexElementMap.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

template <unsigned DIM>
class VertexBasedCellPopulation;

/**
 * After topological rearrangments, e.g. T1-3, node merges, or edge splits due to mitosis,
 * junctional srns must be updated accordingly. For example, after edge split due to mitosis
 * a new SRN must be created that inherits modified model variables.
 * This class deals with Cell SRN update after a topological change occurs to the cell.
 */
template <unsigned DIM>
class VertexBasedPopulationSrn
{
private:
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
    }
    VertexBasedCellPopulation<DIM>* mpCellPopulation;
public:
    /**
     * Default constructor/destructor
     */
    VertexBasedPopulationSrn();
    ~VertexBasedPopulationSrn();

    /**
     * Set the cell population
     * @param p_vertex_population
     */
    void SetVertexCellPopulation(VertexBasedCellPopulation<DIM>* p_vertex_population);

    /**
     * This method iterates over edge operations performed on the mesh and
     * calls on RemapCellSrn() with appropriate arguments
     * @param rElementMap the element map to take into account that some elements are deleted after an edge operation is recorded
     */
    void UpdateSrnAfterBirthOrDeath(VertexElementMap& rElementMap);

    /**
     * Remaps cell SRN. If an edge has not been modified, the old edge SRN is mapped to its junction.
     * Otherwise, edge SRNs are updated according to the operation performed.
     * @param parent_srn_edges Edge SRNs before topology is changed
     * @param pSrnCell_1 SrnCellModel of the cell
     * @param pEdgeChange_1 Contains information about which edge changes occurred
     */
    void RemapCellSrn(std::vector<AbstractSrnModelPtr> parent_srn_edges,
                             SrnCellModel* pSrnCell_1,
                             EdgeRemapInfo* pEdgeChange_1);
};
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexBasedPopulationSrn)
#endif /* VERTEXBASEDPOPULATIONSRN_HPP_ */
