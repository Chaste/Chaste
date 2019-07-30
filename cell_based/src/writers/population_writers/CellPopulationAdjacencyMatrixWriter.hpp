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

#ifndef CELLPOPULATIONADJACENCYMATRIXWRITER_HPP_
#define CELLPOPULATIONADJACENCYMATRIXWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A class written using the visitor pattern for writing the cell population
 * adjacency (i.e. connectivity) matrix to file.
 *
 * The output file is called cellpopulationadjacency.dat by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellPopulationAdjacencyMatrixWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    CellPopulationAdjacencyMatrixWriter();

    /**
     * Visit the population and write the data.
     *
     * @param pCellPopulation a pointer to the population to visit.
     */
    void VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the adjacency matrix (as a long vector).
     *
     * Outputs a line of tab-separated values of the form:
     * [number of cells] [adjacency 1 1] [adjacency 1 2] [adjacency 1 3]...
     *
     * where, in the case of N cells, the (i,j)th entry of the adjacency matrix corresponds to the
     * (1 + i + N*j)th entry in the line (the additional 1 is for the initial entry, which gives N).
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * If cells i and j are adjacent and neither are labelled (as determined by the CellLabel property),
     * then we have [adjacency i j] = 1; if they are adjacent and both are labelled, then we have
     * [adjacency i j] = 2; if they are adjacent and exactly one is labelled, then we have
     * [adjacency i j] = 3; otherwise, if they are no adjacent, then we have [adjacency i i] = 0.
     *
     * By default we have [adjacency i i] = 0, i.e. cells are not considered to be adjacent to
     * themselves.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPopulationAdjacencyMatrixWriter)

#endif /* CELLPOPULATIONADJACENCYMATRIXWRITER_HPP_ */
