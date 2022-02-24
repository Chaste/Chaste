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

#ifndef CELLDIVISIONLOCATIONSWRITER_HPP_
#define CELLDIVISIONLOCATIONSWRITER_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractCellPopulationEventWriter.hpp"
#include "ChasteSerialization.hpp"

/**
 * A writer class to output the time, locations and cell id info of cell divisions to a file.
 *
 * The output file is called divisions.dat by default.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellDivisionLocationsWriter : public AbstractCellPopulationEventWriter<ELEMENT_DIM, SPACE_DIM>
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
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellPopulationEventWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:
    /**
     * Default constructor.
     */
    CellDivisionLocationsWriter();

    /**
     * Visit the population and write the data.
     *
     * Outputs [num divisions] lines of tab-separated values of the form:
     *  [Div time 0] [Div 0 x-pos] [Div 0 y-pos] [Div 0 z-pos] [Parent 0 Age] [Parent 0 ID] [Child 0 ID]
     *  [Div time 1] [Div 1 x-pos] ...
     *
     * where [num divisions] denotes the number of Divisions since the previous output step
     * with other imformation about each division event.
     *
     * This line is appended to the output written by AbstractCellBasedEventWriter, which is nothing.
     *
     * @param pCellPopulation a pointer to the population to visit.
     */
    void VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * Just passes through to VisitAnyPopulation
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * Just passes through to VisitAnyPopulation
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * Just passes through to VisitAnyPopulation
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * Just passes through to VisitAnyPopulation
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * Just passes through to VisitAnyPopulation
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellDivisionLocationsWriter)

#endif /* CELLDIVISIONLOCATIONSWRITER_HPP_ */
