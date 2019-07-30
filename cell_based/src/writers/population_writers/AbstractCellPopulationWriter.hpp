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

#ifndef ABSTRACTCELLPOPULATIONWRITER_HPP_
#define ABSTRACTCELLPOPULATIONWRITER_HPP_

#include "AbstractCellBasedWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class MeshBasedCellPopulation;
template<unsigned SPACE_DIM> class CaBasedCellPopulation;
template<unsigned SPACE_DIM> class NodeBasedCellPopulation;
template<unsigned SPACE_DIM> class PottsBasedCellPopulation;
template<unsigned SPACE_DIM> class VertexBasedCellPopulation;

/**
 * Abstract class for a writer that takes information from an AbstractCellPopulation and writes it to file.
 *
 * The key difference between this class and AbstractCellPopulationCountWriter is that writers inheriting
 * from this class ARE compatible with a RoundRobin loop. For such writers, each process can write its part
 * of the information without knowing anything about the other processes as long as they are not interfering
 * with the writing itself.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractCellPopulationWriter : public AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>
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
        archive & boost::serialization::base_object<AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     * @param rFileName the name of the file to write to.
     */
    AbstractCellPopulationWriter(const std::string& rFileName);

    /**
     * Write the header to file.
     *
     * @param pCellPopulation a pointer to the population to be written.
     */
    virtual void WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;

    /**
     * Visit the population and write the data.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)=0;
};

#endif /*ABSTRACTCELLPOPULATIONWRITER_HPP_*/
