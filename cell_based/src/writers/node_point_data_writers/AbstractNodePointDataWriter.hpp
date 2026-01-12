/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef ABSTRACTNODEPOINTDATAWRITER_HPP_
#define ABSTRACTNODEPOINTDATAWRITER_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "Identifiable.hpp"

#include <string>
#include <vector>

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractOffLatticeCellPopulation;

/**
 * Abstract class for a writer that provides node-wise point data from an AbstractOffLatticeCellPopulation.
 * This is used, for instance, to output node regions to VTK, and is used by cell populations such as
 * ImmersedBoundary and SEM.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractNodePointDataWriter : public Identifiable
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive> void serialize(Archive& archive, const unsigned int version)
    {
        archive & mFieldName;
    }

protected:
    /** The name of the point data field. */
    std::string mFieldName = {};

public:
    /**
     * Default constructor.
     */
    AbstractNodePointDataWriter() = default;

    /**
     * Virtual destructor.
     */
    virtual ~AbstractNodePointDataWriter() = default;

    /**
     * Form a vector of per-node point data values for the given cell population.
     *
     * @param pCellPopulation pointer to the cell population containing the per-node-values
     * @return vector of per-node point data values for the given population
     */
    [[nodiscard]] virtual std::vector<double> GetPointData(AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) const =0;

    /**
     * Get the name of the point data field produced by this writer.
     *
     * @return const reference to the point data field name.
     */
    [[nodiscard]] const std::string& rGetFieldName() const noexcept;
};

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractNodePointDataWriter)

#endif /*ABSTRACTNODEPOINTDATAWRITER_HPP_*/
