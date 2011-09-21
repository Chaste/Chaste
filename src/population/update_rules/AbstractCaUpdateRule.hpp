/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef ABSTRACTCAUPDATERULE_HPP_
#define ABSTRACTCAUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "CaBasedCellPopulation.hpp"
#include "Identifiable.hpp"

/**
 * An abstract CA update rule class.
 */
template<unsigned DIM>
class AbstractCaUpdateRule : public Identifiable
{
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
    }

public:

    /**
     * Default constructor.
     */
    AbstractCaUpdateRule();

    /**
     * Destructor.
     */
    virtual ~AbstractCaUpdateRule();

    /**
     * Calculates the new location for a particular cell.
     *
     * This method must be overridden in concrete classes.
     *
     * @param currentLocationIndex reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     * @param dt timestep of the simulation to calculate probability of movement in current timestep
     */
    virtual unsigned GetNewLocationOfCell(unsigned currentLocationIndex,
                                          CaBasedCellPopulation<DIM>& rCellPopulation,
                                          double dt)=0;

    /**
     * Output update rule to file. Call OutputUpdateRuleParameters() to output
     * all member variables to file.
     *
     * @param rParamsFile a file stream
     */
    void OutputUpdateRuleInfo(out_stream& rParamsFile);

    /**
     * Output update rule parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile a file stream
     */
    virtual void OutputUpdateRuleParameters(out_stream& rParamsFile)=0;
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractCaUpdateRule)

#endif /*ABSTRACTCAUPDATERULE_HPP_*/
