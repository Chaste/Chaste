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

#ifndef ABSTRACTPOTTSUPDATERULE_HPP_
#define ABSTRACTPOTTSUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "Identifiable.hpp"

class PottsBasedCellPopulation; // Circular definition

/**
 * An abstract Potts update rule class, for use in cell-based simulations
 * using the cellular Potts model.
 */
template<unsigned DIM>
class AbstractPottsUpdateRule : public Identifiable
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Serialization of singleton objects must be done with care.
     * Before the object is serialized via a pointer, it *MUST* be
     * serialized directly, or an assertion will trip when a second
     * instance of the class is created on de-serialization.
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
    AbstractPottsUpdateRule();

    /**
     * Destructor.
     */
    virtual ~AbstractPottsUpdateRule();

    /**
	 * Calculate the contribution to the Hamiltonian.
	 *
	 * @param currentNodeIndex The index of the current node/lattice site
	 * @param targetNodeIndex The index of the target node/lattice site
	 * @param rCellPopulation The cell population
	 *
	 * @return The difference in the Hamiltonian with the current configuration and
	 * the configuration with the target node having the same spin as the current node.
	 */
    virtual double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                   unsigned targetNodeIndex,
                                                   PottsBasedCellPopulation& rCellPopulation)=0;

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

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractPottsUpdateRule)

#endif /*ABSTRACTPOTTSUPDATERULE_HPP_*/
