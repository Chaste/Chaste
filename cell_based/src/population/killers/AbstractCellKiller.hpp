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

#ifndef ABSTRACTCELLKILLER_HPP_
#define ABSTRACTCELLKILLER_HPP_

#include "AbstractCellPopulation.hpp"

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "Identifiable.hpp"

/**
 * An abstract cell killer class, for use in cell-based simulations.
 */
template<unsigned SPACE_DIM>
class AbstractCellKiller : public Identifiable
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
        // Archiving of mpCellPopulation is implemented in load_construct_data of subclasses
    }

protected:

    /** The cell population. */
    AbstractCellPopulation<SPACE_DIM>* mpCellPopulation;

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population.
     */
    AbstractCellKiller(AbstractCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Destructor.
     */
    virtual ~AbstractCellKiller();

    /**
     * Pure method which should call StartApoptosis() on any cell
     * which should be about to undergo programmed death, or Kill()
     * on any cell which should die immediately.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void CheckAndLabelCellsForApoptosisOrDeath()=0;

    /**
     * Get a pointer to the cell population.
     *
     * @return A const pointer to the mpCellPopulation
     */
    const AbstractCellPopulation<SPACE_DIM>* GetCellPopulation() const;

    /**
     * Outputs cell killer used in the simulation to file and then calls
     * OutputCellKillerParameters to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerInfo(out_stream& rParamsFile);

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellKillerParameters(out_stream& rParamsFile)=0;
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractCellKiller)

#endif /*ABSTRACTCELLKILLER_HPP_*/
