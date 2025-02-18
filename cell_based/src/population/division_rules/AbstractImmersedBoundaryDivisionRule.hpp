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

#ifndef ABSTRACTIMMERSEDBOUNDARYDIVISIONRULE_HPP_
#define ABSTRACTIMMERSEDBOUNDARYDIVISIONRULE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "ImmersedBoundaryCellPopulation.hpp"

// Forward declaration prevents circular include chain
template <unsigned SPACE_DIM>
class ImmersedBoundaryCellPopulation;

/**
 * An abstract cell division rule for use in immersed boundary simulations.
 *
 * The purpose of this class is to return a vector parallel to the new boundary
 * between the daughter cells.
 */
template <unsigned SPACE_DIM>
class AbstractImmersedBoundaryDivisionRule : public Identifiable
{
private:
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
    }

protected:
    /**
     * Output any parameters associated with the division rule. Currently empty
     * since this class has no member variables. Should be overridden by any
     * child classes that have parameters.
     *
     * @param rParamsFile  The stream of the parameter file
     */
    virtual void OutputCellImmersedBoundaryDivisionRuleParameters(out_stream& rParamsFile);

public:
    /**
     * Default constructor.
     */
    AbstractImmersedBoundaryDivisionRule();

    /**
     * Empty destructor.
     */
    virtual ~AbstractImmersedBoundaryDivisionRule();

    /**
     * Return the vector that will divide the two halves of the existing cell
     * to form the boundary between parent and daughter cell.
     *
     * As this method is pure virtual, it must be overridden in subclasses.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The immersed boundary cell population
     *
     * @return the division vector.
     */
    virtual c_vector<double, SPACE_DIM> CalculateCellDivisionVector(
        CellPtr pParentCell,
        ImmersedBoundaryCellPopulation<SPACE_DIM>& rCellPopulation) = 0;

    /**
     * Output the name of the concrete class and call
     * OutputCellImmersedBoundaryDivisionRuleParameters().
     *
     * @param rParamsFile  The stream of the parameter file
     */
    void OutputCellImmersedBoundaryDivisionRuleInfo(out_stream& rParamsFile);
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractImmersedBoundaryDivisionRule)

#endif /*ABSTRACTIMMERSEDBOUNDARYDIVISIONRULE_HPP_*/
