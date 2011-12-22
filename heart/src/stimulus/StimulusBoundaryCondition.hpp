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


#ifndef STIMULUSBOUNDARYCONDITION_HPP_
#define STIMULUSBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"
#include "AbstractStimulusFunction.hpp"

/**
 * Boundary condition defined by an AbstractStimlus object.
 */
template<unsigned SPACE_DIM>
class StimulusBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{
private:
    /** A pointer to the stimulus that is to be applied as a boundary condition */
    AbstractStimulusFunction* mpStimulus;

public:
    /**
     * Create a new boundary condition object.
     *
     * @param pStimulus Stimulus object defining the parameters of the boundary condition
     */
    StimulusBoundaryCondition(AbstractStimulusFunction* pStimulus);

    /**
     * @param rX The point at which this boundary condition is to be evaluated.
     * @return The constant value given in the constructor.
     */
    double GetValue(const ChastePoint<SPACE_DIM>& rX) const;
};

#endif /*STIMULUSBOUNDARYCONDITION_HPP_*/
