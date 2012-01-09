/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef DIFFUSIONCAUPDATERULE_HPP_
#define DIFFUSIONCAUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCaUpdateRule.hpp"
#include "CaBasedCellPopulation.hpp"

/**
 * An update rule class to model diffusion on the lattice.
 *
 * A diffusion coefficient (defaulting to 1.0) is passed into the constructor.
 * At each time step, a uniform random number r_i is generated for each cell i;
 * if r_i > D*dt, where D denotes the diffusion coefficient and dt denotes the
 * time step, then the cell is moved to a random free neighbouring site (if there
 * are no free neighbouring sites, then the cell is not moved).
 */
template<unsigned DIM>
class DiffusionCaUpdateRule : public AbstractCaUpdateRule<DIM>
{
private:

    /** Diffusion constant. */
    double mDiffusionConstant;

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
        archive & boost::serialization::base_object<AbstractCaUpdateRule<DIM> >(*this);
        archive & mDiffusionConstant;
    }

public:

    /**
     * Constructor.
     *
     * @param diffusionConstant the diffusion constant (defaults to 1.0)
     */
    DiffusionCaUpdateRule(double diffusionConstant=1.0);

    /**
     * Destructor.
     */
    ~DiffusionCaUpdateRule();

    /**
     * Overridden GetNewLocationOfCell() method.
     *
     * This randomly moves the cell in one of the possible directions defined by the cell population method
     * GetNeighbouringNodeIndices().
     *
     * @param currentLocationIndex the current location index corresponding to a cell
     * @param rCellPopulation reference to the cell population
     * @param dt simulation time step, used to calculate the probability of movement
     *
     * @return the new location index of the cell
     */
    unsigned GetNewLocationOfCell(unsigned currentLocationIndex,
                                  CaBasedCellPopulation<DIM>& rCellPopulation,
                                  double dt);

    /** @return mDiffusionConstant */
    double GetDiffusionConstant();

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionCaUpdateRule)

#endif /*DIFFUSIONCAUPDATERULE_HPP_*/
