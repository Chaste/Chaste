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

#ifndef REPULSIONFORCE_HPP_
#define REPULSIONFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"

/**
 * A class for a simple two-body repulsion force law. Designed
 * for use in node-based simulations
 *
 * The force just creates a linear repulsive force between cells
 * with a nonlinear separation less than 2. This force does not
 * take a cell's age or cell cycle phase into account.
 */
template<unsigned DIM>
class RepulsionForce : public GeneralisedLinearSpringForce<DIM>
{
private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM> >(*this);
    }

public :

    /**
     * Constructor.
     */
    RepulsionForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the CellPopulation
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
            AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Outputs force Parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RepulsionForce)

#endif /*REPULSIONFORCE_HPP_*/
