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

#ifndef VERTEXCRYPTBOUNDARYFORCE_HPP_
#define VERTEXCRYPTBOUNDARYFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

/**
 * A boundary force class for use in vertex-based crypt simulations
 * to prevent cells moving below the bottom of the crypt (y=0).
 *
 * The boundary force is taken to be zero for y>0, and proportional to y^2
 * for y<0 (hence is continuously differentiable at y=0). The constant of
 * proportionality is given by the parameter mForceStrength, whose value is
 * set in the constructor.
 */
template<unsigned DIM>
class VertexCryptBoundaryForce  : public AbstractForce<DIM>
{
friend class TestVertexCryptBoundaryForce;

private:

    /** Parameter determining the strength of the force acting on nodes below y=0. */
    double mForceStrength;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mForceStrength;
    }

public:

    /**
     * Constructor.
     *
     * @param forceStrength the force strength
     */
    VertexCryptBoundaryForce(double forceStrength=1.0);

    /**
     * Destructor.
     */
    ~VertexCryptBoundaryForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the boundary force on each node in the vertex-based cell population.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);

    /** @return mForceStrength */
    double GetForceStrength() const;

    /**
     * Outputs force parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexCryptBoundaryForce)

#endif /*VERTEXCRYPTBOUNDARYFORCE_HPP_*/
