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

#ifndef CRYPTPROJECTIONFORCE_HPP_
#define CRYPTPROJECTIONFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A force law for use in crypt projection simulations.
 */
class CryptProjectionForce : public GeneralisedLinearSpringForce<2>
{
    friend class TestCryptProjectionForce;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<2> >(*this);
        archive & mA;
        archive & mB;
        archive & mIncludeWntChemotaxis;
        archive & mWntChemotaxisStrength;
    }

    /**
     * The value of the constant a in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mA;

    /**
     * The value of the constant b in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mB;

    /**
     * Whether to include Wnt-dependent chemotaxis for stem cells.
     */
    bool mIncludeWntChemotaxis;

    /**
     * Strength of Wnt-based chemotactic force.
     */
    double mWntChemotaxisStrength;

    /**
     * Map node indices to 3D locations on the crypt surface.
     */
    std::map<unsigned, c_vector<double, 3> > mNode3dLocationMap;

    /**
     * Fix up the mappings between node indices and 3D locations.
     *
     * @param rCellPopulation the cell population
     */
    void UpdateNode3dLocationMap(AbstractCellPopulation<2>& rCellPopulation);

    /**
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param nodeAGlobalIndex the index of the first node
     * @param nodeBGlobalIndex the index of the second node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double,2> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<2>& rCellPopulation);

public:

    /**
     * Constructor.
     */
    CryptProjectionForce();

    /**
     * Destructor.
     */
    ~CryptProjectionForce();

    /**
     * @return mA.
     */
    double GetA() const;

    /**
     * @return mB.
     */
    double GetB() const;

    /**
     * @return mWntChemotaxisStrength
     */
    double GetWntChemotaxisStrength();

    /**
     * Set mWntChemotaxisStrength.
     *
     * @param wntChemotaxisStrength the new value of mWntChemotaxisStrength
     */
    void SetWntChemotaxisStrength(double wntChemotaxisStrength);

    /**
     * Set mIncludeWntChemotaxis.
     *
     * @param includeWntChemotaxis whether to include Wnt-dependent chemotaxis
     */
    void SetWntChemotaxis(bool includeWntChemotaxis);

    /**
     * Calculates the height of the crypt surface given by
     *     z = f(r) = a*r^b
     * at a point whose 2D position is a distance r from the centre of the cell population.
     * This assumes that the cell population is centred at the origin.
     *
     * @param rNodeLocation
     * @return the z component corresponding to rNodeLocation
     */
    double CalculateCryptSurfaceHeightAtPoint(const c_vector<double,2>& rNodeLocation);

    /**
     * Calculates the derivative df/dr of the crypt surface function z=f(r) at a point
     * whose 2D position is a distance r from the centre of the cell_population, which we assume
     * to be at (0,0).
     *
     * @param rNodeLocation the 2D location of a node
     * @return the gradient
     */
    double CalculateCryptSurfaceDerivativeAtPoint(const c_vector<double,2>& rNodeLocation);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(std::vector<c_vector<double,2> >& rForces,
                              AbstractCellPopulation<2>& rCellPopulation);

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

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CryptProjectionForce)

#endif /*CRYPTPROJECTIONFORCE_HPP_*/
