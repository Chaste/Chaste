/*

Copyright (C) University of Oxford, 2005-2010

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

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

/**
 * A force law for use in crypt projection simulations.
 */
class CryptProjectionForce : public GeneralisedLinearSpringForce<2>
{
    friend class TestForcesNotForRelease;

private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<2> >(*this);
        archive & mA;
        archive & mB;
        archive & mIncludeWntChemotaxis;
    }

    /**
     *  The value of the constant a in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mA;

    /**
     *  The value of the constant b in the definition of the crypt surface
     *      z = f(r) = a*r^b.
     */
    double mB;

    /**
     *  Map node indices to 3D locations on the crypt surface.
     */
    std::map<unsigned, c_vector<double, 3> > mNode3dLocationMap;

    /**
     * Whether to include Wnt-dependent chemotaxis for stem cells.
     */
    bool mIncludeWntChemotaxis;

    /**
     * Fix up the mappings between node indices and 3D locations.
     *
     * @param rTissue the tissue
     */
    void UpdateNode3dLocationMap(AbstractTissue<2>& rTissue);

    /**
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by rCalculateVelocitiesOfEachNode()
     *
     * @param nodeAGlobalIndex the index of the first node
     * @param nodeBGlobalIndex the index of the second node
     * @param rTissue the tissue
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double,2> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractTissue<2>& rTissue);

public :

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
     * Set mIncludeWntChemotaxis.
     *
     * @param includeWntChemotaxis whether to include Wnt-dependent chemotaxis
     */
    void SetWntChemotaxis(bool includeWntChemotaxis);

    /**
     *  Calculates the height of the crypt surface given by
     *      z = f(r) = a*r^b
     *  at a point whose 2D position is a distance r from the centre of the tissue.
     *  This assumes that the tissue is centred at the origin.
     *
     *  @param rNodeLocation
     *
     *  @return the z component corresponding to rNodeLocation
     */
    double CalculateCryptSurfaceHeightAtPoint(const c_vector<double,2>& rNodeLocation);


    /**
     *  Calculates the derivative df/dr of the crypt surface function z=f(r) at a point
     *  whose 2D position is a distance r from the centre of the tissue, which we assume
     *  to be at (0,0).
     *
     *  @param rNodeLocation the 2D location of a node
     *  @return the gradient
     */
    double CalculateCryptSurfaceDerivativeAtPoint(const c_vector<double,2>& rNodeLocation);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double,2> >& rForces,
                              AbstractTissue<2>& rTissue);

};

// Declare identifier for the serializer
#include "TemplatedExport.hpp"
CHASTE_CLASS_EXPORT(CryptProjectionForce)

#endif /*CRYPTPROJECTIONFORCE_HPP_*/
