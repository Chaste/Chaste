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

#ifndef BUSKECOMPRESSIONFORCE_HPP_
#define BUSKECOMPRESSIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "TetrahedralMesh.hpp"

/**
 * A force law employed by Buske et al (2011) in their overlapping spheres
 * model of the intestinal crypt (doi:10.1371/journal.pcbi.1001045).
 *
 * Length is scaled by natural length. \todo does this mean natural radius of a cell? If so at what age? (#1764)
 * Time is in hours.
 * 
 * This class specifically calculates the compression force which forms part of equation (A6) in the Buske paper.
 */
template<unsigned DIM>
class BuskeCompressionForce  : public AbstractForce<DIM>
{
    
    friend class TestForcesNotForRelease;

private:

    /**
     * Get the magnitude of the chemotactic force.
     *
     * @param concentration the local nutrient concentration
     * @param concentrationGradientMagnitude the magnitude of the local nutrient concentration gradient
     */
//    double GetChemotacticForceMagnitude(const double concentration, const double concentrationGradientMagnitude);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }
    
    protected:
    
    /**
     * Compression energy parameter
     *
     * Represented by the parameter K in the model by Buske et al (2011) in
     * their off-lattice model of the intestinal crypt
     * (doi:10.1371/journal.pcbi.1001045).
     *
     * Note: K is the bulk modulus of the spheres.
     */
    double mCompressionEnergyParameter;

public:

    /**
     * Constructor.
     */
    BuskeCompressionForce();

    /**
     * Destructor.
     */
    ~BuskeCompressionForce();

    /**
     * Given a node index, returns the set of neighbouring node indices which
     * lie within a specific interaction distance. (Sum of the two cell radii.)
     *
     * @param index the node index
     * @param rCellPopulation a cell population object
     *
     * @return the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeWithinInteractionDistance(unsigned index, AbstractCellPopulation<DIM>& rCellPopulation);
    
    /**
     * Overridden AddForceContribution() method.
     * 
     * @param rForces a vector of forces on notes
     * @param rCellPopulation a cell population object
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeCompressionForce)

#endif /*BUSKECOMPRESSIONFORCE_HPP_*/
