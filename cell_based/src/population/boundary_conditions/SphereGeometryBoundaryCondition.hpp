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

#ifndef SPHEREGEOMETRYBOUNDARYCONDITION_HPP_
#define SPHEREGEOMETRYBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A spherical cell population boundary condition class, which restricts nodes to lie
 * on the surface of a sphere in the domain. Although the name of this class suggests
 * it is specific to 3D, it is actually also implemented in 2D, for which it is really
 * a circle geometry boundary condition.
 */
template<unsigned DIM>
class SphereGeometryBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /** The centre of the sphere. */
    c_vector<double, DIM> mCentreOfSphere;

    /** The radius of the sphere. */
    double mRadiusOfSphere;

    /** The maximum distance from the surface of the sphere that cells may be. */
    double mMaximumDistance;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
        archive & mMaximumDistance;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param centre the centre of the sphere
     * @param radius the radius of the sphere
     * @param distance the maximum distance from the surface of the sphere that cells may be (defaults to 1e-5)
     */
    SphereGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                    c_vector<double, DIM> centre,
                                    double radius,
                                    double distance=1e-5);

    /**
     * @return mCentreOfSphere.
     */
    const c_vector<double, DIM>& rGetCentreOfSphere() const;

    /**
     * @return mRadiusOfSphere.
     */
    double GetRadiusOfSphere() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::vector< c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SphereGeometryBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SphereGeometryBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SphereGeometryBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> point = t->rGetCentreOfSphere();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << point[i];
    }

    // Archive other member variables
    double radius = t->GetRadiusOfSphere();
    ar << radius;
}

/**
 * De-serialize constructor parameters and initialize a SphereGeometryBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SphereGeometryBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Retrieve c_vectors one component at a time
    c_vector<double, DIM> point;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> point[i];
    }

    // Retrieve other member variables
    double radius;
    ar >> radius;

    // Invoke inplace constructor to initialise instance
    ::new(t)SphereGeometryBoundaryCondition<DIM>(p_cell_population, point, radius);
}
}
} // namespace ...

#endif /*SPHEREGEOMETRYBOUNDARYCONDITION_HPP_*/
