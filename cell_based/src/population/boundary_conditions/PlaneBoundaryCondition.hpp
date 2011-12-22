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

#ifndef PLANEBOUNDARYCONDITION_HPP_
#define PLANEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A plane cell population boundary condition class, which stops nodes moving through
 * a specified plane in the domain. Although the name of this class suggests it is
 * specific to 3D, it is actually also implemented for 1D and 2D, for which it is
 * really a 'point' and 'line' boundary condition respectively.
 */
template<unsigned DIM>
class PlaneBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /**
     * A point on the boundary plane.
     */
    c_vector<double, DIM> mPointOnPlane;

    /**
     * The outward-facing unit normal vector to the boundary plane.
     */
    c_vector<double, DIM> mNormalToPlane;

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
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param point a point on the boundary plane
     * @param normal the outward-facing unit normal vector to the boundary plane
     */
    PlaneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                           c_vector<double, DIM> point,
                           c_vector<double, DIM> normal);

    /**
     * @return mPointOnPlane.
     */
    const c_vector<double, DIM>& rGetPointOnPlane() const;

    /**
     * @return mNormalToPlane.
     */
    const c_vector<double, DIM>& rGetNormalToPlane() const;

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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlaneBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const PlaneBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> point = t->rGetPointOnPlane();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << point[i];
    }
    c_vector<double, DIM> normal = t->rGetNormalToPlane();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << normal[i];
    }
}

/**
 * De-serialize constructor parameters and initialize a PlaneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, PlaneBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> point;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> point[i];
    }
    c_vector<double, DIM> normal;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> normal[i];
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)PlaneBoundaryCondition<DIM>(p_cell_population, point, normal);
}
}
} // namespace ...

#endif /*PLANEBOUNDARYCONDITION_HPP_*/
