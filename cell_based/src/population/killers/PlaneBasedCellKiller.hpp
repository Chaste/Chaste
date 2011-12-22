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

#ifndef PLANEBASEDCELLKILLER_HPP_
#define PLANEBASEDCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A cell killer that kills cells if they are outside the domain.
 * defined by a point, mPointOnPlane, and an outward pointing normal, mNormalToPlane.
 * Works for all CellPopulations.
 */
template<unsigned DIM>
class PlaneBasedCellKiller : public AbstractCellKiller<DIM>
{
private:

    /**
     * A point on the plane which nodes cannot cross.
     */
    c_vector<double, DIM> mPointOnPlane;

    /**
     * The outward pointing unit normal to the boundary plane.
     */
    c_vector<double, DIM> mNormalToPlane;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a cell population
     * @param point point on the plane which nodes cannot cross
     * @param normal the outward pointing unit normal to the boundary plane
     */
    PlaneBasedCellKiller(AbstractCellPopulation<DIM>* pCellPopulation,
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
     * Loops over cells and kills cells outside boundary.
     */
    virtual void TestAndLabelCellsForApoptosisOrDeath();

    /**
     * Overridden OutputCellKillerParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlaneBasedCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a PlaneBasedCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const PlaneBasedCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
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
 * De-serialize constructor parameters and initialise a PlaneBasedCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, PlaneBasedCellKiller<DIM> * t, const unsigned int file_version)
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
    ::new(t)PlaneBasedCellKiller<DIM>(p_cell_population, point, normal);
}
}
} // namespace ...

#endif /*PLANEBASEDCELLKILLER_HPP_*/
