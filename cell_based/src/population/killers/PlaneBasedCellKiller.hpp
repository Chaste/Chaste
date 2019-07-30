/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
    virtual void CheckAndLabelCellsForApoptosisOrDeath();

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
    Archive & ar, const PlaneBasedCellKiller<DIM> * t, const unsigned int file_version)
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
