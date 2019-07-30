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
     * @return #mCentreOfSphere.
     */
    const c_vector<double, DIM>& rGetCentreOfSphere() const;

    /**
     * @return #mRadiusOfSphere.
     */
    double GetRadiusOfSphere() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

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
    Archive & ar, const SphereGeometryBoundaryCondition<DIM>* t, const unsigned int file_version)
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
