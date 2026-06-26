/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef ATTRACTINGPLANEBOUNDARYCONDITION_HPP_
#define ATTRACTINGPLANEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

/**
 * A boundary condition defined by an (infinite) plane that attracts nearby cells onto it.
 *
 * Unlike PlaneBoundaryCondition, which only pushes back nodes that have crossed
 * to the outward side, this condition also pulls in nodes that lie within
 * mAttractionThreshold of the plane on the inward side, snapping them onto it.
 * This is useful for, e.g., holding cells against a moving grip in a simulated
 * tensile test; the plane can be repositioned during a simulation via
 * SetPointOnPlane().
 *
 * Despite the name, the condition is not specific to 3D: in 2D the 'plane' is a
 * line. It is not implemented in 1D.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class AttractingPlaneBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * A point on the boundary plane.
     */
    c_vector<double, SPACE_DIM> mPointOnPlane;

    /**
     * The outward-facing unit normal to the boundary plane. Cells on the side
     * this points towards are considered to have crossed the plane.
     */
    c_vector<double, SPACE_DIM> mNormalToPlane;

    /**
     * Whether to add a small random perturbation to each node as it is snapped
     * onto the plane, to stop nodes from piling up exactly on it.
     */
    bool mUseJiggledNodesOnPlane;

    /**
     * Attraction band width: inward-side nodes within this distance of the
     * plane are pulled onto it (in addition to any node that has crossed it).
     */
    double mAttractionThreshold;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mUseJiggledNodesOnPlane;
        archive & mAttractionThreshold;
    }

public:
    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param point a point on the boundary plane
     * @param normal a vector normal to the plane, facing outward (i.e. towards
     *               the side cells are kept off); it need not be a unit vector,
     *               as it is normalised internally
     */
    AttractingPlaneBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                                     c_vector<double, SPACE_DIM> point, c_vector<double, SPACE_DIM> normal);

    /** @return the attraction band width (#mAttractionThreshold). */
    double GetAttractionThreshold() const;

    /** @return whether nodes are jiggled as they are snapped onto the plane (#mUseJiggledNodesOnPlane). */
    bool GetUseJiggledNodesOnPlane() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Snap onto the plane every node that has crossed to its outward side or
     * that lies within the attraction threshold on the inward side, optionally
     * jiggling each snapped node (see #mUseJiggledNodesOnPlane).
     *
     * @param rOldLocations the node locations before any boundary conditions are applied (unused)
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations) override;

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile) override;

    /** @return the unit normal to the boundary plane (#mNormalToPlane). */
    const c_vector<double, SPACE_DIM>& rGetNormalToPlane() const;

    /** @return the stored point on the boundary plane (#mPointOnPlane). */
    const c_vector<double, SPACE_DIM>& rGetPointOnPlane() const;

    /**
     * Set the width of the inward-side band within which cells are attracted onto
     * the plane (see #mAttractionThreshold).
     *
     * @param attractionThreshold the attraction band width
     */
    void SetAttractionThreshold(double attractionThreshold);

    /**
     * Move the plane to pass through a new point, keeping the same normal. Allows
     * the boundary to be repositioned over the course of a simulation.
     *
     * @param rPoint a point on the boundary plane
     */
    void SetPointOnPlane(const c_vector<double, SPACE_DIM>& rPoint);

    /**
     * Set whether to jiggle nodes as they are snapped onto the plane, which can
     * help stop them overcrowding on it (see #mUseJiggledNodesOnPlane).
     *
     * @param useJiggledNodesOnPlane whether to jiggle the snapped nodes
     */
    void SetUseJiggledNodesOnPlane(bool useJiggledNodesOnPlane);

    /**
     * Overridden VerifyBoundaryCondition() method, called after
     * ImposeBoundaryCondition() to confirm the condition still holds.
     *
     * @return whether every cell centre is on the inward side of the plane
     */
    bool VerifyBoundaryCondition() override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AttractingPlaneBoundaryCondition)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a AttractingPlaneBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void save_construct_data(
        Archive& ar, const AttractingPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
        ar << p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> point = t->rGetPointOnPlane();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << point[i];
        }
        c_vector<double, SPACE_DIM> normal = t->rGetNormalToPlane();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << normal[i];
        }
    }

    /**
     * De-serialize constructor parameters and initialize a AttractingPlaneBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void load_construct_data(
        Archive& ar, AttractingPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
        ar >> p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> point;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> point[i];
        }
        c_vector<double, SPACE_DIM> normal;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> normal[i];
        }

        // Invoke inplace constructor to initialise instance
        ::new (t) AttractingPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, point, normal);
    }
} // namespace serialization
} // namespace boost

#endif // ATTRACTINGPLANEBOUNDARYCONDITION_HPP_
