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

#ifndef OFFLATTICESIMULATION_HPP_
#define OFFLATTICESIMULATION_HPP_

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractForce.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "AbstractNumericalMethod.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

/**
 * Run an off-lattice 2D or 3D cell-based simulation using an off-lattice
 * cell population.
 *
 * In cell-centre-based cell populations, each cell is represented by a
 * single node (corresponding to its centre), and connectivity is defined
 * either by a Delaunay triangulation or a radius of influence. In vertex-
 * based cell populations, each cell is represented by a polytope
 * (corresponding to its membrane) with a variable number of vertices.
 * Alternative cell populations may be defined by the user.
 *
 * The OffLatticeSimulation is constructed with a CellPopulation, which
 * updates the correspondence between each Cell and its spatial representation
 * and handles cell division (governed by the CellCycleModel associated
 * with each cell). Once constructed, one or more Force laws may be passed
 * to the OffLatticeSimulation object, to define the mechanical properties
 * of the CellPopulation. Similarly, one or more CellKillers may be passed
 * to the OffLatticeSimulation object to specify conditions in which Cells
 * may die, and one or more CellPopulationBoundaryConditions to specify
 * regions in space beyond which Cells may not move.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class OffLatticeSimulation : public AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    friend class TestOffLatticeSimulation;
    friend class TestOffLatticeSimulationWithNodeBasedCellPopulation;

    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mForceCollection;
        archive & mBoundaryConditions;
        archive & mpNumericalMethod;
    }

protected:

    /** The mechanics used to determine the new location of the cells, a list of the forces. */
    std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > > mForceCollection;

    /** List of boundary conditions. */
    std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > > mBoundaryConditions;

    /** The numerical method to use in this simulation. Defaults to the explicit forward Euler method. */
    boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > mpNumericalMethod;

    /**
     * Overridden UpdateCellLocationsAndTopology() method.
     *
     * Calculate forces and update node positions.
     */
    virtual void UpdateCellLocationsAndTopology();

    /**
     * Sends nodes back to the positions given in the input map. Used after a failed step
     * when adaptivity is turned on.
     *
     * @param oldNodeLoctions A map linking nodes to their old positions.
     */
    void RevertToOldLocations(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > oldNodeLoctions);

    /**
     * Applies any boundary conditions.
     *
     * @param oldNodeLoctions Mapping between node indices and old node locations
     */
    void ApplyBoundaries(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > oldNodeLoctions);

    /**
     * Overridden SetupSolve() method to clear the forces applied to the nodes.
     */
    virtual void SetupSolve();

    /**
     * Overridden WriteVisualizerSetupFile() method.
     */
    virtual void WriteVisualizerSetupFile();

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation Reference to a cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading
     *     from an archive)
     */
    OffLatticeSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                         bool deleteCellPopulationInDestructor=false,
                         bool initialiseCells=true);

    /**
     * Add a force to be used in this simulation (use this to set the mechanics system).
     *
     * @param pForce pointer to a force law
     */
    void AddForce(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce);

    /**
     * Remove all the forces.
     */
    void RemoveAllForces();

    /**
     * Add a cell population boundary condition to be used in this simulation.
     *
     * @param pBoundaryCondition pointer to a boundary condition
     */
    void AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> >  pBoundaryCondition);

    /**
     * Method to remove all the cell population boundary conditions
     */
    void RemoveAllCellPopulationBoundaryConditions();

    /**
     * Set the numerical method to be used in this simulation (use this to solve the mechanics system).
     *
     * @param pNumericalMethod pointer to a numerical method object
     */
    void SetNumericalMethod(boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > pNumericalMethod);

    /**
     * @return the current numerical method.
     */
    const boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > GetNumericalMethod() const;

    /**
     * Overridden OutputAdditionalSimulationSetup() method.
     *
     * Output any force, boundary condition or numerical method information.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputAdditionalSimulationSetup(out_stream& rParamsFile);

    /**
     * Overridden OutputSimulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationParameters(out_stream& rParamsFile);

    /**
     * Directly access the forces attached to this simulation, to allow their manipulation after archiving.
     *
     * @return mForceCollection the vector of pointers to forces attached to this simulation
     */
    const std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >& rGetForceCollection() const;
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OffLatticeSimulation.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise an OffLatticeSimulation.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance, middle two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise cells.
    ::new(t)OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>(*p_cell_population, true, false);
}
}
} // namespace

#endif /*OFFLATTICESIMULATION_HPP_*/
