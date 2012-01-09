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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTCREATINGANDUSINGANEWCELLPOPULATIONBOUNDARYCONDITIONTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLPOPULATIONBOUNDARYCONDITIONTUTORIAL_HPP_

/*
 * = An example showing how to create and use a new cell population boundary condition =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how to create a new cell population boundary condition
 * class to specify a fixed domain within which cells are constrained to lie, and
 * how to use this in a cell-based simulation.
 *
 * EMPTYLINE
 *
 * == 1. Including header files ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header defines a base class for cell population boundary conditions,
 * from which the new class will inherit. */
#include "AbstractCellPopulationBoundaryCondition.hpp"
/* The remaining header files define classes that will be used in the cell-based
 * simulation test. You will have encountered some these files already in previous
 * cell-based Chaste tutorials. */
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SmartPointers.hpp"

/*
 * EMPTYLINE
 *
 * == Defining the cell population boundary condition class ==
 *
 * As an example, let us consider a boundary condition for a two-dimensional cell-based
 * simulation, in which all cells are constrained to lie within the domain given in
 * Cartesian coordinates by 0 <= y <= 5. To implement this we define a cell population
 * boundary condition class, {{{MyBoundaryCondition}}}, which inherits from
 * {{{AbstractCellPopulationBoundaryCondition}}} and overrides the methods
 * {{{ImposeBoundaryCondition()}}}, {{{VerifyBoundaryCondition()}}} and
 * {{{OutputCellPopulationBoundaryConditionParameters()}}}.
 */
class MyBoundaryCondition : public AbstractCellPopulationBoundaryCondition<2>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<2> >(*this);
    }

public:
    /* The first public method is a default constructor, which calls the base
     * constructor. There is a single input argument, a pointer to a cell population.
     */
    MyBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation)
    {
    }

    /* The second public method overrides {{{ImposeBoundaryCondition()}}}.
     * This method is called during the {{{Solve()}}} method in {{{OffLatticeSimulation}}}
     * at the end of each timestep, just after the position of each node
     * in the cell population has been updated according to its equation of motion.
     * The method iterates over all cells in the population, and moves any cell whose
     * centre has y coordinate less than 0 or greater than 5 back into the domain.
     *
     * Implicit in this method is the assumption that, when a node hits the
     * boundary of the domain, it does so inelastically. This means, for example,
     * that a node hitting the boundary at y=0 has its location moved to y=0. A
     * more physically realistic modelling assumption might be to assume that
     * momentum is conserved in the collision.
     *
     * Also implicit in this method is the assumption that we are using a cell-centre
     * based population. If we were using a vertex-based population then each node
     * would correspond not to a cell centre but to a vertex.
     */
    void ImposeBoundaryCondition(const std::vector< c_vector<double, 2> >& rOldLocations)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);
            double y_coordinate = p_node->rGetLocation()[1];

            if (y_coordinate > 5.0)
            {
                p_node->rGetModifiableLocation()[1] = 5.0;
            }
            else if (y_coordinate < 0.0)
            {
                p_node->rGetModifiableLocation()[1] = 0.0;
            }
        }
    }

    /* The third public method overrides {{{VerifyBoundaryCondition()}}}.
     * This method is called during the {{{Solve()}}} method in {{{OffLatticeSimulation}}}
     * at the end of each timestep, just after {{{ImposeBoundaryCondition()}}}, and checks
     * that each cell in the population now satisfies {{{MyBoundaryCondition}}}.
     */
    bool VerifyBoundaryCondition()
    {
        bool condition_satisfied = true;

        for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, 2> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double y_coordinate = cell_location(1);

            if ((y_coordinate < 0.0) || (y_coordinate > 5.0))
            {
                condition_satisfied = false;
                break;
            }
        }
        return condition_satisfied;
    }

    /* Just as we encountered in [wiki:UserTutorials/CreatingAndUsingANewCellKiller], here we must override
     * a method that outputs any member variables to a specified results file {{{rParamsFile}}}.
     * In our case, there are no parameters, so we simply call the method on the base class.
     * Nonetheless, we still need to override the method, since it is pure virtual in the base
     * class.
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
    {
        AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
    }
};

/* As mentioned in previous cell-based Chaste tutorials, we need to include the next block
 * of code to be able to archive the cell population boundary condition object in a cell-based
 * simulation, and to obtain a unique identifier for our new boundary condition for writing
 * results to file.
 */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyBoundaryCondition)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyBoundaryCondition)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const MyBoundaryCondition * t, const BOOST_PFTO unsigned int file_version)
        {
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, MyBoundaryCondition * t, const unsigned int file_version)
        {
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            ::new(t)MyBoundaryCondition(p_cell_population);
        }
    }
}

/*
 * This completes the code for {{{MyBoundaryCondition}}}. Note that usually this code
 * would be separated out into a separate declaration in a .hpp file and definition
 * in a .cpp file.
 *
 * EMPTYLINE
 *
 * === The Tests ===
 *
 * EMPTYLINE
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingANewCellPopulationBoundaryConditionTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the cell population boundary condition ==
     *
     * EMPTYLINE
     *
     * We now test that our new cell population boundary condition is implemented correctly.
     */
    void TestMyBoundaryCondition() throw(Exception)
    {
        /* We first create a {{{MeshBasedCellPopulation}}} using the helper
         * classes {{{HoneycombMeshGenerator}}} and {{{CellsGenerator}}},
         * as in previous cell-based Chaste tutorials.
         */
        HoneycombMeshGenerator generator(7, 7);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We now use the cell population to construct a cell population boundary condition object.
         */
        MyBoundaryCondition bc(&cell_population);

        /* We start by verifying that some cells do not satisfy the boundary condition:
         */
        bool population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, false);

        std::vector<c_vector<double, 2> > old_node_locations;
        old_node_locations.reserve(p_mesh->GetNumNodes());
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            old_node_locations[node_index] = cell_population.GetNode(node_index)->rGetLocation();
        }

        /* To test that we have implemented the cell population boundary condition correctly,
         * we call the overridden method {{{ImposeBoundaryCondition()}}}...
         */
        bc.ImposeBoundaryCondition(old_node_locations);

        /* ... and check that the cell population does indeed now satisfy the boundary condition:
         */
        population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, true);

        /* The last block of code provides an archiving test for the cell population boundary
         * condition, in a similar way to previous cell-based Chaste tutorials:
         */
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "my_bc.arch";
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            AbstractCellPopulationBoundaryCondition<2>* const p_bc = new MyBoundaryCondition(NULL);
            output_arch << p_bc;
            delete p_bc;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCellPopulationBoundaryCondition<2>* p_bc;
            input_arch >> p_bc;

            delete p_bc;
        }
    }

    /*
     * == Using the boundary condition in a cell-based simulation ==
     *
     * We now provide a test demonstrating how {{{MyBoundaryCondition}}} can be used
     * in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithMyBoundaryCondition() throw(Exception)
    {
        /* Once again we create a {{{MeshBasedCellPopulation}}}. */
        HoneycombMeshGenerator generator(7, 7, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We use the cell population to construct a cell population boundary condition object. */
        MAKE_PTR_ARGS(MyBoundaryCondition, p_bc, (&cell_population));

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple, and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMyBoundaryCondition");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(1.0);

        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        /* We now pass the cell population boundary condition into the cell-based simulation. */
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * When you visualize the results with
     *
     * {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestOffLatticeSimulationWithMyBoundaryCondition/results_from_time_0}}}
     *
     * you should see that cells are restricted to the domain 0 <= y <= 5.
     *
     */
};

#endif /*TESTCREATINGANDUSINGANEWCELLPOPULATIONBOUNDARYCONDITIONTUTORIAL_HPP_*/
