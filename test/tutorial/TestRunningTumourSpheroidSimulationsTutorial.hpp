/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef TESTRUNNINGTUMOURSPHEROIDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGTUMOURSPHEROIDSIMULATIONSTUTORIAL_HPP_
/*
 * = An example showing how to run tumour spheroid simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste is used to run discrete tumour
 * spheroid simulations. Like crypt simulations, tumour spheroid simulations
 * include cell cycle models and force laws to determine how cells divide and
 * move. In tumour spheroid simulations, however, these are also coupled to a
 * system of partial differential equations that determine the concentration
 * of specified nutrients (e.g. oxygen) throughout the tissue. Also, unlike
 * in crypt simulation, the tissue grows substantially as the tissue simulation
 * progresses.
 *
 * In summary, the main differences between this tutorial and the crypt simulation
 * tutorials are
 *
 *  * a PDE is defined, to be used in the simulation, and
 *  * a non-periodic mesh is used.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * The first thing to do is include the following header, which allows us
 * to use certain methods in our test (this header file should be included
 * in any Chaste test):
 */
#include <cxxtest/TestSuite.h>
/* This header file defines a helper class for generating a suitable mesh: */
#include "HoneycombMeshGenerator.hpp"
/* These are the classes that will be used in these tests (note that we use a
 * tissue simulation subclass called {{{TissueSimulationWithNutrients}}}):
 */
#include "TissueSimulationWithNutrients.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "CellwiseNutrientSinkPde.hpp"
/* !PetscSetupAndFinalize.hpp must be included in all tests which use Petsc. This is
 * a suite of data structures and routines that are used in the finite element
 * PDE solvers, which is how we solve the nutrient PDE(s).
 */
#include "PetscSetupAndFinalize.hpp"


/* Next, we define the test class, which inherits from {{{CxxTest::TestSuite}}}.
 */
class TestRunningTumourSpheroidSimulationsTutorial : public CxxTest::TestSuite
{
public:
    void TestSpheroidTutorial() throw(Exception)
    {
        /* This first line can be ignored, it's a macro which just says
         * don't run this test if in parallel. */
        EXIT_IF_PARALLEL; // defined in PetscTools.hpp

        /* The first thing to do, as before, is to set up the start time and
         * reset the parameters. */
        SimulationTime::Instance()->SetStartTime(0.0);
        TissueConfig::Instance()->Reset();
        TissueConfig::Instance()->SetHepaOneParameters();

        /* Now we want to create a ''non-periodic'' 'honeycomb' mesh.
         * We use the honeycomb mesh generator, as before, saying 10 cells wide
         * and 10 cells high. Note that the thickness of the ghost nodes layer is
         * 0, ie no ghost nodes, and the {{{false}}} indicates not cylindrical.
         */
        HoneycombMeshGenerator generator(10, 10, 0, false);
        /* Get the mesh. Note we call {{{GetMesh()}}} rather than {{{GetCyclindricalMesh}}},
         * and that a {{{MutableMesh}}} is returned. */
        MutableMesh<2,2>* p_mesh = generator.GetMesh();


        /* Next, we need to create some cells. Unlike before, we don't just use
         * a {{{CellsGenerator}}} class, but do it manually, in a loop. First,
         * define the cells vector. */
        std::vector<TissueCell> cells;
        /* then loop over the nodes... */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            /*.. then create a cell, and giving a particular cell cycle model
             * - {{{SimpleOxygenBasedCellCycleModel}}}.  The cell cycle model is
             * parameterised by the dimension of the problem. */
            TissueCell cell(STEM, HEALTHY, new SimpleOxygenBasedCellCycleModel(2));

            /* We now define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a !HepaOne cell, and t,,2,, is the basic S+G,,2,,+M phases duration.
             */
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  TissueConfig::Instance()->GetHepaOneCellG1Duration()
                                  + TissueConfig::Instance()->GetSG2MDuration() );
            /* .. then we set the birth time and push the cell back into the vector
             * of cells. */
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        /* Now that we have defined the cells, we can define the Tissue. This time it
         * is just a mesh-based tissue (ie not a {{{MeshBasedTissueWithGhostNodes()}}}.
         * Again, the constructor takes in the mesh and the cells vector. */
        MeshBasedTissue<2> tissue(*p_mesh, cells);

        /* Recall that in the Wnt based crypt simulation, we defined a singleton class
         * which cell-cycles used to get the wnt concentration. Here, we do the same kind
         * of thing, but using the singletom {{{CellwiseData}}} class, which stores the
         * value of the current nutrient concentration, for each cell. We have to
         * tell the {{{CellwiseData}}} object how many nodes and variables per node there
         * are (in this case, 1 variable per node, ie the oxygen concentration), and
         * the tissue.
         */
        CellwiseData<2>::Instance()->SetNumNodesAndVars(p_mesh->GetNumNodes(),1);
        CellwiseData<2>::Instance()->SetTissue(tissue);
        /* Then we have to initialise the oxygen concentration for each node (to 1.0), by
         * calling {{{SetValue}}}. This takes in the concentration, and the node
         * which this concentration is for .*/
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellwiseData<2>::Instance()->SetValue(1.0, p_mesh->GetNode(i));
        }

        /* Next we instantiate an instance of the PDE class which we defined above.
         * This will be passed into the simulator. The !CellwiseNutrientSinkPde is
         * a Pde class which inherits from !AbstractLinearEllipticPde, and represents
         * the PDE: u_xx + u_yy = k(x) u, where k(x) = 0.03 (the coefficient below)
         * if x is in a live cell, and k(x)=0 if x is within a apoptotic cell
         */
        CellwiseNutrientSinkPde<2> pde(tissue, 0.03);

        /* We must now create one or more force laws, which determine the mechanics of
         * the tissue. For this test, we assume that a cell experiences a force from each
         * neighbour that can be represented as a linear overdamped spring. Since this
         * model was first proposed in the context of crypt modelling by Meineke ''et al''
         * (Cell Prolif. 34:253-266, 2001), we call this object a
         * {{{GeneralisedLinearSpringForce}}}. We pass a pointer to this force into a vector.
         * Note that we have called the method {{{UseCutoffPoint}}} on the
         * {{{GeneralisedLinearSpringForce}}} before passing it into the collection of force
         * laws - this modifies the force law so that two neighbouring cells do not impose
         * a force on each other if they are located more than 3 units (=3 cell widths)
         * away from each other. This modification is necessary when no ghost nodes are used,
         * for example to avoid artificially large forces between cells that lie close together
         * on the spheroid boundary.
         */
        GeneralisedLinearSpringForce<2> linear_force;
        linear_force.UseCutoffPoint(3);
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&linear_force);

        /*
         * The simulator object for these problems is
         * {{{TissueSimulationWithNutrients}}}. We pass in the tissue, the
         * mechanics system, and the PDE.
         */
        TissueSimulationWithNutrients<2> simulator(tissue, force_collection, &pde);

        /* As with {{{CryptSimulation2d}}} (which inherits from the same base class
         * as {{{TissueSimulationWithNutrients}}}), we can set the output directory
         * and end time. */
        simulator.SetOutputDirectory("SpheroidTutorial");
        simulator.SetEndTime(10.0);

        /* Solve. */
        simulator.Solve();

        /* Finally, call {{{Destroy()}}} on the singleton classes. The results
         * can be visualised as in the previous test. */
        SimulationTime::Destroy();
        CellwiseData<2>::Destroy();
    }
};
#endif /*TESTRUNNINGTUMOURSPHEROIDSIMULATIONSTUTORIAL_HPP_*/
