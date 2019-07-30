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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTCELLBASEDDEMOTUTORIAL_HPP_
#define TESTCELLBASEDDEMOTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and cell-based simulations in Chaste =
 *
 * == Introduction ==
 *
 * This tutorial is designed to give you a quick introduction to running cell-based
 * simulations in Chaste. Full details are postponed until later tutorials.
 *
 * We begin with a simple monolayer simulation and see how to:
 *   * change the cell-level model;
 *   * how to impose boundaries;
 *   * how to impose periodic conditions;
 *   * how to specify how to remove cells; and
 *   * how to change cell-cycle models.
 *
 * == The test ==
 *
 * We begin by including the necessary header files. These will be described in detail in
 * subsequent cell-based tutorials.
 */
#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "VoronoiDataWriter.hpp"

#include "FakePetscSetup.hpp"

/* Next, we define the test class which inherits from {{{AbstractCellBasedTestSuite}}}.
 * We inherit from {{{AbstractCellBasedTestSuite}}} rather than {{{CxxTest::TestSuite}}} directly because
 * this class sets up and destroys some singleton objects for us. Singletons are objects that we want to exist only
 * once in each simulation and will be covered in detail in later tutorials.
 * Since we are using {{{AbstractCellBasedTestSuite}}} the singleton {{{SimulationTime}}} is initialised to zero at the beginning of the test and destroyed at the end
 * of the test; {{{RandomNumberGenerator}}} is re-seeded with zero at the beginning and destroyed at the end of the test;
 * and {{{CellPropertyRegistry}}} (which stores {{{CellProperties}}}, you learn about these in a later tutorial
 * [wiki:UserTutorials/CreatingAndUsingANewCellProperty]) is cleared at the beginning of the test.
 * This makes for cleaner code.
 */
class TestCellBasedDemoTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
     * == Test 1 - a basic vertex-based simulation ==
     *
     * In the first test, we run a simple vertex-based simulation of an epithelial monolayer.
     * Each cell in the simulation is assigned a simple stochastic cell-cycle model, the cells will divide randomly and never stop proliferating.
     */
    void TestVertexBasedMonolayer()
    {
        /* The first thing we define is a 2D (specified by the <2,2>) mesh which holds the spatial information of the simulation. To do this we use one of a
         * number of {{{MeshGenerators}}}.*/
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* We now generate a collection of cells. We do this by using a {{{CellsGenerator}}} and we specify the proliferative
         * behaviour of the cell by choosing a {{{CellCycleModel}}}, here we choose a {{{UniformG1GenerationalCellCycleModel}}} where
         * each cell is given a division time, drawn from a uniform distribution when it is created.
         * (Note that here we need to use a phase based cell cycle model so that we can use the target area modifiers which are needed by the vertex
         * based simulations).
         * For a vertex simulation
         * we need as may cells as elements in the mesh.*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        /* We now create a {{{CellPopulation}}} object (passing in the mesh and cells) to connect the mesh and the cells together.
         * Here that is a {{{VertexBasedCellPopulation}}} and the dimension is <2>.*/
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /*
         * We now create an {{{OffLatticeSimulation}}} object and pass in the
         * {{{CellPopulation}}}. We also set some options on the simulation
         * like output directory, output multiple (so we don't visualize every
         * timestep), and end time.
         */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo1");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(20.0);

        /* To specify how cells move around, we create a "shared pointer" to a
         * {{{Force}}} object and pass it to the {{{OffLatticeSimulation}}}. This is done using the MAKE_PTR macro as follows.
         */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* A {{{NagaiHondaForce}}} has to be used together with a child class of {{{AbstractTargetAreaModifier}}}.
         * This modifies the target area of individual cells and thus alters the relative forces
         * between neighbouring cells.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * We are checking that we reached the end time of the simulation
         * with the correct number of cells. If different simulation input parameters are being explored
         * the lines should be removed.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

    /*
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/CellBasedDemo1/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     * EMPTYLINE
     *
     * The {{{make_a_movie}}} script can be used to generate a video based on the results of your simulation.
     * To do this, first visualize the results using {{{Visualize2dVertexCells}}} as described above. Click
     * on the box marked "Output" and play through the whole simulation to generate a sequence of {{{.png}}}
     * images, one for each time step. Next, still in the {{{anim}}} folder, do: {{{./make_a_movie}}}.
     * This reads in the {{{.png}}} files and creates a video file called {{{simulation.mpeg}}}.
     *
     * EMPTYLINE
     *
     * Results can also be visualized using Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * EMPTYLINE
     *
     * == Test 2 - basic node-based simulation ==
     *
     * We next show how to modify the previous test to implement a 'node-based' simulation,
     * in which cells are represented by overlapping spheres (actually circles, since we're
     * in 2D).
     */
    void TestNodeBasedMonolayer()
    {
        /* We now need to create a {{{NodesOnlyMesh}}} we do this by first creating a {{{MutableMesh}}}
         * and passing this to a helper method {{{ConstructNodesWithoutMesh}}} along with a interaction cut off length
         * that defines the connectivity in the mesh.
         */
        HoneycombMeshGenerator generator(2, 2); //**Changed**//
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //**Changed**//
        NodesOnlyMesh<2> mesh; //**Changed**//
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5); //**Changed**//

        /* We create the cells as before, only this time we need one cell per node.*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type); //**Changed**//

        /* This time we create a {{{NodeBasedCellPopulation}}} as we are using a {{{NodesOnlyMesh}}}.*/
        NodeBasedCellPopulation<2> cell_population(mesh, cells);//**Changed**//

        /* We create an {{{OffLatticeSimulation}}} object as before, all we change is the output directory
         * and output results more often as a larger default timestep is used for these simulations. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo2"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12); //**Changed**//
        simulator.SetEndTime(20.0);

        /* We use a different {{{Force}}} which is suitable for node based simulations.
         */
        MAKE_PTR(RepulsionForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        /* In all types of simulation you may specify how cells are removed from the simulation by specifying
         * a {{{CellKiller}}}. You create these in the same was as the {{{Force}}} and pass them to the {{{CellBasedSimulation}}}.
         * Note that here the constructor for {{{RandomCellKiller}}} requires some arguments to be passed to it, therefore we use the
         * {{{MAKE_PTR_ARGS}}} macro.
         */
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.01)); //**Changed**//
        simulator.AddCellKiller(p_cell_killer);

        /* Again we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * Again, we are checking that we reached the end time of the simulation
         * with the correct number of cells.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 7u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

    /*
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/CellBasedDemo2/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * EMPTYLINE
     *
     * As described above, the {{{make_a_movie}}} script can be used to generate a video based on the results of your simulation.
     *
     * == Test 3 - basic mesh-based simulation ==
     *
     * We next show how to modify the previous test to implement a 'mesh-based' simulation,
     * in which cells are represented by their centres and a Voronoi tessellation is used to
     * find nearest neighbours.
     */
    void TestMeshBasedMonolayer()
    {
        /* This time we just create a {{{MutableMesh}}} and use that to specify the spatial locations of cells.*/
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();  //**Changed**//

        /* We create the same number of cells as the previous test.*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        /* This time we create a {{{MeshBasedCellPopulation}}} as we are using a {{{MutableMesh}}}.*/
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells); //**Changed**//

        /* To view the results of this and the subsequent mesh based tutorials in Paraview it is necessary to explicitly
        * generate the required .vtu files. This is detailed in the [wiki:UserTutorials/VisualizingWithParaview] tutorial.
        * Note that the results in Paraview may appear different to those in the java based visualizer. This is related
        * to the different methods used to generate voronoi tesselations in each and is resolved through the use of
        * 'ghost nodes', as shown in the next test. */
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        /* We create an {{{OffLatticeSimulation}}} object as before, all we change is the output directory.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo3"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        /* We use a different {{{Force}}} which is suitable for mesh based simulations.*/
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        /* Again we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

    /*
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the
     * previous test, with the results directory changed from {{{CellBasedDemo2}}} to {{{CellBasedDemo3}}}.
     *
     * == Test 4 - basic mesh-based simulation with ghost nodes ==
     *
     * We next show how to modify the previous test to include 'ghost nodes', which do not
     * correspond to cells but are sometimes needed when using a Voronoi tessellation. We
     * will discuss ghost nodes in more detail in subsequent cell-based tutorials.
     */
    void TestMeshBasedMonolayerWithGhostNodes()
    {
        /* This time we just create a {{{MutableMesh}}} and use that to specify the spatial locations of cells.
         * Here we pass an extra argument to the {{{HoneycombMeshGenerator}}} which adds another 2 rows of
         * nodes round the mesh, known as ghost nodes.*/
        HoneycombMeshGenerator generator(2, 2, 2); //**Changed**//
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* We only want to create cells for non ghost nodes. To find these we get them from the {{{HoneycombMeshGenerator}}}
         * using the method {{{GetCellLocationIndices}}}. We also use a different {{{CellCycleModel}}}. Here we use a
         * {{{TysonNovakCellCycleModel}}} which solves a coupled set of ODEs for each cell to calculate when each cell divides. */
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();//**Changed**//
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<TysonNovakCellCycleModel, 2> cells_generator; //**Changed**//
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type); //**Changed**//

        /* This time we create a {{{MeshBasedCellPopulation}}} as we are using a {{{MutableMesh}}} and have ghost nodes.
         * We also need to pass the indices of non ghost nodes as an extra argument.*/
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); //**Changed**//

        /* Again Paraview output is explicitly requested.*/
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        /* We create an {{{OffLatticeSimulation}}} object as before, all we change is the output directory and the end time.
         * The Tyson Novak model is for yeast cells and therefore cells proliferate much more often and so we run the simulation for
         * less time to keep cell numbers relatively small for this demo.
         */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo4"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(2.0); //**Changed**//

        /* We use the same {{{Force}}} as before and run the simulation in the same way.*/
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 2.0, 1e-10);
    }

    /*
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the
     * previous test, with the results directory changed from {{{CellBasedDemo3}}} to {{{CellBasedDemo4}}}.
     *
     * == Test 5 - basic periodic mesh-based simulation ==
     *
     * We next show how to modify the previous test to implement a periodic boundary to the
     * left and right of the domain.
     */
    void TestMeshBasedMonolayerPeriodic()
    {
        /* We now want to impose periodic boundaries on the domain. To do this we create a {{{Cylindrical2dMesh}}}
         * using a {{{CylindricalHoneycombMeshGenerator}}}.*/
        CylindricalHoneycombMeshGenerator generator(5, 2, 2); //**Changed**//
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh(); //**Changed**//

        /* Again we create one cell for each non ghost node. Note that we have changed back to using a {{{UniformG1GenerationalCellCycleModel}}}.*/
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator; //**Changed**//
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);

        /* We use the same {{{CellPopulation}}}, {{{CellBasedSimulation}}} (only changing the output directory and end time) and {{{Force}}} as before and run the simulation.*/
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /* Again Paraview output is explicitly requested.*/
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo5"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0); //**Changed**//

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 29u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

    /*
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the
     * previous test, with the results directory changed from {{{CellBasedDemo4}}} to {{{CellBasedDemo5}}}.
     *
     * == Test 6 - basic periodic mesh-based simulation with obstructions ==
     *
     * We next show how to modify the previous test to include one
     * or more 'obstructions' within the domain.
     */
    void TestMeshBasedMonolayerPeriodicSolidBottomBoundary()
    {
        /* We make the same {{{Mesh}}}, {{{Cells}}}, {{{CellPopulation}}},
         * {{{CellBasedSimulation}}} and forces as before, all we change is the output directory.*/
        CylindricalHoneycombMeshGenerator generator(5, 2, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_stem_type); //**Changed**//

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellBasedDemo6"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(20.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        /* We now want to impose the condition y>0 on the cells. To do this we create a "shared pointer" to a {{{PlaneBoundaryCondition}}}.
         * Much like the {{{RandomCellKiller}}} earlier we pass arguments to the constructor (a point (0,0) on the plane (line in 2D) and an outward pointing normal to the plane (0,-1) ) using the {{{MAKE_PTR_ARGS}}} macro.*/
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        /* Finally we call the {{{Solve}}} method as in all other simulations.*/
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 23u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }
    /*
     * The results may be visualized using {{{Visualize2dCentreCells}}} as described in the
     * previous test, with the results directory changed from {{{CellBasedDemo5}}} to {{{CellBasedDemo6}}}.
     *
     * == Test 7 - basic Potts-based simulation ==
     *
     * In the final test we show how to modify the earlier tests (using off lattice models) to implement a 'Potts-based' simulation,
     * in which cells are represented by collections of sites on a fixed lattice.
     */
    void TestPottsBasedMonolayer()
    {
        /* In common with the off lattice simulations we begin by creating a mesh. Here we use the {{{PottsMeshGenerator}}}
         * class to generate a {{{PottsMesh}}} each element in the mesh is a collection of lattice sites (represented by nodes at their centres).
         * All the connectivity between lattice sites is defined by the {{{PottsMeshGenerator}}},
         * and there are arguments to make the domains periodic.
         */
        PottsMeshGenerator<2> generator(20, 2, 4, 20, 2, 4); //**Changed**//
        PottsMesh<2>* p_mesh = generator.GetMesh(); //**Changed**//

        /* We generate one cell for each element as in vertex based simulations.*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        /* As we have a {{{PottsMesh}}} we use a {{{PottsBasedCellPopulation}}}. Note here we also change the
         * "temperature" of the Potts simulation to make cells more motile.*/
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);//**Changed**//
        cell_population.SetTemperature(1.0);

        /* As a Potts simulation is restricted to a lattice we create a {{{OnSimulation}}} object and pass in the {{{CellPopulation}}} in much the same
         * way as an {{{OffLatticeSimulation}}} in the above examples. We also set some
         * options on the simulation like output directory and end time.
         */
        OnLatticeSimulation<2> simulator(cell_population);//**Changed**//
        simulator.SetOutputDirectory("CellBasedDemo7"); //**Changed**//
        simulator.SetEndTime(20.0);

        /* In order to specify how cells move around we create "shared pointers" to
         * {{{UpdateRule}}} objects and pass them to the {{{OnLatticeSimulation}}}.
         * This is analogous to {{{Forces}}} in earlier examples.
         */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule); //**Changed**//
        simulator.AddUpdateRule(p_volume_constraint_update_rule); //**Changed**//
        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule); //**Changed**//
        simulator.AddUpdateRule(p_surface_area_update_rule); //**Changed**//
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule); //**Changed**//
        simulator.AddUpdateRule(p_adhesion_update_rule); //**Changed**//

        /* We can add {{{CellKillers}}} as before.*/
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.01));
        simulator.AddCellKiller(p_cell_killer);

        /* Again we run the simulation by calling the {{{Solve}}} method.*/
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }
    /*
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/CellBasedDemo7/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     */
};

#endif /*TESTCELLBASEDDEMOTUTORIAL_HPP_*/
