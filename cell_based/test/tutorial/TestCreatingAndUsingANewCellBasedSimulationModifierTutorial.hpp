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

#ifndef TESTCREATINGANDUSINGANEWCELLBASEDSIMULATIONMODIFIERTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLBASEDSIMULATIONMODIFIERTUTORIAL_HPP_


/*
 * = An example showing how to create a new cell-based simulation modifier and use it in a simulation =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial, we show how to create a new cell-based simulation modifier
 * and use this in a cell-based simulation. The simulation modifier class
 * hierarchy is used to implement setup, update and finalise methods in cell-based
 * simulations.
 *
 * == 1. Including header files ==
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary
 * header file and archiving headers.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header defines a base class for cell-based simulation modifiers.
 * Our new modifier class will inherit from this abstract class. */
#include "AbstractCellBasedSimulationModifier.hpp"
/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based
 * Chaste tutorials. */
#include "AbstractForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellsGenerator.hpp"
#include "UniformCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "RepulsionForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * == Defining the cell-based simulation modifier class ==
 *
 * As an example, let us consider a simulation modifier that, at each simulation
 * time step, calculates each cell's height (y coordinate) in a two-dimensional
 * domain and stores it in in the CellData property as "height". This might be
 * used, for example in cell-based simulations where cell behaviour is dictated
 * through some form of positional information along a tissue axis.
 *
 * Note that usually this code would be separated out into a separate declaration
 * in a .hpp file and definition in a .cpp file.
 * Also, while the abstract simulation modifier class is templated over dimensions,
 * for simplicity we hardcode this concrete modifier class to work in 2D only.
 */
class CellHeightTrackingModifier : public AbstractCellBasedSimulationModifier<2,2>
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<2,2> >(*this);
    }

/* The first public method is a default constructor, which simply calls the base
 * constructor. */
public:

    CellHeightTrackingModifier()
        : AbstractCellBasedSimulationModifier<2,2>()
    {}

    /*
     * The next public method is a destructor, which calls the base destructor.
     */
    ~CellHeightTrackingModifier()
    {}

    /*
     * Next, we override the {{{UpdateAtEndOfTimeStep()}}} method, which specifies what
     * to do to the simulation at the end of each time step. In this class, we simply
     * call the method {{{UpdateCellData()}}} on the cell population; this method is
     * defined later in the class definition.
     */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        UpdateCellData(rCellPopulation);
    }

    /*
     * The next overridden method, {{{SetupSolve()}}}, specifies what to do to the
     * simulation before the start of the time loop. In this class, we call
     * {{{UpdateCellData()}}} on the cell population, just as in
     * {{{UpdateAtEndOfTimeStep()}}}. This is needed because otherwise
     * {{{CellData}}} will not have been fully initialised when we enter
     * the main time loop of the simulation.
     */
    void SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, std::string outputDirectory)
    {

        UpdateCellData(rCellPopulation);
    }

    /*
     * Next, we define the {{{UpdateCellData()}}} method itself. This is a helper
     * method that computes the height (y coordinate) of each cell in the population
     * and stores this in the {{{CellData}}} property.
     */
    void UpdateCellData(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        /*
         * We begin by calling {{{Update()}}} on the cell population, which ensures that
         * it is in a coherent state.
         */
        rCellPopulation.Update();

        /*
         * Next, we iterate over the cell population...
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            /*
             * ...find its height...
             */
            double cell_height = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

            /*
             * ...and store this in the {{{CellData}}} item "height".
             */
            cell_iter->GetCellData()->SetItem("height", cell_height);
        }
    }

    /*
     * Finally, we must override the {{{OutputSimulationModifierParameters()}}} method, which
     * outputs to file any parameters that are defined in the class. In this class, there are
     * no such parameters to output, so we simply call the method defined on the direct
     * parent class (in this case, the abstract class).
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
    }
};

/*
 * This concludes the definition of the {{{CellHeightTrackingModifier}}} class.
 *
 * As mentioned in previous cell-based Chaste tutorials, we need to include the next block
 * of code to be able to archive the simulation modifier object in a cell-based simulation,
 * and to obtain a unique identifier for our new class for when writing results to file.
 *
 * The identifiers for this class are defined together here, since we can only have each
 * #include once in this source file. Normally the first include and export would go in
 * the class's header file, and the second #include and export would go in in the .cpp file.
 */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellHeightTrackingModifier)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CellHeightTrackingModifier)

/*
 * EMPTYLINE
 *
 * == The Tests ==
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingANewCellBasedSimulationModifierTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * === Using the modifier in a cell-based simulation ===
     *
     * We conclude with a brief test demonstrating how {{{CellHeightTrackingModifier}}} can be used
     * in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithCellHeightTrackingModifier()
    {
        /*
         * In this case, we choose to create a small {{{NodeBasedCellPopulation}}} comprising 25 cells.
         * We choose a cut-off for mechanical interactions between cells of 1.5 units and add a
         * simple {{{ReplusionForce}}} to the simulation. We use a {{{UniformCellCycleModel}}}
         * to implement some random proliferation in the simulation.
         */
        HoneycombMeshGenerator generator(2, 2, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithCellHeightTrackingModifier");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

        /*
         * Finally, we add a {{{CellHeightTrackingModifier}}} to the simulation.
         */
        MAKE_PTR(CellHeightTrackingModifier, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
};
/*
 * It is most straightforward to visualize the results of this simulation in Paraview.
 * Load the file {{{/tmp/$USER/testoutput/TestOffLatticeSimulationWithCellHeightTrackingModifier/results_from_time_0/results.pvd}}},
 * and add glyphs to represent cells.
 */

#endif /* TESTCREATINGANDUSINGANEWCELLBASEDSIMULATIONMODIFIERTUTORIAL_HPP_ */
