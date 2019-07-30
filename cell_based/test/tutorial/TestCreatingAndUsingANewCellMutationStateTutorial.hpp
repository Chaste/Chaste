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

#ifndef TESTCREATINGANDUSINGANEWCELLMUTATIONSTATETUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLMUTATIONSTATETUTORIAL_HPP_

/*
 * = An example showing how to create a new cell mutation state and use it in a cell-based simulation =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In the tumour spheroid tutorial we noted that a cell mutation state is always required
 * when constructing a cell. In this tutorial, we show how to create a new cell mutation
 * state class, and how this can be used in a cell-based simulation.
 *
 * EMPTYLINE
 *
 * == 1. Including header files ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary
 * header file and archiving headers.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header defines a base class for cell mutation states. Our new
 * cell mutation state will inherit from this abstract class. */
#include "AbstractCellMutationState.hpp"
/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based
 * Chaste tutorials. */
#include "HoneycombMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellsGenerator.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

/*
 * EMPTYLINE
 *
 * == Defining the cell mutation state class ==
 *
 * As an example, let us consider a cell mutation state representing the p53
 * 172R-H gain-of-function mutant, which is equivalent to the common 175R-H
 * human breast cancer mutant; for further details on this mutant, see for
 * example Murphy et al, FASEB J. 14:2291-2302 (2000).
 *
 * Wild-type p53 has been referred to as the "guardian of the genome",
 * responding to DNA damage or checkpoint failure by either arresting cell
 * cycle progression to facilitate DNA repair or initiating an apoptotic
 * pathway to remove damaged cells. Approximately 40% of human breast cancers
 * contain alterations in p53.
 *
 * As we can see, apart from a serialize() method and a constructor, this class
 * does not contain any member variables or methods. This is because generally
 * a cell's mutation state is used, much like a flag, by other classes when
 * determining a cell's behaviour (whether a cell should undergo
 * apoptosis following prolonged stress, for example, or alter its proliferative
 * behaviour).
 */
class P53GainOfFunctionCellMutationState : public AbstractCellMutationState
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    /* The only public method is a default constructor, which just calls the base
     * constructor with a single unsigned parameter. This sets the value of the
     * base class member variable {{{mColour}}}, which can be used by visualization tools
     * to paint cells with this mutation state a distinct colour if required. */
    P53GainOfFunctionCellMutationState()
        : AbstractCellMutationState(5)
    {
    }
};

/* As mentioned in previous cell-based Chaste tutorials, we need to include the next block
 * of code to be able to archive the cell mutation state object in a cell-based
 * simulation, and to obtain a unique identifier for our new cell mutation state for writing
 * results to file.
 */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(P53GainOfFunctionCellMutationState)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(P53GainOfFunctionCellMutationState)

/*
 * This completes the code for {{{P53GainOfFunctionCellMutationState}}}. Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 *
 * EMPTYLINE
 *
 * === The Tests ===
 *
 * EMPTYLINE
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingANewCellMutationStateTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * EMPTYLINE
     *
     * == Testing the cell mutation state ==
     *
     * EMPTYLINE
     *
     * We begin by testing that our new cell mutation state is implemented correctly.
     */
    void TestP53GainOfFunctionCellMutationState()
    {
        /* We begin by testing that some of the base class methods work correctly.
         * We typically use shared pointers to create and access cell mutation states, as
         * follows. This is because it makes sense for all cells that have the same mutation
         * to share a pointer to the same cell mutation state object (although strictly speaking,
         * they are not required to).*/
        MAKE_PTR(P53GainOfFunctionCellMutationState, p_state);

        /* Each cell mutation state has a member variable, {{{mCellCount}}}, which
         * stores the number of cells with this mutation state. In fact, {{{mCellCount}}}
         * is defined in the class {{{AbstractCellProperty}}}, from which
         * {{{AbstractCellMutationState}}} inherits, as well as other cell properties
         * such as {{{CellLabel}}}. We can test whether {{{mCellCount}}} is being
         * updated correctly by our cell mutation state, as follows. */
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        p_state->IncrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 1u);
        p_state->DecrementCellCount();
        TS_ASSERT_EQUALS(p_state->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_state->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");

        /* We can also test that {{{mColour}}} has been set correctly by our constructor, as follows. */
        TS_ASSERT_EQUALS(p_state->GetColour(), 5u);

        /* We can also test whether our cell mutation state is of a given type, as follows. */
        TS_ASSERT_EQUALS(p_state->IsType<WildTypeCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_state->IsType<P53GainOfFunctionCellMutationState>(), true);

        /* We can also test that archiving is implemented correctly for our cell
         * mutation state, as follows (further details on how to implement and
         * test archiving can be found at ChasteGuides/BoostSerialization).  */
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "p53_mutation.arch";

        {
            AbstractCellProperty* const p_const_state = new P53GainOfFunctionCellMutationState();
            p_const_state->IncrementCellCount();

            TS_ASSERT_EQUALS(p_const_state->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(dynamic_cast<AbstractCellMutationState*>(p_const_state)->GetColour(), 5u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_const_state;

            delete p_const_state;
        }

        {
            AbstractCellProperty* p_arch_state;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_arch_state;

            TS_ASSERT_EQUALS(p_arch_state->GetCellCount(), 1u);

            P53GainOfFunctionCellMutationState* p_real_state = dynamic_cast<P53GainOfFunctionCellMutationState*>(p_arch_state);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 5u);

            delete p_arch_state;
        }
    }

    /*
     * EMPTYLINE
     *
     * == Using the cell mutation state in a cell-based simulation ==
     *
     * EMPTYLINE
     *
     * We conclude with a brief test demonstrating how {{{P53GainOfFunctionCellMutationState}}} can be used
     * in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithP53GainOfFunctionCellMutationState()
    {
        /* We use the {{{HoneycombMeshGenerator}}} to create a honeycomb mesh covering a
         * circular domain of given radius, as follows. */
        HoneycombMeshGenerator generator(10, 10);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(5);

        /* We now create a shared pointer to our new cell mutation state, as follows. */
        MAKE_PTR(P53GainOfFunctionCellMutationState, p_state);

        /* Next, we create some cells, as follows. */
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        /* We now assign the mutation to the 11th and 51st cells.*/
        cells[10]->SetMutationState(p_state);
        cells[50]->SetMutationState(p_state);

        /* Now that we have defined the mesh and cells, we can define the cell population. The constructor
         * takes in the mesh and the cells vector. */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to visualize labelled cells we need to use the following command.*/
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple, and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNewMutationState");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * When you visualize the results with
     *
     * {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestOffLatticeSimulationWithNewMutationState/results_from_time_0}}}
     *
     * you should see two cells in black which are the cells with the new mutation. If we want these cells to behave differently we
     * would need to write an new {{{CellCycleModel}}}, {{{CellKiller}}}, {{{Force}}}, or {{{CellPopulationBoundaryCondition}}}
     * which checks for the new mutation.
     *
     */
};

#endif /*TESTCREATINGANDUSINGANEWCELLMUTATIONSTATETUTORIAL_HPP_*/
