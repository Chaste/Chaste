/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef TESTCREATINGANDUSINGNEWCELLBASEDWRITERSTUTORIAL_HPP_
#define TESTCREATINGANDUSINGNEWCELLBASEDWRITERSTUTORIAL_HPP_

/*
 * = An example showing how to create a new cell writer and use it in a cell-based simulation =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * This tutorial assumes you have already read [wiki:UserTutorials/CreatingAndUsingANewCellProperty].
 *
 * EMPTYLINE
 *
 * In the [wiki:UserTutorials/CreatingAndUsingANewCellProperty] we showed how to create a new cell
 * property class, {{{MotileCellProperty}}}, and how this can be used in a cell-based simulation.
 * In this tutorial, we show how to create a new cell writer class, which can be used to output
 * different data from a cell-based simulation. Our example will be a writer class for outputting
 * information about those cells in a simulation that have the {{{MotileCellProperty}}}.
 *
 * == 1. Including header files ==
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header file and
 * archiving headers. For simplicity, we duplicate the code presented in [wiki:UserTutorials/CreatingAndUsingANewCellProperty]
 * that defines the {{{MotileCellProperty}}} class. As before, note that usually this code would be
 * separated out into a declaration in a .hpp file and a definition in a .cpp file.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header defines a base class for cell properties. Our new
 * cell property will inherit from this abstract class. */
#include "AbstractCellProperty.hpp"
/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based
 * Chaste tutorials. */
#include "AbstractForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellMutationStatesWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

class MotileCellProperty : public AbstractCellProperty
{
private:
    unsigned mColour;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    MotileCellProperty(unsigned colour=5)
        : AbstractCellProperty(),
          mColour(colour)
    {
    }

    ~MotileCellProperty()
    {}

    unsigned GetColour() const
    {
        return mColour;
    }
};

/*
 * EMPTYLINE
 *
 * == Defining a cell writer class ==
 *
 * We define a class that writes information about each cell in the population, using
 * the public methods of the population class. We inherit from the base class, {{{AbstractCellWriter}}},
 * which encapsulates how files are written. To specify how this writer should act on each cell we implement
 * the `VisitCell` method.
 *
 * Note that usually this code would be separated out into a declaration in a .hpp file and
 * definition in a .cpp file.
 */
#include "AbstractCellWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellMotilityWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /*
     * The serialize method defines how a cell writer object itself can be written to file.
     * In almost all cases it should just call the base class serializer, using the code below.
     * If the new cell writer class has any data members, they should be serialized in this method.
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /*
     * The constructor method calls the base class constructor, with the name of the output file as
     * a parameter. In this case the filename written to will be "cellmotilityresults.dat".
     */
    CellMotilityWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellmotilityresults.dat")
    {
    }

    /*
     * The implmenetation of the `VisitCell` method defines the data this writer commits to the file.
     * Data can be streamed into the member variable mpOutputStream using the << operator. This data is
     * the written to file.
     * In this example, for each cell `pCell` in the population `pCellPopulation`, we first write its
     * location index, followed by a space, followed by its location co-ordinates (space separated)
     * followed by a 0 or 1 indicating whether the cell has the motile property.
     * The base class handles writing the timestamp and newline for each simulation timestep, so we 
     * don't need to worry about that.
     */
    void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

        c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << cell_location[i] << " ";
        }

        *this->mpOutStream << pCell->HasCellProperty<MotileCellProperty>() << " ";
    }
};

/* 
 * As mentioned in previous cell-based Chaste tutorials, we need to include the next block
 * of code to be able to archive the cell property and writer objects in a cell-based simulation,
 * and to obtain a unique identifier for our new classes for when writing results to file.
 *
 * Identifiers for both classes are defined together here, since we can only have each #include once
 * in this source file.  Normally the first include and export would go in the class' header, and the second
 * include and export in the .cpp file for each respective class.
 */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MotileCellProperty)
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellMotilityWriter)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MotileCellProperty)
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellMotilityWriter)

/*
 * TODO: This completes the code for {{{MotileCellProperty and  {{{CellMotilityWriter}}}. 
 * Note that usually this code would be separated out into a separate declaration in a .hpp 
 * file and definition in a .cpp file.
 *
 * EMPTYLINE
 *
 * == The Test ==
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingNewCellBasedWritersTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * === Using the cell-based writers in a cell-based simulation ===
     *
     * We conclude with a brief test demonstrating how {{{MotileCellProperty}}} and {{{CellMotilityWriter}}}
     * can be used in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithMotileCellPropertyAndWriters() throw(Exception)
    {
        /** We include the next line because HoneycombMeshGenerator, used in this test, is not
         *  yet implemented in parallel. */
        EXIT_IF_PARALLEL;

        /* We use the {{{HoneycombMeshGenerator}}} to create a honeycomb mesh covering a
         * circular domain of given radius, and use this to generate a {{{NodesOnlyMesh}}}
         * as follows. */
        HoneycombMeshGenerator generator(10, 10);
        MutableMesh<2,2>* p_generating_mesh = generator.GetCircularMesh(5);

        NodesOnlyMesh<2> mesh;
        /* We construct the mesh using the generating mesh and a cut-off 1.5 which defines the
         * connectivity in the mesh.
         */
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        /* We now create a shared pointer to our new property, as follows. */
        MAKE_PTR(MotileCellProperty, p_motile);
        /*
         * Also create a shared pointer to a cell label so we can visualize the
         * different cell types. Note that this is also a {{{CellProperty}}}.
         */
        MAKE_PTR(CellLabel, p_label);

        /* Next, we create some cells. We don't use a {{{CellsGenerator}}} as we want to give some cells the new cell property, therefore
         * we create the cells in a loop, as follows.*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            /* For each node we create a cell with our cell-cycle model and the wild-type cell mutation state.
             * We then add the property {{{MotileCellProperty}}} to a random selection of the cells, as follows. */
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();

            CellPropertyCollection collection;
            if (RandomNumberGenerator::Instance()->ranf() < 0.2)
            {
                collection.AddProperty(p_motile);
                collection.AddProperty(p_label);
            }

            CellPtr p_cell(new Cell(p_state, p_model, false, collection));
            p_cell->SetCellProliferativeType(p_diff_type);

            /* Now, we define a random birth time, chosen from [-T,0], where
             * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
             * of a stem cell, and t,,2,, is the basic S+G,,2,,+M phases duration.
             */
            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (p_model->GetStemCellG1Duration()
                                        + p_model->GetSG2MDuration());

            /* Finally, we set the birth time and push the cell back into the vector of cells. */
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Now that we have defined the mesh and cells, we can define the cell population. The constructor
         * takes in the mesh and the cells vector. */
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        /* In order to write cell motility data using our writer, we must add it to the list of writers
         * used by the population, which is done using the templated `AddWriter` method.*/
        cell_population.AddWriter<CellMotilityWriter>();

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple, and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMotileCellPropertyAndWriters");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * When you visualize the results with
     *
     * {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestOffLatticeSimulationWithMotileCellProperty/results_from_time_0}}}
     *
     * you should see a collection of cells with the cell property (labelled dark blue) moving towards the origin.     *
     */
};

#endif /*TESTCREATINGANDUSINGNEWCELLBASEDWRITERSTUTORIAL_HPP_*/
