/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTCREATINGANDUSINGANEWCELLPROPERTYTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWCELLPROPERTYTUTORIAL_HPP_

/*
 * = An example showing how to create a new cell property and use it in a cell-based simulation =
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * This tutorial assumes you have already read [wiki:UserTutorials/CreatingAndUsingANewForce].
 *
 * EMPTYLINE
 *
 * In the  [wiki:UserTutorials/CreatingAndUsingANewCellMutationState] we showed how to create a new cell mutation
 * state class, and how this can be used in a cell-based simulation. As well as
 * mutation states, cells may be given much more general properties, using the cell
 * property class hierarchy. In this tutorial, we show how to create a new cell property
 * class, and how this can be used in a cell-based simulation. We will also use a simple
 * new force to illustrate what you can do with cell properties (and also mutations).
 *
 * == 1. Including header files ==
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary
 * header file and archiving headers.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

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
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

/*
 * == Defining the cell property class ==
 *
 * As an example, let us consider a cell property class that is used to label
 * those cells that are "motile". This cell property could then be used when
 * implementing some form of chemotaxis down an imposed chemoattractant gradient,
 * as occurs for example when macrophages migrate within a tumour towards high
 * concentrations of the vascular endothelial growth factor VEGF; for further
 * details, see for example Owen ''et al.'', J. Theor. Biol.
 * 226: 377-391 (2004).
 *
 * Note that usually this code would be separated out into a separate declaration
 * in a .hpp file and definition in a .cpp file.
 */
class MotileCellProperty : public AbstractCellProperty
{
private:

    /* We define a member variable {{{mColour}}}, which can be used by visualization tools
     * to paint cells with this mutation state a distinct colour if required. */
    unsigned mColour;

    /* The next block of code allows us to archive (save or load) the cell property object
     * in a cell-based simulation. The code consists of a serialize() method, in which we first
     * archive the cell property using the serialization code defined in the base class
     * {{{AbstractCellProperty}}}, then archive the member variable {{{mColour}}}. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    /* The default constructor allows us to specify a value for the member variable {{{mColour}}},
     * or leave it with a default value. */
    MotileCellProperty(unsigned colour=5)
        : AbstractCellProperty(),
          mColour(colour)
    {
    }

    /* We then define a destructor and a get method for the member variable {{{mColour}}}. */
    ~MotileCellProperty()
    {}

    unsigned GetColour() const
    {
        return mColour;
    }
};
/* This completes the code for {{{MotileCellProperty}}}.  Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 */

/*
 * EMPTYLINE
 *
 * == Defining the motive force class ==
 *
 * In order to illustrate the use of cell properties we make a simple force law which
 * causes all cells with the {{{MotileCellProperty}}} to move towards the origin. To do this we
 * create a new force class, {{{MyMotiveForce}}}, which inherits from
 * {{{AbstractForce}}} and overrides the methods {{{AddForceContribution()}}} and
 * {{{OutputForceParameters()}}}.
 *
 * Note that usually this code would be separated out into a separate declaration
 * in a .hpp file and definition in a .cpp file.
 */
class MyMotiveForce : public AbstractForce<2>
{
private:

    /* This force class includes a member variable, {{{mStrength}}}, which
     * defines the strength of the force. This member variable will be set
     * in the constructor.
     */
    double mStrength;

    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the force model object in a cell-based simulation.
     * The code consists of a serialize method, in which we first archive the force
     * using the serialization code defined in the base class {{{AbstractForce}}},
     * then archive the member variable. */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mStrength;
    }

public:
    /* The first public method is a default constructor, which calls the base
     * constructor. There is a single input argument, which defines the strength
     * of the force. We provide a default value of 2.0 for this argument. Inside
     * the method, we add an assertion to make sure that the strength is strictly
     * positive.
     */
    MyMotiveForce(double strength=2.0)
        : AbstractForce<2>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    /* The second public method overrides {{{AddForceContribution()}}}.
     * This method takes in two arguments: a reference to a vector of
     * total forces on nodes in a cell population, which is update to by the
     * force object; and a reference to the cell population itself.
     */
    void AddForceContribution(std::vector<c_vector<double, 2> >& rForces,
                              AbstractCellPopulation<2>& rCellPopulation)
    {
        /* Inside the method, we loop over cells, and add a vector to
         * each node associated with cells with the {{{MotileCellProperty}}}, which is proportional (with constant {{{mStrength}}}) to the negative of the position. Causing
         * cells to move inwards towards the origin. Note that this will currently only work with subclasses of {{{AbstractCentreBasedCellPopulation}}}s as
         * we associate cells with nodes in the force calculation. However, this could easily be modified to make it work for {{{VertexBasedCellPopulation}}}s.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            if (cell_iter->HasCellProperty<MotileCellProperty>())
            {
                unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

                c_vector<double, 2> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                rForces[node_index] -= mStrength *  location;
            }
        }
    }

    /* Just as we encountered in [wiki:UserTutorials/CreatingAndUsingANewCellKiller], here we must override
     * a method that outputs any member variables to a specified results file {{{rParamsFile}}}.
     * In our case, we output the member variable {{{mStrength}}}, then call the method on the base class.
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2>::OutputForceParameters(rParamsFile);
    }
};

/* As mentioned in previous cell-based Chaste tutorials, we need to include the next block
 * of code to be able to archive the cell property and force objects in a cell-based simulation,
 * and to obtain a unique identifier for our new classes for when writing results to file.
 *
 * Identifiers for both classes are defined together here, since we can only have each #include once
 * in this source file.  Normally the first include and export would go in the class' header, and the second
 * include and export in the .cpp file for each respective class.
 */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MotileCellProperty)
CHASTE_CLASS_EXPORT(MyMotiveForce)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MotileCellProperty)
CHASTE_CLASS_EXPORT(MyMotiveForce)

/*
 * This completes the code for {{{MyMotiveForce}}}. Note that usually this code
 * would be separated out into a separate declaration in a .hpp file and definition
 * in a .cpp file.
 *
 * EMPTYLINE
 *
 * == The Tests ==
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingANewCellPropertyTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * === Testing the cell property ===
     *
     * We begin by testing that our new cell property is implemented correctly.
     */
    void TestMotileCellProperty() throw(Exception)
    {
        /* We begin by testing that some of the base class methods work correctly.
         * We typically use shared pointers to create and access a cell property
         * like {{{MotileCellProperty}}}, for which it makes sense for all cells
         * that have the same mutation to share a pointer to the same cell property
         * object (although strictly speaking, they are not required to). Observe that
         * in this case we have provided a value for the member variable {{{mColour}}}
         * in the {{{MotileCellProperty}}} constructor.*/
        MAKE_PTR_ARGS(MotileCellProperty, p_property, (8));

        /* Each cell property has a member variable, {{{mCellCount}}}, which
         * stores the number of cells with this cell property. We can test whether
         * {{{mCellCount}}} is being updated correctly by our cell property, as follows. */
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        p_property->IncrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);
        p_property->DecrementCellCount();
        TS_ASSERT_EQUALS(p_property->GetCellCount(), 0u);
        TS_ASSERT_THROWS_THIS(p_property->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");

        /* We can also test whether our cell property is of a given type, as follows. */
        TS_ASSERT_EQUALS(p_property->IsType<WildTypeCellMutationState>(), false);
        TS_ASSERT_EQUALS(p_property->IsType<MotileCellProperty>(), true);

        /* We can also test that archiving is implemented correctly for our cell
         * property, as follows (further details on how to implement and
         * test archiving can be found at ChasteGuides/BoostSerializationW).  */
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "property.arch";

        {
            AbstractCellProperty* const p_const_property = new MotileCellProperty(7);
            p_const_property->IncrementCellCount();

            TS_ASSERT_EQUALS(p_const_property->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(dynamic_cast<MotileCellProperty*>(p_const_property)->GetColour(), 7u);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_const_property;

            delete p_const_property;
        }

        {
            AbstractCellProperty* p_property;

            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_property;

            TS_ASSERT_EQUALS(p_property->GetCellCount(), 1u);

            MotileCellProperty* p_real_property = dynamic_cast<MotileCellProperty*>(p_property);
            TS_ASSERT(p_real_property != NULL);
            TS_ASSERT_EQUALS(p_real_property->GetColour(), 7u);

            delete p_property;
        }
    }

    /*
     * === Using the cell property in a cell-based simulation ===
     *
     * We conclude with a brief test demonstrating how {{{MotileCellProperty}}} can be used
     * in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithMotileCellProperty() throw(Exception)
    {
        /* We use the {{{HoneycombMeshGenerator}}} to create a honeycomb mesh covering a
         * circular domain of given radius, and use this to generate a {{{NodesOnlyMesh}}}
         * as follows. */
        HoneycombMeshGenerator generator(10, 10);
        MutableMesh<2,2>* p_generating_mesh = generator.GetCircularMesh(5);

        NodesOnlyMesh<2>* p_mesh = new(NodesOnlyMesh<2>);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

        /* We now create a shared pointer to our new property, as follows. */
        MAKE_PTR(MotileCellProperty, p_motile);
        /* Also create a shared pointer to a cell label so we can visualise the different cell types.
         * Note that this is also a {{{CellProperty}}}.
         */
        MAKE_PTR(CellLabel, p_label);

        /* Next, we create some cells. We don't use a Generator as we want to give some cells the new cell property, therefore
         * we create the cells in a loop, as follows.*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            /* For each node we create a cell with our cell-cycle model and the wild-type cell mutation state.
             * We then add the property {{{MotileCellProperty}}} to a random selection of the cells, as follows. */
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPropertyCollection collection;
            if (RandomNumberGenerator::Instance()->ranf() < 0.2)
            {
                collection.AddProperty(p_motile);
                collection.AddProperty(p_label);
            }

            CellPtr p_cell(new Cell(p_state, p_model, false, collection));

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
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        /* To run node based simulations you need to define a cut of length, which
         * defines the connectivity of the nodes by defining a radius of interaction.
         */
        cell_population.SetMechanicsCutOffLength(1.5);

        /* In order to visualise labelled cells we need to use the following command.*/
        cell_population.SetOutputCellMutationStates(true);

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple, and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMotileCellProperty");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        /* Now create a {{{MotlieForce}}} and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(MyMotiveForce, p_motive_force);
        simulator.AddForce(p_motive_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* We conclude by deleting any pointers that we created in the test to avoid memory leaks.*/
        delete p_mesh;
    }
    /*
     * When you visualise the results with
     *
     * {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestOffLatticeSimulationWithMotileCellProperty/results_from_time_0}}}
     *
     * you should see a collection of cells with the cell property (labelled dark blue) moving towards the origin.     *
     */
};

#endif /*TESTCREATINGANDUSINGANEWCELLPROPERTYTUTORIAL_HPP_*/
