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

#ifndef TESTCREATINGANDUSINGANEWSRNMODELTUTORIAL_HPP_
#define TESTCREATINGANDUSINGANEWSRNMODELTUTORIAL_HPP_

/*
 * = An example showing how to create a new subcellular reaction network (SRN) model and use it in a cell-based simulation. =
 *
 * == Introduction ==
 *
 * In the previous cell-based Chaste tutorials, we used existing cell-cycle and SRN models to define how cells
 * proliferate and update and subcellular model. In this tutorial, we show how to create a new SRN model class, and how this
 * can be used in a cell-based simulation.
 *
 * == Including header files ==
 *
 * We begin by including the necessary header files. */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"
/* The next header includes the NEVER_REACHED macro, used in one of the methods below. */
#include "Exception.hpp"

/* The next header defines a base class for ode-based SRN models.
 * Our new SRN model will inherit from this abstract class. */
#include "AbstractOdeSrnModel.hpp"

/* These headers specify the methods to solve the ODE system. */
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/* This header specifies the ODE solvers. */
#include "CellCycleModelOdeSolver.hpp"

/* The following headers are needed for checkpointing. */
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based Chaste
 * tutorials. */
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * == Defining the SRN model and ODE system classes ==
 *
 * As an example, let us consider a SRN model in which we solve a simple ODE
 * dx/dt = -0.25*y
 * dy/dt = x
 * This has exact solution x = A cos 0.5t + B sin 0.5t
 * where A and B are determined by the initial condions.
 *
 * To implement this model we define a new SRN model, {{{MySrnModel}}},
 * which inherits from {{{AbstractOdeSrnModel}}} and
 * contains a {{{MyOdeSystem}}}.
 *
 * Note that usually this code would be separated out into a separate declaration in
 * a .hpp file and definition in a .cpp file.
 */
class MyOdeSystem : public AbstractOdeSystem
{
private:
    friend class boost::serialization::access;
    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the ODE system (and therefore the SRN model) object in a cell-based simulation.
     * The code consists of a serialize method, in which we archive the ODE system
     * using the serialization code defined in the base class
     * {{{AbstractOdeSystem}}}.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:
    MyOdeSystem(std::vector<double> stateVariables=std::vector<double>()) : AbstractOdeSystem(2)
    {
        mpSystemInfo = OdeSystemInformation<MyOdeSystem>::Instance();

        if (stateVariables != std::vector<double>())
        {
            SetStateVariables(stateVariables);
        }
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        rDY[0] = -0.25*rY[1];
        rDY[1] = rY[0];
    }
};

/* As in the ODE tutorials we need to define the ODE system information.
 */
template<>
void OdeSystemInformation<MyOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("x");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}


class MySrnModel : public AbstractOdeSrnModel
{
private:

    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the SRN model object in a cell-based simulation.
     * The code consists of a serialize method, in which we archive the SRN
     * model using the serialization code defined in the base class
     * {{{AbstractOdeSrnModel}}}.
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:
    /**
     * We need to define a protected copy-constructor for use by CreateSrnModel.
     * The only way for external code to create a copy of a SRN model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * Note that the parent SRN model will have had ResetForDivision() called just before CreateSrnModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     */
    MySrnModel(const MySrnModel& rModel)
        : AbstractOdeSrnModel(rModel)
    {
        /*
         * These lines copy the ODE system.
         */
        assert(rModel.GetOdeSystem());
        SetOdeSystem(new MyOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
    }


    /* The first public method is a constructor, which just calls the base
     * constructor.  Note you can include an optional argument to specify the ODE solver.*/
public:

    MySrnModel()
        : AbstractOdeSrnModel(2, boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {

        mpOdeSolver = CellCycleModelOdeSolver<MySrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.1);

        assert(mpOdeSolver->IsSetUp());
    }

    /* The second public method overrides {{{CreateSrnModel()}}}. This is a
     * builder method to create new copies of the SRN model. We call
     * the (protected) copy constructor which creates a copy of the cell cycle model.
     *
     */
    AbstractSrnModel* CreateSrnModel()
    {
        return new MySrnModel(*this);
    }

    /* The third public method overrides {{{Initialise()}}}. */
    void Initialise()
    {
        AbstractOdeSrnModel::Initialise(new MyOdeSystem);
    }

    /* The fourth public method runs the ODEs at each timestep and saves some results to {{{CellData}}}. */
    void SimulateToCurrentTime()
    {
        // run the ODE simulation as needed
        AbstractOdeSrnModel::SimulateToCurrentTime();

        /* this line outputs the ODE system variable to {{{CellData}}}. */
        mpCell->GetCellData()->SetItem("x",mpOdeSystem->rGetStateVariables()[0]);
    }

    /* Finally we define a method to output any parameters in our model, this needs to be included in every SRN model.*/
    void OutputSrnModelParameters(out_stream& rParamsFile)
    {
        // No new parameters to output, so just call method on direct parent class
        AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
    }
};

/* We need to include the next block of code if you want to be able to archive (save or load)
 * the SRN model object in a cell-based simulation. It is also required for writing out
 * the parameters file describing the settings for a simulation - it provides the unique
 * identifier for our new SRN model. Thus every SRN model class must provide this,
 * or you'll get errors when running simulations. */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyOdeSystem)
CHASTE_CLASS_EXPORT(MySrnModel)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(MySrnModel)

/* Since we're defining the new SRN model and ODEs within the test file, we need to include the
 * following stanza as well, to make the code work with newer versions of the Boost libraries.
 * Normally the above export declaration would occur in the SRN model's .hpp file, and
 * the following lines would appear in the .cpp file.  See ChasteGuides/BoostSerialization for
 * more information.
 */
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyOdeSystem)
CHASTE_CLASS_EXPORT(MySrnModel)

/*
 * Need to re-include this after {{{SerializationExportWrapperForCpp.hpp}}}. This is to export the
 * components that would normally be in a seperate cpp file.
 */
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(MySrnModel)

/*
 * This completes the code for {{{MySrnModel}}}. Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 *
 * === The Tests ===
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingANewSrnModelTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * == Testing the SRN model ==
     *
     * We begin by testing that our new cell-cycle model is implemented correctly.
     */
    void TestMySrnModel()
    {
        /* Test that we can construct a {{{MySrnModel}}} object: */
        TS_ASSERT_THROWS_NOTHING(MySrnModel srn_model);

        /* Now we construct and initialise a cell with a {{{MySrnModel}}}.*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
        MySrnModel* p_srn_model = new MySrnModel;
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        /* Now increment time and check the ODE in {{{MySrnModel}}} is solved correctly. */
        double end_time = 10;
        unsigned num_steps = 1000;
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();

            double current_time = SimulationTime::Instance()->GetTime();

            /* Check that the ODE system is solved correctly */
            p_srn_model->SimulateToCurrentTime();

            // Test converged to steady state
            TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("x"), cos(0.5*current_time), 1e-4);
        }

        /* Lastly, we briefly test that archiving of {{{MySrnModel}}} has
         * been implemented correctly. Create an {{{OutputFileHandler}}} and use
         * this to define a filename for the archive.
         */
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "my_srn_model.arch";

        /* Create an output archive. */
        {
            /* Destroy the current instance of {{{SimulationTime}}} and create another instance.
             * Set the start time, end time and number of time steps. */
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 4);

            /* Create a cell with associated srn and cell-cycle model. */
            UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
            AbstractSrnModel* p_srn_model = new MySrnModel;
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            /* Move forward two time steps. */
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();
            /* Solve the SRN. */
            p_srn_model->SimulateToCurrentTime();

            double current_time = 1.5;
            TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("x"), cos(0.5*current_time), 1e-4);

            /* Now archive the cell-cycle model through its cell. */
            CellPtr const p_const_cell = p_cell;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_const_cell;
        }

        /* Now create an input archive. Begin by again destroying the current
         * instance of {{{SimulationTime}}} and creating another instance. Set
         * the start time, end time and number of time steps. note that this is
         * overwritten when you load the archive.
         */
        {
            SimulationTime::Destroy();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.0, 1e-4);

            /* Create a pointer to a cell. */
            CellPtr p_cell;

            /* Create an input archive and restore the cell from the archive. */
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell;

            /* Test that the state of the ODES has been restored correctly. */
            double current_time = 1.5;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), current_time, 1e-4);
            TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("x"), cos(0.5*current_time), 1e-4);

            /* Move forward two more time steps. */
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();
            /* Solve the SRN. */
            p_cell->GetSrnModel()->SimulateToCurrentTime();

            /* Check it's moved on OK */
            current_time = 3.0;
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), current_time, 1e-4);
            TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("x"), cos(0.5*current_time), 1e-4);
        }
    }

    /*
     * == Using the SRN model in a cell-based simulation ==
     *
     * We conclude with a brief test demonstrating how {{{MySrnModel}}} can be used
     * in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithMySrnModel()
    {
        /* We use the honeycomb vertex mesh generator to create a vertex mesh.
         */
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* Next, we create some cells. First, define the cells vector. */
        std::vector<CellPtr> cells;
        /* We must create a shared_ptr to a {{{CellMutationState}}} with which to bestow the cells.
         * We make use of the macro MAKE_PTR to do this: the first argument is the class and
         * the second argument is the name of the shared_ptr. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        /* Then we loop over the nodes. */
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            /* For each node we create a cell with our SRN model and simple Stochastic cell cycle model. */
            UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
            MySrnModel* p_srn_model = new MySrnModel;

            /* We choose to initialise the concentrations to random levels in each cell. */
            std::vector<double> initial_conditions;
            initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
            initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_stem_type);


            /* Now, we define a random birth time, chosen from [-T,0], where
             * T is the typical cell cycle duration
             */
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_cell_cycle_model->GetAverageStemCellCycleTime();
            /* We then set the birth time and push the cell back into the vector of cells. */
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Now that we have defined the mesh and cells, we can define the cell population, forces, areas modifier, and simulation
         * in the same way as the other tutorials. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMySrnModel");
        simulator.SetEndTime(10.0);
        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /* Finally to run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information
     *
     * Load the file {{{/tmp/$USER/testoutput/TestOffLatticeSimulationWithMySrnModel/results_from_time_0/results.pvd}}},
     * and color by {{{x}}}.
     */
};

#endif /*TESTCREATINGANDUSINGANEWSRNMODELTUTORIAL_HPP_*/
