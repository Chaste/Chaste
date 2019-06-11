//
// Created by twin on 14/01/19.
//

#ifndef TESTRUNNINGMULTIODECELLS_HPP_
#define TESTRUNNINGMULTIODECELLS_HPP_

/*
 * = An example showing how to run Delta/Notch simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a growing cell monolayer culture
 * into which a simple model of Delta/Notch signalling is incorporated. This model was developed
 * by Collier et al. ("Pattern formation by lateral inhibition with feedback: a mathematical
 * model of delta-notch intercellular signalling", J. Theor. Biol. 183:429-446) and comprises
 * two ODEs to describe the evolution in concentrations of Delta and Notch in each cell. The ODE
 * for Notch includes a reaction term that depends on the mean Delta concentration among neighbouring
 * cells. Thus in this simulation each cell needs to be able to access information about its
 * neighbours. We use the {{{CellData}}} class to facilitate this, and introduce a subclass
 * of {{{OffLatticeSimulation}}} called {{{DeltaNotchOffLatticeSimulation}}} to handle the updating
 * of {{{CellData}}} at each time step as cell neighbours change.
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous tutorials, we begin by including the necessary header files. We have
 * encountered these files already. Recall that often, either {{{CheckpointArchiveTypes.hpp}}}
 * or {{{CellBasedSimulationArchiver.hpp}}} must be included the first Chaste header.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"


/*
 * The next header file defines a simple subcellular reaction network model that includes the functionality
 * for solving each cell's Delta/Notch signalling ODE system at each time step, using information about neighbouring
 * cells through the {{{CellData}}} class.
 */
#include "DeltaNotchEdgeSrnModel.hpp"

#include "CellEdgeSrnModel.hpp"

#include "CellEdgeDeltaNotchTrackingModifier.hpp"



/* Having included all the necessary header files, we proceed by defining the test class.
 */
class TestVertexCellEdgeDeltaNotchODESimulation : public AbstractCellBasedTestSuite
{
public:

    void TestDeltaNotchEdgeSrnCorrectBehaviour()
    {
        TS_ASSERT_THROWS_NOTHING(DeltaNotchEdgeSrnModel srn_model);

        // Create cell edge srn with four edges
        auto p_cell_edge_srn_model = new CellEdgeSrnModel();
        for(int i = 0 ; i < 4; i++)
        {
            boost::shared_ptr<DeltaNotchEdgeSrnModel> p_delta_notch_edge_srn_model(new DeltaNotchEdgeSrnModel());

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            starter_conditions.push_back(0.5);
            p_delta_notch_edge_srn_model->SetInitialConditions(starter_conditions);
            p_cell_edge_srn_model->AddEdgeSrn(p_delta_notch_edge_srn_model);

        }

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_cell_edge_srn_model, false, CellPropertyCollection()));
        p_cell->SetCellProliferativeType(p_diff_type);
        std::vector<double> p_mean_delta = {1.0, 1.0, 1.0, 1.0};
        p_cell->GetCellEdgeData()->SetItem("mean delta", p_mean_delta);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        // Now updated to initial conditions
        for(unsigned i = 0; i < p_cell_edge_srn_model->GetNumEdgeSrn(); i++)
        {
            auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(p_cell_edge_srn_model->GetEdgeSrn(i));

            TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(), 0.5, 1e-4);
            TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(), 0.5, 1e-4);
        }


        // Now update the SRN
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        double end_time = 10.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();
            p_cell_edge_srn_model->SimulateToCurrentTime();
        }


        // Test converged to steady state
        for(unsigned i = 0; i < p_cell_edge_srn_model->GetNumEdgeSrn(); i++)
        {
            auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(p_cell_edge_srn_model->GetEdgeSrn(i));

            TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetNotch(), 0.9900, 1e-4);
            TS_ASSERT_DELTA(p_delta_notch_edge_srn_model->GetDelta(), 0.0101, 1e-4);
        }

    }

    void TestDeltaNotchEdgeSrnCreateCopy()
    {
        int numEdges = 4;

        auto p_cell_edge_srn_model = new CellEdgeSrnModel();
        for(int i = 0 ; i < numEdges; i++)
        {
            boost::shared_ptr<DeltaNotchEdgeSrnModel> p_delta_notch_edge_srn_model(new DeltaNotchEdgeSrnModel());


            // Set ODE system
            std::vector<double> state_variables;
            state_variables.push_back(2.0);
            state_variables.push_back(3.0);
            p_delta_notch_edge_srn_model->SetOdeSystem(new DeltaNotchOdeSystem(state_variables));
            p_delta_notch_edge_srn_model->SetInitialConditions(state_variables);
            p_cell_edge_srn_model->AddEdgeSrn(p_delta_notch_edge_srn_model);

        }



        // Create a copy
        CellEdgeSrnModel* p_cell_edge_srn_model2 = static_cast<CellEdgeSrnModel*> (p_cell_edge_srn_model->CreateSrnModel());

        for(int i = 0 ; i < numEdges; i++)
        {
            auto p_delta_notch_edge_srn_model = boost::static_pointer_cast<DeltaNotchEdgeSrnModel>(p_cell_edge_srn_model2->GetEdgeSrn(i));
            // Check correct initializations
            TS_ASSERT_EQUALS(p_delta_notch_edge_srn_model->GetNotch(), 2.0);
            TS_ASSERT_EQUALS(p_delta_notch_edge_srn_model->GetDelta(), 3.0);

        }





        // Destroy models
        delete p_cell_edge_srn_model;
        delete p_cell_edge_srn_model2;
    }

    void TestArchiveDeltaNotchSrnModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_edge_srn.arch";

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractSrnModel* p_srn_model = new DeltaNotchSrnModel;

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            // We must create a cell to be able to initialise the cell srn model's ODE system
            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->GetCellData()->SetItem("mean delta", 10.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            p_cell->SetBirthTime(0.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Read mean Delta from CellData
            static_cast<DeltaNotchSrnModel*>(p_srn_model)->UpdateDeltaNotch();
            TS_ASSERT_DELTA(static_cast<DeltaNotchSrnModel*>(p_srn_model)->GetMeanNeighbouringDelta(), 10.0, 1e-12);

            output_arch << p_srn_model;

            // Note that here, deletion of the cell-cycle model and srn is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractSrnModel* p_srn_model;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_srn_model;

            TS_ASSERT_DELTA(static_cast<DeltaNotchSrnModel*>(p_srn_model)->GetMeanNeighbouringDelta(), 10.0, 1e-12);

            delete p_srn_model;
        }
    }




    /*
     * EMPTYLINE
     *
     * == Test 1: a vertex-based monolayer with Delta/Notch signalling ==
     *
     * EMPTYLINE
     *
     * In the first test, we demonstrate how to simulate a monolayer that incorporates
     * Delta/Notch signalling, using a vertex-based approach.
     */
    void TestRunningMultiODECells()
    {
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;

        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(2, 1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* We then create some cells, each with a cell-cycle model, {{{UniformG1GenerationalCellCycleModel}}} and a subcellular reaction network model
         * {{{DeltaNotchEdgeSrnModel}}}, which
         * incorporates a Delta/Notch ODE system, here we use the hard coded initial conditions of 1.0 and 1.0.
         * In this example we choose to make each cell differentiated,
         * so that no cell division occurs. */
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {

            /* Initalise cell cycle */
            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);


            /* Initialise edge based SRN */
            auto p_element = p_mesh->GetElement(elem_index);

            auto p_cell_edge_srn_model = new CellEdgeSrnModel();

            /* We choose to initialise the concentrations to random levels */
            auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
            auto notch_concentration = RandomNumberGenerator::Instance()->ranf();

            double total_edge_length = 0.0;
            for(unsigned i = 0 ; i < p_element->GetNumEdges() ; i ++)
                total_edge_length += p_element->GetEdge(i)->rGetLength();


            /* Gets the edges of the element and create an srn for each edge */
            for(unsigned i = 0 ; i < p_element->GetNumEdges() ; i ++)
            {
                auto p_elem_edge = p_element->GetEdge(i);

                auto p_edge_length = p_elem_edge->rGetLength();
                std::vector<double> initial_conditions;

                /* Initial concentration of delta and notch vary depending on the edge length */
                initial_conditions.push_back( p_edge_length/total_edge_length * delta_concentration);
                initial_conditions.push_back( p_edge_length/total_edge_length * notch_concentration);

                boost::shared_ptr<AbstractOdeSrnModel> p_srn_model(new DeltaNotchEdgeSrnModel());
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrn(p_srn_model);
            }

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. We can make the simulation run for longer to see more patterning by increasing
         * the end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexCellEdgeDeltaNotchODESimulation");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10.0);

        /* Then, we define the modifier class, which automatically updates the values of Delta and Notch within
         * the cells in {{{CellData}}} and passes it to the simulation.*/
        MAKE_PTR(CellEdgeDeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* This modifier assigns target areas to each cell, which are required by the {{{NagaiHondaForce}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        simulator.Solve();
    }


};


#endif //TESTRUNNINGMULTIODECELLS_HPP_
