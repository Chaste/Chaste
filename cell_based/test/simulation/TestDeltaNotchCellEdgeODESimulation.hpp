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
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "AlwaysDivideCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"


/*
 * The next header file defines a simple subcellular reaction network model that includes the functionality
 * for solving each cell's Delta/Notch signalling ODE system at each time step, using information about neighbouring
 * cells through the {{{CellEdgeData}}} class.
 */
#include "DeltaNotchEdgeSrnModel.hpp"
#include "SrnCellModel.hpp"
#include "CellEdgeDeltaNotchTrackingModifier.hpp"



/* Having included all the necessary header files, we proceed by defining the test class.
 */
class TestDeltaNotchCellEdgeODESimulation : public AbstractCellBasedTestSuite
{
public:



    /*
     * EMPTYLINE
     *
     * == A running simulation of a vertex-based monolayer with Delta/Notch signalling with srn cell edge representation ==
     *
     * EMPTYLINE
     *
     * Test showing how to put together a vertex-based monolayer with Delta/Notch signalling with srn cell edge representation.
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

            auto p_cell_edge_srn_model = new SrnCellModel();

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

                MAKE_PTR(DeltaNotchEdgeSrnModel, p_srn_model);
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
        simulator.SetOutputDirectory("TestDeltaNotchCellEdgeODESimulation");
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
