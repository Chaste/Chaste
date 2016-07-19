/*

Copyright (c) 2005-2016, University of Oxford.
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

// Needed for the test environment
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Includes from trunk
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CellVolumesWriter.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DiagonalVertexBasedDivisionRule.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"

// Includes from Immersed Boundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"
#include "ImmersedBoundaryCellCellInteractionForce.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

///\todo Vary the cell population geometry across tests
class TestImmersedBoundaryCellPopulation : public AbstractCellBasedTestSuite
{
public:

    void TestGetAndSetMethods() throw(Exception)
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // Test that GetDampingConstant() returns the correct value
        for (unsigned node_index=0; node_index<cell_population.GetNumNodes(); node_index++)
        {
            TS_ASSERT_DELTA(cell_population.GetDampingConstant(node_index), 0.0, 1e-6);
        }

        // Test that GetIntrinsicSpacing() returns the correct value
        TS_ASSERT_DELTA(cell_population.GetIntrinsicSpacing(), 0.01, 1e-6);

        // Test GetInteractionDistance() and SetInteractionDistance() work correctly
        TS_ASSERT_DELTA(cell_population.GetInteractionDistance(), 0.0139, 1e-4);

        cell_population.SetInteractionDistance(0.1234);
        TS_ASSERT_DELTA(cell_population.GetInteractionDistance(), 0.1234, 1e-6);

        // Test DoesPopulationHaveActiveSources() and SetIfPopulationHasActiveSources() work correctly
        TS_ASSERT_EQUALS(cell_population.DoesPopulationHaveActiveSources(), false);

        cell_population.SetIfPopulationHasActiveSources(true);
        TS_ASSERT_EQUALS(cell_population.DoesPopulationHaveActiveSources(), true);
    }

    void TestMeshMethods() throw(Exception)
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // Test that GetNumElements() and GetNumNodes() return the correct values
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 6u); ///\todo Is this 6, not 5, because of the basement lamina?
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 583u);

        // Test that rGetMesh() returns the mesh correctly
        ImmersedBoundaryMesh<2, 2>& r_mesh = cell_population.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 583u);
        TS_ASSERT_EQUALS(r_mesh.GetNumGridPtsX(), 256u);
        TS_ASSERT_EQUALS(r_mesh.GetNumGridPtsY(), 256u);
        TS_ASSERT_DELTA(r_mesh.GetCharacteristicNodeSpacing(), 0.0115, 1e-4);
        TS_ASSERT_DELTA(r_mesh.GetSpacingRatio(), 2.9610, 1e-4);

        // Test that GetElement() returns an element correctly
        ImmersedBoundaryElement<2, 2>* p_element_0 = cell_population.GetElement(0);
        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 83u);

        // Note: we cannot yet test the corner nodes or average node spacing of this element;
        //       these are set in the ImmersedBoundarMesh method DivideElement()

        // Test that GetElementCorrespondingToCell() returns an element correctly
        CellPtr p_cell_0 = *(cell_population.Begin());
        ImmersedBoundaryElement<2, 2>* p_element_0_again = cell_population.GetElementCorrespondingToCell(p_cell_0);
        TS_ASSERT_EQUALS(p_element_0_again, p_element_0);

        // Test that GetNode() returns a node correctly
        Node<2>* p_node_3 = cell_population.GetNode(3);
        TS_ASSERT_EQUALS(p_node_3->IsBoundaryNode(), true);
        TS_ASSERT_DELTA(p_node_3->rGetLocation()[0], 0.0421, 1e-4);
        TS_ASSERT_DELTA(p_node_3->rGetLocation()[1], 0.2910, 1e-4);

        // Test GetLocationOfCellCentre() returns the correct location
        c_vector<double, 2> cell_0_centre = cell_population.GetLocationOfCellCentre(p_cell_0);
        TS_ASSERT_DELTA(cell_0_centre[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(cell_0_centre[1], 0.0, 1e-6);

        ///\todo Test GetNeighbouringNodeIndices() and GetNeighbouringLocationIndices() and implement these methods

        // Test GetWidth() returns the correct values
        TS_ASSERT_DELTA(cell_population.GetWidth(0), 0.9891, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetWidth(1), 0.4371, 1e-4);

        // Test GetVolumeOfCell() returns the correct value
        CellPtr p_cell_1 = *(++(cell_population.Begin()));
        double cell_1_volume = cell_population.GetVolumeOfCell(p_cell_1);
        TS_ASSERT_DELTA(cell_1_volume, 0.0769, 1e-4);

        ///\todo With this mesh generator, does cell 0 correspond to the basement lamina?

        // Test SetNode() works correctly
        c_vector<double, 2> new_location = cell_population.GetNode(0)->rGetLocation();
        new_location[0] += 1e-2;
        new_location[1] += 1e-2;
        ChastePoint<2> new_location_point(new_location);
        cell_population.SetNode(0, new_location_point);

        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[0], new_location[0], 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[1], new_location[1], 1e-12);
    }

    ///\todo Test AddNode(), UpdateNodeLocations(), AddCell(), RemoveDeadCells(), IsCellAssociatedWithADeletedLocation() and Update()

    void TestVertexBasedDivisionRuleMethods() throw (Exception)
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        CellPtr p_cell_0 = *(cell_population.Begin());

        // Test GetVertexBasedDivisionRule() correctly returns a ShortAxisVertexBasedDivisionRule
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_rule = cell_population.GetVertexBasedDivisionRule();
        TS_ASSERT(p_rule != NULL);
        ///\todo The following lines throw an error if uncommented, since CalculateCellDivisionVector() requires a VertexBasedCellPopulation
//        c_vector<double, 2> division_vector = p_rule->CalculateCellDivisionVector(p_cell_0, cell_population);
//        TS_ASSERT_DELTA(division_vector[0], 1.0, 1e-6);
//        TS_ASSERT_DELTA(division_vector[1], 1.0, 1e-6);

        // Set the division rule for our population to be the diagonal division rule
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_rule_to_set(new DiagonalVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_rule_to_set);

        // Test GetVertexBasedDivisionRule() now correctly returns a DiagonalVertexBasedDivisionRule
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_new_rule = cell_population.GetVertexBasedDivisionRule();
        TS_ASSERT(p_new_rule != NULL);
        ///\todo The following lines throw an error if uncommented, since CalculateCellDivisionVector() requires a VertexBasedCellPopulation
//        c_vector<double, 2> new_division_vector = p_new_rule->CalculateCellDivisionVector(p_cell_0, cell_population);
//        TS_ASSERT_DELTA(new_division_vector[0], 1.0, 1e-9);
//        TS_ASSERT_DELTA(new_division_vector[1], 1.0, 1e-9);
    }

    void TestCalculateCellDivisionVector() throw (Exception)
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        CellPtr p_cell_0 = *(cell_population.Begin());

        // Test CalculateCellDivisionVector() returns the correct vector
        c_vector<double, 2> division_vector = cell_population.CalculateCellDivisionVector(p_cell_0);
        TS_ASSERT_DELTA(division_vector[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(division_vector[1], 0.0, 1e-6);
    }

    ///\todo Check output files by eye too
    void TestWritersWithImmersedBoundaryCellPopulation() throw (Exception)
    {
        // Set up SimulationTime (needed if VTK is used)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        boost::shared_ptr<AbstractCellProperty> p_stem(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_transit(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_diff(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_wildtype(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        std::vector<CellPtr> cells;
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            UniformlyDistributedCellCycleModel* p_model = new UniformlyDistributedCellCycleModel();

            CellPtr p_cell(new Cell(p_wildtype, p_model));
            if (elem_index%3 == 0)
            {
                p_cell->SetCellProliferativeType(p_stem);
            }
            else if (elem_index%3 == 1)
            {
                p_cell->SetCellProliferativeType(p_transit);
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff);
            }

            double birth_time = 0.0 - elem_index;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.InitialiseCells();
        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "ImmersedBoundaryCellPopulation-2");

        // Allocate some cells to have a different cell mutation state, cell label or apoptotic cell property
        cell_population.GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_apc1(cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(cell_population.GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(cell_population.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(cell_population.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(1)->SetMutationState(p_apc1);
        cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);;
        cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);;
        cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_apoptotic_state);;
        cell_population.SetCellAncestorsToLocationIndices();

        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellLocationIndexWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        // Coverage of writing CellData to VTK
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("var0", 0.0);
            cell_iter->GetCellData()->SetItem("var1", 3.0);
        }

        std::string output_directory = "TestWritersWithImmersedBoundaryCellPopulation";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteResultsToFiles(output_directory);

        SimulationTime::Instance()->IncrementTimeOneStep();
        cell_population.Update();
        cell_population.WriteResultsToFiles(output_directory);

        cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizelements", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.vizelements").CompareFiles();

        FileComparison(results_dir + "cellages.dat", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/cellages.dat").CompareFiles();
        FileComparison(results_dir + "results.vizancestors", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.vizancestors").CompareFiles();
        FileComparison(results_dir + "loggedcell.dat", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/loggedcell.dat").CompareFiles();
        FileComparison(results_dir + "results.vizlabels", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.vizlabels").CompareFiles();
        FileComparison(results_dir + "results.vizlocationindices", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.vizlocationindices").CompareFiles();
        FileComparison(results_dir + "results.vizmutationstates", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.vizmutationstates").CompareFiles();
        FileComparison(results_dir + "results.vizcellphases", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.vizcellphases").CompareFiles();
        FileComparison(results_dir + "results.vizcelltypes", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.vizcelltypes").CompareFiles();
        FileComparison(results_dir + "cellareas.dat", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/cellareas.dat").CompareFiles();

        FileComparison(results_dir + "cellmutationstates.dat", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/cellmutationstates.dat").CompareFiles();
        FileComparison(results_dir + "cellcyclephases.dat", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/cellcyclephases.dat").CompareFiles();
        FileComparison(results_dir + "celltypes.dat", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/celltypes.dat").CompareFiles();

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 3u);
        TS_ASSERT_EQUALS(cell_mutation_states[1], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[2], 1u);
        TS_ASSERT_EQUALS(cell_mutation_states[3], 1u);

        // Test the GetCellProliferativeTypeCount() function
        std::vector<unsigned> cell_types = cell_population.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 4u);
        TS_ASSERT_EQUALS(cell_types[0], 2u);
        TS_ASSERT_EQUALS(cell_types[1], 2u);
        TS_ASSERT_EQUALS(cell_types[2], 2u);
        TS_ASSERT_EQUALS(cell_types[3], 0u);

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        FileComparison( results_dir + "results.parameters", "projects/ImmersedBoundary/test/data/TestWritersWithImmersedBoundaryCellPopulation/results.parameters").CompareFiles();

#ifdef CHASTE_VTK
        // Test that VTK writer has produced some files

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_1.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());

        // PVD file
        FileFinder vtk_file3(results_dir + "results.pvd", RelativeTo::Absolute);
        TS_ASSERT(vtk_file3.Exists());
 #endif //CHASTE_VTK
    }

    ///\todo Test archiving?
};
