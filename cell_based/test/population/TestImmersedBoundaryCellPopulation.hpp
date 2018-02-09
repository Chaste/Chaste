/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYCELLPOPULATION_HPP_
#define TESTIMMERSEDBOUNDARYCELLPOPULATION_HPP_

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
#include "FixedVertexBasedDivisionRule.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"

// Includes from Immersed Boundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryLinearInteractionForce.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

#include "Debug.hpp"

///\todo Vary the cell population geometry across tests
class TestImmersedBoundaryCellPopulation : public AbstractCellBasedTestSuite
{
public:

    void TestGetAndSetMethods()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
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

        TS_ASSERT_EQUALS(cell_population.GetReMeshFrequency(), UINT_MAX);
        cell_population.SetReMeshFrequency(5u);
        TS_ASSERT_EQUALS(cell_population.GetReMeshFrequency(), 5u);
    }

    void TestMeshMethods()
    {
        // Create an immersed boundary cell population object
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // Test that GetNumElements() and GetNumNodes() return the correct values
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 5u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 583u);

        // Test that rGetMesh() returns the mesh correctly
        ImmersedBoundaryMesh<2, 2>& r_mesh = cell_population.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 583u);
        TS_ASSERT_EQUALS(r_mesh.GetNumGridPtsX(), 128u);
        TS_ASSERT_EQUALS(r_mesh.GetNumGridPtsY(), 128u);
        TS_ASSERT_DELTA(r_mesh.GetCharacteristicNodeSpacing(), 0.0115, 1e-4);
        TS_ASSERT_DELTA(r_mesh.GetSpacingRatio(), 1.4805, 1e-4);

        // Test that GetElement() returns an element correctly
        ImmersedBoundaryElement<2, 2>* p_element_0 = cell_population.GetElement(0);
        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 100u);

        // Test that GetElement() returns the lamina correctly
        ImmersedBoundaryElement<1, 2>* p_lamina_0 = cell_population.GetLamina(0);
        TS_ASSERT_EQUALS(p_lamina_0->GetNumNodes(), 83u);

        // Test that GetElementCorrespondingToCell() returns an element correctly
        CellPtr p_cell_0 = *(cell_population.Begin());
        ImmersedBoundaryElement<2, 2>* p_element_0_again = cell_population.GetElementCorrespondingToCell(p_cell_0);
        TS_ASSERT_EQUALS(p_element_0_again, p_element_0);

        // Test that GetNode() returns a node correctly
        Node<2>* p_node_3 = cell_population.GetNode(3);
        TS_ASSERT_EQUALS(p_node_3->IsBoundaryNode(), true);
        TS_ASSERT_DELTA(p_node_3->rGetLocation()[0], 0.0421, 1e-4);
        TS_ASSERT_DELTA(p_node_3->rGetLocation()[1], 0.2800, 1e-4);

        // Test GetLocationOfCellCentre() returns the correct location
        c_vector<double, 2> cell_0_centre = cell_population.GetLocationOfCellCentre(p_cell_0);
        TS_ASSERT_DELTA(cell_0_centre[0], 0.1950, 1e-6);
        TS_ASSERT_DELTA(cell_0_centre[1], 0.505641, 1e-6);

        ///\todo Test GetNeighbouringNodeIndices() and GetNeighbouringLocationIndices() and implement these methods

        // Test GetWidth() returns the correct values
        TS_ASSERT_DELTA(cell_population.GetWidth(0), 0.9891, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetWidth(1), 0.4481, 1e-4);

        // Test GetVolumeOfCell() returns the correct value
        CellPtr p_cell_1 = *(++(cell_population.Begin()));
        double cell_1_volume = cell_population.GetVolumeOfCell(p_cell_1);
        TS_ASSERT_DELTA(cell_1_volume, 0.0774, 1e-4);

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

    ///\todo Check output files by eye too
    void TestWritersWithImmersedBoundaryCellPopulation()
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
            UniformCellCycleModel* p_model = new UniformCellCycleModel();

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

        std::string output_directory = "TestImmersedBoundaryPopulationWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteResultsToFiles(output_directory);

        SimulationTime::Instance()->IncrementTimeOneStep();
        cell_population.Update();
        cell_population.WriteResultsToFiles(output_directory);

        cell_population.CloseWritersFiles();

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizelements").CompareFiles();
//\todo: figure out what these tests are actually doing...
//        FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/cellages.dat").CompareFiles();
//        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizancestors").CompareFiles();
//        FileComparison(results_dir + "loggedcell.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/loggedcell.dat").CompareFiles();
//        FileComparison(results_dir + "results.vizlabels", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizlabels").CompareFiles();
//        FileComparison(results_dir + "results.vizlocationindices", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizlocationindices").CompareFiles();
//        FileComparison(results_dir + "results.vizmutationstates", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizmutationstates").CompareFiles();
//        FileComparison(results_dir + "results.vizcellphases", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizcellphases").CompareFiles();
//        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.vizcelltypes").CompareFiles();
//        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/cellareas.dat").CompareFiles();
//
//        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/cellmutationstates.dat").CompareFiles();
//        FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/cellcyclephases.dat").CompareFiles();
//        FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/celltypes.dat").CompareFiles();

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        FileComparison( results_dir + "results.parameters", "cell_based/test/data/TestImmersedBoundaryPopulationWriters/results.parameters").CompareFiles();

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

#endif /*TESTIMMERSEDBOUNDARYCELLPOPULATION_HPP_*/
