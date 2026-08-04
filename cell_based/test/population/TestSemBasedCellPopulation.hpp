/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTSEMBASEDCELLPOPULATION_HPP_
#define TESTSEMBASEDCELLPOPULATION_HPP_


#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "SemBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "CellDivisionLocationsWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "OutputFileHandler.hpp"

#include "AbstractCellBasedTestSuite.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 // Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif // CHASTE_VTK

class TestSemBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    static std::vector<Node<2>*> CreateNodes()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1u, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(2u, false, 0.05, 0.0));
        return nodes;
    }

    static std::vector<SemElement<2>*> CreateElements(std::vector<Node<2>*>& rNodes)
    {
        std::vector<SemElement<2>*> elements;

        std::vector<Node<2>*> element_0_nodes;
        element_0_nodes.push_back(rNodes[0]);
        element_0_nodes.push_back(rNodes[1]);
        elements.push_back(new SemElement<2>(0u, element_0_nodes));

        std::vector<Node<2>*> element_1_nodes;
        element_1_nodes.push_back(rNodes[1]);
        element_1_nodes.push_back(rNodes[2]);
        elements.push_back(new SemElement<2>(1u, element_1_nodes));

        return elements;
    }

    struct TwoElementSemMesh
    {
        std::vector<Node<2>*> nodes;
        std::vector<SemElement<2>*> elements;
        SemMesh<2> mesh;

        TwoElementSemMesh()
            : nodes(CreateNodes()),
              elements(CreateElements(nodes)),
              mesh(nodes, elements)
        {
        }
    };

    static std::vector<CellPtr> CreateCells(unsigned numCells)
    {
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, numCells);
        return cells;
    }

public:

    void TestConstructorRespectsLocationIndices()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        CellPtr p_cell_0 = cells[0];
        CellPtr p_cell_1 = cells[1];
        std::vector<unsigned> location_indices;
        location_indices.push_back(1u);
        location_indices.push_back(0u);

        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells, false, true, location_indices);

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell_0), 1u);
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell_1), 0u);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(1u), p_cell_0);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(0u), p_cell_1);
    }

    void TestConstructorRejectsInvalidLocationIndices()
    {
        {
            TwoElementSemMesh fixture;
            std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
            std::vector<unsigned> location_indices;
            location_indices.push_back(0u);
            location_indices.push_back(0u);

            TS_ASSERT_THROWS_THIS(SemBasedCellPopulation<2> cell_population(fixture.mesh, cells, false, true, location_indices),
                                  "A SemElement location index is assigned to more than one cell");
        }

        {
            TwoElementSemMesh fixture;
            std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
            std::vector<unsigned> location_indices;
            location_indices.push_back(0u);
            location_indices.push_back(2u);

            TS_ASSERT_THROWS_THIS(SemBasedCellPopulation<2> cell_population(fixture.mesh, cells, false, true, location_indices),
                                  "A supplied location index does not correspond to a SemElement");
        }
    }

    void TestConstructorRejectsWrongNumberOfCells()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(1u);

        TS_ASSERT_THROWS_THIS(SemBasedCellPopulation<2> cell_population(fixture.mesh, cells),
                              "There must be precisely one CellPtr for each SemElement");
    }

    void TestBasicAccessorsWidthAndDefaultTimeStep()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        TS_ASSERT_EQUALS(&(cell_population.rGetMesh()), &(fixture.mesh));
        TS_ASSERT_EQUALS(cell_population.GetElement(0u), fixture.mesh.GetElement(0u));
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(cell_population.GetNode(2u), fixture.mesh.GetNode(2u));
        TS_ASSERT_DELTA(cell_population.GetWidth(0u), 0.05, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetWidth(1u), 1.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetDefaultTimeStep(), 0.002, 1e-12);

        ChastePoint<2> new_location(0.25, 0.25);
        cell_population.SetNode(2u, new_location);
        TS_ASSERT_DELTA(cell_population.GetNode(2u)->rGetLocation()[0], 0.25, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(2u)->rGetLocation()[1], 0.25, 1e-12);
    }

    void TestUpdatePopulatesNeighbouringNodeAndLocationIndices()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        c_vector<double, 4> domain;
        domain[0] = -0.1;
        domain[1] =  0.2;
        domain[2] = -0.1;
        domain[3] =  1.1;
        fixture.mesh.SetUpBoxCollection(0.2, domain);

        cell_population.Update(false);

        std::set<unsigned> neighbouring_node_indices = cell_population.GetNeighbouringNodeIndices(0u);
        TS_ASSERT_EQUALS(neighbouring_node_indices.size(), 1u);
        TS_ASSERT_EQUALS(neighbouring_node_indices.count(2u), 1u);

        std::set<unsigned> neighbouring_location_indices = cell_population.GetNeighbouringLocationIndices(cell_population.GetCellUsingLocationIndex(0u));
        TS_ASSERT_EQUALS(neighbouring_location_indices.size(), 1u);
        TS_ASSERT_EQUALS(neighbouring_location_indices.count(1u), 1u);
    }

    /**
     * The neighbour searches examine each node pair from both ends, and which end a given node or
     * element lands on depends only on the order the pair happens to be stored in. The queries in
     * TestUpdatePopulatesNeighbouringNodeAndLocationIndices land on one end; these are the reverse
     * queries over the same single pair, which land on the other.
     */
    void TestNeighbourQueriesFromBothEndsOfANodePair()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        c_vector<double, 4> domain;
        domain[0] = -0.1;
        domain[1] =  0.2;
        domain[2] = -0.1;
        domain[3] =  1.1;
        fixture.mesh.SetUpBoxCollection(0.2, domain);

        cell_population.Update(false);

        std::set<unsigned> neighbours_of_node_2 = cell_population.GetNeighbouringNodeIndices(2u);
        TS_ASSERT_EQUALS(neighbours_of_node_2.size(), 1u);
        TS_ASSERT_EQUALS(neighbours_of_node_2.count(0u), 1u);

        std::set<unsigned> neighbours_of_element_1
            = cell_population.GetNeighbouringLocationIndices(cell_population.GetCellUsingLocationIndex(1u));
        TS_ASSERT_EQUALS(neighbours_of_element_1.size(), 1u);
        TS_ASSERT_EQUALS(neighbours_of_element_1.count(0u), 1u);
    }

    /**
     * AddNode and GetElementCorrespondingToCell are thin pass-throughs to the mesh, but both are
     * part of the population's public interface.
     */
    void TestAddNodeAndGetElementCorrespondingToCell()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 3u);

        // The mesh takes ownership of the node, so it must not be deleted here
        Node<2>* p_new_node = new Node<2>(3u, false, 0.5, 0.5);
        TS_ASSERT_EQUALS(cell_population.AddNode(p_new_node), 3u);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(cell_population.GetNode(3u), p_new_node);

        CellPtr p_cell_1 = cell_population.GetCellUsingLocationIndex(1u);
        TS_ASSERT_EQUALS(cell_population.GetElementCorrespondingToCell(p_cell_1), fixture.mesh.GetElement(1u));
    }

    /**
     * The growing-domain PDE modifiers build their finite element mesh from the population, which
     * a SEM population cannot supply. Both hooks must say so rather than returning a sentinel: the
     * caller dereferences the returned mesh immediately without a null check, so returning nullptr
     * segfaults instead of reporting the problem.
     *
     * The box-domain PDE modifiers generate their own mesh and call neither of these.
     */
    void TestGrowingDomainPdeModifierHooksAreRejected()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        TS_ASSERT_THROWS_CONTAINS(cell_population.GetTetrahedralMeshForPdeModifier(),
                                  "Currently can't solve PDEs on a SemMesh");

        std::string item = "unused";
        TS_ASSERT_THROWS_CONTAINS(cell_population.GetCellDataItemAtPdeNode(0u, item, false, 0.0),
                                  "Currently can't solve PDEs on a SemMesh");
    }

    /**
     * The population accepts count and event writers, both of which are deliberate no-ops for SEM
     * because element division and removal events are not tracked.
     */
    void TestAcceptsPopulationCountAndEventWriters()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        boost::shared_ptr<AbstractCellPopulationCountWriter<2, 2> > p_count_writer(
            new CellMutationStatesCountWriter<2, 2>());
        TS_ASSERT_THROWS_NOTHING(cell_population.AcceptPopulationCountWriter(p_count_writer));

        boost::shared_ptr<AbstractCellPopulationEventWriter<2, 2> > p_event_writer(
            new CellDivisionLocationsWriter<2, 2>());
        TS_ASSERT_THROWS_NOTHING(cell_population.AcceptPopulationEventWriter(p_event_writer));
    }

    /**
     * Once every element containing a node has been deleted there is no cell left to take a damping
     * constant from, so the population refuses rather than returning a zero that would divide
     * through into an infinite velocity.
     */
    void TestDampingConstantThrowsWhenNoLiveElementContainsTheNode()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        // Removing both cells deletes both elements and unregisters every node from them
        cell_population.GetCellUsingLocationIndex(0u)->Kill();
        cell_population.GetCellUsingLocationIndex(1u)->Kill();
        TS_ASSERT_EQUALS(cell_population.RemoveDeadCells(), 2u);

        TS_ASSERT_THROWS_CONTAINS(cell_population.GetDampingConstant(0u),
            "Node 0 is not contained in any live SEM elements");
    }

    /**
     * Validate() tolerates an element that was already marked as deleted before the population was
     * built, which is what allows a population to be reconstructed around a mesh carrying dead
     * cells. The cell mapped to that element is simply not a real cell: the population iterator
     * skips it, and it is not counted among the real cells.
     *
     * This exercises the only branch of Validate() that is reachable at all. Every throwing check
     * in it is pre-empted by the constructor's own validation or by that same iterator skip.
     */
    void TestValidateToleratesAnAlreadyDeletedElement()
    {
        TwoElementSemMesh fixture;
        fixture.mesh.GetElement(1u)->MarkAsDeleted();

        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());

        // The constructor takes the cells by reference and empties the vector, so it can only be
        // called once. Validate() runs inside it, and failing it would throw out of this line.
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        // The deleted element still occupies its index, but its cell is not a real one
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);
        TS_ASSERT(cell_population.GetElement(1u)->IsDeleted());
        TS_ASSERT(cell_population.IsCellAssociatedWithADeletedLocation(cell_population.GetCellUsingLocationIndex(1u)));
    }

    void TestDampingConstantUsesContainingElementCells()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);
        cell_population.SetDampingConstantNormal(2.0);
        cell_population.SetDampingConstantMutant(8.0);

        boost::shared_ptr<AbstractCellProperty> p_apc1(cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        cell_population.GetCellUsingLocationIndex(1u)->SetMutationState(p_apc1);

        TS_ASSERT_DELTA(cell_population.GetDampingConstant(0u), 2.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetDampingConstant(2u), 8.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetDampingConstant(1u), 5.0, 1e-12);
    }

    void TestRemoveDeadCellsDeletesElementsAndMappings()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        CellPtr p_dead_cell = cell_population.GetCellUsingLocationIndex(1u);
        p_dead_cell->Kill();

        TS_ASSERT(!cell_population.GetElement(1u)->IsDeleted());
        TS_ASSERT(!cell_population.IsCellAssociatedWithADeletedLocation(p_dead_cell));

        TS_ASSERT_EQUALS(cell_population.RemoveDeadCells(), 1u);
        TS_ASSERT(cell_population.GetElement(1u)->IsDeleted());
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);
        TS_ASSERT_EQUALS(fixture.mesh.GetNode(2u)->rGetContainingElementIndices().count(1u), 0u);
    }

    void TestAddCellIsExplicitlyUnsupported()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);
        std::vector<CellPtr> new_cells = CreateCells(1u);

        TS_ASSERT_THROWS_THIS(cell_population.AddCell(new_cells[0], cell_population.GetCellUsingLocationIndex(0u)),
                              "SemBasedCellPopulation does not support AddCell() because SEM element division is not implemented");
    }

    void TestCheckForStepSizeException()
    {
        TwoElementSemMesh fixture;
        std::vector<CellPtr> cells = CreateCells(fixture.mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(fixture.mesh, cells);

        fixture.mesh.SetMaximumInteractionDistance(0.5);

        // A displacement smaller than the interaction distance must not throw
        c_vector<double, 2> small_displacement;
        small_displacement[0] = 0.1;
        small_displacement[1] = 0.1;
        TS_ASSERT_THROWS_NOTHING(cell_population.CheckForStepSizeException(0u, small_displacement, 0.01));

        // A displacement larger than the interaction distance must throw a StepSizeException,
        // signalling to the numerical method that the timestep is too large
        c_vector<double, 2> large_displacement;
        large_displacement[0] = 0.6;
        large_displacement[1] = 0.6;
        TS_ASSERT_THROWS_ANYTHING(cell_population.CheckForStepSizeException(0u, large_displacement, 0.01));
    }

    void TestWriteVtkResultsWithPerElementData()
    {
#ifdef CHASTE_VTK
        // Two disjoint triangular elements, offset along x
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0u, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1u, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2u, false, 0.5, 1.0));
        nodes.push_back(new Node<2>(3u, false, 5.0, 0.0));
        nodes.push_back(new Node<2>(4u, false, 6.0, 0.0));
        nodes.push_back(new Node<2>(5u, false, 5.5, 1.0));

        std::vector<Node<2>*> element_0_nodes(nodes.begin(), nodes.begin() + 3);
        std::vector<Node<2>*> element_1_nodes(nodes.begin() + 3, nodes.end());
        std::vector<SemElement<2>*> elements;
        elements.push_back(new SemElement<2>(0u, element_0_nodes));
        elements.push_back(new SemElement<2>(1u, element_1_nodes));
        SemMesh<2> mesh(nodes, elements);

        // One VTK cell (a point-cloud cell) per element, so cell data maps one-to-one to elements
        mesh.SetOutputElementSurfacesToVtk(false);

        std::vector<CellPtr> cells = CreateCells(mesh.GetNumElements());
        SemBasedCellPopulation<2> cell_population(mesh, cells);

        // A per-cell CellData item, and the cell volume (via CellVolumesWriter)
        for (unsigned elem_index = 0u; elem_index < mesh.GetNumElements(); ++elem_index)
        {
            cell_population.GetCellUsingLocationIndex(elem_index)->GetCellData()->SetItem("my data", 100.0 + elem_index);
        }
        cell_population.AddCellWriter<CellVolumesWriter>();

        // The VTK writer stamps files with the elapsed step count, so set up the time stepper
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        std::string output_directory = "TestSemBasedCellPopulationVtk";
        OutputFileHandler output_file_handler(output_directory, false);
        cell_population.OpenWritersFiles(output_file_handler);
        cell_population.WriteVtkResultsToFile(output_directory);
        cell_population.CloseWritersFiles();

        // Read the VTK output back and check the per-element cell-data arrays
        std::string results_file = OutputFileHandler::GetChasteTestOutputDirectory()
                                   + output_directory + "/results_0.vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridReader> p_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        p_reader->SetFileName(results_file.c_str());
        p_reader->Update();
        vtkUnstructuredGrid* p_grid = p_reader->GetOutput();

        TS_ASSERT_EQUALS(p_grid->GetNumberOfCells(), 2);

        vtkDataArray* p_element_index = p_grid->GetCellData()->GetArray("SemElementIndex");
        vtkDataArray* p_my_data = p_grid->GetCellData()->GetArray("my data");
        vtkDataArray* p_volumes = p_grid->GetCellData()->GetArray("Cell volumes");
        TS_ASSERT(p_my_data != nullptr);
        TS_ASSERT(p_volumes != nullptr);

        for (unsigned cell = 0u; cell < 2u; ++cell)
        {
            unsigned elem_index = static_cast<unsigned>(p_element_index->GetTuple1(cell));
            TS_ASSERT_DELTA(p_my_data->GetTuple1(cell), 100.0 + elem_index, 1e-12);

            // The written volume must equal GetVolumeOfCell() for the corresponding cell
            double expected_volume = cell_population.GetVolumeOfCell(cell_population.GetCellUsingLocationIndex(elem_index));
            TS_ASSERT_LESS_THAN(0.0, expected_volume);
            TS_ASSERT_DELTA(p_volumes->GetTuple1(cell), expected_volume, 1e-6);
        }
#endif // CHASTE_VTK
    }
};

#endif /*TESTSEMBASEDCELLPOPULATION_HPP_*/
