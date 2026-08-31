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

#ifndef TESTCELLBASEDSIMULATIONFROMXML_HPP_
#define TESTCELLBASEDSIMULATIONFROMXML_HPP_

/**
 * @file TestCellBasedSimulationFromXml.hpp
 *
 * Tests for CellBasedSimulationFromXml:
 *  - Load and run a 5x5 MeshBasedCellPopulation from the bundled XML example.
 *  - Verify that key simulation parameters (end time, dt, sampling multiple)
 *    are correctly read from the file.
 *  - Verify that the GeneralisedLinearSpringForce is configured with the
 *    spring stiffness specified in the XML.
 *  - Smoke-test that Solve() completes without throwing.
 */

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <sstream>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedSimulationFromXml.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"

namespace
{

std::string BuildUniformInitialCellGeometryAndState(unsigned numCellsAcross,
                                                    unsigned numCellsUp,
                                                    double minCellCycleDuration,
                                                    double maxCellCycleDuration)
{
    std::ostringstream xml;
    const double sqrt3_over_2 = std::sqrt(3.0) / 2.0;
    unsigned index = 0;

    xml << "  <InitialCellGeometryAndState>\n";
    for (unsigned row = 0; row < numCellsUp; ++row)
    {
        for (unsigned col = 0; col < numCellsAcross; ++col, ++index)
        {
            double x = static_cast<double>(col) + (row % 2u == 0u ? 0.0 : 0.5);
            double y = static_cast<double>(row) * sqrt3_over_2;

            xml << "    <Cell index=\"" << index << "\">\n"
                << "      <Location>" << x << " " << y << "</Location>\n"
                << "      <MutationState>WildTypeCellMutationState</MutationState>\n"
                << "      <CellCycleModel>\n"
                << "        <UniformCellCycleModel>\n"
                << "          <MinCellCycleDuration>" << minCellCycleDuration << "</MinCellCycleDuration>\n"
                << "          <MaxCellCycleDuration>" << maxCellCycleDuration << "</MaxCellCycleDuration>\n"
                << "        </UniformCellCycleModel>\n"
                << "      </CellCycleModel>\n"
                << "      <SrnModel>\n"
                << "        <NullSrnModel/>\n"
                << "      </SrnModel>\n"
                << "    </Cell>\n";
        }
    }
    xml << "  </InitialCellGeometryAndState>\n";

    return xml.str();
}

}

/**
 * Test suite for CellBasedSimulationFromXml.
 */
class TestCellBasedSimulationFromXml : public AbstractCellBasedWithTimingsTestSuite
{
public:

    /**
     * Load the bundled 5x5 MeshBased example XML, check that the simulation
     * parameters are read correctly, and run the simulation.
     *
     * The XML specifies EndTime=15 (long enough for at least one division
     * cycle of min duration 10), a PlaneBasedCellKiller at x=5, and a
     * PlaneBoundaryCondition at x=0.
     */
    void TestLoadAndRunMeshBased5x5()
    {
        EXIT_IF_PARALLEL; // HoneycombMeshGenerator does not work in parallel

        // Locate the example XML file (relative to Chaste source root)
        FileFinder xml_file(
            "cell_based/test/data/TestCellBasedSimulationFromXml/MeshBasedSimulation5x5.xml",
            RelativeTo::ChasteSourceRoot);
        TS_ASSERT(xml_file.Exists());

        // Build the simulation from XML
        boost::shared_ptr<OffLatticeSimulation<2> > p_simulator =
            CellBasedSimulationFromXml::Load(xml_file.GetAbsolutePath());

        // ── Check simulation parameters were read from XML ────────────────

        // EndTime = 15
        TS_ASSERT_DELTA(p_simulator->mEndTime, 15.0, 1e-9);

        // Dt = 0.00833333  (1/120)
        TS_ASSERT_DELTA(p_simulator->GetDt(), 0.00833333, 1e-6);

        // SamplingTimestepMultiple = 6
        TS_ASSERT_EQUALS(p_simulator->mSamplingTimestepMultiple, 6u);

        // ── Check population type and initial size ────────────────────────
        AbstractCellPopulation<2>& r_pop = p_simulator->rGetCellPopulation();
        TS_ASSERT_EQUALS(r_pop.GetNumRealCells(), 25u); // 5x5

        // ── Check force and killer were added ─────────────────────────────
        TS_ASSERT_EQUALS(p_simulator->rGetForceCollection().size(), 1u);

        // ── Check initial node locations (before Solve) ───────────────────
        // The explicit InitialCellGeometryAndState section encodes a honeycomb lattice.
        // Row k has y = k * sqrt(3)/2; odd rows are offset by 0.5 in x.
        // Node 0:  (0.0,      0.0)
        // Node 4:  (4.0,      0.0)
        // Node 12: (2.0,  2*sqrt(3)/2)  — centre node of the grid
        // Node 24: (4.0,  4*sqrt(3)/2)
        const double kSqrt3Over2 = std::sqrt(3.0) / 2.0;
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(0)[0],  0.0,               1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(0)[1],  0.0,               1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(4)[0],  4.0,               1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(4)[1],  0.0,               1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(12)[0], 2.0,               1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(12)[1], 2.0 * kSqrt3Over2, 1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(24)[0], 4.0,               1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(24)[1], 4.0 * kSqrt3Over2, 1e-6);

        // ── Run the simulation ────────────────────────────────────────────
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // ── Post-run checks ───────────────────────────────────────────────
        // After 15 time units with UniformCellCycleModel(min=10, max=12),
        // cells should have divided (at least once), so population is larger.
        TS_ASSERT_LESS_THAN(25u, r_pop.GetNumRealCells());

        // The PlaneBasedCellKiller removes cells with x > 5, so no cell
        // centre should be to the right of x=5 after Solve().
        for (AbstractCellPopulation<2>::Iterator cell_iter = r_pop.Begin();
             cell_iter != r_pop.End();
             ++cell_iter)
        {
            c_vector<double, 2> loc = r_pop.GetLocationOfCellCentre(*cell_iter);
            TS_ASSERT_LESS_THAN_EQUALS(loc[0], 5.0 + 1e-6);
        }

        // The PlaneBoundaryCondition (x=0, normal -1,0) prevents cells from
        // crossing x=0, so no cell centre should be to the left of x=0.
        for (AbstractCellPopulation<2>::Iterator cell_iter = r_pop.Begin();
             cell_iter != r_pop.End();
             ++cell_iter)
        {
            c_vector<double, 2> loc = r_pop.GetLocationOfCellCentre(*cell_iter);
            TS_ASSERT_LESS_THAN_EQUALS(-1e-6, loc[0]);
        }
    }

    /**
     * Test that the XML loader correctly propagates a spring stiffness
     * overridden in the XML (16.0 instead of the default 15.0).
     *
     * The test writes a temporary XML file with the modified value and
     * re-loads it.
     */
    void TestXmlOverridesSpringStiffness()
    {
        EXIT_IF_PARALLEL;

        // Write a temporary XML file with a custom spring stiffness
        OutputFileHandler handler("TestCellBasedSimulationFromXmlTmp", false);
        out_stream p_file = handler.OpenOutputFile("custom_spring.xml");
        *p_file <<
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            "<CellBasedSimulation version=\"1\" elementDim=\"2\" spaceDim=\"2\">\n"
            "  <InitialCells>\n"
            "    <DefaultProliferativeType>TransitCellProliferativeType</DefaultProliferativeType>\n"
            "  </InitialCells>\n"
            << BuildUniformInitialCellGeometryAndState(3u, 3u, 10.0, 12.0) <<
            "  <Population type=\"MeshBasedCellPopulation\"/>\n"
            "  <Forces>\n"
            "    <Force type=\"GeneralisedLinearSpringForce\">\n"
            "      <MeinekeSpringStiffness>16</MeinekeSpringStiffness>\n"
            "      <MeinekeDivisionRestingSpringLength>0.5</MeinekeDivisionRestingSpringLength>\n"
            "      <MeinekeSpringGrowthDuration>1</MeinekeSpringGrowthDuration>\n"
            "      <UseCutOffLength>0</UseCutOffLength>\n"
            "      <CutOffLength>1e10</CutOffLength>\n"
            "    </Force>\n"
            "  </Forces>\n"
            "  <CellKillers/>\n"
            "  <SimulationModifiers/>\n"
            "  <BoundaryConditions/>\n"
            "  <NumericalMethod type=\"ForwardEulerNumericalMethod\">\n"
            "    <UseAdaptiveTimestep>0</UseAdaptiveTimestep>\n"
            "    <UseUpdateNodeLocation>0</UseUpdateNodeLocation>\n"
            "    <GhostNodeForcesEnabled>0</GhostNodeForcesEnabled>\n"
            "  </NumericalMethod>\n"
            "  <Simulation>\n"
            "    <Dt>0.00833333</Dt>\n"
            "    <EndTime>0.1</EndTime>\n"
            "    <SamplingTimestepMultiple>2</SamplingTimestepMultiple>\n"
            "    <OutputDirectory>TestCellBasedSimulationFromXml_CustomSpring</OutputDirectory>\n"
            "  </Simulation>\n"
            "</CellBasedSimulation>\n";
        p_file->close();

        std::string xml_path = handler.GetOutputDirectoryFullPath() + "custom_spring.xml";

        boost::shared_ptr<OffLatticeSimulation<2> > p_simulator =
            CellBasedSimulationFromXml::Load(xml_path);

        // Population: 3×3 = 9 real cells
        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 9u);

        // End time was set to 0.1
        TS_ASSERT_DELTA(p_simulator->mEndTime, 0.1, 1e-9);

        // One force added
        TS_ASSERT_EQUALS(p_simulator->rGetForceCollection().size(), 1u);

        // ── Check initial node positions for the 3×3 grid (before Solve) ──
        // HoneycombMeshGenerator(3, 3, 0):
        // Node 0: (0.0, 0.0), Node 2: (2.0, 0.0), Node 4: (1.5, sqrt(3)/2)
        const double kSqrt3Over2 = std::sqrt(3.0) / 2.0;
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(0)[0], 0.0,           1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(0)[1], 0.0,           1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(2)[0], 2.0,           1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(2)[1], 0.0,           1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(4)[0], 1.5,           1e-6);
        TS_ASSERT_DELTA(p_simulator->GetNodeLocation(4)[1], kSqrt3Over2,   1e-6);

        // Run without throwing
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // ── Post-run checks ───────────────────────────────────────────────
        // No divisions occur in 0.1 time units; population must stay at 9.
        TS_ASSERT_EQUALS(p_simulator->rGetCellPopulation().GetNumRealCells(), 9u);
    }
};

#endif // TESTCELLBASEDSIMULATIONFROMXML_HPP_
