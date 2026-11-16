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

#ifndef TESTCELLPOPULATIONEXTENTWRITER_HPP_
#define TESTCELLPOPULATIONEXTENTWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "CellPopulationExtentWriter.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "ChasteCuboid.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "OutputFileHandler.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemMesh.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemSphericalElementMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestCellPopulationExtentWriter : public AbstractCellBasedTestSuite
{
private:

    /**
     * Run the writer over a population and read back the numbers it wrote, one vector of values
     * per line. The first entry of each line is the time stamp, so the per-axis values follow it.
     *
     * @param rDirectory the output directory to write into
     * @param pCellPopulation the population to visit
     * @param numVisits how many time stamps to write
     * @return the parsed contents of the file, one vector per line
     */
    template <unsigned DIM, class POPULATION_TYPE>
    std::vector<std::vector<double> > WriteAndReadBack(const std::string& rDirectory,
                                                       POPULATION_TYPE* pCellPopulation,
                                                       unsigned numVisits = 1u)
    {
        OutputFileHandler handler(rDirectory, false);
        CellPopulationExtentWriter<DIM, DIM> writer;

        writer.OpenOutputFile(handler);
        for (unsigned visit = 0; visit < numVisits; ++visit)
        {
            writer.WriteTimeStamp();
            writer.Visit(pCellPopulation);
            writer.WriteNewline();
        }
        writer.CloseFile();

        std::vector<std::vector<double> > lines;
        std::ifstream file(handler.GetOutputDirectoryFullPath() + "cellpopulationextent.dat");
        std::string line;
        while (std::getline(file, line))
        {
            std::istringstream stream(line);
            std::vector<double> values;
            double value;
            while (stream >> value)
            {
                values.push_back(value);
            }
            if (!values.empty())
            {
                lines.push_back(values);
            }
        }
        return lines;
    }

    /**
     * @param numCells the number of cells to generate
     * @return that many cells, each with no cell cycle model
     */
    template <unsigned DIM>
    std::vector<CellPtr> CreateCells(unsigned numCells)
    {
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, DIM> cells_generator;
        cells_generator.GenerateBasicRandom(cells, numCells);
        return cells;
    }

public:

    /**
     * A 3x3x3 grid of unit scale puts nine nodes at each of three evenly spaced positions along
     * every axis. The displacements from the centroid are then -s, 0 and +s in equal numbers, so
     * the root-mean-square displacement is s*sqrt(2/3) exactly, on all three axes.
     */
    void TestWritesTheRootMeanSquareDisplacementOfAKnownGrid()
    {
        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<3>("TestCellPopulationExtentWriter", &population);

        TS_ASSERT_EQUALS(lines.size(), 1u);
        TS_ASSERT_EQUALS(lines[0].size(), 4u);

        // The time stamp comes first, then one value per axis
        TS_ASSERT_DELTA(lines[0][0], 0.0, 1e-9);

        const double expected = (1.0 / 3.0) * sqrt(2.0 / 3.0);
        TS_ASSERT_DELTA(lines[0][1], expected, 1e-6);
        TS_ASSERT_DELTA(lines[0][2], expected, 1e-6);
        TS_ASSERT_DELTA(lines[0][3], expected, 1e-6);
    }

    /**
     * The reason this writer measures a radius of gyration rather than a bounding box. A ball
     * carved from a close-packed lattice is isotropic, but the lattice reaches further along some
     * axes than others, so its bounding box is markedly anisotropic. Averaging over every node
     * rather than taking the two most extreme ones recovers the isotropy.
     */
    void TestIsInsensitiveToLatticeOrientationWhereABoundingBoxIsNot()
    {
        SemSphericalElementMeshGenerator<3> generator(200, 0.25);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        // The bounding box of this very same mesh is far from isotropic
        const ChasteCuboid<3> bounding_box = p_mesh->CalculateBoundingBox();
        double min_width = DBL_MAX;
        double max_width = 0.0;
        for (unsigned dim = 0; dim < 3; ++dim)
        {
            min_width = std::min(min_width, bounding_box.GetWidth(dim));
            max_width = std::max(max_width, bounding_box.GetWidth(dim));
        }
        TS_ASSERT_LESS_THAN(1.15, max_width / min_width);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<3>("TestCellPopulationExtentWriterIsotropy", &population);

        double min_value = DBL_MAX;
        double max_value = 0.0;
        for (unsigned dim = 0; dim < 3; ++dim)
        {
            min_value = std::min(min_value, lines[0][dim + 1]);
            max_value = std::max(max_value, lines[0][dim + 1]);
        }
        // Strictly more isotropic than the bounding box, and close to isotropic in absolute terms
        TS_ASSERT_LESS_THAN(max_value / min_value, max_width / min_width);
        TS_ASSERT_LESS_THAN(max_value / min_value, 1.10);

        // A uniformly filled ball of radius R has a per-axis radius of gyration of R/sqrt(5)
        for (unsigned dim = 0; dim < 3; ++dim)
        {
            TS_ASSERT_DELTA(lines[0][dim + 1], 0.25 / sqrt(5.0), 0.02);
        }
    }

    /**
     * A uniform deformation must scale the reported value in exactly the same proportion, since
     * that is what makes the strain measurable from it.
     */
    void TestScalesWithAUniformDeformation()
    {
        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        const std::vector<std::vector<double> > before
            = WriteAndReadBack<3>("TestCellPopulationExtentWriterBefore", &population);

        // Stretch by a half along z and squash by a fifth along x, leaving y alone
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            p_mesh->GetNode(i)->rGetModifiableLocation()[0] *= 0.8;
            p_mesh->GetNode(i)->rGetModifiableLocation()[2] *= 1.5;
        }

        const std::vector<std::vector<double> > after
            = WriteAndReadBack<3>("TestCellPopulationExtentWriterAfter", &population);

        // The file carries the default stream precision of about six significant figures
        TS_ASSERT_DELTA(after[0][1] / before[0][1], 0.8, 1e-5);
        TS_ASSERT_DELTA(after[0][2] / before[0][2], 1.0, 1e-5);
        TS_ASSERT_DELTA(after[0][3] / before[0][3], 1.5, 1e-5);
    }

    /**
     * Nothing about this writer is specific to subcellular element populations.
     */
    void TestWorksWithANodeBasedPopulation()
    {
        // Eight cells at the corners of a unit cube: the centroid is its centre and every node is
        // half a unit from it along each axis, so the RMS displacement is exactly a half
        std::vector<Node<3>*> nodes;
        unsigned index = 0u;
        for (unsigned k = 0; k < 2; ++k)
        {
            for (unsigned j = 0; j < 2; ++j)
            {
                for (unsigned i = 0; i < 2; ++i)
                {
                    nodes.push_back(new Node<3>(index++, false, static_cast<double>(i),
                                                static_cast<double>(j), static_cast<double>(k)));
                }
            }
        }

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        std::vector<CellPtr> cells = CreateCells<3>(mesh.GetNumNodes());
        NodeBasedCellPopulation<3> population(mesh, cells);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<3>("TestCellPopulationExtentWriterNodeBased", &population);

        TS_ASSERT_DELTA(lines[0][1], 0.5, 1e-9);
        TS_ASSERT_DELTA(lines[0][2], 0.5, 1e-9);
        TS_ASSERT_DELTA(lines[0][3], 0.5, 1e-9);

        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }

    /**
     * The writer is used once per output timestep over a whole simulation, so successive visits
     * must accumulate as separate lines in the one file.
     */
    void TestSuccessiveVisitsAppendLines()
    {
        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 1.0);
        auto p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells = CreateCells<3>(p_mesh->GetNumElements());
        SemBasedCellPopulation<3> population(*p_mesh, cells);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<3>("TestCellPopulationExtentWriterAppend", &population, 3u);

        TS_ASSERT_EQUALS(lines.size(), 3u);
        for (unsigned line = 0; line < 3u; ++line)
        {
            TS_ASSERT_EQUALS(lines[line].size(), 4u);
            TS_ASSERT_DELTA(lines[line][1], lines[0][1], 1e-9);
        }
    }

    /**
     * The writer offers a Visit() for every population type, and each defers to the same
     * measurement, so each overload needs exercising to show the dispatch reaches it.
     *
     * A CA population sits on a 4x4 lattice of unit spacing, so the displacements from the centroid
     * are +/-1.5 and +/-0.5, each four times over, on both axes. The mean square is therefore
     * (2*1.5^2 + 2*0.5^2)/4 = 1.25 and the value written is sqrt(1.25) exactly.
     */
    void TestWorksWithACaBasedPopulation()
    {
        PottsMeshGenerator<2> generator(4, 0, 0, 4, 0, 0);
        boost::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells = CreateCells<2>(4u);
        std::vector<unsigned> location_indices = { 7u, 11u, 12u, 13u };
        CaBasedCellPopulation<2> population(*p_mesh, cells, location_indices);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<2>("TestCellPopulationExtentWriterCaBased", &population);

        TS_ASSERT_EQUALS(lines.size(), 1u);
        TS_ASSERT_EQUALS(lines[0].size(), 3u);
        TS_ASSERT_DELTA(lines[0][1], sqrt(1.25), 1e-5);
        TS_ASSERT_DELTA(lines[0][2], sqrt(1.25), 1e-5);
    }

    /**
     * A Potts population over the same 4x4 lattice gives the same value: the measurement depends
     * only on where the nodes are, not on how they are grouped into elements.
     */
    void TestWorksWithAPottsBasedPopulation()
    {
        PottsMeshGenerator<2> generator(4, 1, 2, 4, 1, 2);
        boost::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells = CreateCells<2>(p_mesh->GetNumElements());
        PottsBasedCellPopulation<2> population(*p_mesh, cells);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<2>("TestCellPopulationExtentWriterPottsBased", &population);

        TS_ASSERT_EQUALS(lines[0].size(), 3u);
        TS_ASSERT_DELTA(lines[0][1], sqrt(1.25), 1e-5);
        TS_ASSERT_DELTA(lines[0][2], sqrt(1.25), 1e-5);
    }

    /**
     * Mesh-based populations are the one case where the writer's ELEMENT_DIM and SPACE_DIM can
     * differ, so this also pins the <2,2> instantiation of the mesh-based overload.
     */
    void TestWorksWithAMeshBasedPopulation()
    {
        HoneycombMeshGenerator generator(5, 5, 0);
        boost::shared_ptr<MutableMesh<2, 2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells = CreateCells<2>(p_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> population(*p_mesh, cells);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<2>("TestCellPopulationExtentWriterMeshBased", &population);

        TS_ASSERT_EQUALS(lines[0].size(), 3u);

        // A 5x5 honeycomb spans a few cell diameters, so both axes must report a positive extent
        TS_ASSERT_LESS_THAN(0.0, lines[0][1]);
        TS_ASSERT_LESS_THAN(0.0, lines[0][2]);

        // The generated patch is wider than it is tall, and the measure must reflect that
        TS_ASSERT_LESS_THAN(lines[0][2], lines[0][1]);
    }

    void TestWorksWithAVertexBasedPopulation()
    {
        HoneycombVertexMeshGenerator generator(4, 6);
        boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells = CreateCells<2>(p_mesh->GetNumElements());
        VertexBasedCellPopulation<2> population(*p_mesh, cells);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<2>("TestCellPopulationExtentWriterVertexBased", &population);

        TS_ASSERT_EQUALS(lines[0].size(), 3u);
        TS_ASSERT_LESS_THAN(0.0, lines[0][1]);
        TS_ASSERT_LESS_THAN(0.0, lines[0][2]);
    }

    void TestWorksWithAnImmersedBoundaryPopulation()
    {
        ImmersedBoundaryPalisadeMeshGenerator generator(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2, 2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells = CreateCells<2>(p_mesh->GetNumElements());
        ImmersedBoundaryCellPopulation<2> population(*p_mesh, cells);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<2>("TestCellPopulationExtentWriterImmersedBoundary", &population);

        TS_ASSERT_EQUALS(lines[0].size(), 3u);
        TS_ASSERT_LESS_THAN(0.0, lines[0][1]);
        TS_ASSERT_LESS_THAN(0.0, lines[0][2]);

        // An immersed boundary mesh lives in the unit square, so no axis can exceed it
        TS_ASSERT_LESS_THAN(lines[0][1], 1.0);
        TS_ASSERT_LESS_THAN(lines[0][2], 1.0);
    }

    /**
     * A population with no nodes has no centroid to measure from, so the writer reports zero on
     * every axis rather than dividing by a node count of zero.
     */
    void TestWritesZerosForAPopulationWithNoNodes()
    {
        // The single-argument constructor takes ownership of the mesh, so it must be heap allocated
        SemMesh<3>* p_mesh = new SemMesh<3>();
        SemBasedCellPopulation<3> population(*p_mesh);

        TS_ASSERT_EQUALS(population.GetNumNodes(), 0u);

        const std::vector<std::vector<double> > lines
            = WriteAndReadBack<3>("TestCellPopulationExtentWriterEmpty", &population);

        TS_ASSERT_EQUALS(lines.size(), 1u);
        TS_ASSERT_EQUALS(lines[0].size(), 4u);
        TS_ASSERT_DELTA(lines[0][1], 0.0, 1e-12);
        TS_ASSERT_DELTA(lines[0][2], 0.0, 1e-12);
        TS_ASSERT_DELTA(lines[0][3], 0.0, 1e-12);
    }

    void TestArchiving()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellPopulationExtentWriter.arch";

        {
            AbstractCellBasedWriter<3, 3>* const p_writer = new CellPopulationExtentWriter<3, 3>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_writer;

            delete p_writer;
        }

        {
            AbstractCellBasedWriter<3, 3>* p_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_writer_2;

            CellPopulationExtentWriter<3, 3>* p_loaded
                = dynamic_cast<CellPopulationExtentWriter<3, 3>*>(p_writer_2);
            TS_ASSERT(p_loaded != nullptr);

            delete p_writer_2;
        }
    }
};

#endif /*TESTCELLPOPULATIONEXTENTWRITER_HPP_*/
