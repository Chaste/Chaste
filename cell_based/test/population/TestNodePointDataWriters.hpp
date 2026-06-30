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

#ifndef TESTNOTEPOINTDATAWRITERS_HPP_
#define TESTNOTEPOINTDATAWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "ArchiveOpener.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SemMultiElementMeshGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "SemMesh.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemSingleElementMeshGenerator.hpp"

// Node point data writers
#include "BoundaryNodePointDataWriter.hpp"
#include "ElementIdNodePointDataWriter.hpp"
#include "NodeRegionPointDataWriter.hpp"

// Boost
#include <boost/shared_ptr.hpp>

#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"


class TestCellWriters : public AbstractCellBasedTestSuite
{
public:

    void TestBoundaryNodePointDataWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a regular vertex mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        const boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        const boost::shared_ptr<TransitCellProliferativeType> p_transit_type;
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and test the point data writer object
        const BoundaryNodePointDataWriter<2, 2> point_data_writer;
        TS_ASSERT_EQUALS(point_data_writer.rGetFieldName(), "Boundary Node");

        const std::vector<double> boundary_node_info = point_data_writer.GetPointData(&cell_population);
        TS_ASSERT_DELTA(boundary_node_info[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[6], 0.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[7], 1.0, 1e-12);
    }

    void TestBoundaryNodePointDataWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BoundaryNodePointDataWriter.arch";
        // serialise
        {
            std::shared_ptr<AbstractNodePointDataWriter<2,2>> p_cell_writer = std::make_shared<BoundaryNodePointDataWriter<2,2>>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_cell_writer;
        }
        // deserialize
        {
            std::shared_ptr<AbstractNodePointDataWriter<2,2>> p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;
        }
    }

    void TestElementIdNodePointDataWriter()
    {
        EXIT_IF_PARALLEL;

        SemMultiElementMeshGenerator<2> generator({ 2, 2 }, {1, 2}, 0.5);
        auto p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create and test the point data writer object
        const ElementIdNodePointDataWriter<2, 2> point_data_writer;
        TS_ASSERT_EQUALS(point_data_writer.rGetFieldName(), "Element ID");

        const std::vector<double> boundary_node_info = point_data_writer.GetPointData(&cell_population);
        TS_ASSERT_DELTA(boundary_node_info[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[1], 0.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[2], 0.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[3], 0.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[4], 1.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[5], 1.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[6], 1.0, 1e-12);
        TS_ASSERT_DELTA(boundary_node_info[7], 1.0, 1e-12);
    }

    void TestElementIdNodePointDataWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ElementIdNodePointDataWriter.arch";
        // serialise
        {
            std::shared_ptr<AbstractNodePointDataWriter<2,2>> p_cell_writer = std::make_shared<ElementIdNodePointDataWriter<2,2>>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_cell_writer;
        }
        // deserialize
        {
            std::shared_ptr<AbstractNodePointDataWriter<2,2>> p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;
        }
    }

    void TestNodeRegionPointDataWriter()
    {
        EXIT_IF_PARALLEL;

        SemSingleElementMeshGenerator<3> generator({ 3, 3, 3 }, 0.5);
        auto p_mesh = generator.GetMesh();

        // Set all but the very corner nodes to region 1, others to 3
        {
            c_vector<double, 3> centroid = p_mesh->GetCentroidOfElement(0u);
            auto p_elem_0 = p_mesh->GetElement(0);
            const double threshold = sqrt(3) / 6.0 - 1e-6;  // A tiny bit less than dist from corner to centroid
            for (unsigned i = 0; i < p_elem_0->GetNumNodes(); ++i)
            {
                auto p_node = p_elem_0->GetNode(i);
                if (norm_2(p_mesh->GetVectorFromAtoB(centroid, p_node->rGetLocation())) < threshold)
                {
                    p_node->SetRegion(1u);
                }
                else
                {
                    p_node->SetRegion(3u);
                }
            }
        }

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Create and test the point data writer object
        const NodeRegionPointDataWriter<3, 3> point_data_writer;
        TS_ASSERT_EQUALS(point_data_writer.rGetFieldName(), "Node Region");

        const std::vector<double> node_regions = point_data_writer.GetPointData(&cell_population);
        const std::vector<double> correct = {
            3.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 1.0, 3.0,
            1.0, 1.0, 1.0, 3.0, 1.0, 3.0
        };
        for (unsigned i = 0; i < node_regions.size(); ++i)
        {
            TS_ASSERT_DELTA(node_regions[i], correct[i], 1e-12);
        }
    }

    void TestNodeRegionPointDataWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "NodeRegionPointDataWriter.arch";
        // serialise
        {
            std::shared_ptr<AbstractNodePointDataWriter<2,2>> p_cell_writer = std::make_shared<NodeRegionPointDataWriter<2,2>>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_cell_writer;
        }
        // deserialize
        {
            std::shared_ptr<AbstractNodePointDataWriter<2,2>> p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;
        }
    }
};

#endif /*TESTNOTEPOINTDATAWRITERS_HPP_*/
