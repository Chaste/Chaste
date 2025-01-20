/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYFORCES_HPP_
#define TESTIMMERSEDBOUNDARYFORCES_HPP_

// Needed for test framework
#include "AbstractCellBasedTestSuite.hpp"

#include <memory>

#include "CellsGenerator.hpp"
#include "CellLabel.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "ImmersedBoundaryElement.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "NoCellCycleModel.hpp"
#include "Node.hpp"
#include "SmartPointers.hpp"
#include "UblasCustomFunctions.hpp"

#include <boost/pointer_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

// Immersed boundary forces tested in this test suite
#include "ImmersedBoundaryLinearInteractionForce.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryMorseInteractionForce.hpp"
#include "ImmersedBoundaryMorseMembraneForce.hpp"
#include "ImmersedBoundaryKinematicFeedbackForce.hpp"
#include "ImmersedBoundaryLinearDifferentialAdhesionForce.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryForces : public AbstractCellBasedTestSuite
{
public:

    void TestImmersedBoundaryLinearInteractionForceMethods()
    {
        // No additive noise
        {
            // Create a small 1x1 mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(1, 1, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.2, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryLinearInteractionForce<2> force;
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }

        // With additive noise
        {
            // Create a small mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(3, 3, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            std::vector<std::pair<Node<2>*, Node<2>*>> pairs;
            pairs.push_back(std::pair<Node<2>*, Node<2>*>(p_mesh->GetNode(0), p_mesh->GetNode(1)));
            p_mesh->GetNode(0)->ClearAppliedForce();
            p_mesh->GetNode(1)->ClearAppliedForce();

            ImmersedBoundaryLinearInteractionForce<2> force;
            force.SetAdditiveNormalNoise(true);
            for (auto iter = p_mesh->GetNodeIteratorBegin(); iter != p_mesh->GetNodeIteratorEnd(); ++iter)
            {
                iter->ClearAppliedForce();
            }
            force.AddImmersedBoundaryForceContribution(pairs, cell_population);
        }

        // With noise & laminas
        {
            // Create a small mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(2, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);
            auto p_mesh = &mesh;
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            for (auto node_iter = p_mesh->GetNodeIteratorBegin();
                 node_iter != p_mesh->GetNodeIteratorEnd();
                 ++node_iter)
            {
                node_iter->ClearAppliedForce();
            }
            std::vector<std::pair<Node<2>*, Node<2>*>> pairs;
            pairs.push_back(std::pair<Node<2>*, Node<2>*>(p_mesh->GetNode(0), p_mesh->GetNode(1)));

            ImmersedBoundaryLinearInteractionForce<2> force;
            force.SetAdditiveNormalNoise(true);
            force.AddImmersedBoundaryForceContribution(pairs, cell_population);
        }
    }

    void TestArchivingOfImmersedBoundaryLinearInteractionForce()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundaryLinearInteractionForce.arch";

        {
            ImmersedBoundaryLinearInteractionForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetSpringConst(1.23);
            force.SetRestLength(2.34);
            force.SetLaminaSpringConstMult(3.45);
            force.SetLaminaRestLengthMult(4.56);
            force.SetAdditiveNormalNoise(true);
            force.SetNormalNoiseMean(5.67);
            force.SetNormalNoiseStdDev(6.78);

            // Serialize via pointer to most abstract class possible
            AbstractImmersedBoundaryForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractImmersedBoundaryForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;
            ImmersedBoundaryLinearInteractionForce<2>* p_derived_force = static_cast<ImmersedBoundaryLinearInteractionForce<2>*>(p_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_derived_force->GetSpringConst(), 1.23, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetRestLength(), 2.34, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaSpringConstMult(), 3.45, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaRestLengthMult(), 4.56, 1e-6);
            TS_ASSERT(p_derived_force->GetAdditiveNormalNoise());
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseMean(), 5.67, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseStdDev(), 6.78, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestImmersedBoundaryMorseInteractionForceMethods()
    {
        {
            // Create a small 3x3 mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(3, 3, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.102, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryMorseInteractionForce<2> force;
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }

        { // Test with laminas
            // Create a small mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            for (auto p_node : nodes)
            {
                p_node->ClearAppliedForce();
            }

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(3, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(4, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(5, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);
            auto p_mesh = &mesh;

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.7, 0.7);
            Node<2> node1(0, true, 0.702, 0.7);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);
            node1.SetRegion(LAMINA_REGION);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryMorseInteractionForce<2> force;
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }

        { // Coverage for add normal noise
            // Create a small 3x3 mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(3, 3, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            auto nodes = p_mesh->rGetNodes();
            for (auto& p_node : nodes)
            {
                p_node->ClearAppliedForce();
            }

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.102, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryMorseInteractionForce<2> force;
            force.SetAdditiveNormalNoise(true);
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }
    }

    void TestArchivingOfImmersedBoundaryMorseInteractionForce()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundaryLinearInteractionForce.arch";

        {
            ImmersedBoundaryMorseInteractionForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetWellDepth(1.23);
            force.SetRestLength(2.34);
            force.SetLaminaWellDepthMult(3.45);
            force.SetLaminaRestLengthMult(4.56);
            force.SetWellWidth(5.67);
            force.SetAdditiveNormalNoise(true);
            force.SetNormalNoiseMean(6.78);
            force.SetNormalNoiseStdDev(7.89);

            // Serialize via pointer to most abstract class possible
            AbstractImmersedBoundaryForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractImmersedBoundaryForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;
            ImmersedBoundaryMorseInteractionForce<2>* p_derived_force = static_cast<ImmersedBoundaryMorseInteractionForce<2>*>(p_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_derived_force->GetWellDepth(), 1.23, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetRestLength(), 2.34, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaWellDepthMult(), 3.45, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaRestLengthMult(), 4.56, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetWellWidth(), 5.67, 1e-6);
            TS_ASSERT(p_derived_force->GetAdditiveNormalNoise());
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseMean(), 6.78, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseStdDev(), 7.89, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestImmersedBoundaryLinearMembraneForce()
    {
        {
            // Create a small 3x3 mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(3, 3, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.102, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryLinearMembraneForce<2> force;
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }
        { // Coverage with noise
            // Create a small 3x3 mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(3, 3, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.102, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryLinearMembraneForce<2> force;
            force.SetAdditiveNormalNoise(true);
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }
    }

    void TestArchivingOfImmersedBoundaryLinearMembraneForce()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundaryLinearMembraneForce.arch";

        {
            ImmersedBoundaryLinearMembraneForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetElementSpringConst(1.23);
            force.SetElementRestLength(2.34);
            force.SetLaminaSpringConst(3.45);
            force.SetLaminaRestLength(4.56);
            force.SetAdditiveNormalNoise(true);
            force.SetNormalNoiseMean(5.67);
            force.SetNormalNoiseStdDev(6.78);

            // Serialize via pointer to most abstract class possible
            AbstractImmersedBoundaryForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractImmersedBoundaryForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;
            ImmersedBoundaryLinearMembraneForce<2>* p_derived_force = static_cast<ImmersedBoundaryLinearMembraneForce<2>*>(p_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_derived_force->GetElementSpringConst(), 1.23, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetElementRestLength(), 2.34, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaSpringConst(), 3.45, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaRestLength(), 4.56, 1e-6);
            TS_ASSERT(p_derived_force->GetAdditiveNormalNoise());
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseMean(), 5.67, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseStdDev(), 6.78, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestImmersedBoundaryMorseMembraneForce()
    {
        {
            // Create a small 3x3 mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(3, 3, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.102, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryMorseMembraneForce<2> force;
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }

        { // With noise
            // Create a small 3x3 mesh
            ImmersedBoundaryHoneycombMeshGenerator gen(3, 3, 3, 0.1, 0.3);
            ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
            p_mesh->SetNumGridPtsXAndY(32);

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.102, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryMorseMembraneForce<2> force;
            force.SetAdditiveNormalNoise(true);
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }

        { // With laminas

            // Create a small mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));
            elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes));

            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(2, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);
            auto p_mesh = &mesh;

            // Create a minimal cell population
            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.SetInteractionDistance(0.01);

            // Create two nodes and put them in a vector of pairs
            Node<2> node0(0, true, 0.1, 0.1);
            Node<2> node1(0, true, 0.102, 0.1);
            std::vector<std::pair<Node<2>*, Node<2>*> > node_pair;
            node_pair.push_back(std::pair<Node<2>*, Node<2>*>(&node0, &node1));

            // Put the nodes in different elements so force calculation is triggered
            node0.AddElement(0);
            node1.AddElement(1);

            node0.ClearAppliedForce();
            node1.ClearAppliedForce();

            ImmersedBoundaryMorseMembraneForce<2> force;
            force.AddImmersedBoundaryForceContribution(node_pair, cell_population);
        }
    }

    void TestArchivingOfImmersedBoundaryMorseMembraneForce()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundaryMorseMembraneForce.arch";

        {
            ImmersedBoundaryMorseMembraneForce<2> force;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            force.SetElementWellDepth(1.23);
            force.SetElementRestLength(2.34);
            force.SetLaminaWellDepth(3.45);
            force.SetLaminaRestLength(4.56);
            force.SetWellWidth(5.67);
            force.SetAdditiveNormalNoise(true);
            force.SetNormalNoiseMean(6.78);
            force.SetNormalNoiseStdDev(7.89);

            // Serialize via pointer to most abstract class possible
            AbstractImmersedBoundaryForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractImmersedBoundaryForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;
            ImmersedBoundaryMorseMembraneForce<2>* p_derived_force = static_cast<ImmersedBoundaryMorseMembraneForce<2>*>(p_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_derived_force->GetElementWellDepth(), 1.23, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetElementRestLength(), 2.34, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaWellDepth(), 3.45, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaRestLength(), 4.56, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetLaminaRestLength(), 4.56, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetWellWidth(), 5.67, 1e-6);
            TS_ASSERT(p_derived_force->GetAdditiveNormalNoise());
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseMean(), 6.78, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseStdDev(), 7.89, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestImmersedBoundaryKinematicFeedbackForce()
    {
        // Test member variables
        {
            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            // Set member variables
            p_force->SetSpringConst(1.23);
            TS_ASSERT_DELTA(p_force->GetSpringConst(), 1.23, 1e-6);
        }

        // Test CalculateRelativeVelocityComponent() helper method
        {
            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            // Two nodes previous positions (0,0) and (1,0), new positions (0,0) and (1,1).
            // Unit perp should be (0,1), and velocity in that direction is just 1/dt
            c_vector<double, 2> unit_perp;
            c_vector<double, 2> previous_disp = Create_c_vector(1.0, 0.0);
            c_vector<double, 2> current_disp = Create_c_vector(1.0, 1.0);
            double vel_comp = p_force->CalculateRelativeVelocityComponent(previous_disp, current_disp, unit_perp);

            TS_ASSERT_DELTA(vel_comp, 1.0 / dt, 1e-6);
            TS_ASSERT_DELTA(unit_perp[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(unit_perp[1], 1.0, 1e-6);

            //Two nodes previous positions (2,2) and (0,0), new positions (2,0) and (0,0).
            // Unit perp should be root2*(0.5,-0.5), and velocity in that direction is just -root(2)/dt
            previous_disp = Create_c_vector(-2.0, -2.0);
            current_disp = Create_c_vector(-2.0, 0.0);
            vel_comp = p_force->CalculateRelativeVelocityComponent(previous_disp, current_disp, unit_perp);

            TS_ASSERT_DELTA(vel_comp, -sqrt(2.0) / dt, 1e-6);
            TS_ASSERT_DELTA(unit_perp[0], 0.5 * sqrt(2.0), 1e-6);
            TS_ASSERT_DELTA(unit_perp[1], -0.5 * sqrt(2.0), 1e-6);

            // Same as previous example, but new positions (2,4) and (0,0).
            // Exactly the same, but velocity component is minus what it was before.
            previous_disp = Create_c_vector(-2.0, -2.0);
            current_disp = Create_c_vector(-2.0, -4.0);
            vel_comp = p_force->CalculateRelativeVelocityComponent(previous_disp, current_disp, unit_perp);

            TS_ASSERT_DELTA(vel_comp, sqrt(2.0) / dt, 1e-6);
            TS_ASSERT_DELTA(unit_perp[0], 0.5 * sqrt(2.0), 1e-6);
            TS_ASSERT_DELTA(unit_perp[1], -0.5 * sqrt(2.0), 1e-6);
        }

        // Test correct force is added
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            double force_factor = mesh.GetAverageNodeSpacingOfElement(0u, false) / population.GetIntrinsicSpacing();
            p_force->SetSpringConst(1.0 / force_factor);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            p_force->mPreviousLocations = {Create_c_vector(0.5, 0.5),
                                           Create_c_vector(0.6, 0.6),
                                           Create_c_vector(0.6, 0.5), //unmoved
                                           Create_c_vector(0.6, 0.5),
                                           Create_c_vector(0.7, 0.7),
                                           Create_c_vector(0.7, 0.7)}; //unmoved

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};

            for (auto&& p_node : mesh.rGetNodes())
            {
                p_node->ClearAppliedForce();
            }

            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);

            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[1], -1.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[0], -0.5, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[1], 0.5, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[0], 0.0, 1e-6); //not involved
            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[1], 0.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[0], 0.0, 1e-6); //opposite to node 0
            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[1], 1.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[0], 0.5, 1e-6); //opposite to node 1
            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[1], -0.5, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[0], 0.0, 1e-6); //not involved
            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[1], 0.0, 1e-6);

            // If we do the exact same thing again, but with the pairs swapped, the forces should be identical
            for (auto&& p_node : mesh.rGetNodes())
            {
                p_node->ClearAppliedForce();
            }

            // Need to reset the previous locations, as they will have been overwritten in the previous calculation
            p_force->mPreviousLocations = {Create_c_vector(0.5, 0.5),
                                           Create_c_vector(0.6, 0.6),
                                           Create_c_vector(0.6, 0.5), //unmoved
                                           Create_c_vector(0.6, 0.5),
                                           Create_c_vector(0.7, 0.7),
                                           Create_c_vector(0.7, 0.7)}; //unmoved

            std::vector<std::pair<Node<2>*, Node<2>*>> swapped_node_pairs = {std::make_pair(nodes[3], nodes[0]),
                                                                             std::make_pair(nodes[4], nodes[1]),
                                                                             std::make_pair(nodes[5], nodes[2])};

            p_force->AddImmersedBoundaryForceContribution(swapped_node_pairs, population);

            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[1], -1.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[0], -0.5, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[1], 0.5, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[0], 0.0, 1e-6); //not involved
            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[1], 0.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[0], 0.0, 1e-6); //opposite to node 0
            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[1], 1.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[0], 0.5, 1e-6); //opposite to node 1
            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[1], -0.5, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[0], 0.0, 1e-6); //not involved
            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[1], 0.0, 1e-6);
        }

        // Coverage for lamina regions
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.back()->SetRegion(LAMINA_REGION);
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            double force_factor = mesh.GetAverageNodeSpacingOfElement(0u, false) / population.GetIntrinsicSpacing();
            p_force->SetSpringConst(1.0 / force_factor);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            p_force->mPreviousLocations = {Create_c_vector(0.5, 0.5),
                                           Create_c_vector(0.6, 0.6),
                                           Create_c_vector(0.6, 0.5), //unmoved
                                           Create_c_vector(0.6, 0.5),
                                           Create_c_vector(0.7, 0.7),
                                           Create_c_vector(0.7, 0.7)}; //unmoved

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};
            for (auto&& p_node : mesh.rGetNodes())
            {
                p_node->ClearAppliedForce();
            }

            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);
        }

        // Coverage for nodes in same element
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            double force_factor = mesh.GetAverageNodeSpacingOfElement(0u, false) / population.GetIntrinsicSpacing();
            p_force->SetSpringConst(1.0 / force_factor);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            p_force->mPreviousLocations = {Create_c_vector(0.5, 0.5),
                                           Create_c_vector(0.6, 0.6),
                                           Create_c_vector(0.6, 0.5), //unmoved
                                           Create_c_vector(0.6, 0.5),
                                           Create_c_vector(0.7, 0.7),
                                           Create_c_vector(0.7, 0.7)}; //unmoved

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[1]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};
            for (auto&& p_node : mesh.rGetNodes())
            {
                p_node->ClearAppliedForce();
            }

            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);
        }

        // Coverage for adding normal noise
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.back()->SetRegion(LAMINA_REGION);
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            double force_factor = mesh.GetAverageNodeSpacingOfElement(0u, false) / population.GetIntrinsicSpacing();
            p_force->SetSpringConst(1.0 / force_factor);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            p_force->mPreviousLocations = {Create_c_vector(0.5, 0.5),
                                           Create_c_vector(0.6, 0.6),
                                           Create_c_vector(0.6, 0.5), //unmoved
                                           Create_c_vector(0.6, 0.5),
                                           Create_c_vector(0.7, 0.7),
                                           Create_c_vector(0.7, 0.7)}; //unmoved

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};
            for (auto&& p_node : mesh.rGetNodes())
            {
                p_node->ClearAppliedForce();
            }

            p_force->SetAdditiveNormalNoise(true);
            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);
        }

        // Coverage for number of nodes in population changing
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.back()->SetRegion(LAMINA_REGION);
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            // Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            double force_factor = mesh.GetAverageNodeSpacingOfElement(0u, false) / population.GetIntrinsicSpacing();
            p_force->SetSpringConst(1.0 / force_factor);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            p_force->mPreviousLocations = {Create_c_vector(0.5, 0.5),
                                           Create_c_vector(0.6, 0.6),
                                           Create_c_vector(0.6, 0.5), //unmoved
                                           Create_c_vector(0.6, 0.5),
                                           Create_c_vector(0.7, 0.7),
                                           Create_c_vector(0.7, 0.7)}; //unmoved

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};
            for (auto&& p_node : mesh.rGetNodes())
            {
                p_node->ClearAppliedForce();
            }

            Node<2>* p_new_node = new Node<2>(6u, Create_c_vector(0.1, 0.1));
            p_new_node->ClearAppliedForce();
            population.AddNode(p_new_node);
            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);
        }

        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.back()->SetRegion(LAMINA_REGION);
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            // Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            double force_factor = mesh.GetAverageNodeSpacingOfElement(0u, false) / population.GetIntrinsicSpacing();
            p_force->SetSpringConst(1.0 / force_factor);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            p_force->mPreviousLocations = {Create_c_vector(0.5, 0.5),
                                           Create_c_vector(0.6, 0.6),
                                           Create_c_vector(0.6, 0.5), //unmoved
                                           Create_c_vector(0.6, 0.5),
                                           Create_c_vector(0.7, 0.7),
                                           Create_c_vector(0.7, 0.7)}; //unmoved

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};
            for (auto&& p_node : mesh.rGetNodes())
            {
                p_node->ClearAppliedForce();
            }

            Node<2>* p_new_node = new Node<2>(6u, Create_c_vector(0.1, 0.1));
            p_new_node->ClearAppliedForce();
            population.AddNode(p_new_node);
            p_force->UpdatePreviousLocations(population);
        }
    }

    void TestArchivingOfImmersedBoundaryKinematicFeedbackForce()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundaryKinematicFeedbackForce.arch";

        {
            auto p_force = boost::make_shared<ImmersedBoundaryKinematicFeedbackForce<2>>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            p_force->SetSpringConst(1.23);
            p_force->SetAdditiveNormalNoise(true);
            p_force->SetNormalNoiseMean(3.45);
            p_force->SetNormalNoiseStdDev(4.56);

            // Serialize via pointer to the base class
            auto p_base = boost::static_pointer_cast<AbstractImmersedBoundaryForce<2>>(p_force);
            output_arch << p_base;
        }

        {
            boost::shared_ptr<AbstractImmersedBoundaryForce<2>> p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;
            auto p_derived_force = boost::dynamic_pointer_cast<ImmersedBoundaryKinematicFeedbackForce<2>>(p_force);

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(p_derived_force->GetSpringConst(), 1.23, 1e-6);
            TS_ASSERT(p_derived_force->GetAdditiveNormalNoise());
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseMean(), 3.45, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseStdDev(), 4.56, 1e-6);
        }
    }

    void TestImmersedBoundaryLinearDifferentialAdhesionForce()
    {
        // Test member variables
        {
            auto p_force = std::make_unique<ImmersedBoundaryLinearDifferentialAdhesionForce<2>>();

            // Set member variables
            p_force->SetLabelledCellToLabelledCellSpringConst(1.23);
            p_force->SetLabelledCellToCellSpringConst(2.34);
            p_force->SetCellToCellSpringConst(3.45);
            p_force->SetRestLength(4.56);

            TS_ASSERT_DELTA(p_force->GetLabelledCellToLabelledCellSpringConst(), 1.23, 1e-6);
            TS_ASSERT_DELTA(p_force->GetLabelledCellToCellSpringConst(), 2.34, 1e-6);
            TS_ASSERT_DELTA(p_force->GetCellToCellSpringConst(), 3.45, 1e-6);
            TS_ASSERT_DELTA(p_force->GetRestLength(), 4.56, 1e-6);
        }

        // Test correct force is added
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryLinearDifferentialAdhesionForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            p_force->SetLabelledCellToLabelledCellSpringConst(1.0);
            p_force->SetLabelledCellToCellSpringConst(1.0);
            p_force->SetCellToCellSpringConst(1.0);
            p_force->SetRestLength(0.01);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);


            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};

            for (auto&& node : mesh.rGetNodes())
            {
                node->ClearAppliedForce();
            }

            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);

            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[0], 0.333333, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[1], 0.333333, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[0], -0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[1], 0.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[0], 0.6291, 1e-4);
            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[1], 1.2582, 1e-4);

            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[0], -0.333333, 1e-6); //opposite to node 0
            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[1], -0.333333, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[0], 0.0, 1e-6); //opposite to node 1
            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[1], -0.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[0], -0.6291, 1e-4); //not involved
            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[1], -1.2582, 1e-4);

            // If we do the exact same thing again, but with the pairs swapped, the forces should be identical
            for (auto&& node : mesh.rGetNodes())
            {
                node->ClearAppliedForce();
            }

            std::vector<std::pair<Node<2>*, Node<2>*>> swapped_node_pairs = {std::make_pair(nodes[3], nodes[0]),
                                                                             std::make_pair(nodes[4], nodes[1]),
                                                                             std::make_pair(nodes[5], nodes[2])};

            p_force->AddImmersedBoundaryForceContribution(swapped_node_pairs, population);

            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[0], 0.333333, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(0u)->rGetAppliedForce()[1], 0.333333, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[0], -0.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1u)->rGetAppliedForce()[1], 0.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[0], 0.6291, 1e-4); //not involved
            TS_ASSERT_DELTA(mesh.GetNode(2u)->rGetAppliedForce()[1], 1.2582, 1e-4);

            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[0], -0.333333, 1e-6); //opposite to node 0
            TS_ASSERT_DELTA(mesh.GetNode(3u)->rGetAppliedForce()[1], -0.333333, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[0], 0.0, 1e-6); //opposite to node 1
            TS_ASSERT_DELTA(mesh.GetNode(4u)->rGetAppliedForce()[1], -0.0, 1e-6);

            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[0], -0.6291, 1e-4); //not involved
            TS_ASSERT_DELTA(mesh.GetNode(5u)->rGetAppliedForce()[1], -1.2582, 1e-4);
        }

        // Coverage for laminas
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(2, nodes));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, lams, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);

            auto p_force = std::make_shared<ImmersedBoundaryLinearDifferentialAdhesionForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            p_force->SetLabelledCellToLabelledCellSpringConst(1.0);
            p_force->SetLabelledCellToCellSpringConst(1.0);
            p_force->SetCellToCellSpringConst(1.0);
            p_force->SetRestLength(0.01);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};

            TS_ASSERT_THROWS_CONTAINS(p_force->AddImmersedBoundaryForceContribution(node_pairs, population), "presence of lamina");
        }

        // Coverage for labelled cells
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);
            MAKE_PTR(CellLabel, p_label);
            population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);

            auto p_force = std::make_shared<ImmersedBoundaryLinearDifferentialAdhesionForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            p_force->SetLabelledCellToLabelledCellSpringConst(1.0);
            p_force->SetLabelledCellToCellSpringConst(1.0);
            p_force->SetCellToCellSpringConst(1.0);
            p_force->SetRestLength(0.01);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);

            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};

            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);

        }

        // Coverage for two labelled cells
        {
            // Create a minimal cell population
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, Create_c_vector(0.5, 0.5)));
            nodes.push_back(new Node<2>(1u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(2u, Create_c_vector(0.6, 0.5))); //unused
            nodes.push_back(new Node<2>(3u, Create_c_vector(0.6, 0.6)));
            nodes.push_back(new Node<2>(4u, Create_c_vector(0.7, 0.6)));
            nodes.push_back(new Node<2>(5u, Create_c_vector(0.7, 0.7))); //unused

            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0u, {nodes[0], nodes[1], nodes[2]}));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1u, {nodes[3], nodes[4], nodes[5]}));

            ImmersedBoundaryMesh<2, 2> mesh(nodes, elements, {}, 8u, 8u);

            //Create cells
            std::vector<CellPtr> cells;
            auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
            CellsGenerator<NoCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);

            // Create cell population
            ImmersedBoundaryCellPopulation<2> population(mesh, cells);
            population.SetInteractionDistance(10.0);
            MAKE_PTR(CellLabel, p_label);
            population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
            population.GetCellUsingLocationIndex(1)->AddCellProperty(p_label);

            auto p_force = std::make_shared<ImmersedBoundaryLinearDifferentialAdhesionForce<2>>();

            // Fudge so that the force on nodes is what would be "expected" if we did not have to correct for
            // immersed boundary interpolation onto fluid grid
            p_force->SetLabelledCellToLabelledCellSpringConst(1.0);
            p_force->SetLabelledCellToCellSpringConst(1.0);
            p_force->SetCellToCellSpringConst(1.0);
            p_force->SetRestLength(0.01);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10u);
            double dt = SimulationTime::Instance()->GetTimeStep();
            TS_ASSERT_DELTA(dt, 0.1, 1e-6);


            std::vector<std::pair<Node<2>*, Node<2>*>> node_pairs = {std::make_pair(nodes[0], nodes[3]),
                                                                     std::make_pair(nodes[1], nodes[4]),
                                                                     std::make_pair(nodes[2], nodes[5])};

            p_force->AddImmersedBoundaryForceContribution(node_pairs, population);
        }
    }

    void TestArchivingOfImmersedBoundaryLinearDifferentialAdhesionForce()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundaryLinearDifferentialAdhesionForce.arch";

        {
            auto p_force = boost::make_shared<ImmersedBoundaryLinearDifferentialAdhesionForce<2>>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            p_force->SetAdditiveNormalNoise(true);
            p_force->SetNormalNoiseMean(3.45);
            p_force->SetNormalNoiseStdDev(4.56);

            // Serialize via pointer to the base class
            auto p_base = boost::static_pointer_cast<AbstractImmersedBoundaryForce<2>>(p_force);
            output_arch << p_base;
        }

        {
            boost::shared_ptr<AbstractImmersedBoundaryForce<2>> p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;
            auto p_derived_force = boost::dynamic_pointer_cast<ImmersedBoundaryLinearDifferentialAdhesionForce<2>>(p_force);

            // Check member variables have been correctly archived
            TS_ASSERT(p_derived_force->GetAdditiveNormalNoise());
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseMean(), 3.45, 1e-6);
            TS_ASSERT_DELTA(p_derived_force->GetNormalNoiseStdDev(), 4.56, 1e-6);
        }
    }

    void TestImmersedBoundaryForceOutputParameters()
    {
        EXIT_IF_PARALLEL;
        std::string output_directory = "TestForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with ImmersedBoundaryLinearInteractionForce
        {
            ImmersedBoundaryLinearInteractionForce<2> cell_cell_force;
            cell_cell_force.SetSpringConst(1.23);
            cell_cell_force.SetRestLength(2.34);
            cell_cell_force.SetLaminaSpringConstMult(3.45);
            cell_cell_force.SetLaminaRestLengthMult(4.56);
            cell_cell_force.SetAdditiveNormalNoise(true);
            cell_cell_force.SetNormalNoiseMean(5.67);
            cell_cell_force.SetNormalNoiseStdDev(6.78);

            TS_ASSERT_EQUALS(cell_cell_force.GetIdentifier(), "ImmersedBoundaryLinearInteractionForce-2");

            out_stream cell_cell_force_parameter_file = output_file_handler.OpenOutputFile("ib_linear_int.parameters");
            cell_cell_force.OutputImmersedBoundaryForceParameters(cell_cell_force_parameter_file);
            cell_cell_force_parameter_file->close();

            {
                FileFinder generated_file = output_file_handler.FindFile("ib_linear_int.parameters");
                FileFinder reference_file("cell_based/test/data/TestForces/ib_linear_int.parameters",
                                          RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated_file, reference_file);
                TS_ASSERT(comparer.CompareFiles());
            }
        }

        // Test with ImmersedBoundaryMorseInteractionForce
        {
            ImmersedBoundaryMorseInteractionForce<2> cell_cell_force;
            cell_cell_force.SetWellDepth(1.23);
            cell_cell_force.SetRestLength(2.34);
            cell_cell_force.SetLaminaWellDepthMult(3.45);
            cell_cell_force.SetLaminaRestLengthMult(4.56);
            cell_cell_force.SetWellWidth(5.67);
            cell_cell_force.SetAdditiveNormalNoise(true);
            cell_cell_force.SetNormalNoiseMean(6.78);
            cell_cell_force.SetNormalNoiseStdDev(7.89);

            TS_ASSERT_EQUALS(cell_cell_force.GetIdentifier(), "ImmersedBoundaryMorseInteractionForce-2");

            out_stream cell_cell_force_parameter_file = output_file_handler.OpenOutputFile("ib_morse_int.parameters");
            cell_cell_force.OutputImmersedBoundaryForceParameters(cell_cell_force_parameter_file);
            cell_cell_force_parameter_file->close();

            {
                FileFinder generated_file = output_file_handler.FindFile("ib_morse_int.parameters");
                FileFinder reference_file("cell_based/test/data/TestForces/ib_morse_int.parameters",
                                          RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated_file, reference_file);
                TS_ASSERT(comparer.CompareFiles());
            }
        }

        // Test with ImmersedBoundaryLinearMembraneForce
        {
            ImmersedBoundaryLinearMembraneForce<2> membrane_force;
            membrane_force.SetElementSpringConst(5.67);
            membrane_force.SetElementRestLength(6.78);
            membrane_force.SetLaminaSpringConst(7.89);
            membrane_force.SetLaminaRestLength(8.91);
            membrane_force.SetAdditiveNormalNoise(true);
            membrane_force.SetNormalNoiseMean(1.23);
            membrane_force.SetNormalNoiseStdDev(2.34);

            TS_ASSERT_EQUALS(membrane_force.GetIdentifier(), "ImmersedBoundaryLinearMembraneForce-2");

            out_stream membrane_force_parameter_file = output_file_handler.OpenOutputFile("ib_linear_mem.parameters");
            membrane_force.OutputImmersedBoundaryForceParameters(membrane_force_parameter_file);
            membrane_force_parameter_file->close();

            {
                FileFinder generated_file = output_file_handler.FindFile("ib_linear_mem.parameters");
                FileFinder reference_file("cell_based/test/data/TestForces/ib_linear_mem.parameters",
                                          RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated_file, reference_file);
                TS_ASSERT(comparer.CompareFiles());
            }
        }

        // Test with ImmersedBoundaryMorseMembraneForce
        {
            ImmersedBoundaryMorseMembraneForce<2> membrane_force;
            membrane_force.SetElementWellDepth(1.23);
            membrane_force.SetElementRestLength(2.34);
            membrane_force.SetLaminaWellDepth(3.45);
            membrane_force.SetLaminaRestLength(4.56);
            membrane_force.SetWellWidth(5.67);
            membrane_force.SetAdditiveNormalNoise(true);
            membrane_force.SetNormalNoiseMean(6.78);
            membrane_force.SetNormalNoiseStdDev(7.89);

            TS_ASSERT_EQUALS(membrane_force.GetIdentifier(), "ImmersedBoundaryMorseMembraneForce-2");

            out_stream membrane_force_parameter_file = output_file_handler.OpenOutputFile("ib_morse_mem.parameters");
            membrane_force.OutputImmersedBoundaryForceParameters(membrane_force_parameter_file);
            membrane_force_parameter_file->close();

            {
                FileFinder generated_file = output_file_handler.FindFile("ib_morse_mem.parameters");
                FileFinder reference_file("cell_based/test/data/TestForces/ib_morse_mem.parameters",
                                          RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated_file, reference_file);
                TS_ASSERT(comparer.CompareFiles());
            }
        }

        // Test with ImmersedBoundaryDifferentialAdhesion
        {
            ImmersedBoundaryLinearDifferentialAdhesionForce<2> cell_cell_force;
            cell_cell_force.SetRestLength(2.34);
            cell_cell_force.SetAdditiveNormalNoise(true);
            cell_cell_force.SetNormalNoiseMean(5.67);
            cell_cell_force.SetNormalNoiseStdDev(6.78);

            TS_ASSERT_EQUALS(cell_cell_force.GetIdentifier(), "ImmersedBoundaryLinearDifferentialAdhesionForce-2");

            out_stream cell_cell_force_parameter_file = output_file_handler.OpenOutputFile("ib_diff_ad.parameters");
            cell_cell_force.OutputImmersedBoundaryForceParameters(cell_cell_force_parameter_file);
            cell_cell_force_parameter_file->close();

            {
                FileFinder generated_file = output_file_handler.FindFile("ib_diff_ad.parameters");
                FileFinder reference_file("cell_based/test/data/TestForces/ib_diff_ad.parameters",
                                          RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated_file, reference_file);
                TS_ASSERT(comparer.CompareFiles());
            }
        }

        // Test with ImmersedBoundaryKinematicFeedbackForce
        {
            ImmersedBoundaryKinematicFeedbackForce<2> kinematic_cell_force;
            kinematic_cell_force.SetAdditiveNormalNoise(true);
            kinematic_cell_force.SetNormalNoiseMean(5.67);
            kinematic_cell_force.SetNormalNoiseStdDev(6.78);

            TS_ASSERT_EQUALS(kinematic_cell_force.GetIdentifier(), "ImmersedBoundaryKinematicFeedbackForce-2");

            out_stream kinematic_cell_force_parameter_file = output_file_handler.OpenOutputFile("ib_kinematic.parameters");
            kinematic_cell_force.OutputImmersedBoundaryForceParameters(kinematic_cell_force_parameter_file);
            kinematic_cell_force_parameter_file->close();

            {
                FileFinder generated_file = output_file_handler.FindFile("ib_kinematic.parameters");
                FileFinder reference_file("cell_based/test/data/TestForces/ib_kinematic.parameters",
                                          RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated_file, reference_file);
                TS_ASSERT(comparer.CompareFiles());
            }
        }
    }
};

#endif /*TESTIMMERSEDBOUNDARYFORCES_HPP_*/
