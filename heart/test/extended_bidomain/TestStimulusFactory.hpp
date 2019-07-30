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


#ifndef TESTSTIMULUSFACTORY_HPP_
#define TESTSTIMULUSFACTORY_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <iostream>
#include <vector>

#include "AbstractStimulusFactory.hpp"
#include "HeartConfig.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractStimulusFunction.hpp"
#include "SimpleStimulus.hpp"
#include "ElectrodesStimulusFactory.hpp"
#include "ExtendedBidomainProblem.hpp"

#include "ElectrodesStimulusFactory.hpp"
#include "ChasteNodesList.hpp"
#include "ChasteCuboid.hpp"
#include "AbstractChasteRegion.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "ZeroNetChargeElectrodes.hpp"
#include "PetscSetupAndFinalize.hpp"

class MyStimulusFactory: public AbstractStimulusFactory<3>
{

public:
    MyStimulusFactory() : AbstractStimulusFactory<3>()
    {
    }

    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];

        boost::shared_ptr<SimpleStimulus> p_stimulus;
        if ((x < 0.0005) && (y < 0.0005) && (z < 0.0005))
        {
            p_stimulus.reset ( new SimpleStimulus(-428000, 1.0, 0.1) );
        }
        else
        {
            p_stimulus.reset ( new SimpleStimulus(0.0, 0.5, 0.1) );

        }
        return p_stimulus;
    }
};

class TestStimulusFactory : public CxxTest::TestSuite
{
public:

    void TestDefaultImplementation()
    {
        HeartConfig::Instance()->Reset();

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2,2,2);

        AbstractStimulusFactory<3> my_factory;
        my_factory.SetMesh(&mesh);

        for (AbstractMesh<3,3>::NodeIterator iter=mesh.GetNodeIteratorBegin();
             iter != mesh.GetNodeIteratorEnd();
             ++iter)
        {
            boost::shared_ptr<AbstractStimulusFunction> stimulus = my_factory.CreateStimulusForNode(&(*iter));

            TS_ASSERT_EQUALS(stimulus->GetStimulus(0), 0);
            TS_ASSERT_EQUALS(stimulus->GetStimulus(1), 0);
            TS_ASSERT_EQUALS(stimulus->GetStimulus(2), 0);
            TS_ASSERT_EQUALS(stimulus->GetStimulus(258), 0);
        }
    }

    void TestOneFactory()
    {
        HeartConfig::Instance()->Reset();

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2,2,2);

        MyStimulusFactory my_factory;

        //cover getting the mesh when it was not set
        TS_ASSERT_THROWS_THIS(my_factory.GetMesh(),
        "The mesh object has not been set in the stimulus factory");

        my_factory.SetMesh(&mesh);
        TS_ASSERT_EQUALS(my_factory.GetNumberOfCells(), 27u);

        //probe node 0 (in the stimulated corner)
        boost::shared_ptr<AbstractStimulusFunction> stimulus = my_factory.CreateStimulusForNode(mesh.GetNode(0));

        TS_ASSERT_EQUALS(stimulus->GetStimulus(0), 0);//before it starts
        TS_ASSERT_EQUALS(stimulus->GetStimulus(1.01), -428000);//during the stimulus
        TS_ASSERT_EQUALS(stimulus->GetStimulus(2), 0);//after it finished

        //probe node 15 (not in the stimulated corner)
        boost::shared_ptr<AbstractStimulusFunction> stimulus_zero = my_factory.CreateStimulusForNode(mesh.GetNode(15));

        TS_ASSERT_EQUALS(stimulus_zero->GetStimulus(0), 0);
        TS_ASSERT_EQUALS(stimulus_zero->GetStimulus(1.01), 0);
        TS_ASSERT_EQUALS(stimulus_zero->GetStimulus(2), 0);
    }

    void TestComputeContributiontoRHS()
    {
        HeartConfig::Instance()->Reset();

        /**
         * In this test the mesh NEEDS to be a Distributed mesh.
         * The loop over nodes in the ComputeElectrodeTotalFlux method only loops over the nodes we own (using a NodeIterator).
         * If we run in parallel and uses a Tetrahedral mesh, each proc owns the whole mesh, hence the final flux computation
         * will be N times the one in sequential (where N is number of procesors).
         * In a normal simulation, this is not a problem, because the computed fluxes are used to compute their ratio,
         * hence if both electrodes fluxes are mulitplied by N it doesn't really matter.
         * Here, instead, we are testing the ComputeElectrodeTotalFlux method specifically,
         * so we want consistent numbers to compare against on any number of procs.
         *
         */

        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(30);

        std::vector<std::pair<AbstractChasteRegion<1>*, AbstractChasteRegion<1>*> > electrode_pairs;
        std::vector<double> magnitudes;
        std::vector<double> durations;
        std::vector<double> periods;
        std::vector<double> starts;
        std::vector<double> ends;

        ChastePoint<1> LowerPoint_1(-0.01);
        ChastePoint<1> UpperPoint_1(0.15);//1 node here

        ChastePoint<1> LowerPoint_2(28.5);
        ChastePoint<1> UpperPoint_2(31.0);//3 nodes here

        ChasteCuboid<1>* p_cuboid_1 = new ChasteCuboid<1>(LowerPoint_1, UpperPoint_1);
        ChasteCuboid<1>* p_cuboid_2 = new ChasteCuboid<1>(LowerPoint_2, UpperPoint_2);
        std::pair<AbstractChasteRegion<1>*, AbstractChasteRegion<1>* > first_pair;
        first_pair.first = p_cuboid_1;
        first_pair.second = p_cuboid_2;

        electrode_pairs.push_back(first_pair);
        magnitudes.push_back(75.0);
        durations.push_back(15.0);
        periods.push_back(20000);
        starts.push_back(0);
        ends.push_back(DBL_MAX);

        ElectrodesStimulusFactory<1> electrodes(electrode_pairs, magnitudes, durations, periods, starts, ends);
        electrodes.SetMesh(&mesh);
        TS_ASSERT_DELTA(electrodes.ComputeElectrodeTotalFlux(p_cuboid_1, magnitudes[0]), 62.4999, 1e-3);
        TS_ASSERT_DELTA(electrodes.ComputeElectrodeTotalFlux(p_cuboid_2, magnitudes[0]), 62.4999*3, 1e-3);

        delete p_cuboid_1;
        delete p_cuboid_2;
    }

    /**
     * Convenience heartconfig setups for the tests below
     */
    void SetupParameters()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01,0.1,1.0);
        HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");//block jacobi is slow

        HeartConfig::Instance()->SetSimulationDuration(0.1);  //ms.

        /** Output visualization options, we ask for meshalyzer and cmgui **/
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(false);
    }


    void TestRegularCubeChasteNodesList()
    {
        SetupParameters();

        std::string directory = "RegularCube";
        std::string filename = "extended3d";

        /**Where to utput stuff**/
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix(filename);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> first_cell_factory(0.0);
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> second_cell_factory(0.0);

        ExtendedBidomainProblem<3> extended_problem( &first_cell_factory , &second_cell_factory);
        extended_problem.SetMesh(&mesh);

        std::vector<std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>*> > electrode_pairs;
        std::vector<double> magnitudes;
        std::vector<double> durations;
        std::vector<double> periods;
        std::vector<double> starts;
        std::vector<double> ends;

        std::vector<Node<3>* > arary_of_nodes_of_electrode_1;
        std::vector<Node<3>* > arary_of_nodes_of_electrode_2;

        std::vector<unsigned> indices;
        indices.push_back(0u);//first electrode, near origin
        indices.push_back(4u);//first electrode
        indices.push_back(21u);//first electrode
        indices.push_back(6u);//second electrope, opposite corner
        indices.push_back(30u);//second electrode

        if (mesh.rGetNodePermutation().size() > 0)
        {
            for (unsigned i = 0; i < indices.size(); i++)
            {
                indices[i] = mesh.rGetNodePermutation()[ indices[i] ];
            }
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(indices[0]))
        {
            arary_of_nodes_of_electrode_1.push_back(mesh.GetNode(indices[0]));
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(indices[1]))
        {
            arary_of_nodes_of_electrode_1.push_back(mesh.GetNode(indices[1]));
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(indices[2]))
        {
            arary_of_nodes_of_electrode_1.push_back(mesh.GetNode(indices[2]));
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(indices[3]))
        {
            arary_of_nodes_of_electrode_2.push_back(mesh.GetNode(indices[3]));
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(indices[4]))
        {
            arary_of_nodes_of_electrode_2.push_back(mesh.GetNode(indices[4]));
        }

        ChasteNodesList<3>* p_node_list_1 = new ChasteNodesList<3>(arary_of_nodes_of_electrode_1);
        ChasteNodesList<3>* p_node_list_2 = new ChasteNodesList<3>(arary_of_nodes_of_electrode_2);
        std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>* > first_pair;
        first_pair.first = p_node_list_1;
        first_pair.second = p_node_list_2;

        electrode_pairs.push_back(first_pair);
        magnitudes.push_back(-91800.0);
        durations.push_back(15.0);
        periods.push_back(20000);
        starts.push_back(0);
        ends.push_back(DBL_MAX);

        ElectrodesStimulusFactory<3> electrodes(electrode_pairs, magnitudes, durations, periods, starts, ends);
        extended_problem.SetExtracellularStimulusFactory(&electrodes);
        extended_problem.Initialise();
        extended_problem.Solve();

        delete p_node_list_1;
        delete p_node_list_2;
    }

    void TestIrregularCubeChastecuboid()
    {
        SetupParameters();

        std::string directory = "IrregularCube";
        std::string filename = "extended3d";

        /**Where to utput stuff**/
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix(filename);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> first_cell_factory(0.0);
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> second_cell_factory(0.0);

        std::vector<std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>*> > electrode_pairs;
        std::vector<double> magnitudes;
        std::vector<double> durations;
        std::vector<double> periods;
        std::vector<double> starts;
        std::vector<double> ends;

        ChastePoint<3> LowerPoint_1(-0.01, -0.01, -0.01);
        ChastePoint<3> UpperPoint_1(0.013, 0.013, 0.013);

        ChastePoint<3> LowerPoint_2(0.02474, -0.01, 0.0248);
        ChastePoint<3> UpperPoint_2(0.0251, 0.013, 0.0289);

        ChasteCuboid<3>* p_cuboid_1 = new ChasteCuboid<3>(LowerPoint_1, UpperPoint_1);
        ChasteCuboid<3>* p_cuboid_2 = new ChasteCuboid<3>(LowerPoint_2, UpperPoint_2);
        std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>* > first_pair;
        first_pair.first = p_cuboid_1;
        first_pair.second = p_cuboid_2;

        electrode_pairs.push_back(first_pair);
        magnitudes.push_back(-91800.0);
        durations.push_back(15.0);
        periods.push_back(20000);
        starts.push_back(0);
        ends.push_back(DBL_MAX);

        ElectrodesStimulusFactory<3> electrodes(electrode_pairs, magnitudes, durations, periods, starts, ends);

        //test this in a real problem.
        ExtendedBidomainProblem<3> extended_problem( &first_cell_factory , &second_cell_factory);
        extended_problem.SetMesh(&mesh);
        extended_problem.SetExtracellularStimulusFactory(&electrodes);
        TS_ASSERT_THROWS_NOTHING(extended_problem.Initialise());

        delete p_cuboid_1;
        delete p_cuboid_2;
    }

    void TestRegularCubeZeroNetChargeElectrodes()
    {
        SetupParameters();

        std::string directory = "IrregularCube";
        std::string filename = "extended3d";

        /**Where to utput stuff**/
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix(filename);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> first_cell_factory(0.0);
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> second_cell_factory(0.0);

        ExtendedBidomainProblem<3> extended_problem( &first_cell_factory , &second_cell_factory);
        extended_problem.SetMesh(&mesh);

        std::vector<std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>*> > electrode_pairs;
        std::vector<double> magnitudes;
        std::vector<double> durations;
        std::vector<double> periods;
        std::vector<double> starts;
        std::vector<double> ends;

        ChastePoint<3> LowerPoint_1(-0.01, -0.01, -0.01);
        ChastePoint<3> UpperPoint_1(0.1, 2.0, 0.1);

        ChastePoint<3> LowerPoint_2(0.9, -0.1, 0.9);
        ChastePoint<3> UpperPoint_2(1.1, 2.0, 1.1);

        ChasteCuboid<3>* p_cuboid_1 = new ChasteCuboid<3>(LowerPoint_1, UpperPoint_1);
        ChasteCuboid<3>* p_cuboid_2 = new ChasteCuboid<3>(LowerPoint_2, UpperPoint_2);
        std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>* > first_pair;
        first_pair.first = p_cuboid_1;
        first_pair.second = p_cuboid_2;

        electrode_pairs.push_back(first_pair);
        magnitudes.push_back(-91800.0);
        durations.push_back(15.0);
        periods.push_back(20000);
        starts.push_back(0);
        ends.push_back(DBL_MAX);

        ZeroNetChargeElectrodes<3> electrodes(electrode_pairs, magnitudes, durations, periods, starts, ends);
        extended_problem.SetExtracellularStimulusFactory(&electrodes);

        extended_problem.Initialise();
        extended_problem.Solve();

        delete p_cuboid_1;
        delete p_cuboid_2;
    }

    void TestGroundingSecondElectrode()
    {
        SetupParameters();

        //this one needs to be a bit longer so that we are sure the second electrode stays at zero...
        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms.

        std::string directory = "GroundedElectrode";
        std::string filename = "extended3d";

        /**Where to utput stuff**/
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix(filename);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> first_cell_factory(0.0);
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> second_cell_factory(0.0);

        ExtendedBidomainProblem<3> extended_problem( &first_cell_factory , &second_cell_factory);
        extended_problem.SetMesh(&mesh);

        std::vector<std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>*> > electrode_pairs;
        std::vector<double> magnitudes;
        std::vector<double> durations;
        std::vector<double> periods;
        std::vector<double> starts;
        std::vector<double> ends;

        ChastePoint<3> LowerPoint_1(-0.01, -0.01, -0.01);
        ChastePoint<3> UpperPoint_1(0.1, 2.0, 0.1);

        ChastePoint<3> LowerPoint_2(0.9, -0.1, 0.9);
        ChastePoint<3> UpperPoint_2(1.1, 2.0, 1.1);

        ChasteCuboid<3>* p_cuboid_1 = new ChasteCuboid<3>(LowerPoint_1, UpperPoint_1);
        ChasteCuboid<3>* p_cuboid_2 = new ChasteCuboid<3>(LowerPoint_2, UpperPoint_2);
        std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>* > first_pair;
        first_pair.first = p_cuboid_1;
        first_pair.second = p_cuboid_2;
        electrode_pairs.push_back(first_pair);
        magnitudes.push_back(-9180.0);
        durations.push_back(15.0);
        periods.push_back(20000);
        starts.push_back(0);
        ends.push_back(DBL_MAX);

        ElectrodesStimulusFactory<3> electrodes(electrode_pairs, magnitudes, durations, periods, starts, ends);

        electrodes.SetMesh(&mesh);
        electrodes.GroundSecondElectrode(true);

        extended_problem.SetExtracellularStimulusFactory(&electrodes);

        extended_problem.Initialise();
        extended_problem.Solve();

        std::vector<AbstractChasteRegion<3>* > grounded_regions = electrodes.GetRegionsToBeGrounded();
        TS_ASSERT_LESS_THAN(0u, grounded_regions.size());

        unsigned probe_index_1 = 5u;
        unsigned probe_index_2 = 6u;

        // In parallel work out the new indices
        if (mesh.rGetNodePermutation().size() > 0)
        {
            probe_index_1 =  mesh.rGetNodePermutation()[probe_index_1];
            probe_index_2 =  mesh.rGetNodePermutation()[probe_index_2];
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(probe_index_1))
        {
            TS_ASSERT(grounded_regions[0]->DoesContain(mesh.GetNode(probe_index_1)->GetPoint()));
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(probe_index_2))
        {
            TS_ASSERT(grounded_regions[0]->DoesContain(mesh.GetNode(probe_index_2)->GetPoint()));
        }

        //now check that the second electrode was actually grounded
        Hdf5DataReader reader_extended(directory, filename);

        for (unsigned node_index = 0; node_index < mesh.GetNumNodes(); node_index++)
        {
            if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(node_index))
            {
                if (grounded_regions[0]->DoesContain(mesh.GetNode(node_index)->GetPoint()))
                {
                    unsigned node = mesh.GetNode(node_index)->GetIndex();
                    std::vector<double> phi_e_extended = reader_extended.GetVariableOverTime("Phi_e", node);
                    for (unsigned j = 0; j <  phi_e_extended.size(); j++)
                    {
                        TS_ASSERT_EQUALS(phi_e_extended[j], 0);//should be zero at all times
                    }
                }
            }
        }
        delete p_cuboid_1;
        delete p_cuboid_2;
    }

    void TestRegularCubeIntersecting()
    {
        DistributedTetrahedralMesh<3,3> mesh;
        double width = 2.0;
        double height=2.0;
        double depth=2.0;
        double inter_node_space = 0.5;
        mesh.ConstructRegularSlabMesh(inter_node_space, width, height, depth);

        std::vector<std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>*> > electrode_pairs;
        std::vector<double> magnitudes;
        std::vector<double> durations;
        std::vector<double> periods;
        std::vector<double> starts;
        std::vector<double> ends;

        //first electrodes, x=0 plane
        ChastePoint<3> LowerPoint_1(-1, -1, -1);
        ChastePoint<3> UpperPoint_1(0.15, 3.0, 3.0);

        //second electrode, z=0 plane -> should intersect the first one at x=0,z=0, any y
        ChastePoint<3> LowerPoint_2(-1, -1, -1);
        ChastePoint<3> UpperPoint_2(3.0, 3.0, 0.15);

        ChasteCuboid<3>* p_cuboid_1 = new ChasteCuboid<3>(LowerPoint_1, UpperPoint_1);
        ChasteCuboid<3>* p_cuboid_2 = new ChasteCuboid<3>(LowerPoint_2, UpperPoint_2);
        std::pair<AbstractChasteRegion<3>*, AbstractChasteRegion<3>* > first_pair;
        first_pair.first = p_cuboid_1;
        first_pair.second = p_cuboid_2;
        electrode_pairs.push_back(first_pair);
        magnitudes.push_back(-91800.0);
        durations.push_back(15.0);
        periods.push_back(20000);
        starts.push_back(0);
        ends.push_back(DBL_MAX);

        ElectrodesStimulusFactory<3> electrodes(electrode_pairs, magnitudes, durations, periods, starts, ends);
        electrodes.SetMesh(&mesh);
        TS_ASSERT_THROWS_THIS( electrodes.SetCompatibleExtracellularStimulus(),
                        "Two or more electrodes intersect with each other");

        //covers exception in constructor.
        periods.push_back(25000);
        TS_ASSERT_THROWS_THIS(ElectrodesStimulusFactory<3> electrodes_2(electrode_pairs, magnitudes, durations, periods, starts, ends),
                                "Vector of electrode pairs and vector of stimulation paremeters must have the same size");

        delete p_cuboid_1;
        delete p_cuboid_2;

    }
};

#endif /*TESTSTIMULUSFACTORY_HPP_*/
