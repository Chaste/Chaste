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


#ifndef TESTEXTENDEDBIDOMAINTISSUE_HPP_
#define TESTEXTENDEDBIDOMAINTISSUE_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <iostream>
#include <vector>

#include "ArchiveOpener.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"
#include "ExtendedBidomainTissue.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractStimulusFunction.hpp"
#include "DistributedVectorFactory.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "TetrahedralMesh.hpp"
#include "AbstractStimulusFactory.hpp"
#include <petsc.h>
#include "PetscSetupAndFinalize.hpp"

class StimulatedCellFactory: public AbstractCardiacCellFactory<3>
{
private:
        boost::shared_ptr<SimpleStimulus> mpStimulus;
        //static const double magnitude = -105.0*1400;
public:
    StimulatedCellFactory() : AbstractCardiacCellFactory<3>(),
            mpStimulus ( new SimpleStimulus(-105.0*1400, 1.0))/*amplitude, duration (ms)*/
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        CellLuoRudy1991FromCellML* first_cell;
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];

        if ((x < 0.005) && (y < 0.005) && (z < 0.005))
        {
            first_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            first_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }

        return first_cell;
    }
};

class UnStimulatedCellFactory: public AbstractCardiacCellFactory<3>
{

public:
    UnStimulatedCellFactory() : AbstractCardiacCellFactory<3>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        CellLuoRudy1991FromCellML* second_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        return second_cell;
    }
};

class ExtracellularStimulusFactory: public AbstractStimulusFactory<3>
{

public:
    ExtracellularStimulusFactory() : AbstractStimulusFactory<3>()
    {
    }

    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];

        boost::shared_ptr<SimpleStimulus> p_stimulus;
        if ((x < 0.005) && (y < 0.005) && (z < 0.005))
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

/**this is just for testing in 2D*/
class StimulusFactory2D: public AbstractStimulusFactory<2>
{
public:
    StimulusFactory2D() : AbstractStimulusFactory<2>()
    {}
    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<2>* pNode)
    {
        boost::shared_ptr<SimpleStimulus> p_stimulus ( new SimpleStimulus(-428000, 1.0, 0.1) );
        return p_stimulus;
    }
};

class SimpleConductivityModifier : public AbstractConductivityModifier<2,2>
{
private:
    c_matrix<double,2,2> mTensor;

public:
    SimpleConductivityModifier()
        : AbstractConductivityModifier<2,2>(),
          mTensor(zero_matrix<double>(2,2))
          // Conductivity tensors are diagonal, so we can only need to do the diagonal if we initialise our tensor with the above line
    {
    }

    c_matrix<double,2,2>& rCalculateModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,2,2>& rOriginalConductivity, unsigned domainIndex)
    {
        // Increase by factor of two for element 1
        if (elementIndex == 1)
        {
            for ( unsigned dim=0; dim<2; dim++)
            {
                mTensor(dim,dim) = 2.0*rOriginalConductivity(dim,dim);
            }
        }
        else
        {
            for ( unsigned dim=0; dim<2; dim++)
            {
                mTensor(dim,dim) = rOriginalConductivity(dim,dim);
            }
        }
        return mTensor;
    }
};
class TestExtendedBidomainTissue : public CxxTest::TestSuite
{
public:

    void TestExtendedBidomainTissueSolveCellSystems( void )
    {
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2,2,2);

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;
        ExtracellularStimulusFactory extracellular_stimulus_factory;

        stimulated_cell_factory.SetMesh(&mesh);
        unstimulated_cell_factory.SetMesh(&mesh);
        extracellular_stimulus_factory.SetMesh(&mesh);

        ExtendedBidomainTissue<3>  extended_bidomain_tissue( &stimulated_cell_factory,  &unstimulated_cell_factory, &extracellular_stimulus_factory);

        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;
        double initial_phie=0.0;

        // initial condition;
        DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
        Vec extended_vec = p_factory->CreateVec(3);
        DistributedVector extended_bidomain_potentials = p_factory->CreateDistributedVector(extended_vec);

        DistributedVector::Stripe voltage_first_cell(extended_bidomain_potentials,0);
        DistributedVector::Stripe voltage_second_cell(extended_bidomain_potentials,1);
        DistributedVector::Stripe phi_e(extended_bidomain_potentials,2);

        for (DistributedVector::Iterator index=extended_bidomain_potentials.Begin();
             index != extended_bidomain_potentials.End();
             ++index)
        {
            voltage_first_cell[index] = initial_voltage;
            voltage_second_cell[index] = initial_voltage;
            phi_e[index] = initial_phie;
        }

        extended_bidomain_potentials.Restore();
        extended_bidomain_tissue.SolveCellSystems(extended_vec, 0.0, 0.5);

        //Check ionic currents against hardcoded values
        TS_ASSERT_DELTA(extended_bidomain_tissue.rGetIionicCacheReplicated()[0], 0.0034, 1e-4);
        TS_ASSERT_DELTA(extended_bidomain_tissue.rGetIionicCacheReplicatedSecondCell()[0], 0.0034, 1e-4);

        // Check that the first cell and extracellular stimulus have the right stimulus value at node 0 and 1
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularStimulusCacheReplicated()[0], -105.0*1400);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularStimulusCacheReplicated()[0], -428000);

        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularStimulusCacheReplicated()[1], 0);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularStimulusCacheReplicated()[1], 0);

        // Second cell is unstimulated throughout
        for (unsigned node_index = mesh.GetDistributedVectorFactory()->GetLow();
             node_index < mesh.GetDistributedVectorFactory()->GetHigh();
             node_index++)
        {
            TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularStimulusCacheReplicatedSecondCell()[node_index], 0);
        }

        PetscTools::Destroy(extended_vec);
    }

    /**This test checks heterogeneous conductivities*/
    void TestExtendedTissueHeterogeneous3D()
    {
        HeartConfig::Instance()->Reset();
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<ChasteCuboid<3> > heterogeneity_area;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;

        //first cuboid include element 0
        ChastePoint<3> cornerA(-1, -1, 0);
        ChastePoint<3> cornerB(0.1, 0.2, 0.2);
        ChasteCuboid<3> cuboid_1(cornerA, cornerB);
        heterogeneity_area.push_back(cuboid_1);

        //second cuboid include element 4
        ChastePoint<3> cornerC(0.11, 0.0, 0);
        ChastePoint<3> cornerD(0.2, 0.11, 0.2);
        ChasteCuboid<3> cuboid_2(cornerC, cornerD);

        heterogeneity_area.push_back(cuboid_2);

        //within the first area
        intra_conductivities.push_back( Create_c_vector(1.0, 2.0, 3.0) );
        extra_conductivities.push_back( Create_c_vector(51.0, 52.0, 53.0) );

        //within the second area
        intra_conductivities.push_back( Create_c_vector(11.0, 22.0, 33.0) );
        extra_conductivities.push_back( Create_c_vector(151.0, 152.0, 153.0) );

        HeartConfig::Instance()->SetConductivityHeterogeneities(heterogeneity_area, intra_conductivities, extra_conductivities);

        //elsewhere
        double isotropic_intra_conductivity=15.0;
        double isotropic_extra_conductivity=65.0;

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(isotropic_intra_conductivity, isotropic_intra_conductivity, isotropic_intra_conductivity));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(isotropic_extra_conductivity, isotropic_extra_conductivity, isotropic_extra_conductivity));

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;
        ExtracellularStimulusFactory extracellular_stimulus_factory;

        stimulated_cell_factory.SetMesh(&mesh);
        unstimulated_cell_factory.SetMesh(&mesh);
        extracellular_stimulus_factory.SetMesh(&mesh);

        ExtendedBidomainTissue<3>  extended_bidomain_tissue( &stimulated_cell_factory,  &unstimulated_cell_factory, &extracellular_stimulus_factory);

        extended_bidomain_tissue.SetIntracellularConductivitiesSecondCell(Create_c_vector(isotropic_intra_conductivity, isotropic_intra_conductivity, isotropic_intra_conductivity));
        extended_bidomain_tissue.CreateIntracellularConductivityTensorSecondCell();

        //first cell
//        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(0u)(0,0),1.0);//within first cuboid
        //Line above commented due to curious problem with IntelProduction interprocedural optimisation
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(4u)(0,0),11.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(4u)(1,1),22.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(8u)(0,0),15.0);//elsewhere, e.g. element 8

        //second cell, should be the same as first cell
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(0u)(0,0),1.0);//within first cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(4u)(0,0),11.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(4u)(1,1),22.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(8u)(0,0),15.0);//elsewhere, e.g. element 8

        //sigma_e
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(0u)(0,0),51.0);//within first cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(4u)(0,0),151.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(4u)(1,1),152.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(8u)(0,0),65.0);//elsewhere, e.g. element 8

    }

    void TestExtendedTissueHeterogeneousConductivities2D()
    {
        HeartConfig::Instance()->Reset();
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        //HeartConfig setup needs to be in 3D anyway (hardcoded in HeartConfig).
        std::vector<ChasteCuboid<3> > heterogeneity_area;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;

        //first cuboid includes element 0
        ChastePoint<3> cornerA(-1, -1,-1);
        ChastePoint<3> cornerB(0.48, 2.0, 0.48);
        ChasteCuboid<3> cuboid_1(cornerA, cornerB);

        heterogeneity_area.push_back(cuboid_1);

        //second cuboid includes element 2
        ChastePoint<3> cornerC(0.52, -1, 0.52);
        ChastePoint<3> cornerD(2, 2, 2);
        ChasteCuboid<3> cuboid_2(cornerC, cornerD);

        heterogeneity_area.push_back(cuboid_2);

        //within the first area
        intra_conductivities.push_back( Create_c_vector(1.0, 2.0, 3.0) );
        extra_conductivities.push_back( Create_c_vector(51.0, 52.0, 53.0) );

        //within the second area
        intra_conductivities.push_back( Create_c_vector(11.0, 22.0, 33.0) );
        extra_conductivities.push_back( Create_c_vector(151.0, 152.0, 153.0) );

        HeartConfig::Instance()->SetConductivityHeterogeneities(heterogeneity_area, intra_conductivities, extra_conductivities);

        //elsewhere
        double isotropic_intra_conductivity=15.0;
        double isotropic_extra_conductivity=65.0;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(isotropic_intra_conductivity, isotropic_intra_conductivity, isotropic_intra_conductivity));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(isotropic_extra_conductivity, isotropic_extra_conductivity, isotropic_extra_conductivity));

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory_1;
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory_2;
        StimulusFactory2D extracellular_stimulus_factory;

        cell_factory_1.SetMesh(&mesh);
        cell_factory_2.SetMesh(&mesh);
        extracellular_stimulus_factory.SetMesh(&mesh);

        //2D tissue
        ExtendedBidomainTissue<2>  extended_bidomain_tissue( &cell_factory_1,  &cell_factory_2, &extracellular_stimulus_factory);

        // Do conductivity modifier here too (for coverage)
        SimpleConductivityModifier conductivity_modifier;
        extended_bidomain_tissue.SetConductivityModifier( &conductivity_modifier );

        extended_bidomain_tissue.SetIntracellularConductivitiesSecondCell(Create_c_vector(isotropic_intra_conductivity, isotropic_intra_conductivity));
        extended_bidomain_tissue.CreateIntracellularConductivityTensorSecondCell();

        //first cell
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(0u)(0,0),1.0);//within first cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(1u)(0,0),30.0);//within no cuboid (modified from 15 to 30!)
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(2u)(1,1),22.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensor(3u)(0,0),15.0);//within no cuboid

        //second cell, should be the same as first cell
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(0u)(0,0),1.0);//within first cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(1u)(0,0),30.0);//within no cuboid (modified from 15 to 30!)
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(2u)(1,1),22.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetIntracellularConductivityTensorSecondCell(3u)(0,0),15.0);//within no cuboid

        //sigma_e
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(0u)(0,0),51.0);//within first cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(1u)(0,0),130.0);//within no cuboid (modified from 65 to 130!)
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(2u)(1,1),152.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetExtracellularConductivityTensor(3u)(0,0),65.0);//within no cuboid
    }

    void TestExtendedTissueHeterogeneousGgap3D()
    {
        HeartConfig::Instance()->Reset();
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<boost::shared_ptr<AbstractChasteRegion<3> > > heterogeneity_areas;
        std::vector<double> Ggap_values;

        //first cuboid include node 0
        ChastePoint<3> cornerA(-1, -1, 0);
        ChastePoint<3> cornerB(0.001, 0.001, 0.001);
        boost::shared_ptr<ChasteCuboid<3> > p_cuboid_1(new ChasteCuboid<3>(cornerA, cornerB));
        heterogeneity_areas.push_back(p_cuboid_1);

        //second cuboid include node 6
        ChastePoint<3> cornerC(0.199, 0.199, 0.199);
        ChastePoint<3> cornerD(0.25, 0.25, 0.25);
        boost::shared_ptr<ChasteCuboid<3> > p_cuboid_2(new ChasteCuboid<3>(cornerC, cornerD));
        ChasteCuboid<3> cuboid_2(cornerC, cornerD);

        heterogeneity_areas.push_back(p_cuboid_2);

        //within the first area
        Ggap_values.push_back(143.0);
        //within the second area
        Ggap_values.push_back(9143.0);
        //elsewhere
        double isotropic_ggap=586.0;

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;
        ExtracellularStimulusFactory extracellular_stimulus_factory;

        stimulated_cell_factory.SetMesh(&mesh);
        unstimulated_cell_factory.SetMesh(&mesh);
        extracellular_stimulus_factory.SetMesh(&mesh);

        ExtendedBidomainTissue<3>  extended_bidomain_tissue( &stimulated_cell_factory,  &unstimulated_cell_factory, &extracellular_stimulus_factory);

        extended_bidomain_tissue.SetGGap(isotropic_ggap);//this is what the problem class does first

        extended_bidomain_tissue.SetGgapHeterogeneities(heterogeneity_areas, Ggap_values);
        extended_bidomain_tissue.CreateGGapConductivities();

        unsigned probe_node_1 = 0u;
        unsigned probe_node_2 = 6u;
        unsigned probe_node_3 = 5u;

        const std::vector<unsigned>& r_permutation = mesh.rGetNodePermutation();
        if (!r_permutation.empty())
        {
            probe_node_1 = r_permutation[0u];//within first cuboid
            probe_node_2 = r_permutation[6u];//within second cuboid
            probe_node_3 = r_permutation[5u];//elsewhere
        }

        Vec vector =  mesh.GetDistributedVectorFactory()->CreateVec();
        DistributedVector dist_solution = mesh.GetDistributedVectorFactory()->CreateDistributedVector(vector);
        for (DistributedVector::Iterator index = dist_solution.Begin();
             index != dist_solution.End();
             ++index)
        {
            extended_bidomain_tissue.UpdateAdditionalCaches(index.Global, index.Local, 2.0);
        }
        extended_bidomain_tissue.ReplicateAdditionalCaches();

        //Ggap
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetGgapCacheReplicated()[probe_node_1],143.0);//within first cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetGgapCacheReplicated()[probe_node_2],9143.0);//within second cuboid
        TS_ASSERT_EQUALS(extended_bidomain_tissue.rGetGgapCacheReplicated()[probe_node_3],586.0);//elsewhere

        PetscTools::Destroy(vector);
    }

    void TestExtendedBidomainTissueParameters()
    {
        HeartConfig::Instance()->Reset();

        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2,2,2);

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;
        ExtracellularStimulusFactory extracellular_stimulus_factory;

        stimulated_cell_factory.SetMesh(&mesh);
        unstimulated_cell_factory.SetMesh(&mesh);
        extracellular_stimulus_factory.SetMesh(&mesh);

        ExtendedBidomainTissue<3>  extended_bidomain_tissue( &stimulated_cell_factory,  &unstimulated_cell_factory, &extracellular_stimulus_factory);

        //cover the set and get method for the flag about the extracellular stimulus
        TS_ASSERT_EQUALS(extended_bidomain_tissue.HasTheUserSuppliedExtracellularStimulus(), false);
        extended_bidomain_tissue.SetUserSuppliedExtracellularStimulus(true);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.HasTheUserSuppliedExtracellularStimulus(), true);

        //set some values
        extended_bidomain_tissue.SetAmFirstCell(1.0);
        extended_bidomain_tissue.SetAmSecondCell(2.0);
        extended_bidomain_tissue.SetAmGap(3.0);
        extended_bidomain_tissue.SetCmFirstCell(4.0);
        extended_bidomain_tissue.SetCmSecondCell(5.0);
        extended_bidomain_tissue.SetGGap(6.0);

        //and pick them up
        TS_ASSERT_EQUALS(extended_bidomain_tissue.GetAmFirstCell(), 1.0);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.GetAmSecondCell(), 2.0);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.GetAmGap(), 3.0);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.GetCmFirstCell(), 4.0);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.GetCmSecondCell(), 5.0);
        TS_ASSERT_EQUALS(extended_bidomain_tissue.GetGGap(), 6.0);
    }

    /**
     * Tests archiving of the tissue object.
     * It creates one, changes the default values of some member variables and saves.
     * Then it tries to load from the archive and checks that the member variables are with the right values.
     */
    void TestSaveAndLoadExtendedBidomainTissue()
    {
        HeartConfig::Instance()->Reset();
        // Archive settings
        FileFinder archive_dir("extended_tissue_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "extended_bidomain_tissue.arch";

        bool cache_replication_saved = false;
        double saved_printing_timestep = 2.0;
        double default_printing_timestep = HeartConfig::Instance()->GetPrintingTimeStep();

        c_matrix<double, 3, 3> intra_tensor_before_archiving;
        c_matrix<double, 3, 3> intra_tensor_second_cell_before_archiving;
        c_matrix<double, 3, 3> extra_tensor_before_archiving;

        //creation and save
        {
            // This call is required to set the appropriate conductivity media and to make sure that HeartConfig
            // knows the mesh filename despite we use our own mesh reader.
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_136_elements");

            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
            DistributedTetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            UnStimulatedCellFactory first_cell;
            StimulatedCellFactory second_cell;
            ExtracellularStimulusFactory extra_factory;
            first_cell.SetMesh(&mesh);
            second_cell.SetMesh(&mesh);
            extra_factory.SetMesh(&mesh);
            ExtendedBidomainTissue<3> extended_tissue( &first_cell, &second_cell , &extra_factory);

            //set a value different from default for the conductivities of the second cell
            extended_tissue.SetIntracellularConductivitiesSecondCell(Create_c_vector(25.0,26.0,27.0));
            //this is normally done by the problem class, but in this test we do it manually
            extended_tissue.CreateIntracellularConductivityTensorSecondCell();
            extended_tissue.SetCacheReplication(cache_replication_saved); // Not the default to check it is archived...

            //shuffle default values to check if they get archived properly
            extended_tissue.SetAmFirstCell(11.0);
            extended_tissue.SetAmSecondCell(22.0);
            extended_tissue.SetAmGap(33.0);
            extended_tissue.SetCmFirstCell(44.0);
            extended_tissue.SetCmSecondCell(55.0);
            extended_tissue.SetGGap(66.0);

            //again, away from default value to check for archiving
            extended_tissue.SetUserSuppliedExtracellularStimulus(true);

            //set some heterogeneities in Ggap
            std::vector<boost::shared_ptr<AbstractChasteRegion<3> > > heterogeneity_areas;
            std::vector<double> Ggap_values;
            ChastePoint<3> cornerA(-1, -1, 0);
            ChastePoint<3> cornerB(0.001, 0.001, 0.001);
            boost::shared_ptr<ChasteCuboid<3> > p_cuboid_1(new ChasteCuboid<3>(cornerA, cornerB));
            heterogeneity_areas.push_back(p_cuboid_1);
            //within the first area
            Ggap_values.push_back(143.0);
            extended_tissue.SetGgapHeterogeneities(heterogeneity_areas, Ggap_values);
            extended_tissue.CreateGGapConductivities();

            // Some checks to make sure HeartConfig is being saved and loaded by this too.
            HeartConfig::Instance()->SetPrintingTimeStep(saved_printing_timestep);
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);

            intra_tensor_before_archiving = extended_tissue.rGetIntracellularConductivityTensor(0);
            intra_tensor_second_cell_before_archiving = extended_tissue.rGetIntracellularConductivityTensorSecondCell(0);
            extra_tensor_before_archiving = extended_tissue.rGetExtracellularConductivityTensor(0);

            // Save
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<3>* const p_archive_bidomain_tissue = &extended_tissue;
            (*p_arch) << p_archive_bidomain_tissue;

            HeartConfig::Reset();
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), default_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep);
        }
        //load
        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<3>* p_abstract_tissue;
            (*p_arch) >> p_abstract_tissue;
            assert(p_abstract_tissue!=NULL);

            //dynamic cast so we are able to test specific variables of ExtendedBidomainTissue
            ExtendedBidomainTissue<3>* p_extended_tissue = dynamic_cast<ExtendedBidomainTissue<3>*>(p_abstract_tissue);
            assert(p_extended_tissue != NULL);

            const c_matrix<double, 3, 3>& intra_tensor_after_archiving = p_extended_tissue->rGetIntracellularConductivityTensor(0);
            const c_matrix<double, 3, 3>& intra_tensor_second_cell_after_archiving = p_extended_tissue->rGetIntracellularConductivityTensorSecondCell(0);
            const c_matrix<double, 3, 3>& extra_tensor_after_archiving = p_extended_tissue->rGetExtracellularConductivityTensor(0);

            //check before archiving = after archiving
            for (unsigned i=0; i<3; i++)
            {
                for (unsigned j=0; j<3; j++)
                {
                    TS_ASSERT_DELTA(intra_tensor_before_archiving(i,j), intra_tensor_after_archiving(i,j), 1e-9);
                    TS_ASSERT_DELTA(intra_tensor_second_cell_before_archiving(i,j), intra_tensor_second_cell_after_archiving(i,j), 1e-9);
                    TS_ASSERT_DELTA(extra_tensor_before_archiving(i,j), extra_tensor_after_archiving(i,j), 1e-9);
                }
            }

            //check that the member variable mIntracellularConductivitiesSecondCell was archived properly
            TS_ASSERT_EQUALS(p_extended_tissue->GetIntracellularConductivitiesSecondCell()(0),25.0);
            TS_ASSERT_EQUALS(p_extended_tissue->GetIntracellularConductivitiesSecondCell()(1),26.0);
            TS_ASSERT_EQUALS(p_extended_tissue->GetIntracellularConductivitiesSecondCell()(2),27.0);

            //check that we get the same values from the archive which are different from the default
            TS_ASSERT_EQUALS(p_extended_tissue->GetAmFirstCell(), 11.0);
            TS_ASSERT_EQUALS(p_extended_tissue->GetAmSecondCell(), 22.0);
            TS_ASSERT_EQUALS(p_extended_tissue->GetAmGap(), 33.0);
            TS_ASSERT_EQUALS(p_extended_tissue->GetCmFirstCell(), 44.0);
            TS_ASSERT_EQUALS(p_extended_tissue->GetCmSecondCell(), 55.0);
            TS_ASSERT_EQUALS(p_extended_tissue->GetGGap(), 66.0);

            // We shouldn't need to re-build the mesh, but we use it to check that the new tissue has the same mesh
            // Also, when testing in parallel, we use it to get the vector factory to loop over the nodes we own.
            // this is because  p_extended_tissue->pGetMesh()->GetDistributedVectorFactory() doesn't compile (discards qualifier stuff caused by use of const).
            TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
            DistributedTetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), p_extended_tissue->pGetMesh()->GetNumNodes());//note: this is allowed because GetNumNodes has const in the signature

            //check archiving of stimulus for first cell at some random times (it is unstimulated everywhere at all times)
            for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
            {
                if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(i))
                {
                    TS_ASSERT_EQUALS(p_extended_tissue->GetCardiacCell(i)->GetIntracellularStimulus(0.0), 0.0);
                    TS_ASSERT_EQUALS(p_extended_tissue->GetCardiacCell(i)->GetIntracellularStimulus(0.1), 0.0);
                    TS_ASSERT_EQUALS(p_extended_tissue->GetCardiacCell(i)->GetIntracellularStimulus(2.5), 0.0);
                    TS_ASSERT_EQUALS(p_extended_tissue->GetCardiacCell(i)->GetIntracellularStimulus(28.9), 0.0);
                }
            }

            //for second cell and other stuff, we probe nodes 0 and 1.
            unsigned node_0 = 0u;
            unsigned node_1 = 1u;
            //If the test is run in parallel, we need to work out the new indices
            const std::vector<unsigned>& r_permutation = mesh.rGetNodePermutation();
            if (!r_permutation.empty())
            {
                node_0 = r_permutation[0u];
                node_1 = r_permutation[1u];
            }

            //second cell is stimulated in the corner (node 0) from time 0 to 1. check it gets all this after loading
            if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(node_0))
            {
                TS_ASSERT_EQUALS(p_extended_tissue->GetCardiacSecondCell(node_0)->GetIntracellularStimulus(0.5), -105.0*1400);
                TS_ASSERT_EQUALS(p_extended_tissue->GetCardiacSecondCell(node_0)->GetIntracellularStimulus(2.5), 0.0);

                //find local index of (the new) node_0, it should be in the heterogeneity region
                unsigned ownership_range_low = mesh.GetDistributedVectorFactory()->GetLow();
                unsigned  local_index = node_0 - ownership_range_low;
                //std::cout<<local_index<<std::endl;
                TS_ASSERT_EQUALS(p_extended_tissue->rGetGapsDistributed()[local_index],143.0);//g_gap value inside heterogeneity region
            }

            //node 0 has extracellular stimulus (1 ms from 0.1)
            if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(node_0))
            {
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_0)->GetStimulus(0.0), 0);
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_0)->GetStimulus(0.5), -428000);
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_0)->GetStimulus(1.0), -428000);
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_0)->GetStimulus(1.15), 0);
            }

            //node 1 doesn't
            if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(node_1))
            {
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_1)->GetStimulus(0.0), 0);
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_1)->GetStimulus(0.5), 0);
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_1)->GetStimulus(1.0), 0);
                TS_ASSERT_EQUALS(p_extended_tissue->GetExtracellularStimulus(node_1)->GetStimulus(1.15), 0);

                //find local index of (the new) node_1, it should NOT be in the heterogeneity region
                unsigned ownership_range_low = mesh.GetDistributedVectorFactory()->GetLow();
                unsigned  local_index = node_1 - ownership_range_low;
                TS_ASSERT_EQUALS(p_extended_tissue->rGetGapsDistributed()[local_index],66.0);//standard g_gap value, outside heterogeneity region
            }

            //check the archiving of the flag (it would be false by default, but we set it to true before archiving)
            TS_ASSERT_EQUALS(p_extended_tissue->HasTheUserSuppliedExtracellularStimulus(),true);

            TS_ASSERT_EQUALS(cache_replication_saved, p_extended_tissue->GetDoCacheReplication());
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep); // Test we are testing something in case default changes

            delete p_extended_tissue;
        }
    }
};

#endif /*TESTEXTENDEDBIDOMAINTISSUE_HPP_*/
