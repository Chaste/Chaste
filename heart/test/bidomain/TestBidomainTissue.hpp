/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef TESTBIDOMAINTISSUE_HPP_
#define TESTBIDOMAINTISSUE_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <iostream>
#include <vector>

#include "ArchiveOpener.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"
#include "MonodomainTissue.hpp"
#include "BidomainTissue.hpp"
#include "OdeSolution.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "OrthotropicConductivityTensors.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include <petsc.h>


// cell factory for creating 2 cells with both intra and extracellular stimuli
template <unsigned PROBLEM_DIM=1>
class MyCardiacCellFactory : public AbstractCardiacCellFactory<PROBLEM_DIM>
{
private:
    boost::shared_ptr<AbstractStimulusFunction> mpStimulus;
public:

    MyCardiacCellFactory()
        : AbstractCardiacCellFactory<PROBLEM_DIM>(),
          mpStimulus(new SimpleStimulus(-80.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (node==0)
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }

    ~MyCardiacCellFactory(void)
    {
    }
};



class SimpleConductivityModifier : public AbstractConductivityModifier<1,1>
{
private:
    c_matrix<double,1,1> mTensor;

public:
    SimpleConductivityModifier()
        : AbstractConductivityModifier<1,1>()
    {
    }

    c_matrix<double,1,1>& rGetModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,1,1>& rOriginalConductivity)
    {
        mTensor(0,0) = (elementIndex+2.0)*rOriginalConductivity(0,0); //so conductivity on element 0 gets scaled by 2, and by 3 on element 1
        return mTensor;
    }
};


class TestBidomainTissue : public CxxTest::TestSuite
{
public:

    void TestBidomainTissueSolveCellSystems( void )
    {
        // This call is required to set the appropriate conductivity media and to make sure that
        // HeartConfig knows the mesh filename despite we use our own mesh reader.
        HeartConfig::Instance()->SetMeshFileName("linear_mesh", cp::media_type::NoFibreOrientation);

        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(1);

        double big_time_step = 0.5;
        MyCardiacCellFactory<> cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );
        BidomainTissue<1>     bidomain_tissue( &cell_factory );

        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;

        // initial condition;
        DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
        Vec monodomain_vec = p_factory->CreateVec();
        DistributedVector monodomain_voltage = p_factory->CreateDistributedVector(monodomain_vec);
        Vec bidomain_vec = p_factory->CreateVec(2);
        DistributedVector bidomain_ic = p_factory->CreateDistributedVector(bidomain_vec);
        DistributedVector::Stripe bidomain_voltage(bidomain_ic,0);

        for (DistributedVector::Iterator index=monodomain_voltage.Begin();
             index != monodomain_voltage.End();
             ++index)
        {
            monodomain_voltage[index] = initial_voltage;
            bidomain_voltage[index] = initial_voltage;
        }

        monodomain_voltage.Restore();
        bidomain_ic.Restore();

        monodomain_tissue.SolveCellSystems(monodomain_vec, 0, big_time_step);
        bidomain_tissue.SolveCellSystems(bidomain_vec, 0, big_time_step);


        // Check that both the monodomain and bidomain tissue have the same ionic cache
        for (unsigned node_index = mesh.GetDistributedVectorFactory()->GetLow();
             node_index < mesh.GetDistributedVectorFactory()->GetHigh();
             node_index++)
        {
            TS_ASSERT_EQUALS(monodomain_tissue.rGetIionicCacheReplicated()[node_index], bidomain_tissue.rGetIionicCacheReplicated()[node_index]);
        }

        // Check that the bidomain tissue has the right intracellular stimulus at node 0 and 1
        TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularStimulusCacheReplicated()[0], -80);
        TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularStimulusCacheReplicated()[1], 0);

        VecDestroy(monodomain_vec);
        VecDestroy(bidomain_vec);
    }


    void TestBidomainTissueWithHeterogeneousConductivitiesDistributed() throw (Exception)
    {
        HeartConfig::Instance()->Reset();

        // This call is required to set the appropriate conductivity media and to make sure that
        // HeartConfig knows the mesh filename despite we use our own mesh reader.
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_2mm_12_elements", cp::media_type::NoFibreOrientation);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Check that if we're in parallel no single process owns every element (to ensure that the conductivities
        // really are distributed).
        if (PetscTools::IsParallel())
        {
            TS_ASSERT_DIFFERS( mesh.GetNumElements(), mesh.GetNumLocalElements() );
        }

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

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML,3> cell_factory_for_het;
        cell_factory_for_het.SetMesh(&mesh);

        //CreateIntracellularConductivityTensor called in the constructor
        BidomainTissue<3> bidomain_tissue( &cell_factory_for_het );

        if (mesh.CalculateDesignatedOwnershipOfElement(0u))
        {
             TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(0u)(0,0),1.0);//within first cuboid
             TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(0u)(0,0),51.0);//within first cuboid
        }

        if (mesh.CalculateDesignatedOwnershipOfElement(4u))
        {
            TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(4u)(0,0),11.0);//within second cuboid
            TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(4u)(1,1),22.0);//within second cuboid
            TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(4u)(0,0),151.0);//within second cuboid
            TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(4u)(1,1),152.0);//within second cuboid
        }

        if (mesh.CalculateDesignatedOwnershipOfElement(8u))
        {
            TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(8u)(0,0),15.0);//elsewhere, e.g. element 8
            TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(8u)(0,0),65.0);//elsewhere, e.g. element 8
        }


    }

    void TestBidomainTissueWithHeterogeneousConductivitiesEllipsoid() throw (Exception)
    {
        HeartConfig::Instance()->Reset();

        // This call is required to set the appropriate conductivity media and to make sure that
        // HeartConfig knows the mesh filename despite we use our own mesh reader.
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/cube_2mm_12_elements", cp::media_type::NoFibreOrientation);

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<ChasteEllipsoid<3> > heterogeneity_area;
        std::vector< c_vector<double,3> > intra_conductivities;
        std::vector< c_vector<double,3> > extra_conductivities;

        //first small ellipsoid including element 0 centroid
        ChastePoint<3> centre_1(0.025, 0.075, 0.05);
        ChastePoint<3> radii_1(0.1, 0.1, 0.1);
        ChasteEllipsoid<3> ellipsoid_1(centre_1, radii_1);
        heterogeneity_area.push_back(ellipsoid_1);

        //second small ellipsoid including element 4 centroid
        ChastePoint<3> centre_2(0.175, 0.025, 0.05);
        ChastePoint<3> radii_2(0.1, 0.1, 0.1);
        ChasteEllipsoid<3> ellipsoid_2(centre_2, radii_2);

        heterogeneity_area.push_back(ellipsoid_2);

        //within the first area
        intra_conductivities.push_back( Create_c_vector(1.0, 2.0, 3.0) );
        extra_conductivities.push_back( Create_c_vector(51.0, 52.0, 53.0) );

        //within the second area
        intra_conductivities.push_back( Create_c_vector(11.0, 22.0, 33.0) );
        extra_conductivities.push_back( Create_c_vector(151.0, 152.0, 153.0) );

        HeartConfig::Instance()->SetConductivityHeterogeneitiesEllipsoid(heterogeneity_area, intra_conductivities, extra_conductivities);


        //elsewhere
        double isotropic_intra_conductivity=15.0;
        double isotropic_extra_conductivity=65.0;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(isotropic_intra_conductivity, isotropic_intra_conductivity, isotropic_intra_conductivity));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(isotropic_extra_conductivity, isotropic_extra_conductivity, isotropic_extra_conductivity));

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML,3> cell_factory_for_het;
        cell_factory_for_het.SetMesh(&mesh);

        //CreateIntracellularConductivityTensor called in the constructor
        BidomainTissue<3> bidomain_tissue( &cell_factory_for_het );

        TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(0u)(0,0),1.0);//within first ellipsoid
        TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(4u)(0,0),11.0);//within second ellipsoid
        TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(4u)(1,1),22.0);//within second ellipsoid
        TS_ASSERT_EQUALS(bidomain_tissue.rGetIntracellularConductivityTensor(8u)(0,0),15.0);//elsewhere, e.g. element 8

        TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(0u)(0,0),51.0);//within first ellipsoid
        TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(4u)(0,0),151.0);//within second ellipsoid
        TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(4u)(1,1),152.0);//within second ellipsoid
        TS_ASSERT_EQUALS(bidomain_tissue.rGetExtracellularConductivityTensor(8u)(0,0),65.0);//elsewhere, e.g. element 8
    }

    void TestGetConductivityAndConductivityModifier() throw(Exception)
    {
        HeartConfig::Instance()->Reset();
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, ie 2 node mesh

        MyCardiacCellFactory<1> cell_factory;
        cell_factory.SetMesh(&mesh);

        BidomainTissue<1> bidomain_tissue( &cell_factory );

        double orig_conductivity_0 = bidomain_tissue.rGetExtracellularConductivityTensor(0)(0,0);
        double orig_conductivity_1 = bidomain_tissue.rGetExtracellularConductivityTensor(1)(0,0);
        TS_ASSERT_DELTA(orig_conductivity_0, 7.0, 1e-9); // hard-coded using default
        TS_ASSERT_DELTA(orig_conductivity_1, 7.0, 1e-9); // hard-coded using default

        SimpleConductivityModifier modifier;
        bidomain_tissue.SetConductivityModifier(&modifier);

        TS_ASSERT_DELTA(bidomain_tissue.rGetExtracellularConductivityTensor(0)(0,0), 2*orig_conductivity_0, 1e-9);
        TS_ASSERT_DELTA(bidomain_tissue.rGetExtracellularConductivityTensor(1)(0,0), 3*orig_conductivity_1, 1e-9);
    }


    void TestSaveAndLoadCardiacPDE()
    {
        HeartConfig::Instance()->Reset();
        // Archive settings
        FileFinder archive_dir("tissue_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "bidomain_tissue.arch";

        bool cache_replication_saved = false;
        double saved_printing_timestep = 2.0;
        double default_printing_timestep = HeartConfig::Instance()->GetPrintingTimeStep();

        std::vector<cp::media_type> media_types;
        media_types.push_back(cp::media_type::Orthotropic);
        media_types.push_back(cp::media_type::Axisymmetric);
        media_types.push_back(cp::media_type::NoFibreOrientation);

        for (std::vector<cp::media_type>::iterator it = media_types.begin();
             it != media_types.end();
             ++it)
        {
            c_matrix<double, 3, 3> intra_tensor_before_archiving;
            c_matrix<double, 3, 3> extra_tensor_before_archiving;
            {
                // This call is required to set the appropriate conductivity media and to make sure that HeartConfig
                // knows the mesh filename despite we use our own mesh reader.
                HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_Single_tetrahedron_element", *it);

                TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
                TetrahedralMesh<3,3> mesh;
                mesh.ConstructFromMeshReader(mesh_reader);

                MyCardiacCellFactory<3> cell_factory;
                cell_factory.SetMesh(&mesh);

                BidomainTissue<3> bidomain_tissue( &cell_factory );
                bidomain_tissue.SetCacheReplication(cache_replication_saved); // Not the default to check it is archived...

                // Some checks to make sure HeartConfig is being saved and loaded by this too.
                HeartConfig::Instance()->SetPrintingTimeStep(saved_printing_timestep);
                TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);

                intra_tensor_before_archiving = bidomain_tissue.rGetIntracellularConductivityTensor(0);
                extra_tensor_before_archiving = bidomain_tissue.rGetExtracellularConductivityTensor(0);

                // Save
                ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
                boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

                AbstractCardiacTissue<3>* const p_archive_bidomain_tissue = &bidomain_tissue;
                (*p_arch) << p_archive_bidomain_tissue;

                HeartConfig::Reset();
                TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), default_printing_timestep, 1e-9);
                TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep);
            }

            {
                ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
                boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

                AbstractCardiacTissue<3>* p_bidomain_tissue;
                (*p_arch) >> p_bidomain_tissue;

                assert(p_bidomain_tissue!=NULL);

                const c_matrix<double, 3, 3>& intra_tensor_after_archiving = p_bidomain_tissue->rGetIntracellularConductivityTensor(0);
                const c_matrix<double, 3, 3>& extra_tensor_after_archiving = dynamic_cast<BidomainTissue<3>*>(p_bidomain_tissue)->rGetExtracellularConductivityTensor(0); //Naughty Gary using dynamic cast, but only for testing...

                for(unsigned i=0; i<3; i++)
                {
                    for(unsigned j=0; j<3; j++)
                    {
                        TS_ASSERT_DELTA(intra_tensor_before_archiving(i,j), intra_tensor_after_archiving(i,j), 1e-9);
                        TS_ASSERT_DELTA(extra_tensor_before_archiving(i,j), extra_tensor_after_archiving(i,j), 1e-9);
                    }
                }

                TS_ASSERT_EQUALS(cache_replication_saved, p_bidomain_tissue->GetDoCacheReplication());
                TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);
                TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep); // Test we are testing something in case default changes

                delete p_bidomain_tissue;
            }
        }
    }
};

#endif /*TESTBIDOMAINTISSUE_HPP_*/
