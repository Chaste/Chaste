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


#ifndef _TESTMONODOMAINTISSUE_HPP_
#define _TESTMONODOMAINTISSUE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "CardiacSimulationArchiver.hpp"

#include <vector>

#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"
#include "MonodomainTissue.hpp"
#include "OdeSolution.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "DistributedVector.hpp"
#include "PetscTools.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "UblasCustomFunctions.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "ArchiveOpener.hpp"
#include "DiFrancescoNoble1985.hpp"
#include "MonodomainProblem.hpp"

#include "PetscSetupAndFinalize.hpp"

class MyCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:

    MyCardiacCellFactory()
        : AbstractCardiacCellFactory<1>(),
          mpStimulus(new SimpleStimulus(-80.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        unsigned node_index = pNode->GetIndex();
        if (node_index==0)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    boost::shared_ptr<SimpleStimulus> GetStimulus()
    {
        return mpStimulus;
    }
};

class PurkinjeCellFactory : public AbstractPurkinjeCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    PurkinjeCellFactory()
        : AbstractPurkinjeCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-6000.0, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        ChastePoint<2> location = pNode->GetPoint();

        if (fabs(location[0])<1e-6)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    AbstractCardiacCell* CreatePurkinjeCellForTissueNode(Node<2>* pNode,
                                                         AbstractCardiacCellInterface* pCardiacCell)
    {
        return new CellDiFrancescoNoble1985FromCellML(mpSolver, mpZeroStimulus);
    }
};

class TestMonodomainTissue : public CxxTest::TestSuite
{
public:
    void TestMonodomainTissueBasic()
    {
        if (PetscTools::GetNumProcs() > 2u)
        {
            // There are only 2 nodes in this simulation
            TS_TRACE("This test is not suitable for more than 2 processes.");
            return;
        }
        HeartConfig::Instance()->Reset();
        unsigned num_nodes=2;
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, i.e. A mesh with 2 nodes
        assert(mesh.GetNumNodes()==num_nodes);

        double start_time = 0;
        double big_time_step = 0.5;

        boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        // Coverage
        TS_ASSERT_THROWS_THIS(cell_factory.FillInCellularTransmuralAreas(),
              "To get here you have probably asked for Epi/Mid/Endo CellularHeterogeneities in your HeartConfig "
              "options or configuration .xml file, to use this you will need to provide a method"
              " `FillInCellularTransmuralAreas()` in your cell factory to override this one.");

        // Stimulus function to use at node 0. Node 1 is not stimulated.
        boost::shared_ptr<SimpleStimulus> p_stimulus = cell_factory.GetStimulus();
        boost::shared_ptr<ZeroStimulus> p_zero_stim(new ZeroStimulus);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );

        // check the purkinje cells vector is empty
        TS_ASSERT(!monodomain_tissue.HasPurkinje());
        TS_ASSERT_THROWS_ANYTHING(monodomain_tissue.rGetPurkinjeCellsDistributed().size());
        TS_ASSERT_THROWS_ANYTHING(monodomain_tissue.rGetPurkinjeIionicCacheReplicated().GetSize());

        // voltage that gets passed in solving ode
        double initial_voltage = -83.853;

        // initial condition;
        Vec voltage = PetscTools::CreateAndSetVec(num_nodes, initial_voltage);

        // Solve 1 (PDE) timestep using MonodomainTissue
        monodomain_tissue.SolveCellSystems(voltage, start_time, start_time+big_time_step);

        // Check results by solving ODE systems directly
        // Check node 0
        double value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[0];
        CellLuoRudy1991FromCellML ode_system_stimulated(p_solver, p_stimulus);
        ode_system_stimulated.ComputeExceptVoltage(start_time, start_time + big_time_step);
        double value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 0.000001);

        // shouldn't be different when called again as reset not yet been called
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[0];
        TS_ASSERT_DELTA(value_tissue, value_ode, 0.000001);

        // Check node 1
        CellLuoRudy1991FromCellML ode_system_not_stim(p_solver, p_zero_stim);
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[1];
        ode_system_not_stim.ComputeExceptVoltage(start_time, start_time + big_time_step);
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 0.000001);

        // Reset the voltage vector from ODE systems
        DistributedVector dist_voltage = mesh.GetDistributedVectorFactory()->CreateDistributedVector(voltage);
        for (DistributedVector::Iterator index = dist_voltage.Begin();
             index != dist_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                dist_voltage[index] = ode_system_stimulated.rGetStateVariables()[0];
            }
            if (index.Global==1)
            {
                dist_voltage[index] = ode_system_not_stim.rGetStateVariables()[0];
            }
        }
        dist_voltage.Restore();

        // Use MonodomainTissue to solve a second (PDE) time step
        monodomain_tissue.SolveCellSystems(voltage, start_time, start_time+big_time_step);
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[0];

        // Check node 0 by solving ODE system directly
        ode_system_stimulated.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_ode = ode_system_stimulated.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 1e-10);

        // Check node 1 by solving ODE system directly
        ode_system_not_stim.ComputeExceptVoltage( start_time + big_time_step, start_time + 2*big_time_step );
        value_tissue = monodomain_tissue.rGetIionicCacheReplicated()[1];
        value_ode = ode_system_not_stim.GetIIonic();
        TS_ASSERT_DELTA(value_tissue, value_ode, 1e-10);

        TS_ASSERT_THROWS_THIS(monodomain_tissue.rGetExtracellularConductivityTensor(0),
                              "Monodomain tissues do not have extracellular conductivity tensors.");
        TS_ASSERT_THROWS_THIS(monodomain_tissue.rGetIntracellularConductivityTensor(1),
                              "Conductivity tensor requested for element with global_index=1, but there are only 1 elements in the mesh.");

        PetscTools::Destroy(voltage);
    }

    void TestMonodomainTissueGetCardiacCell()
    {
        if (PetscTools::GetNumProcs() > 2u)
        {
            // There are only 2 nodes in this simulation
            TS_TRACE("This test is not suitable for more than 2 processes.");
            return;
        }
        HeartConfig::Instance()->Reset();
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, i.e. 2 nodes in the mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0))
        {
            AbstractCardiacCellInterface* cell = monodomain_tissue.GetCardiacCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),-80,1e-10);
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1))
        {
            AbstractCardiacCellInterface* cell = monodomain_tissue.GetCardiacCell(1);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }
    }

    void TestSolveCellSystemsInclUpdateVoltage()
    {
        if (PetscTools::GetNumProcs() > 2u)
        {
            // There are only 2 nodes in this simulation
            TS_TRACE("This test is not suitable for more than 2 processes.");
            return;
        }
        HeartConfig::Instance()->Reset();
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, i.e. a 2 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory );

        Vec voltage = PetscTools::CreateAndSetVec(2, -81.4354); // something that isn't resting potential
        monodomain_tissue.SolveCellSystems(voltage, 0, 1, false); // solve for 1ms without updating the voltage

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0))
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -81.4354, 1e-3);
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1))
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -81.4354, 1e-3);
        }

        Vec voltage2 = PetscTools::CreateAndSetVec(2, -75);
        monodomain_tissue.SolveCellSystems(voltage2, 1, 2, true); // solve another ms, using this new voltage, but now updating the voltage too
        ReplicatableVector voltage2_repl(voltage2); // should have changed following solve

        // check the new voltage in the cell is NEAR -75 (otherwise the passed in voltage wasn't used, but
        // NOT EXACTLY -75, ie that the voltage was solved for.
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0))
        {
            // check has been updated
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75);
            // check near -75
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75, 2.0); // within 2mV
            // check the passed in voltage was updated
            TS_ASSERT_DELTA(voltage2_repl[0], monodomain_tissue.GetCardiacCell(0)->GetVoltage(), 1e-10);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1))
        {
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[1], monodomain_tissue.GetCardiacCell(1)->GetVoltage(), 1e-10);
        }

        PetscTools::Destroy(voltage);
        PetscTools::Destroy(voltage2);
    }

    void TestNodeExchange()
    {
        HeartConfig::Instance()->Reset();

        HeartConfig::Instance()->Reset();
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(0.1, 1.0); // [0,1] with h=0.1, ie 11 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory, true );

        if (PetscTools::GetNumProcs() == 1)
        {
            TS_ASSERT_EQUALS(mesh.GetNumHaloNodes(), 0u);
        }
        else
        {
            if (PetscTools::AmMaster() || PetscTools::AmTopMost())
            {
                TS_ASSERT_EQUALS(mesh.GetNumHaloNodes(), 1u);
            }
            else
            {
                TS_ASSERT_EQUALS(mesh.GetNumHaloNodes(), 2u);
            }
        }

        for (DistributedTetrahedralMesh<1,1>::HaloNodeIterator it=mesh.GetHaloNodeIteratorBegin();
             it != mesh.GetHaloNodeIteratorEnd();
             ++it)
        {
            AbstractCardiacCellInterface* cell = monodomain_tissue.GetCardiacCellOrHaloCell( (*it)->GetIndex() );
            TS_ASSERT_DELTA(cell->GetStimulus(0.001),0,1e-10);
        }

        if (PetscTools::AmMaster())
        {
            // Master owns node 0
            AbstractCardiacCellInterface* cell = monodomain_tissue.GetCardiacCellOrHaloCell(0);
            TS_ASSERT_DELTA(cell->GetStimulus(0.001), -80.0, 1e-10);
        }
        else
        {
            // Zero is not halo owned by any process (unless we have a lot of them).
            TS_ASSERT_THROWS_CONTAINS(monodomain_tissue.GetCardiacCellOrHaloCell(0),
                                      "Requested node/halo 0 does not belong to processor ");
        }
    }

    void TestSolveCellSystemsInclUpdateVoltageWithNodeExchange()
    {
        if (PetscTools::GetNumProcs() > 2u)
        {
            // There are only 2 nodes in this simulation
            TS_TRACE("This test is not suitable for more than 2 processes.");
            return;
        }
        HeartConfig::Instance()->Reset();
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0); // [0,1] with h=1.0, i.e. 2 node mesh

        MyCardiacCellFactory cell_factory;
        cell_factory.SetMesh(&mesh);

        MonodomainTissue<1> monodomain_tissue( &cell_factory, true );

        Vec voltage = PetscTools::CreateAndSetVec(2, -81.4354); // something that isn't resting potential
        monodomain_tissue.SolveCellSystems(voltage, 0, 1, false); // solve for 1ms without updating the voltage

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0)) // Is node 0 locally owned?
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -81.4354, 1e-3);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), -81.4354, 1e-3);
        }

        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1)) // Is node 1 locally owned?
        {
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), -81.4354, 1e-3);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -81.4354, 1e-3);
        }

        Vec voltage2 = PetscTools::CreateAndSetVec(2, -75);
        monodomain_tissue.SolveCellSystems(voltage2, 1, 2, true); // solve another ms, using this new voltage, but now updating the voltage too
        ReplicatableVector voltage2_repl(voltage2); // should have changed following solve

        // check the new voltage in the cell is NEAR -75 (otherwise the passed in voltage wasn't used, but
        // NOT EXACTLY -75, ie that the voltage was solved for.
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0)) // Is node 0 locally owned?
        {
            // check has been updated
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75);
            // check near -75
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(0)->GetVoltage(), -75, 2.0); // within 2mV
            // check the passed in voltage was updated
            TS_ASSERT_DELTA(voltage2_repl[0], monodomain_tissue.GetCardiacCell(0)->GetVoltage(), 1e-10);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1)) // Is node 1 locally owned?
        {
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCell(1)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[1], monodomain_tissue.GetCardiacCell(1)->GetVoltage(), 1e-10);
        }

        // now check the new voltages have been communicated
        // check the new voltage in the cell is NEAR -75 (otherwise the passed in voltage wasn't used, but
        // NOT EXACTLY -75, ie that the voltage was solved for.
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0)) // Is node 0 locally owned?
        {
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[1], monodomain_tissue.GetCardiacCellOrHaloCell(1)->GetVoltage(), 1e-10);
        }
        if (mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(1)) // Is node 1 locally owned?
        {
            TS_ASSERT_DIFFERS(monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), -75);
            TS_ASSERT_DELTA(monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), -75, 2.0); // within 2mV
            TS_ASSERT_DELTA(voltage2_repl[0], monodomain_tissue.GetCardiacCellOrHaloCell(0)->GetVoltage(), 1e-10);
        }

        PetscTools::Destroy(voltage);
        PetscTools::Destroy(voltage2);
    }

    void TestSaveAndLoadCardiacTissue()
    {
        HeartConfig::Instance()->Reset();
        // Archive settings
        FileFinder archive_dir("monodomain_tissue_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "monodomain_tissue.arch";

        bool cache_replication_saved = false;
        double saved_printing_timestep = 2.0;
        double default_printing_timestep = HeartConfig::Instance()->GetPrintingTimeStep();

        // Info about the first cell on this process (if any)
        bool has_cell = false;
        unsigned cell_v_index = (unsigned)(-1);
        double cell_v = DBL_MAX;

        c_matrix<double, 1, 1> tensor_before_archiving;
        {
            TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            MyCardiacCellFactory cell_factory;
            cell_factory.SetMesh(&mesh);

            MonodomainTissue<1> monodomain_tissue( &cell_factory );
            monodomain_tissue.SetCacheReplication(cache_replication_saved); // Not the default to check it is archived...

            tensor_before_archiving = monodomain_tissue.rGetIntracellularConductivityTensor(1);

            // Get some info about the first cell on this process (if any)
            const std::vector<AbstractCardiacCellInterface*>& r_cells = monodomain_tissue.rGetCellsDistributed();
            has_cell = !r_cells.empty();
            if (has_cell)
            {
                cell_v_index = r_cells[0]->GetVoltageIndex();
                cell_v = r_cells[0]->GetVoltage();
            }

            // Some checks to make sure HeartConfig is being saved and loaded by this too.
            HeartConfig::Instance()->SetPrintingTimeStep(saved_printing_timestep);
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);

            // Save
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<1>* const p_archive_monodomain_tissue = &monodomain_tissue;
            (*p_arch) << p_archive_monodomain_tissue;

            HeartConfig::Reset();
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), default_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep);
        }

        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<1>* p_monodomain_tissue;
            (*p_arch) >> p_monodomain_tissue;

            // Test rGetIntracellularConductivityTensor
            const c_matrix<double, 1, 1>& tensor_after_archiving = p_monodomain_tissue->rGetIntracellularConductivityTensor(1);
            TS_ASSERT_DELTA(tensor_before_archiving(0,0), tensor_after_archiving(0,0), 1e-9);

            TS_ASSERT_EQUALS(cache_replication_saved, p_monodomain_tissue->GetDoCacheReplication());
            TS_ASSERT_DELTA(HeartConfig::Instance()->GetPrintingTimeStep(), saved_printing_timestep, 1e-9);
            TS_ASSERT_DIFFERS(saved_printing_timestep, default_printing_timestep); // Test we are testing something in case default changes

            // Test cardiac cells have also been archived
            const std::vector<AbstractCardiacCellInterface*>& r_cells = p_monodomain_tissue->rGetCellsDistributed();
            TS_ASSERT_EQUALS(has_cell, !r_cells.empty());
            if (has_cell)
            {
                TS_ASSERT_EQUALS(cell_v_index, r_cells[0]->GetVoltageIndex());
                TS_ASSERT_EQUALS(cell_v, r_cells[0]->GetVoltage());
            }

            delete p_monodomain_tissue;
        }
    }

    void TestMonodomainTissueUsingPurkinjeCellFactory()
    {
        HeartConfig::Instance()->Reset();

        TrianglesMeshReader<2,2> reader("mesh/test/data/mixed_dimension_meshes/2D_0_to_1mm_200_elements");
        MixedDimensionMesh<2,2> mixed_mesh;
        mixed_mesh.ConstructFromMeshReader(reader);

        PurkinjeCellFactory cell_factory;
        cell_factory.SetMesh(&mixed_mesh);

        MonodomainTissue<2> tissue( &cell_factory );

        TS_ASSERT(tissue.HasPurkinje());
        TS_ASSERT_EQUALS(tissue.rGetPurkinjeCellsDistributed().size(), tissue.rGetCellsDistributed().size());
        TS_ASSERT_EQUALS(tissue.rGetPurkinjeIionicCacheReplicated().GetSize(),
                         tissue.rGetIionicCacheReplicated().GetSize());

        for (AbstractTetrahedralMesh<2,2>::NodeIterator current_node = mixed_mesh.GetNodeIteratorBegin();
             current_node != mixed_mesh.GetNodeIteratorEnd();
             ++current_node)
        {
            unsigned global_index = current_node->GetIndex();
            AbstractCardiacCellInterface* p_purkinje_cell = tissue.GetPurkinjeCell(global_index);
            double y = current_node->rGetLocation()[1];

            // Cable nodes are on y=0.05 (we don't test by index because indices may be permuted in parallel).
            if (fabs(y-0.05) < 1e-8)
            {
                TS_ASSERT(dynamic_cast<CellDiFrancescoNoble1985FromCellML*>(p_purkinje_cell) != NULL);
            }
            else
            {
                TS_ASSERT(dynamic_cast<FakeBathCell*>(p_purkinje_cell) != NULL);
            }

            TS_ASSERT_EQUALS(tissue.rGetPurkinjeCellsDistributed()[global_index-mixed_mesh.GetDistributedVectorFactory()->GetLow()],
                             p_purkinje_cell);
        }

        // Test archiving too
        FileFinder archive_dir("monodomain_tissue_purkinje_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "monodomain_tissue_purkinje.arch";

        {
            // Save to archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Make sure at least one Purkinje cell has a non-initial-condition state variable to compare
            if (mixed_mesh.GetDistributedVectorFactory()->IsGlobalIndexLocal(0u))
            {
                tissue.GetPurkinjeCell(0u)->SetVoltage(1234.5);
            }

            AbstractCardiacTissue<2>* const p_archive_tissue = &tissue;
            (*p_arch) << p_archive_tissue;
        }

        {
            // Load from archive and compare
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacTissue<2>* p_tissue;
            (*p_arch) >> p_tissue;

            TS_ASSERT(p_tissue->HasPurkinje());
            TS_ASSERT_EQUALS(p_tissue->rGetPurkinjeCellsDistributed().size(), tissue.rGetPurkinjeCellsDistributed().size());
            TS_ASSERT_EQUALS(p_tissue->rGetPurkinjeIionicCacheReplicated().GetSize(), tissue.rGetPurkinjeIionicCacheReplicated().GetSize());

            for (AbstractTetrahedralMesh<2,2>::NodeIterator current_node = p_tissue->mpMesh->GetNodeIteratorBegin();
                 current_node != p_tissue->mpMesh->GetNodeIteratorEnd();
                 ++current_node)
            {
                unsigned global_index = current_node->GetIndex();
                AbstractCardiacCellInterface* p_purkinje_cell = p_tissue->GetPurkinjeCell(global_index);
                double y = current_node->rGetLocation()[1];

                // cable nodes are on y=0.05 (we don't test by index because indices may be permuted in parallel).
                if (fabs(y-0.05) < 1e-8)
                {
                    TS_ASSERT(dynamic_cast<CellDiFrancescoNoble1985FromCellML*>(p_purkinje_cell) != NULL);
                }
                else
                {
                    TS_ASSERT(dynamic_cast<FakeBathCell*>(p_purkinje_cell) != NULL);
                }
                TS_ASSERT_EQUALS(p_purkinje_cell->GetVoltage(),
                                 tissue.GetPurkinjeCell(global_index)->GetVoltage());

                TS_ASSERT_EQUALS(p_tissue->rGetPurkinjeCellsDistributed()[global_index-p_tissue->mpMesh->GetDistributedVectorFactory()->GetLow()],
                                 p_purkinje_cell);
            }

            delete p_tissue;
        }

        const std::string migration_archive_dir("TestMonodomainTissue/purkinje_migration_archive");
        {
            // Save via MonodomainProblem so we can migrate
            // Note: from Chaste release 3.1 onward we no longer support Boost 1.33.
            // The earliest version of Boost supported in 1.34

            // Run the test with b=_hostconfig,boost=1-34_5 to save
            /*
               scons b=_hostconfig,boost=1-34_5 ts=heart/test/monodomain/TestMonodomainTissue.hpp
             *
             */
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mixed_mesh);
            monodomain_problem.Initialise();
            TS_ASSERT(monodomain_problem.GetMonodomainTissue()->HasPurkinje());
            CardiacSimulationArchiver<MonodomainProblem<2> >::Save(monodomain_problem, migration_archive_dir);
        }
        TS_ASSERT_EQUALS(tissue.rGetPurkinjeIionicCacheReplicated().GetSize(), 121u);
    }

    //Failure of this test may mean that the archive from the previous needs to be regenerated
    void TestArchiveMigration()
    { // Load from 5-process archive and compare
        FileFinder saved_archive_dir("heart/test/data/checkpoints/purkinje_migration_archive",
                                     RelativeTo::ChasteSourceRoot);
        MonodomainProblem<2>* p_problem = CardiacSimulationArchiver<MonodomainProblem<2> >::Load(saved_archive_dir);
        MonodomainTissue<2>* p_tissue = p_problem->GetMonodomainTissue();

        TS_ASSERT(p_tissue->HasPurkinje());
        TS_ASSERT_EQUALS(p_tissue->rGetPurkinjeIionicCacheReplicated().GetSize(), 121u);

        delete p_problem;
    }
};

#endif //_TESTMONODOMAINTISSUE_HPP_
