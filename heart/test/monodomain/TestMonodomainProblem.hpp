/*

Copyright (C) University of Oxford, 2005-2012

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


#ifndef _TESTMONODOMAINPROBLEM_HPP_
#define _TESTMONODOMAINPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "FaberRudy2000.hpp"
#include "Hdf5DataReader.hpp"
#include "ReplicatableVector.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ArchiveOpener.hpp"
#include "CmguiMeshWriter.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "NumericFileComparison.hpp"
#include "VtkMeshReader.hpp"

#include <sys/stat.h> // For chmod()

/*
 *  This cell factory introduces a stimulus in the very centre of mesh/test/data/2D_0_to_1mm_400_elements.
 *  This is node 60 at (0.5, 0.5)
 */
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    unsigned mFoundMiddlePoint;

public:
    PointStimulus2dCellFactory()
        : AbstractCardiacCellFactory<2>(),
          mpStimulus(new SimpleStimulus(-6000.0, 0.5)),
          mFoundMiddlePoint(0)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        ChastePoint<2> location = GetMesh()->GetNode(node)->GetPoint();

        if (fabs(location[0]-0.05)<1e-6 && fabs(location[1]-0.05)<1e-6)
        {
            mFoundMiddlePoint++;
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }

    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        unsigned found_middle_point_reduced;
        MPI_Allreduce(&mFoundMiddlePoint, &found_middle_point_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        assert(found_middle_point_reduced == 1); // Only 1 cell should be stimulated
    }
};



class TestMonodomainProblem : public CxxTest::TestSuite
{
private:
    std::vector<double> mVoltageReplicated1d2ms;///<Used to test differences between tests

public:
    void setUp()
    {
        HeartConfig::Reset();
    }

    // This is a test on a very small mesh (where there may be more processes than there are nodes)
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblemSimplestMesh1D() throw(Exception)
    {
        DistributedTetrahedralMesh<1,1> mesh;
        //TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_1_element");
        mesh.ConstructLinearMesh(1);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        if (!PetscTools::IsSequential())
        {
            if (PetscTools::GetMyRank() < 2)
            {
                TS_ASSERT_EQUALS(mesh.GetNumLocalNodes(),1u);
                TS_ASSERT_EQUALS(mesh.GetNumLocalElements(),1u);
            }
            else
            {
                TS_ASSERT_EQUALS(mesh.GetNumLocalNodes(),0u);
                TS_ASSERT_EQUALS(mesh.GetNumLocalElements(),0u);
            }
        }
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1dSimplest");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
        //HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_1_element");
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );
        monodomain_problem.SetMesh(&mesh);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(10.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());

        double atol = 5e-3;

        TS_ASSERT_DELTA(voltage_replicated[0], -52.1698, atol);
        TS_ASSERT_DELTA(voltage_replicated[1], -83.8381, atol);
    }

    // Solve on a 1D string of cells, 1mm long with a space step of 0.1mm.
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblem1D() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
        // Set extra variable (to cover extension of HDF5 with named variables) in a later test
        std::vector<std::string> output_variables;
        output_variables.push_back("cytosolic_calcium_concentration");
        HeartConfig::Instance()->SetOutputVariables( output_variables );

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());

        double atol = 5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

        for (unsigned index=0; index<voltage_replicated.GetSize(); index++)
        {
            mVoltageReplicated1d2ms.push_back(voltage_replicated[index]);
        }
        //How close is our "standard" answer?
        atol = 5e-3;
        TS_ASSERT_DELTA(mVoltageReplicated1d2ms[1], 20.7710232, atol);
        TS_ASSERT_DELTA(mVoltageReplicated1d2ms[3], 21.5319692, atol);
        TS_ASSERT_DELTA(mVoltageReplicated1d2ms[5], 22.9280817, atol);
        TS_ASSERT_DELTA(mVoltageReplicated1d2ms[7], 24.0611303, atol);
        TS_ASSERT_DELTA(mVoltageReplicated1d2ms[9], -0.770330519, atol);
        TS_ASSERT_DELTA(mVoltageReplicated1d2ms[10], -19.2234919, atol);

        // coverage
        monodomain_problem.GetTissue();

        // check a progress report exists
        TS_ASSERT_EQUALS(system(("ls " + OutputFileHandler::GetChasteTestOutputDirectory() + "MonoProblem1d/").c_str()), 0);
    }


    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblem1DWithRelativeTolerance() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-9);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1dRelTol");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());

        /// \todo: If we request "relative" tolerance we shouldn't testing in an "absolute" manner
        double atol = 2e-5;
        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

        for (unsigned index=0; index<voltage_replicated.GetSize(); index++)
        {
            TS_ASSERT_DELTA(voltage_replicated[index], mVoltageReplicated1d2ms[index], 5e-3);
        }

    }

    // Same as TestMonodomainProblem1D, except the 1D mesh is embedded in 3D space.
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblem1Din3D() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_in_3D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1din3d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1din3d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1, 3> cell_factory;
        MonodomainProblem<1,3> monodomain_problem( &cell_factory );
        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol = 5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

        for (unsigned index=0; index<voltage_replicated.GetSize(); index++)
        {
            TS_ASSERT_DELTA(voltage_replicated[index], mVoltageReplicated1d2ms[index],  1e-12);
        }
        // cover get pde
        monodomain_problem.GetTissue();

        // check a progress report exists
        TS_ASSERT_EQUALS(system(("ls " + OutputFileHandler::GetChasteTestOutputDirectory() + "MonoProblem1din3d/").c_str()), 0);
    }

    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblem1DWithAbsoluteTolerance() throw (Exception)
    {
        double atol = 1e-4;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetAbsoluteTolerance(), 2e-4);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(atol/20.0);
        ///\todo this is dependent on the number of processes used
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetAbsoluteTolerance(), 5e-6);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1dAbsTol");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);
        for (unsigned index=0; index<voltage_replicated.GetSize(); index++)
        {
            TS_ASSERT_DELTA(voltage_replicated[index], mVoltageReplicated1d2ms[index],  5e-3);
        }

    }

    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating the left
    // edge.
    // Should behave like the 1D case, extrapolated.
    // See also TestMonodomainSlab.hpp (nightly test) for the 3D version.
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblem2DWithEdgeStimulus() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem2dWithEdgeStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2dWithEdgeStimulus");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;

        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        // Coverage
        TS_ASSERT(!monodomain_problem.GetHasBath());

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);

        //Since we are going to compare voltages that may be owned by
        //various processes it makes sense to replicate the data.
        Vec voltage=monodomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);
        /*
         * Test the top right node against the right one in the 1D case,
         * comparing voltage, and then test all the nodes on the right hand
         * side of the square against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=0.0;

        need_initialisation = true;

        // Test the RHS of the mesh
        AbstractTetrahedralMesh<2,2>& r_mesh=monodomain_problem.rGetMesh();
        for (AbstractTetrahedralMesh<2,2>::NodeIterator iter=r_mesh.GetNodeIteratorBegin();
             iter != r_mesh.GetNodeIteratorEnd();
             ++iter)
        {
            unsigned i=(*iter).GetIndex();
            if ((*iter).GetPoint()[0] == 0.1)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.1cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[i];
                    need_initialisation = false;
                }
                else
                {
                    // Tests the final voltages for all the RHS edge nodes
                    // are close to each other.
                    // This works as we are using the 'criss-cross' mesh,
                    // the voltages would vary more with a mesh with all the
                    // triangles aligned in the same direction.
                    TS_ASSERT_DELTA(voltage_replicated[i], probe_voltage, 7e-4);
                }


                // Check against 1d case - THIS TEST HAS BEEN REMOVED AS THE MESH
                // IS FINER THAN THE 1D MESH SO WE DONT EXPECT THE RESULTS TO BE THE SAME
                // TS_ASSERT_DELTA(p_voltage_array[i], -35.1363, 35*0.1);

                // test the RHS edge voltages
                // hardcoded result that looks accurate - this is a test to see
                // that nothing has changed
                // assumes endtime = 2ms
                TS_ASSERT_DELTA(voltage_replicated[i], -59.6488, 5e-4);
            }
        }

    }


    // Solve on a 2D 1mm by 1mm mesh (space step = 0.1mm), stimulating in the
    // very centre of the mesh.
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainProblem2DWithPointStimulusInTheVeryCentreOfTheMesh() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.3); //ms - needs to be 1.3 ms to pass test
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem2dWithPointStimulus");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2dWithPointStimulus");

        // To time the solve
        time_t start,end;

        time (&start);

        PointStimulus2dCellFactory cell_factory;

        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        monodomain_problem.Solve();

        // To time the solve
        time (&end);
        //double dif;
        //dif = difftime (end,start);
        //printf ("\nSolve took %.2lf seconds. \n", dif );

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars(monodomain_problem);

        /*
         * Test that corners are 'equal', and centres of sides.
         * Irregularities in which way the triangles are oriented make
         * this rather difficult, especially since the edges are sampled
         * during the upstroke.
         */
        DistributedVector voltage = monodomain_problem.GetSolutionDistributedVector();

        // corners -> 0, 10, 110, 120
        // hardcoded result to check against
        // assumes endtime = 1.3
        unsigned corners_checked=0;
        for (DistributedVector::Iterator node_index = voltage.Begin();
             node_index!= voltage.End();
             ++node_index)
        {
            ChastePoint<2> location = monodomain_problem.rGetMesh().GetNode(node_index.Global)->GetPoint();

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 0
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }

            if (fabs(location[0]-0.1)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 10
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 110
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.0)<1e-6) // Corner 120
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], -34.3481, 1e-3);
                corners_checked++;
            }
        }

        unsigned corners_checked_reduced;
        MPI_Allreduce(&corners_checked, &corners_checked_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT(corners_checked_reduced==4);


        // centre of edges -> 5, 55, 65, 115
        // hardcoded result to check against
        // assumes endtime = 1.3
        unsigned edges_checked=0;
        for (DistributedVector::Iterator node_index = voltage.Begin();
             node_index!= voltage.End();
             ++node_index)
        {
            ChastePoint<2> location = monodomain_problem.rGetMesh().GetNode(node_index.Global)->GetPoint();

            if (fabs(location[0]-0.05)<1e-6 && fabs(location[1]-0.0)<1e-6) // Node 5
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }

            if (fabs(location[0]-0.0)<1e-6 && fabs(location[1]-0.05)<1e-6) // Node 55
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }

            if (fabs(location[0]-0.1)<1e-6 && fabs(location[1]-0.05)<1e-6) // Node 65
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }

            if (fabs(location[0]-0.05)<1e-6 && fabs(location[1]-0.1)<1e-6) // Node 115
            {
                TS_ASSERT_DELTA(voltage[node_index.Global], 34.6692, 1e-3);
                edges_checked++;
            }
        }

        unsigned edges_checked_reduced;
        MPI_Allreduce(&edges_checked, &edges_checked_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        TS_ASSERT(edges_checked_reduced==4);
    }


    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestMonodomainProblemPrintsOnlyAtRequestedTimes() throw(Exception)
    {
        HeartConfig::Instance()->SetPrintingTimeStep(0.1);
        HeartConfig::Instance()->SetSimulationDuration(0.3); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1dRequestedTimes");
        HeartConfig::Instance()->SetOutputFilenamePrefix("mono_testPrintTimes");

        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1>* p_monodomain_problem = new MonodomainProblem<1>( &cell_factory );

        p_monodomain_problem->Initialise();
        p_monodomain_problem->SetWriteInfo();

        p_monodomain_problem->Solve();


        // read data entries for the time file and check correct
        //Hdf5DataReader data_reader1("MonoProblem1d", "mono_testPrintTimes");
        Hdf5DataReader data_reader1= p_monodomain_problem->GetDataReader();
        delete p_monodomain_problem;

        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), 4u);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);


    }


    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainWithMeshInMemoryToMeshalyzer() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(0.1);  //ms
        HeartConfig::Instance()->SetOutputDirectory("Monodomain2d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain2d");

        //set the postprocessing information we want
        std::vector<unsigned> origin_nodes;
        origin_nodes.push_back(0u);
        HeartConfig::Instance()->SetConductionVelocityMaps(origin_nodes);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;

        ///////////////////////////////////////////////////////////////////
        // monodomain
        ///////////////////////////////////////////////////////////////////
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(0.01, 0.1, 0.1);
        monodomain_problem.SetMesh(&mesh);
        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        //Clean previous output
        OutputFileHandler handler1("Monodomain2d", true); //This one makes sure that we leave a signature file in Monodomain2d so that we can clean up later
        OutputFileHandler handler("Monodomain2d/output", true);
        PetscTools::Barrier();

        // Checking that the following files don't exist in the output directory before calling Solve()
        unsigned num_files=6;
        std::string test_file_names[6]={"monodomain2d_mesh.pts", "monodomain2d_mesh.tri", "monodomain2d_V.dat",
              "ChasteParameters.xml", "ChasteDefaults.xml", "ConductionVelocityFromNode0.dat"};
        for (unsigned i=0; i<num_files; i++)
        {
            std::string compare_command = "cmp -s ";
            compare_command += handler.GetOutputDirectoryFullPath()+"/"+test_file_names[i];
            compare_command += " ";
            compare_command += "heart/test/data/Monodomain2d/";
            compare_command += test_file_names[i];
            TS_ASSERT_EQUALS(system(compare_command.c_str()), 512);//Not there
        }

        // now solve
        monodomain_problem.Solve();

        PetscTools::Barrier();
        // Compare output files
        for (unsigned i=0; i<num_files; i++)
        {
            if(test_file_names[i] == "monodomain2d_V.dat")
            {
                /*
                 * Since we started using bjacobi as the default preconditioner, parallel and sequential tests
                 * may return different results (always accurate to the tolerance requested). "diff" is unable
                 * to take this consideration into account.
                 *
                 * We will test that the file exists though.
                 */
                std::ifstream vm_file;
                std::string command = handler.GetOutputDirectoryFullPath()+"/"+test_file_names[i];
                vm_file.open(command.c_str());
                TS_ASSERT(vm_file.is_open());
                vm_file.close();
            }
            else
            {
                std::string compare_command = "diff --ignore-matching-lines=\"<cp:ChasteParameters\" --ignore-matching-lines=\"# \" ";
                compare_command += handler.GetOutputDirectoryFullPath()+"/"+test_file_names[i];
                compare_command += " ";
                compare_command += "heart/test/data/Monodomain2d/";
                compare_command += test_file_names[i];

                //Compare the new test file with one from the repository
                unsigned diff_result=system(compare_command.c_str());

                /* In case of failure, compare with known working alternative.
                 * This is for backward compatibility with XSD 2-3, in which case the
                 * above failure will be a little bit noisy (diff output goes to the screen)
                 * but we can live with that if the XSD 2-3 users can.
                 */

                if ((diff_result != 0) && (i==3|| i==4))
                {
                    compare_command += "_alt";
                    TS_ASSERT_EQUALS(system(compare_command.c_str()), 0);
                }
                else
                {
                    TS_ASSERT_EQUALS(diff_result, 0U);
                }
            }
        }
    }

    void TestMonodomainProblemCreates1DGeometry()
    {
        HeartConfig::Instance()->SetSpaceDimension(1);
        HeartConfig::Instance()->SetFibreLength(1, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(5.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("MonodomainCreatesGeometry");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory(-600 * 5000);

        MonodomainProblem<1> monodomain_problem( &cell_factory );
        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol = 5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 13.9682, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 13.9149, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 13.8216, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 13.7275, atol);
    }

    void TestMonodomainProblemCreates2DGeometry()
    {
        HeartConfig::Instance()->SetSpaceDimension(2);
        HeartConfig::Instance()->SetSheetDimensions(0.3, 0.3, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(4.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("MonodomainCreatesGeometry");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain2d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-600 * 5000);

        MonodomainProblem<2> monodomain_problem( &cell_factory );
        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol = 5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 17.5728, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 17.4562, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 17.2559, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 17.0440, atol);
    }

    void TestMonodomainProblemCreates3DGeometry()
    {
        HeartConfig::Instance()->SetSpaceDimension(3);
        HeartConfig::Instance()->SetSlabDimensions(0.1, 0.1, 0.1, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("MonodomainCreatesGeometry");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain3d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-600 * 5000);

        MonodomainProblem<3> monodomain_problem( &cell_factory );

        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);
        HeartConfig::Instance()->SetVisualizerOutputPrecision(6);
        monodomain_problem.Initialise();

        monodomain_problem.Solve();

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol = 5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 31.9227, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 32.1385, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 32.7569, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 33.7192, atol);


#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        if (PetscTools::IsParallel())
        {

            std::string filepath = OutputFileHandler::GetChasteTestOutputDirectory() + "MonodomainCreatesGeometry/vtk_output/";
            std::string basename = filepath + "monodomain3d";
            std::stringstream vtu_path;
            vtu_path << basename << "_" << PetscTools::GetMyRank() << ".vtu";

            FileFinder pvtk_file(basename + ".pvtu", RelativeTo::Absolute);
            TS_ASSERT(pvtk_file.Exists());
            FileFinder vtk_file(vtu_path.str(), RelativeTo::Absolute);
            TS_ASSERT(vtk_file.Exists());
            FileFinder param_file(filepath  + "ChasteParameters.xml", RelativeTo::Absolute);
            TS_ASSERT(param_file.Exists());
            FileFinder defaults_file(filepath  + "ChasteDefaults.xml", RelativeTo::Absolute);
            TS_ASSERT(defaults_file.Exists());
            FileFinder info_file(basename  + "_times.info", RelativeTo::Absolute);
            TS_ASSERT(info_file.Exists());
         }
#endif //CHASTE_VTK


        // We check that the cmgui files generated by the convert command in the problem class are OK
        // We compare against mesh files and one data file that are known to be visualized correctly in Cmgui.

        // The files written in parallel are different from the ones written in sequential because of the different node
        // numbering, therefore we test only the sequential case.
        // Note that the outputs of sequential and parallel simulations look the same when loaded with cmgui.
        // There are also minor rounding differences at the last decimal figure between sequential and parallel.
        EXIT_IF_PARALLEL

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "MonodomainCreatesGeometry/cmgui_output/";
        //the mesh files...

        std::string node_file1 = results_dir + "/monodomain3d.exnode";
        std::string node_file2 = "heart/test/data/CmguiData/monodomain/monodomain3dValid.exnode";
        std::string elem_file1 = results_dir + "/monodomain3d.exelem";
        std::string elem_file2 = "heart/test/data/CmguiData/monodomain/monodomain3dValid.exelem";

        bool comparison_result = CmguiMeshWriter<3,3>::CompareCmguiFiles(node_file1,node_file2);
        TS_ASSERT(comparison_result);
        comparison_result = CmguiMeshWriter<3,3>::CompareCmguiFiles(elem_file1,elem_file2);
        TS_ASSERT(comparison_result);
        //TS_ASSERT(CompareCmguiFiles(node_file1,node_file2));
        //TS_ASSERT(CompareCmguiFiles(elem_file1,elem_file2));

        //TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/monodomain3d.exelem heart/test/data/CmguiData/monodomain/monodomain3dValid.exelem").c_str()), 0);
        //TS_ASSERT_EQUALS(system(("cmp " + results_dir + "/monodomain3d.exnode heart/test/data/CmguiData/monodomain/monodomain3dValid.exnode").c_str()), 0);
        //...and one data file as example
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/monodomain3d_43.exnode heart/test/data/CmguiData/monodomain/monodomain3dValidData.exnode").c_str()), 0);
        //Info file
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + results_dir + "/monodomain3d_times.info heart/test/data/CmguiData/monodomain/monodomain3dValidData_times.info").c_str()), 0);
        //HeartConfig XML
        std::string filename_param = results_dir + "ChasteParameters.xml";
        std::ifstream file_param(filename_param.c_str());
        TS_ASSERT(file_param.is_open());
        file_param.close();
        std::string filename_default = results_dir + "ChasteDefaults.xml";
        std::ifstream file_default(filename_default.c_str());
        TS_ASSERT(file_default.is_open());
        file_default.close();

#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "MonodomainCreatesGeometry/vtk_output/";

        VtkMeshReader<3,3> mesh_reader(results_dir + "monodomain3d.vtu");
        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 1331U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 6000U);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
        TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( next_node[0] , 0.01 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );

        //Last time step and midway time step for V_m
        std::vector<double> v_at_last, v_at_100;
        mesh_reader.GetPointData( "V_000100", v_at_100);
        TS_ASSERT_DELTA( v_at_100[0],    47.9573, 1e-3 );
        TS_ASSERT_DELTA( v_at_100[665],  26.6333, 1e-3 );
        TS_ASSERT_DELTA( v_at_100[1330], -55.0584, 1e-3 );
        mesh_reader.GetPointData( "V_000200", v_at_last);
        TS_ASSERT_DELTA( v_at_last[0],    31.8997, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[665],  32.6873, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[1330], 34.5176, 1e-3 );

        //HeartConfig XML
        filename_param = results_dir + "ChasteParameters.xml";
        std::ifstream file_param2(filename_param.c_str());
        TS_ASSERT(file_param2.is_open());
        file_param2.close();
        filename_default = results_dir + "ChasteDefaults.xml";
        std::ifstream file_default2(filename_default.c_str());
        TS_ASSERT(file_default2.is_open());
        file_default2.close();
#endif //CHASTE_VTK
    }

    // Test the functionality for outputing the values of requested cell state variables
    void TestMonodomainProblemPrintsMultipleVariables() throw (Exception)
    {
        // Set configuration file
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/MultipleVariablesMonodomain.xml");

        // Set up problem
        PlaneStimulusCellFactory<CellFaberRudy2000FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        // Solve
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        // Get a reference to a reader object for the simulation results
        Hdf5DataReader data_reader1=monodomain_problem.GetDataReader();
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        // Check there is information about 101 timesteps (0, 0.01, 0.02, ...)
        TS_ASSERT_EQUALS( times.size(), 11u);
        TS_ASSERT_DELTA( times[0], 0.0, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.01, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.02, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.03, 1e-12);

        // There should be 101 values per variable and node.
        std::vector<double> node_5_v = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5_v.size(), 11u);

        std::vector<double> node_5_cai = data_reader1.GetVariableOverTime("cytosolic_calcium_concentration", 5);
        TS_ASSERT_EQUALS( node_5_cai.size(), 11U);

        std::vector<double> node_5_nai = data_reader1.GetVariableOverTime("ionic_concentrations__Nai", 5);
        TS_ASSERT_EQUALS( node_5_nai.size(), 11U);

        std::vector<double> node_5_ki = data_reader1.GetVariableOverTime("ionic_concentrations__Ki", 5);
        TS_ASSERT_EQUALS( node_5_ki.size(), 11U);
    }

    void TestMonodomainProblemExceptions() throw (Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms

        TS_ASSERT_THROWS_THIS(MonodomainProblem<1> monodomain_problem( NULL ),
                "AbstractCardiacProblem: Please supply a cell factory pointer to your cardiac problem constructor.");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        // Throws because we've not called initialise
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),"Cardiac tissue is null, Initialise() probably hasn\'t been called");

        // Throws because mesh filename is unset
        TS_ASSERT_THROWS_THIS(monodomain_problem.Initialise(),
                "No mesh given: define it in XML parameters file or call SetMesh()\n"
                "No Mesh provided (neither default nor user defined)");

        // Throws because initialise hasn't been called
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),"Cardiac tissue is null, Initialise() probably hasn\'t been called");

        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("");
        HeartConfig::Instance()->SetOutputFilenamePrefix("");

        monodomain_problem.Initialise();

        //Throws because the HDF5 slab isn't on the disk
        TS_ASSERT_THROWS_THIS(monodomain_problem.GetDataReader(),"Data reader invalid as data writer cannot be initialised");

        // throw because end time is negative
        HeartConfig::Instance()->SetSimulationDuration(-1.0); //ms
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),"End time should be in the future");
        HeartConfig::Instance()->SetSimulationDuration( 1.0); //ms

        // throws because output dir and filename are both ""
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),
                "Either explicitly specify not to print output (call PrintOutput(false)) or "
                "specify the output directory and filename prefix");

        // Throws because can't open the results file
        std::string directory("UnwriteableFolder");
        HeartConfig::Instance()->SetOutputDirectory(directory);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        OutputFileHandler handler(directory, false);
        if (PetscTools::AmMaster())
        {
            chmod(handler.GetOutputDirectoryFullPath().c_str(), 0555);
        }
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),
                              "Hdf5DataWriter could not create " + handler.GetOutputDirectoryFullPath() + "results.h5");
        if (PetscTools::AmMaster())
        {
            chmod(handler.GetOutputDirectoryFullPath().c_str(), 0755);
        }
    }

    /**
     * Not a very thorough test yet - just checks we can load a problem, simulate it, and
     * get expected results.
     *
     * This test relies on the h5 file generated in TestMonodomainProblem1D. Always run after!
     */
    void TestArchiving() throw(Exception)
    {
        // Based on TestMonodomainProblem1D()
        FileFinder archive_dir("monodomain_problem_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "monodomain_problem.arch";

        // Values to test against after load
        unsigned num_cells;

        // Save
        {
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("MonoProblemArchive");
            HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);
            // Set extra variable (to cover extension of HDF5 with named variables) - in a later test
            std::vector<std::string> output_variables;
            output_variables.push_back("cytosolic_calcium_concentration");
            HeartConfig::Instance()->SetOutputVariables( output_variables );

            PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();
            HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
            monodomain_problem.Solve();

            num_cells = monodomain_problem.GetTissue()->rGetCellsDistributed().size();

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacProblem<1,1,1>* const p_monodomain_problem = &monodomain_problem;
            (*p_arch) & p_monodomain_problem;
        }

        // Load
        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacProblem<1,1,1> *p_monodomain_problem;
            (*p_arch) >> p_monodomain_problem;

            // Check values
            TS_ASSERT_EQUALS(p_monodomain_problem->GetTissue()->rGetCellsDistributed().size(),
                             num_cells);

            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            p_monodomain_problem->Solve();

            // test whether voltages and gating variables are in correct ranges
            CheckMonoLr91Vars<1>(static_cast<MonodomainProblem<1,1>&>(*p_monodomain_problem));

            // check some voltages
            ReplicatableVector voltage_replicated(p_monodomain_problem->GetSolution());
            double atol = 5e-3;
            TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
            TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
            TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
            TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
            TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
            TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);

            for (unsigned index=0; index<voltage_replicated.GetSize(); index++)
            {
                //Shouldn't differ from the original run at all
                TS_ASSERT_DELTA(voltage_replicated[index], mVoltageReplicated1d2ms[index],  5e-11);
            }
            // check a progress report (or something) exists - in a noisy way
            TS_ASSERT_EQUALS(system(("ls " + OutputFileHandler::GetChasteTestOutputDirectory() + "MonoProblemArchive/").c_str()), 0);

            // check output file contains results for the whole simulation
            TS_ASSERT(CompareFilesViaHdf5DataReader("MonoProblemArchive", "MonodomainLR91_1d", true,
                                                    "MonoProblem1d", "MonodomainLR91_1d", true));

            // Free memory
            delete p_monodomain_problem;
        }
    }

    /**
     * Same as TestMonodomainProblem1D, except run the simulation in 2 halves: first to 1ms,
     * then continue to 2ms.
     *
     * This test relies on the h5 file generated in TestMonodomainProblem1D. Always run after!
     */
    void TestMonodomainProblem1dInTwoHalves() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem1dInTwoHalves");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        // Set extra variable (to cover extension of HDF5 with named variables)
        std::vector<std::string> output_variables;
        output_variables.push_back("cytosolic_calcium_concentration");
        HeartConfig::Instance()->SetOutputVariables( output_variables );

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
        monodomain_problem.Solve();
        ReplicatableVector voltage_replicated_midway(monodomain_problem.GetSolution());
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        monodomain_problem.Solve();

        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // check some voltages
        ReplicatableVector voltage_replicated(monodomain_problem.GetSolution());
        double atol = 5e-3;

        TS_ASSERT_DELTA(voltage_replicated[1], 20.7710232, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 21.5319692, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 22.9280817, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 24.0611303, atol);
        TS_ASSERT_DELTA(voltage_replicated[9], -0.770330519, atol);
        TS_ASSERT_DELTA(voltage_replicated[10], -19.2234919, atol);
        for (unsigned index=0; index<voltage_replicated.GetSize(); index++)
        {
            TS_ASSERT_DELTA(voltage_replicated[index], mVoltageReplicated1d2ms[index], 5e-11);
        }

        // check output file contains results for the whole simulation and agree with normal test
        TS_ASSERT(CompareFilesViaHdf5DataReader("MonoProblem1dInTwoHalves", "MonodomainLR91_1d", true,
                                                "MonoProblem1d", "MonodomainLR91_1d", true));
    }


    void TestMonodomain2dOriginalPermutationInParallel() throw(Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005, 0.0005));
        HeartConfig::Instance()->SetSimulationDuration(0.5); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoProblem2dOriginalPermutation");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_2d");
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory;

        // using the criss-cross mesh so wave propagates properly
        MonodomainProblem<2> monodomain_problem( &cell_factory );

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        HeartConfig::Instance()->SetVisualizeWithVtk(true);
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);

        monodomain_problem.Solve();

        // The following lines check that the output is always in the same permutation
        // order, regardless of whether it has been permuted internally.

        // In sequential mode, no permutation is applied
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), false);
        }
        else
        {
            // In parallel the mesh in memory has been permuted. Will check the output has been unpermuted.
            TS_ASSERT_EQUALS(HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(), true);
        }

        OutputFileHandler handler("MonoProblem2dOriginalPermutation/", false);
        //Note that without the "SetOutputUsingOriginalNodeOrdering" above the following would fail
        //since METIS partitioning will have changed the permutation

        /*
         * Meshalyzer format
         */
        //Mesh
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + handler.GetOutputDirectoryFullPath()
                                + "/output/MonodomainLR91_2d_mesh.pts   heart/test/data/MonoProblem2dOriginalPermutation/MonodomainLR91_2d_mesh.pts").c_str() ), 0);
        //Transmembrane
        std::string file1=handler.GetOutputDirectoryFullPath()+ "/output/MonodomainLR91_2d_V.dat";
        std::string file2="heart/test/data/MonoProblem2dOriginalPermutation/MonodomainLR91_2d_V.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-3)); //This can be quite flexible since the permutation differences will be quite large

        /*
         * Cmgui format
         */
        //Mesh
        TS_ASSERT_EQUALS(system(("diff -a -I \"Created by Chaste\" " + handler.GetOutputDirectoryFullPath()
                                + "/cmgui_output/MonodomainLR91_2d.exnode  heart/test/data/MonoProblem2dOriginalPermutation/MonodomainLR91_2d.exnode").c_str() ), 0);
        //Transmembrane
        file1=handler.GetOutputDirectoryFullPath()+ "/cmgui_output/MonodomainLR91_2d_50.exnode";
        file2="heart/test/data/MonoProblem2dOriginalPermutation/MonodomainLR91_2d_50.exnode";
        NumericFileComparison comp_cmgui(file1, file2);
        TS_ASSERT(comp_cmgui.CompareFiles(1e-3)); //This can be quite flexible since the permutation differences will be quite large

        /*
         * Vtk Format
         */
#ifdef CHASTE_VTK
        // Read in a VTKUnstructuredGrid as a mesh

        VtkMeshReader<2,2> mesh_reader(handler.GetOutputDirectoryFullPath()+"vtk_output/MonodomainLR91_2d.vtu");
        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 221U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 400U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 400U);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
        TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( next_node[0] , 0.01 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );

        //Last time step V_m
        std::vector<double> v_at_last;
        mesh_reader.GetPointData( "V_000050", v_at_last);
        TS_ASSERT_DELTA( v_at_last[0],    13.146, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[110],  13.146, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[220], -83.855, 1e-3 );

#endif // CHASTE_VTK

    }

};

#endif //_TESTMONODOMAINPROBLEM_HPP_
