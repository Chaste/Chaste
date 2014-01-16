/*

Copyright (c) 2005-2014, University of Oxford.
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


#ifndef _TESTMONODOMAINPROBLEM_HPP_
#define _TESTMONODOMAINPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "MonodomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "AbstractCvodeCell.hpp"
#include "LuoRudy1991.hpp"
#include "LuoRudy1991Cvode.hpp"
//#include "LuoRudy1991BackwardEuler.hpp"
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
#include "FileComparison.hpp"
#include "Warnings.hpp"
#include "ChasteSyscalls.hpp"

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

    void FinaliseCellCreation(std::vector<AbstractCardiacCellInterface* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        unsigned found_middle_point_reduced;
        MPI_Allreduce(&mFoundMiddlePoint, &found_middle_point_reduced, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        assert(found_middle_point_reduced == 1); // Only 1 cell should be stimulated
    }
};

#ifdef CHASTE_CVODE
/*
 * Cell factory for TestOutputDoesNotDependOnPrintTimestep. Returns CVODE cells
 */
class Cvode1dCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    Cvode1dCellFactory() : AbstractCardiacCellFactory<1>(),
        mpStimulus(new SimpleStimulus(-70000.0, 1.0, 0.0))
        {
        }

//    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
//    {
//        AbstractCardiacCell* p_cell;
//        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
//
//        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
//
//        if (x<0.3)
//        {
//            p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(p_solver, mpStimulus);
//        }
//        else
//        {
//            p_cell = new CellLuoRudy1991FromCellMLBackwardEuler(p_solver, mpZeroStimulus);
//        }
//
//
//        return p_cell;
//    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        AbstractCvodeCell* p_cell;
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;

        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];

        if (x<0.3)
        {
            p_cell = new CellLuoRudy1991FromCellMLCvode(p_empty_solver, mpStimulus);
        }
        else
        {
            p_cell = new CellLuoRudy1991FromCellMLCvode(p_empty_solver, mpZeroStimulus);
        }

        p_cell->SetMinimalReset(true);
        p_cell->SetTolerances(1e-5,1e-7);

        return p_cell;
    }
};
#endif // CHASTE_CVODE

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
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
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

        // Checking we only get warnings when no cells are allocated to processors with rank 2 or more
        if (PetscTools::GetMyRank() < 2)
        {
            // Rank is 0 or 1
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        }
        else
        {
            // Rank is 2 or more
            // There are only 2 nodes in this simulation so processors with rank 2 or more will have no nodes (and will through a warning)
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
            std::stringstream str1;
            str1 << "No cells were assigned to process " << PetscTools::GetMyRank() << " in AbstractCardiacTissue constructor. Advice: Make total number of processors no greater than number of nodes in the mesh";
            TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(), str1.str() );
        }
        Warnings::QuietDestroy();
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
        /*
         * HOW_TO_TAG Cardiac/Output
         * Calculating and outputting ionic currents ('derived quantities') in a tissue simulation using
         * class:HeartConfig - see also [wiki:ChasteGuides/CodeGenerationFromCellML#Derivedquantities this page].
         */
        // This is how to output an additional state variable
        output_variables.push_back("cytosolic_calcium_concentration");
        // or indeed an additional parameter or derived quantity. In the CellML
        // the variable 'potassium_currents' (sum of potassium currents) is annotated as a
        // derived quantity to be calculated.
        output_variables.push_back("potassium_currents");
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
        OutputFileHandler handler("MonoProblem1d", false);
        TS_ASSERT(handler.FindFile("progress_status.txt").Exists());
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
        OutputFileHandler handler("MonoProblem1din3d", false);
        TS_ASSERT(handler.FindFile("progress_status.txt").Exists());

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
        unsigned num_files = 5;
        std::string test_file_names[5] = {"monodomain2d_mesh.pts", "monodomain2d_mesh.tri", "monodomain2d_V.dat",
                                          "ChasteParameters.xml", "ConductionVelocityFromNode0.dat"};
        for (unsigned i=0; i<num_files; i++)
        {
            TS_ASSERT(!handler.FindFile(test_file_names[i]).Exists());
        }

        // now solve
        monodomain_problem.Solve();

        PetscTools::Barrier();
        // Compare output files
        for (unsigned i=0; i<num_files; i++)
        {
            if (test_file_names[i] == "monodomain2d_V.dat")
            {
                /**
                 * Since we started using bjacobi as the default preconditioner, parallel and sequential tests
                 * may return different results (always accurate to the tolerance requested).
                 * A straight FileComparison is unable to take this consideration into account.
                 * \todo can we use a NumericFileComparison then?
                 *
                 * We will test that the file exists though.
                 */
                TS_ASSERT(handler.FindFile(test_file_names[i]).Exists());
            }
            else
            {
                std::string file_1 = handler.GetOutputDirectoryFullPath()+"/"+test_file_names[i];
                std::string file_2 = "heart/test/data/Monodomain2d/" + test_file_names[i];

                /* In parallel, we want result (below) to be consistent across processes.  For this to happen no process
                 * should abort the file check.  This is why we pretend not to be called collectively.
                 */
                FileComparison comparer(file_1, file_2, false /*Pretend that it's not called collectively*/);
                comparer.SetIgnoreLinesBeginningWith("<cp:ChasteParameters");

                bool result = comparer.CompareFiles(0,false); // Don't throw a TS_ASSERT

                if (!result && (i == 3))
                {
                    // For "ChasteParameters.xml" there is an alternative file
                    // This is because different versions of XSD produce rounded/full precision floating point numbers
                    FileComparison comparer2(file_1, file_2 + "_alt");
                    comparer2.SetIgnoreLinesBeginningWith("<cp:ChasteParameters");
                    TS_ASSERT(comparer2.CompareFiles());
                }
                else
                {
                    TS_ASSERT(comparer.CompareFiles());
                }
            }
        }
    }

    void TestMonodomainProblemCreates1DGeometry()
    {
        HeartConfig::Instance()->SetSpaceDimension(1);
        HeartConfig::Instance()->SetFibreLength(1, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(5.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("MonodomainCreates1dGeometry");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain1d");

        // Switch on 1D VTK output
        ///\todo #2468 deadlocks in parallel
        //HeartConfig::Instance()->SetVisualizeWithVtk(true);
        HeartConfig::Instance()->SetVisualizeWithParallelVtk(true);

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
#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar

        std::string filepath = OutputFileHandler::GetChasteTestOutputDirectory() + "MonodomainCreates1dGeometry/vtk_output/";
        std::string basename = filepath + "monodomain1d";
        if (PetscTools::IsParallel())
        {
            FileFinder pvtk_file(basename + ".pvtu", RelativeTo::Absolute);
            TS_ASSERT(pvtk_file.Exists());
        }
        else
        {
            FileFinder vtu_file(basename + ".vtu", RelativeTo::Absolute);
            TS_ASSERT(vtu_file.Exists());
        }
#endif //CHASTE_VTK
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

        TS_ASSERT_DELTA(voltage_replicated[1], 31.8958, atol);
        TS_ASSERT_DELTA(voltage_replicated[3], 32.1109, atol);
        TS_ASSERT_DELTA(voltage_replicated[5], 32.7315, atol);
        TS_ASSERT_DELTA(voltage_replicated[7], 33.7011, atol);


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
            FileFinder info_file(basename  + "_times.info", RelativeTo::Absolute);
            TS_ASSERT(info_file.Exists());
         }
#else
        std::cout << "This test ran, but did not test VTK-dependent functions as VTK visualization is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
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

        FileComparison comparer(node_file1,node_file2);
        TS_ASSERT(comparer.CompareFiles());

        FileComparison comparer2(elem_file1,elem_file2);
        TS_ASSERT(comparer2.CompareFiles());

        //...and one data file as example
        FileComparison comparer3(results_dir + "/monodomain3d_43.exnode","heart/test/data/CmguiData/monodomain/monodomain3dValidData.exnode");
        TS_ASSERT(comparer3.CompareFiles());

        //Info file
        FileComparison comparer4(results_dir + "/monodomain3d_times.info","heart/test/data/CmguiData/monodomain/monodomain3dValidData_times.info");
        TS_ASSERT(comparer4.CompareFiles());

        //HeartConfig XML
        TS_ASSERT(FileFinder(results_dir + "ChasteParameters.xml").Exists());

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
        TS_ASSERT_DELTA( v_at_100[0],    48.0637, 1e-3 );
        TS_ASSERT_DELTA( v_at_100[665],  26.5404, 1e-3 );
        TS_ASSERT_DELTA( v_at_100[1330], -55.3058, 1e-3 );
        mesh_reader.GetPointData( "V_000200", v_at_last);
        TS_ASSERT_DELTA( v_at_last[0],    31.8730, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[665],  32.6531, 1e-3 );
        TS_ASSERT_DELTA( v_at_last[1330], 34.5152, 1e-3 );

        //HeartConfig XML
        TS_ASSERT(FileFinder(results_dir + "ChasteParameters.xml").Exists());
#else
        std::cout << "This test ran, but did not test VTK-dependent functions as VTK visualization is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
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
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),
                              "Cardiac tissue is null, Initialise() probably hasn\'t been called");

        // Throws because mesh filename is unset
        TS_ASSERT_THROWS_CONTAINS(monodomain_problem.Initialise(),
                "No mesh given: define it in XML parameters file or call SetMesh()\n"
                "No XML element Simulation/Mesh found in parameters when calling");

        // Throws because initialise hasn't been called
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),
                              "Cardiac tissue is null, Initialise() probably hasn\'t been called");

        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("");
        HeartConfig::Instance()->SetOutputFilenamePrefix("");

        monodomain_problem.Initialise();

        //Throws because the HDF5 slab isn't on the disk
        TS_ASSERT_THROWS_THIS(monodomain_problem.GetDataReader(),
                              "Data reader invalid as data writer cannot be initialised");

        // throw because end time is negative
        HeartConfig::Instance()->SetSimulationDuration(-1.0); //ms
        TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),
                              "End time should be in the future");
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
            chmod(handler.GetOutputDirectoryFullPath().c_str(), CHASTE_READ_EXECUTE);
        }
        H5E_BEGIN_TRY //Suppress HDF5 error in this test
        {
            TS_ASSERT_THROWS_THIS(monodomain_problem.Solve(),
                                  "Hdf5DataWriter could not create " + handler.GetOutputDirectoryFullPath() + "results.h5 , H5Fcreate error code = -1");
        }
        H5E_END_TRY;
        if (PetscTools::AmMaster())
        {
            chmod(handler.GetOutputDirectoryFullPath().c_str(), CHASTE_READ_WRITE_EXECUTE);
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
            // Set extra variables (to cover extension of HDF5 with named variables) - in a later test
            std::vector<std::string> output_variables;
            output_variables.push_back("cytosolic_calcium_concentration");
            output_variables.push_back("potassium_currents");
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
            // check a progress report exists
            OutputFileHandler handler("MonoProblemArchive", false);
            TS_ASSERT(handler.FindFile("progress_status.txt").Exists());

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

        // Set extra variables (to cover extension of HDF5 with named variables)
        std::vector<std::string> output_variables;
        output_variables.push_back("cytosolic_calcium_concentration");
        output_variables.push_back("potassium_currents");
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
        FileComparison comparer(handler.GetOutputDirectoryFullPath() + "/output/MonodomainLR91_2d_mesh.pts",
                       "heart/test/data/MonoProblem2dOriginalPermutation/MonodomainLR91_2d_mesh.pts");
        TS_ASSERT(comparer.CompareFiles());

        //Transmembrane
        std::string file1=handler.GetOutputDirectoryFullPath()+ "/output/MonodomainLR91_2d_V.dat";
        std::string file2="heart/test/data/MonoProblem2dOriginalPermutation/MonodomainLR91_2d_V.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(1e-3)); //This can be quite flexible since the permutation differences will be quite large

        /*
         * Cmgui format
         */
        //Mesh
        FileComparison comparer2(handler.GetOutputDirectoryFullPath() + "/cmgui_output/MonodomainLR91_2d.exnode",
                       "heart/test/data/MonoProblem2dOriginalPermutation/MonodomainLR91_2d.exnode");
        TS_ASSERT(comparer2.CompareFiles());

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
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 40U);

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

#else
        std::cout << "This test ran, but did not test VTK-dependent functions as VTK visualization is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK

    }


    void DontTestOutputDoesNotDependOnPrintTimestep() throw(Exception)
    {
#ifdef CHASTE_CVODE

        PetscOptionsSetValue("-ksp_monitor", "");

        const double mesh_spacing = 0.01;
        const double x_size = 1.0;
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructRegularSlabMesh(mesh_spacing, x_size);
        const unsigned max_node_index = mesh.GetNumNodes() - 1;
        Cvode1dCellFactory cell_factory;
        HeartConfig::Instance()->SetSimulationDuration(10.0); //ms
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(2.0));

        // Test two values of print timestep
        const unsigned num_print_steps_to_test = 2u;
        c_vector<double,num_print_steps_to_test> print_steps = Create_c_vector(0.1, 0.01);

        // Always use the same ODE and PDE timestep of 0.01.
        const double ode_and_pde_steps = 0.01; //ms

        for (unsigned i=0; i<num_print_steps_to_test; ++i)
        {
            std::stringstream str_stream;
            str_stream << print_steps[i];
            HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d_" + str_stream.str());
            HeartConfig::Instance()->SetOutputDirectory("TestCvodePrintTimestepDependence" + str_stream.str());
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_and_pde_steps, ode_and_pde_steps, print_steps[i]);

            MonodomainProblem<1> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.SetWriteInfo();
            monodomain_problem.Initialise();
            monodomain_problem.Solve();

            for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
            {
                AbstractCardiacCellInterface* p_cell = monodomain_problem.GetTissue()->GetCardiacCell(node_index);

                // This checks the maximum timestep in each Cvode cell has been set correctly.
                TS_ASSERT_EQUALS(static_cast<AbstractCvodeCell*>(p_cell)->GetTimestep(),
                                 ode_and_pde_steps);
            }
        }

        // Read results in and compare
        std::vector<double> V_to_compare;
        for (unsigned i=0; i<num_print_steps_to_test; ++i)
        {
            std::stringstream str_stream;
            str_stream << print_steps[i];
            Hdf5DataReader simulation_data("TestCvodePrintTimestepDependence" + str_stream.str(),
                                           "MonodomainLR91_1d_" + str_stream.str());
            std::vector<double> V_over_time = simulation_data.GetVariableOverTime("V", max_node_index);
            V_to_compare.push_back(V_over_time.back()); // Voltage at final time.
        }
        TS_ASSERT_DELTA(V_to_compare[0], V_to_compare[1], 1e-12);
#else
        std::cout << "Chaste is not configured to use CVODE on this machine, check your hostconfig settings if required.\n";
#endif // CHASTE_CVODE
    }


};

#endif //_TESTMONODOMAINPROBLEM_HPP_
