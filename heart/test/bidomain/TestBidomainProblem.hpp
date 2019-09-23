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


#ifndef TESTBIDOMAINPROBLEM_HPP_
#define TESTBIDOMAINPROBLEM_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>

#include "LuoRudy1991.hpp"
#include "FaberRudy2000.hpp"
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "Hdf5DataReader.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "HeartEventHandler.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ArchiveOpener.hpp"
#include "CmguiMeshWriter.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "MemfemMeshReader.hpp"
#include "BidomainSolver.hpp"
#include "CompareHdf5ResultsFiles.hpp"
#include "NumericFileComparison.hpp"
#include "Electrodes.hpp"
#include "SimpleBathProblemSetup.hpp"
#include "FileComparison.hpp"
#include "SingleTraceOutputModifier.hpp"

#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkVersion.h>
#endif

/**
 * This stimulus causes an exception to be thrown when normally a simple stimulus would activate.
 *
 * Used for testing exception handling in TestBidomain3dErrorHandling
 */
class ExceptionStimulus : public AbstractStimulusFunction
{
private:
    /** Duration of initial stimulus, typically in milliseconds */
    double mDuration;
    /** The time at which the stimulus starts, typically in milliseconds */
    double mTimeOfStimulus;

public:

    /**
     * Constructor.
     *
     * @param duration  Duration of initial stimulus milliseconds
     * @param timeOfStimulus  The time at which the stimulus starts (defaults to 0.0) milliseconds
     */
    ExceptionStimulus(double duration, double timeOfStimulus)
    {
        mDuration = duration;
        mTimeOfStimulus = timeOfStimulus;

        mDuration += (mDuration+mTimeOfStimulus)*DBL_EPSILON;
    }

    /**
     * @param time  time at which to return the stimulus
     * @return zero or throws an Exception at the time specified by constructor.
     */
    double GetStimulus(double time)
    {
        if (mTimeOfStimulus < time && time <= mDuration+mTimeOfStimulus)
        {
            EXCEPTION("Stimulus has caused a bad thing to happen in a cell model.");
        }
        return 0.0;
    }
};

/**
 * ExceptionStimulusCellFactory provides all cells with an ExceptionStimulus.
 *
 * Used for testing exception handling in TestBidomain3dErrorHandling
 */
template<class CELL, unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class ExceptionStimulusCellFactory : public AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>
{
protected:
    /** The stimulus to apply at stimulated nodes */
    boost::shared_ptr<ExceptionStimulus> mpStimulus;

public:
    /**
     * Constructor
     * @param stimulusDuration  The duration of time the exception should be applied (defaults to 1ms).
     * @param startTime When the exception should be introduced (defaults to 0.0).
     */
    ExceptionStimulusCellFactory(double stimulusDuration=1, double startTime = 0.0)
        : AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>()
    {
        mpStimulus.reset(new ExceptionStimulus(stimulusDuration, startTime));
    }

    /**
     * @param node  The global index of a node
     * @return  A cardiac cell which corresponds to this node.
     */
    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<SPACE_DIM>* pNode)
    {
        return new CELL(this->mpSolver, mpStimulus);
    }
};

class DelayedTotalStimCellFactory : public AbstractCardiacCellFactory<1>
{
private:
    // define a new stimulus
    boost::shared_ptr<SimpleStimulus> mpIntraStimulus;

public:
    DelayedTotalStimCellFactory(double mag)
        : AbstractCardiacCellFactory<1>(),
          mpIntraStimulus(new SimpleStimulus(  mag, 0.1, 0.1))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpIntraStimulus);
    }
};

class TestBidomainProblem : public CxxTest::TestSuite
{
private:
    std::vector<double> mSolutionReplicated1d2ms;///<Used to test differences between tests

public:
    void tearDown()
    {
        HeartConfig::Reset();
    }


    void TestBidomainErrorHandling()
    {
        if (PetscTools::GetNumProcs() > 2u)
        {
            // There are only 2 nodes in this simulation
            TS_TRACE("This test is not suitable for more than 2 processes.");
            return;
        }
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_1_element");
        HeartConfig::Instance()->SetOutputDirectory("BidomainWithErrors");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain");

        // Check the linear system can be solved to a low tolerance (in particular, checks the null space
        // stuff was implemented correctly
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-14);

        // We put in a massive stimulus, and alter the printing timestep from default, to cover
        // some memory handling in the case of exceptions being thrown (#1222).
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

        // 100x bigger stimulus than normal.
        ExceptionStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory(1, 0.02);
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );
        bidomain_problem.Initialise();

        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(),
           "Stimulus has caused a bad thing to happen in a cell model.");
    }

    //This test reproduces the problem size of heart/test/data/xml/bidomain1d_small.xml
    void TestBidomainSmallestMesh()
    {
        if (PetscTools::GetNumProcs() > 3u)
        {
            TS_TRACE("This test is not suitable for more than 3 processes.");
            return;
        }
       HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetFibreLength(0.1, 0.05);
        HeartConfig::Instance()->SetSpaceDimension(1u);
        HeartConfig::Instance()->SetOutputDirectory("Bidomain1dSmall");
        HeartConfig::Instance()->SetOutputFilenamePrefix("small");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.1, 0.1);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }


    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestBidomainDg01DPinned()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("bidomainDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        std::vector<unsigned> pinned_nodes;

        // check throws if the fixed node num isn't valid
        pinned_nodes.push_back(1000);
        bidomain_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);
        TS_ASSERT_THROWS_THIS( bidomain_problem.Solve(), "Fixed node number must be less than total number nodes" );

        /* HOW_TO_TAG Cardiac/Problem definition
         * Fix phi_e at particular nodes (note: this is not required)
         */
        // Pin extracellular potential of node 100 to 0
        pinned_nodes.clear();
        pinned_nodes.push_back(100);
        bidomain_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);

        try
        {
            bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            TS_FAIL(e.GetMessage());
        }

        DistributedVector striped_voltage = bidomain_problem.GetSolutionDistributedVector();
        DistributedVector::Stripe voltage(striped_voltage,0);

        for (DistributedVector::Iterator index = striped_voltage.Begin();
             index != striped_voltage.End();
             ++index)
        {
            // assuming LR model has Ena = 54.4 and Ek = -77
            double Ena   =  54.4;   // mV
            double Ek    = -77.0;   // mV

            TS_ASSERT_LESS_THAN_EQUALS( voltage[index], Ena +  30);
            TS_ASSERT_LESS_THAN_EQUALS(-voltage[index] + (Ek-30), 0);

            std::vector<double> r_ode_vars = bidomain_problem.GetBidomainTissue()->GetCardiacCell(index.Global)->GetStdVecStateVariables();
            for (int j=0; j<8; j++)
            {
                // if not voltage or calcium ion conc, test whether between 0 and 1
                if ((j!=0) && (j!=7))
                {
                    TS_ASSERT_LESS_THAN_EQUALS( r_ode_vars[j], 1.0);
                    TS_ASSERT_LESS_THAN_EQUALS(-r_ode_vars[j], 0.0);
                }
            }

            // wave shouldn't have reached the second half of the mesh so
            // these should all be near the resting potential

            if (index.Global>50)
            {
                TS_ASSERT_DELTA(voltage[index], -83.85, 0.1);
            }

            // final voltages for nodes 0 to 5 produced with ksp_rtol=1e-9
            double test_values[6]={31.0335, 28.9214, 20.0279, -3.92649, -57.9395, -79.7754};

            for (unsigned node=0; node<=5; node++)
            {
                if (index.Global == node)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(voltage[index], test_values[node], 7e-3);
                    //With ksp_rtol set to 1e-6 the starting value may lead to changes of more that 1e-3 in final answer
                }
            }
        }
        DistributedVector::Stripe extracellular_potential(striped_voltage,1);
        if (striped_voltage.IsGlobalIndexLocal(100))
        {
            TS_ASSERT_DELTA(extracellular_potential[100], 0.0, 1e-6);
        }
    }


    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestBidomainDg01DAveragePhiEOverDifferentRows()
    {
        HeartEventHandler::Disable();

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("bidomainDg01d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        // Final values to test against have been produced with ksp_rtol=1e-9
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-8);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);


        // Apply the constraint 'Average phi = 0' to nodes 0, 50, 100.
        for (unsigned node=0; node<bidomain_problem.rGetMesh().GetNumNodes(); node+=50)
        {
            bidomain_problem.SetNodeForAverageOfPhiZeroed(node);

            try
            {
                bidomain_problem.Solve();
            }
            catch (Exception e)
            {
                TS_FAIL(e.GetMessage());
            }

            DistributedVector striped_voltage = bidomain_problem.GetSolutionDistributedVector();
            DistributedVector::Stripe voltage(striped_voltage,0);
            DistributedVector::Stripe phi_e(striped_voltage,1);

            for (DistributedVector::Iterator index = striped_voltage.Begin();
                 index != striped_voltage.End();
                 ++index)
            {
                // assuming LR model has Ena = 54.4 and Ek = -77
                double Ena   =  54.4;   // mV
                double Ek    = -77.0;   // mV

                TS_ASSERT_LESS_THAN_EQUALS( voltage[index], Ena +  30);
                TS_ASSERT_LESS_THAN_EQUALS(-voltage[index] + (Ek-30), 0);

                std::vector<double> r_ode_vars = bidomain_problem.GetBidomainTissue()->GetCardiacCell(index.Global)->GetStdVecStateVariables();
                for (int j=0; j<8; j++)
                {
                    // if not voltage or calcium ion conc, test whether between 0 and 1
                    if ((j!=0) && (j!=7))
                    {
                        TS_ASSERT_LESS_THAN_EQUALS( r_ode_vars[j], 1.0);
                        TS_ASSERT_LESS_THAN_EQUALS(-r_ode_vars[j], 0.0);
                    }
                }

                // wave shouldn't have reached the second half of the mesh so
                // these should all be near the resting potential

                if (index.Global>50)
                {
                    TS_ASSERT_DELTA(voltage[index], -83.85, 0.1);
                }

                // final voltages for nodes 0 to 5 produced with ksp_rtol=1e-9
                double voltage_test_values[6]={31.0335, 28.9214, 20.0279, -3.92649, -57.9395, -79.7754};

                if (index.Global<6)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(voltage[index], voltage_test_values[index.Global], 7e-3);
                }

                // final extracellular potencials for nodes 0 to 5 produced with ksp_rtol=1e-9
                double phi_e_test_values[6]={-55.2567, -54.2006, -49.7538, -37.7767, -10.7701, 0.148278};

                if (index.Global<6)
                {
                    // test against hardcoded value to check nothing has changed
                    TS_ASSERT_DELTA(phi_e[index], phi_e_test_values[index.Global], 7e-3);
                }
            }

            // check mean of extracellular potential is 0
            double local_phi_e=0.0;
            double total_phi_e=0.0;

            for (DistributedVector::Iterator index = striped_voltage.Begin();
                 index != striped_voltage.End();
                 ++index)
            {
                local_phi_e += phi_e[index];
            }

            int ierr = MPI_Allreduce(&local_phi_e, &total_phi_e, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
            TS_ASSERT_EQUALS(ierr, (int)MPI_SUCCESS);

            TS_ASSERT_DELTA(total_phi_e, 0, 1e-4);

            bidomain_problem.Initialise();
        }

        // Coverage of the exception in the solver itself
        BoundaryConditionsContainer<1,1,2> container;

        BidomainSolver<1,1> bidomain_solver(false,
                                            &bidomain_problem.rGetMesh(),
                                            bidomain_problem.GetBidomainTissue(),
                                            &container);

        TS_ASSERT_THROWS_THIS(bidomain_solver.SetRowForAverageOfPhiZeroed(0),
                "Row for applying the constraint \'Average of phi_e = zero\' should be odd in C++ like indexing");

        HeartEventHandler::Enable();
    }

    /*
     * The monodomain equations are obtained by taking the limit of the bidomain
     * equations as sigma_e tends to infinity (corresponding to the extracellular
     * space being grounded). Therefore, if we set sigma_e very large (relative to
     * sigma_i) in a bidomain simulation it should agree with a monodomain
     * simulation with the same parameters.
     *
     * NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
     * surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
     * the equations have been divided through by the surface-area-to-volume ratio.
     * (Historical reasons...)
     *
     */
    void TestCompareBidomainProblemWithMonodomain()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("Monodomain1dVersusBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("monodomain1d");

        Vec monodomain_results;

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;

        // To avoid an issue with the Event handler only one simulation should be
        // in existance at a time: therefore monodomain simulation is defined in a block
        {
            ///////////////////////////////////////////////////////////////////
            // monodomain
            ///////////////////////////////////////////////////////////////////
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();

            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);

            // now solve
            monodomain_problem.Solve();

            VecDuplicate(monodomain_problem.GetSolution(), &monodomain_results);
            VecCopy(monodomain_problem.GetSolution(), monodomain_results);
        }


        ///////////////////////////////////////////////////////////////////
        // bidomain
        ///////////////////////////////////////////////////////////////////

        // keep the intra conductivity to be the same as monodomain
        // and the extra conductivity to be very large in comparison
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(1));
        HeartConfig::Instance()->SetOutputDirectory("Bidomain1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain1d");

        BidomainProblem<1> bidomain_problem( &cell_factory );

        // Not really needed for coverage, but it's worth checking that logic works properly in both scenarios.
        TS_ASSERT(!bidomain_problem.GetHasBath());

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        // now solve
        bidomain_problem.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector dist_bidomain_voltage = bidomain_problem.GetSolutionDistributedVector();
        DistributedVector monodomain_voltage = dist_bidomain_voltage.GetFactory()->CreateDistributedVector(monodomain_results);
        DistributedVector::Stripe bidomain_voltage(dist_bidomain_voltage, 0);
        DistributedVector::Stripe extracellular_potential(dist_bidomain_voltage, 1);

        for (DistributedVector::Iterator index = monodomain_voltage.Begin();
             index != monodomain_voltage.End();
             ++index)
        {
            if (index.Global==0)
            {
                TS_ASSERT_LESS_THAN(0, monodomain_voltage[index]);
            }
            // the mono and bidomains should agree closely
            TS_ASSERT_DELTA(monodomain_voltage[index], bidomain_voltage[index], 0.4);

            // the extracellular potential should be uniform
            TS_ASSERT_DELTA(extracellular_potential[index], 0, 0.06);
        }

        PetscTools::Destroy(monodomain_results);
    }


    ///////////////////////////////////////////////////////////////////
    // Solve a simple simulation and check the output was only
    // printed out at the correct times
    ///////////////////////////////////////////////////////////////////
    void TestBidomainProblemPrintsOnlyAtRequestedTimesAndOnlyRequestedNodes()
    {
        HeartEventHandler::Disable();

        HeartConfig::Instance()->SetPrintingTimeStep(0.1);
        HeartConfig::Instance()->SetPdeTimeStep(0.01);
        HeartConfig::Instance()->SetSimulationDuration(0.3);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("Bidomain1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_testPrintTimes");

        // run testing PrintingTimeSteps
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1>* p_bidomain_problem = new BidomainProblem<1>( &cell_factory );


        /* HOW_TO_TAG Cardiac/Output
         * Only output data for particular nodes
         */
        //Restrict the number of nodes
        std::vector<unsigned> nodes_to_be_output;
        nodes_to_be_output.push_back(0);
        nodes_to_be_output.push_back(5);
        nodes_to_be_output.push_back(10);
        p_bidomain_problem->SetOutputNodes(nodes_to_be_output);

        /* HOW_TO_TAG Cardiac/Output
         * Output data using a light-weight output modifier.  This can be used in addition to regular HDF5
         * output or can replace it.
         */
         // Duplicate this with a single trace at node 5
        boost::shared_ptr<SingleTraceOutputModifier> trace_5(new SingleTraceOutputModifier("trace_5.txt", 5, 0.1));
        p_bidomain_problem->AddOutputModifier(trace_5);

        // for coverage:
        p_bidomain_problem->SetWriteInfo();

        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();

        // Extra trace data
        OutputFileHandler handler("Bidomain1d", false);
        NumericFileComparison comp(handler.GetOutputDirectoryFullPath()+ "trace_5.txt", "heart/test/data/Bidomain1d/trace_5.txt");
        TS_ASSERT(comp.CompareFiles(5e-4));

        // read data entries for the time file and check correct
        Hdf5DataReader data_reader1=p_bidomain_problem->GetDataReader();
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.10, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.20, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.30, 1e-12);

        //Get back node over all times
        std::vector<double> node_0 = data_reader1.GetVariableOverTime("V", 0);
        TS_ASSERT_EQUALS( node_0.size(), 4U);
        TS_ASSERT_DELTA( node_0[0], -83.853, 1e-10);
        TS_ASSERT_DELTA( node_0[1], -83.8354, 5e-4);
        TS_ASSERT_DELTA( node_0[2], -83.8266, 3e-4);
        TS_ASSERT_DELTA( node_0[3], -83.8200, 3e-4);
        std::vector<double> node_5 = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5.size(), 4U);
        std::vector<double> node_10 = data_reader1.GetVariableOverTime("V", 10);
        TS_ASSERT_EQUALS( node_10.size(), 4U);

        //Can't read back this node as it wasn't written
        TS_ASSERT_THROWS_THIS( data_reader1.GetVariableOverTime("V", 1),
                "The incomplete dataset 'Data' does not contain info of node 1");

        delete p_bidomain_problem;

        p_bidomain_problem = new BidomainProblem<1>( &cell_factory );

        // Now check that we can turn off output printing
        // Output should be the same as above: printing every 10th time step
        // because even though we set to print every time step...
        HeartConfig::Instance()->SetPrintingTimeStep(1);
        // ...we have output turned off
        p_bidomain_problem->PrintOutput(false);

        p_bidomain_problem->Initialise();
        p_bidomain_problem->Solve();

        Hdf5DataReader data_reader3=p_bidomain_problem->GetDataReader();
        times = data_reader3.GetUnlimitedDimensionValues();

        TS_ASSERT_EQUALS( times.size(), (unsigned) 4);
        TS_ASSERT_DELTA( times[0], 0.00,  1e-12);
        TS_ASSERT_DELTA( times[1], 0.10,  1e-12);
        TS_ASSERT_DELTA( times[2], 0.20,  1e-12);
        TS_ASSERT_DELTA( times[3], 0.30,  1e-12);

        delete p_bidomain_problem;
        HeartEventHandler::Enable();
    }

    void TestBidomainFallsOverProducesOutput()
    {
#ifndef NDEBUG //Note that this test relies on the debug VerifyStateVariables() method throwing
        HeartConfig::Instance()->SetSimulationDuration(0.3);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainFallsOver");
        HeartConfig::Instance()->SetOutputFilenamePrefix("res");
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer();

        // Something happens at 0.1ms

        // DelayedTotalStimCellFactory bidomain_cell_factory(-6e5); //Normal stimulus
        DelayedTotalStimCellFactory bidomain_cell_factory(-6e6); //Takes sodium out of range
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.Initialise();

        TS_ASSERT_THROWS_CONTAINS(bidomain_problem.Solve(),
               "State variable fast_sodium_current_m_gate__m has gone out of range. Check numerical parameters, for example time and space stepsizes");

        // Make sure that there's time for the files to be written
        // (most files are only written by the master)
        PetscTools::Barrier();

        //Test for regular output
        Hdf5DataReader data_reader=bidomain_problem.GetDataReader();
        std::vector<double> times = data_reader.GetUnlimitedDimensionValues();
        //TS_ASSERT_EQUALS( times.size(), 31U);// For normal stimulation
        TS_ASSERT_EQUALS( times.size(), 31U);  // For fully allocated hdf5
        //TS_ASSERT_EQUALS( times.size(), 21U);// For over stimulation

        TS_ASSERT_DELTA( times[1], 0.01, 1e-12);
        //TS_ASSERT_DELTA( times.back(), 0.30, 1e-12);//For normal stimulation
        //TS_ASSERT_DELTA( times.back(), 0.20, 1e-12);//For over stimulation
        TS_ASSERT_DELTA( times.back(), 0, 1e-12);//For hdf5 writer allocated up to 30 slots
        TS_ASSERT_DELTA( times[20], 0.20, 1e-12);//For hdf5 writer allocated up to 30 slots
        TS_ASSERT_DELTA( times[21],    0, 1e-12);//For hdf5 writer allocated up to 30 slots

        // Test for post-processed output (and don't wipe the directory!)
        OutputFileHandler handler("BidomainFallsOver/output", false);

        std::string files[6] = {"res_mesh.pts","res_mesh.cnnx","ChasteParameters.xml",
                                "res_Phi_e.dat","res_V.dat","res_times.info"};

        for (unsigned i=0; i<6; i++)
        {
            TS_ASSERT(handler.FindFile(files[i]).Exists());
        }
#endif //NDEBUG Note that this test relies on the debug VerifyStateVariables() method throwing
    }


    void TestBidomainProblemExceptions()
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        // Throws because we've not called initialise
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "Cardiac tissue is null, Initialise() probably hasn\'t been called");

        // Throws because mesh filename is unset
        TS_ASSERT_THROWS_CONTAINS(bidomain_problem.Initialise(),
                "No mesh given: define it in XML parameters file or call SetMesh()\n"
                "No XML element Simulation/Mesh found in parameters when calling");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        TS_ASSERT_THROWS_NOTHING(bidomain_problem.Initialise());

        // Negative simulation duration
        HeartConfig::Instance()->SetSimulationDuration(-1.0);  //ms
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "End time should be in the future");
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms

        // set output data to avoid their exceptions (which is covered in TestMonoDg0Assembler
        HeartConfig::Instance()->SetOutputDirectory("temp");
        HeartConfig::Instance()->SetOutputFilenamePrefix("temp");

        // Exception caused by relative tolerance and no clamping
        HeartConfig::Instance()->SetUseRelativeTolerance(2e-3);
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "Bidomain external voltage is not bounded in this simulation - use KSP *absolute* tolerance");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-3);

        // Throws (in AbstractCardiacProblem) as dt does not divide end time
        HeartConfig::Instance()->SetPrintingTimeStep(0.15);
        HeartConfig::Instance()->SetPdeTimeStep(0.15);
        TS_ASSERT_THROWS_THIS( bidomain_problem.Solve(),"PDE timestep does not seem to divide end time - check parameters" );
        HeartConfig::Instance()->SetPdeTimeStep(0.01);
        HeartConfig::Instance()->SetPrintingTimeStep(0.01);

        //Throws because the node number is slightly bigger than the number of nodes in the mesh
        std::vector<unsigned> too_large;
        too_large.push_back(4358743);
        bidomain_problem.SetFixedExtracellularPotentialNodes(too_large);
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "Fixed node number must be less than total number nodes");

        //Explicitly reset the counters so the next test in the test suite doesn't find on
        HeartEventHandler::Reset();
    }


    void TestCompareOrthotropicWithAxisymmetricBidomain()
    {
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/box_shaped_heart/box_heart", cp::media_type::Orthotropic);
        HeartConfig::Instance()->SetOutputDirectory("OrthotropicBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("ortho3d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory;

        ///////////////////////////////////////////////////////////////////
        // orthotropic
        ///////////////////////////////////////////////////////////////////

        BidomainProblem<3> orthotropic_bido( &cell_factory );

        orthotropic_bido.Initialise();
        orthotropic_bido.Solve();

        ///////////////////////////////////////////////////////////////////
        // axisymmetric
        ///////////////////////////////////////////////////////////////////
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/box_shaped_heart/box_heart", cp::media_type::Axisymmetric);
        HeartConfig::Instance()->SetOutputDirectory("AxisymmetricBidomain");
        HeartConfig::Instance()->SetOutputFilenamePrefix("axi3d");

        BidomainProblem<3> axisymmetric_bido( &cell_factory);
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        axisymmetric_bido.Initialise();
        axisymmetric_bido.Solve();

        ///////////////////////////////////////////////////////////////////
        // compare
        ///////////////////////////////////////////////////////////////////
        DistributedVector orthotropic_solution = orthotropic_bido.GetSolutionDistributedVector();
        DistributedVector axisymmetric_solution = axisymmetric_bido.GetSolutionDistributedVector();

        DistributedVector::Stripe ortho_voltage(orthotropic_solution, 0);
        DistributedVector::Stripe axi_voltage(axisymmetric_solution, 0);

        DistributedVector::Stripe ortho_ex_pot(orthotropic_solution, 1);
        DistributedVector::Stripe axi_ex_pot(axisymmetric_solution, 1);

        for (DistributedVector::Iterator index = orthotropic_solution.Begin();
             index != orthotropic_solution.End();
             ++index)
        {
            TS_ASSERT_DELTA(ortho_voltage[index], axi_voltage[index], 1e-7);
            TS_ASSERT_DELTA(ortho_ex_pot[index], axi_ex_pot[index], 1e-7);
        }

        // We check that the cmgui files generated by the convert command in the problem class are OK
        // We compare against mesh files and one data file that are known to be visualized correctly in Cmgui.

        // The files written in parallel are different from the ones written in sequential because of the different node
        // renumbering, therefore we test only the sequential case.
        // Note that the outputs of sequential and parallel simulations look the same when loaded with cmgui.
        // There are also minor rounding differences at the last decimal figure between sequential and parallel.
        EXIT_IF_PARALLEL

        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "AxisymmetricBidomain/cmgui_output/";

         //the mesh files...
        std::string node_file1 = results_dir + "/axi3d.exnode";
        std::string node_file2 = "heart/test/data/CmguiData/bidomain/bidomain3dValid.exnode";
        std::string elem_file1 = results_dir + "/axi3d.exelem";
        std::string elem_file2 = "heart/test/data/CmguiData/bidomain/bidomain3dValid.exelem";

        FileComparison comparer(node_file1,node_file2);
        TS_ASSERT(comparer.CompareFiles());

        FileComparison comparer2(elem_file1,elem_file2);
        TS_ASSERT(comparer2.CompareFiles());

        //...and a couple of data files as examples
        FileComparison comparer3(results_dir + "/axi3d_61.exnode",
                "heart/test/data/CmguiData/bidomain/bidomain3dValidData61.exnode");
        TS_ASSERT(comparer3.CompareFiles());

        FileComparison comparer4(results_dir + "/axi3d_100.exnode",
                "heart/test/data/CmguiData/bidomain/bidomain3dValidData100.exnode");
        TS_ASSERT(comparer4.CompareFiles());

        //info file
        FileComparison comparer5(results_dir + "/axi3d_times.info",
                "heart/test/data/CmguiData/bidomain/axi3d_times.info");
        TS_ASSERT(comparer5.CompareFiles());

        //HeartConfig XML
        TS_ASSERT(FileFinder(results_dir + "ChasteParameters.xml").Exists());

#ifdef CHASTE_VTK
// Requires  "sudo aptitude install libvtk5-dev" or similar
        results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "AxisymmetricBidomain/vtk_output/";

        VtkMeshReader<3,3> mesh_reader(results_dir + "axi3d.vtu");
        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 50U);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 139U);

        std::vector<double> first_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( first_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( first_node[1] , 0.0, 1e-6 );
        TS_ASSERT_DELTA( first_node[2] , 0.0 , 1e-6 );

        std::vector<double> next_node = mesh_reader.GetNextNode();
        TS_ASSERT_DELTA( next_node[0] , 0.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[1] , 5.0 , 1e-6 );
        TS_ASSERT_DELTA( next_node[2] , 0.0 , 1e-6 );

        //V_m and phi_e samples
        std::vector<double> v_at_last, v_at_50, phi_at_50;
        mesh_reader.GetPointData( "V_000050", v_at_50);
        TS_ASSERT_DELTA( v_at_50[0],  -83.6547, 1e-3 );
        TS_ASSERT_DELTA( v_at_50[25], -83.8549, 1e-3 );
        TS_ASSERT_DELTA( v_at_50[49], -83.8557, 1e-3 );
        mesh_reader.GetPointData( "Phi_e_000050", phi_at_50);
        TS_ASSERT_DELTA( phi_at_50[0],  -0.0332, 1e-3 );
        TS_ASSERT_DELTA( phi_at_50[25],  0.0064, 1e-3 );
        TS_ASSERT_DELTA( phi_at_50[49],  0.0062, 1e-3 );
        mesh_reader.GetPointData( "V_000100", v_at_last);
        TS_ASSERT_DELTA( v_at_last[0],   -83.6863, 1e-3 );

        //HeartConfig XML
        TS_ASSERT(FileFinder(results_dir + "ChasteParameters.xml").Exists());
#else
        std::cout << "This test ran, but did not test VTK-dependent functions as VTK visualization is not enabled." << std::endl;
        std::cout << "If required please install and alter your hostconfig settings to switch on chaste support." << std::endl;
#endif //CHASTE_VTK
     }

    // Test the functionality for outputting the values of requested cell state variables
    void TestBidomainProblemPrintsMultipleVariables()
    {
        // Get the singleton in a clean state
        HeartConfig::Instance()->Reset();

        // Set configuration file
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/MultipleVariablesBidomain.xml");

        // Override the variables we are interested in writing.
        std::vector<std::string> output_variables;
        output_variables.push_back("calcium_dynamics__Ca_NSR");
        output_variables.push_back("ionic_concentrations__Nai");
        output_variables.push_back("fast_sodium_current_j_gate__j");
        output_variables.push_back("ionic_concentrations__Ki");

        HeartConfig::Instance()->SetOutputVariables( output_variables );

        // Set up problem
        PlaneStimulusCellFactory<CellFaberRudy2000FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        // Solve
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // Get a reference to a reader object for the simulation results
        Hdf5DataReader data_reader1 = bidomain_problem.GetDataReader();
        std::vector<double> times = data_reader1.GetUnlimitedDimensionValues();

        // Check there is information about 11 timesteps (0, 0.01, 0.02, ...)
        unsigned num_steps = 11u;
        TS_ASSERT_EQUALS( times.size(), num_steps);
        TS_ASSERT_DELTA( times[0], 0.0, 1e-12);
        TS_ASSERT_DELTA( times[1], 0.01, 1e-12);
        TS_ASSERT_DELTA( times[2], 0.02, 1e-12);
        TS_ASSERT_DELTA( times[3], 0.03, 1e-12);

        // There should be 11 values per variable and node.
        std::vector<double> node_5_v = data_reader1.GetVariableOverTime("V", 5);
        TS_ASSERT_EQUALS( node_5_v.size(), num_steps);

        std::vector<double> node_5_phi = data_reader1.GetVariableOverTime("Phi_e", 5);
        TS_ASSERT_EQUALS( node_5_phi.size(), num_steps);

        for (unsigned i=0; i<output_variables.size(); i++)
        {
            unsigned global_index = 2+i*2;
            std::vector<double> values = data_reader1.GetVariableOverTime(output_variables[i], global_index);
            TS_ASSERT_EQUALS( values.size(), num_steps);

            // Check the last values match the cells' state
            if (bidomain_problem.rGetMesh().GetDistributedVectorFactory()->IsGlobalIndexLocal(global_index))
            {
                AbstractCardiacCellInterface* p_cell = bidomain_problem.GetTissue()->GetCardiacCell(global_index);
                TS_ASSERT_DELTA(values.back(), p_cell->GetAnyVariable(output_variables[i],0), 1e-12);
            }
        }
    }


    /* HOW_TO_TAG Cardiac/Output
     * Output all cell model state variables for the cell model used in a particular simulation
     */
    void TestBidomainProblemPrintsAllStateVariables()
    {
        // Get the singleton in a clean state
        HeartConfig::Instance()->Reset();

        // Set configuration file
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/MultipleVariablesBidomain.xml");

        //Override the output directory so that this test goes into a fresh place
        HeartConfig::Instance()->SetOutputDirectory("AllStateVariablesBidomain");

        {
            //MultipleVariablesBidomain.xml uses FaberRudy.  Interrogate a model to get all the state variables
            boost::shared_ptr<ZeroStimulus> p_stimulus(new ZeroStimulus());
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            CellFaberRudy2000FromCellML temporary_fr2000_ode_system(p_solver, p_stimulus);
            TS_ASSERT_EQUALS(temporary_fr2000_ode_system.rGetStateVariableNames().size(), 25u);


            // Copy all the state variable names from the temporary cell into HeartConfig so that they are all printed
            HeartConfig::Instance()->SetOutputVariables( temporary_fr2000_ode_system.rGetStateVariableNames() );
        }

        // Set up problem
        PlaneStimulusCellFactory<CellFaberRudy2000FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        // Solve
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // Get a reference to a reader object for the simulation results
        Hdf5DataReader data_reader1 = bidomain_problem.GetDataReader();


        // Get an obscure state variable
        std::vector<double> obscure_at_5 = data_reader1.GetVariableOverTime("slow_delayed_rectifier_potassium_current_xs2_gate__xs2", 5);
        TS_ASSERT_EQUALS( obscure_at_5.size(), 11u);
        for (unsigned time_step=0; time_step<11; time_step++)
        {
            TS_ASSERT_DELTA( obscure_at_5[time_step], 0.0044, 1e-3); //Does not change over time at this node
        }
    }

    /*
     * Simple bidomain simulation to test against in the archiving tests below
     *
     * NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
     * surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
     * the equations have been divided through by the surface-area-to-volume ratio.
     * (Historical reasons...)
     */
    void TestSimpleBidomain1D()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1dInOneGo");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // Check some voltages
        ReplicatableVector solution_replicated(bidomain_problem.GetSolution());

        double atol=5e-3;

        TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
        TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
        TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
        TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
        TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
        TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
        TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);

        // Save solution for later comparison
        for (unsigned index=0; index<solution_replicated.GetSize(); index++)
        {
            mSolutionReplicated1d2ms.push_back(solution_replicated[index]);
        }
    }

    /*
     * Simple bidomain simulation with a straight permutation applied.
     * HOW_TO_TAG Cardiac/Output
     * Use the `SingleTraceOutputModifier` to output based on a global index
     * (index AFTER any permutation has been applied)
     *
     * NOTE: This test uses NON-PHYSIOLOGICAL parameters values.
     */
    void TestPermutedBidomain1D()
    {

        TetrahedralMesh<1,1> mesh;
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1mm_10_elements");
        mesh.ConstructFromMeshReader(reader);

        std::vector<unsigned> rotation_perm;

        unsigned number_nodes=11;
        for (unsigned index=0; index<(unsigned)number_nodes; index++)
        {
            rotation_perm.push_back( (index + 3) % number_nodes); // 3, 4, ... 0, 1, 2
        }

        //Rotate the permutation
        mesh.PermuteNodes(rotation_perm);

        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        HeartConfig::Instance()->SetSimulationDuration(0.5); //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainUnpermuted1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        // Single trace at node 5
        boost::shared_ptr<SingleTraceOutputModifier> trace_5(new SingleTraceOutputModifier("trace_5.txt", 5u, 0.1));
        bidomain_problem.AddOutputModifier(trace_5);

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // Extra trace data
        OutputFileHandler handler1("BidomainUnpermuted1d", false);
        NumericFileComparison comp1(handler1.GetOutputDirectoryFullPath()+ "trace_5.txt", "heart/test/data/BidomainUnpermuted1d/trace_5.txt");
        TS_ASSERT(comp1.CompareFiles(5e-4));

        //Can't read in the final mesh since it's a 1d example...
        OutputFileHandler handler("BidomainUnpermuted1d/output", false);

        //Mesh
        FileComparison comparer(handler.GetOutputDirectoryFullPath() + "/BidomainLR91_1d_mesh.pts",
                "heart/test/data/BidomainUnpermuted1d/BidomainLR91_1d_mesh.pts");
        TS_ASSERT(comparer.CompareFiles());

        //Transmembrane
        std::string file1=handler.GetOutputDirectoryFullPath()+ "/BidomainLR91_1d_V.dat";
        std::string file2="heart/test/data/BidomainUnpermuted1d/BidomainLR91_1d_V.dat";
        NumericFileComparison comp(file1, file2);
        TS_ASSERT(comp.CompareFiles(3e-3)); //This can be quite flexible since the permutation differences will be quite large
    }

    /*
     * This test is almost identical to TestSimpleBidomain1D
     * and relies on that test generating a h5 file to check against.
     *
     * NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
     * surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
     * the equations have been divided through by the surface-area-to-volume ratio.
     * (Historical reasons...)
     */
    void TestBidomainProblemInTwoHalves()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1dInTwoHalves");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSimulationDuration(1.0);
        bidomain_problem.Solve();
        TS_ASSERT_DELTA(bidomain_problem.GetCurrentTime(), 1.0, 1e-12);

        HeartConfig::Instance()->SetSimulationDuration(2.0);
        bidomain_problem.Solve();
        TS_ASSERT_DELTA(bidomain_problem.GetCurrentTime(), 2.0, 1e-12);

        // Check some voltages
        ReplicatableVector solution_replicated(bidomain_problem.GetSolution());

        double atol=5e-3;

        TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
        TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
        TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
        TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
        TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
        TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
        TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);

        // Compare against saved solution
        for (unsigned index=0; index<solution_replicated.GetSize(); index++)
        {
            TS_ASSERT_DELTA(solution_replicated[index], mSolutionReplicated1d2ms[index], 5e-11);
        }

        // Check output file contains results for the whole simulation and agree with normal test
        TS_ASSERT(CompareFilesViaHdf5DataReader("BidomainSimple1dInTwoHalves", "BidomainLR91_1d", true,
                                                "BidomainSimple1dInOneGo", "BidomainLR91_1d", true));

        // Test that we can keep solving even if the results have been deleted (i.e. by creating a new
        // .h5 file when we realize that there isn't one to extend)
        HeartConfig::Instance()->SetOutputDirectory("BidomainSimple1dInTwoHalves_2");
        OutputFileHandler file_handler("BidomainSimple1dInTwoHalves_2", true);
        FileFinder h5_file = file_handler.FindFile("BidomainLR91_1d.h5");
        TS_ASSERT(!h5_file.Exists());
        HeartConfig::Instance()->SetSimulationDuration(3.0);
        TS_ASSERT_THROWS_NOTHING( bidomain_problem.Solve() );
        TS_ASSERT(h5_file.Exists());
    }

    void TestBidomainProblemWithWriterCache()
    {
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.01);
        HeartConfig::Instance()->SetSimulationDuration(1.0);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainWithWriterCache");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d_with_cache");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        bidomain_problem.SetUseHdf5DataWriterCache(true); // cache on

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        // Doesn't really test that the cache was used (since this test would pass with cache turned off too)...
        TS_ASSERT(CompareFilesViaHdf5DataReader("BidomainWithWriterCache", "BidomainLR91_1d_with_cache", true,
                                                "heart/test/data/BidomainWithWriterCache", "BidomainLR91_1d_with_cache", false,
                                                4e-4));
    }

    void TestBidomainProblemWithWriterCacheIncomplete()
    {
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
        HeartConfig::Instance()->SetOutputDirectory("BidomainWithWriterCacheIncomplete");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d_with_cache_incomplete");
        HeartConfig::Instance()->SetSimulationDuration(1.0);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        bidomain_problem.SetUseHdf5DataWriterCache(true); // cache on

        std::vector<unsigned> nodes_to_be_output;
        nodes_to_be_output.push_back(0);
        nodes_to_be_output.push_back(5);
        nodes_to_be_output.push_back(10);
        bidomain_problem.SetOutputNodes(nodes_to_be_output);

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        TS_ASSERT(CompareFilesViaHdf5DataReader("BidomainWithWriterCacheIncomplete", "BidomainLR91_1d_with_cache_incomplete", true,
                                                "heart/test/data/BidomainWithWriterCache", "BidomainLR91_1d_with_cache_incomplete", false,
                                                4e-4));
    }

    /* Disabled this one for now as it passes in an ugly way */
    void louieTestBidomainProblemWithWriterCacheExtraVarsException()
    {
        HeartConfig::Instance()->SetParametersFile("heart/test/data/xml/MultipleVariablesBidomain.xml");

        std::vector<std::string> output_variables;
        output_variables.push_back("calcium_dynamics__Ca_NSR");
        output_variables.push_back("ionic_concentrations__Nai");
        output_variables.push_back("fast_sodium_current_j_gate__j");
        output_variables.push_back("ionic_concentrations__Ki");
        HeartConfig::Instance()->SetOutputVariables( output_variables );

        PlaneStimulusCellFactory<CellFaberRudy2000FromCellML, 1> cell_factory;
        BidomainProblem<1> bidomain_problem( &cell_factory );
        bidomain_problem.SetUseHdf5DataWriterCache();

        bidomain_problem.Initialise();
        TS_ASSERT_THROWS_THIS(bidomain_problem.Solve(), "Cached writes must write all variables at once.");
    }

    /**
     * Not a very thorough test yet - just checks we can load a problem, simulate it, and
     * get expected results.
     *
     * This test relies on the h5 file generated in TestSimpleBidomain1D. Always run after!
     */
    void TestArchiving()
    {
        FileFinder archive_dir("bidomain_problem_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "bidomain_problem.arch";

        // Values to test against after load
        unsigned num_cells;

        // Save
        {
            HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.0005));
            HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1mm_10_elements");
            HeartConfig::Instance()->SetOutputDirectory("BiProblemArchive");
            HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91_1d");
            HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
            HeartConfig::Instance()->SetCapacitance(1.0);
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

            PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
            BidomainProblem<1> bidomain_problem( &cell_factory );

            bidomain_problem.Initialise();
            HeartConfig::Instance()->SetSimulationDuration(1.0); //ms
            bidomain_problem.Solve();

            num_cells = bidomain_problem.GetTissue()->rGetCellsDistributed().size();

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacProblem<1,1,2>* const p_bidomain_problem = &bidomain_problem;
            (*p_arch) & p_bidomain_problem;
        }

        // Load
        {
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            AbstractCardiacProblem<1,1,2> *p_bidomain_problem;
            (*p_arch) >> p_bidomain_problem;

            // Check values
            TS_ASSERT_EQUALS(p_bidomain_problem->GetTissue()->rGetCellsDistributed().size(),
                             num_cells);

            HeartConfig::Instance()->SetSimulationDuration(2.0); //ms
            p_bidomain_problem->Solve();

            // check some voltages
            ReplicatableVector solution_replicated(p_bidomain_problem->GetSolution());
            double atol=5e-3;
            TS_ASSERT_DELTA(solution_replicated[1], -16.4861, atol);
            TS_ASSERT_DELTA(solution_replicated[2], 22.8117, atol);
            TS_ASSERT_DELTA(solution_replicated[3], -16.4893, atol);
            TS_ASSERT_DELTA(solution_replicated[5], -16.5617, atol);
            TS_ASSERT_DELTA(solution_replicated[7], -16.6761, atol);
            TS_ASSERT_DELTA(solution_replicated[9], -16.8344, atol);
            TS_ASSERT_DELTA(solution_replicated[10], 25.3148, atol);

            for (unsigned index=0; index<solution_replicated.GetSize(); index++)
            {
                //Shouldn't differ from the original run at all
                TS_ASSERT_DELTA(solution_replicated[index], mSolutionReplicated1d2ms[index],  5e-11);
            }
            // check output file contains results for the whole simulation
            TS_ASSERT(CompareFilesViaHdf5DataReader("BiProblemArchive", "BidomainLR91_1d", true,
                                                    "BidomainSimple1dInOneGo", "BidomainLR91_1d", true));

            // Free memory
            delete p_bidomain_problem;
        }
    }
};

#endif /*TESTBIDOMAINPROBLEM_HPP_*/
