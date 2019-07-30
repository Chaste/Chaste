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

#ifndef TEST1DMONODOMAINSHANNONCVODEBENCHMARKS_HPP_
#define TEST1DMONODOMAINSHANNONCVODEBENCHMARKS_HPP_

#include <cxxtest/TestSuite.h>
#include <sstream>
#include <iostream>

#include "TetrahedralMesh.hpp"
#include "MonodomainProblem.hpp"
#include "RegularStimulus.hpp"
#include "Shannon2004.hpp"
#include "Shannon2004Cvode.hpp"
//#include "Shannon2004BackwardEuler.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "HeartConfig.hpp"
#include "CvodeAdaptor.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "Timer.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * Cell Factory defining cells for 1d chain.
 *
 * This gives a shared default solver to all the cells, as standard.
 */
template <class CELL_MODEL>
class ShannonCardiacCellFactory : public AbstractCardiacCellFactory<1>
{
private:

    boost::shared_ptr<RegularStimulus> mpStimulus;

public:
    ShannonCardiacCellFactory()
    : AbstractCardiacCellFactory<1>(),
      mpStimulus(new RegularStimulus(-250000,5,2000,1))
      {
      }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        CELL_MODEL *p_cell;

        if (pNode->GetPoint()[0] == 0.0)
        {
            p_cell = new CELL_MODEL(this->mpSolver,
                                    mpStimulus);
        }
        else
        {
            p_cell = new CELL_MODEL(this->mpSolver,
                                    this->mpZeroStimulus);
        }

        return p_cell;
    }
};

#ifdef CHASTE_CVODE

/**
 *
 * Cell Factory defining Shannon cells for 1d chain.
 *
 * In this version each node has its own Cvode adaptor solver
 */
class ShannonCvodeAdaptorCellFactory : public AbstractCardiacCellFactory<1>
{
private:

    boost::shared_ptr<RegularStimulus> mpStimulus;

public:
    ShannonCvodeAdaptorCellFactory()
    : AbstractCardiacCellFactory<1>(),
      mpStimulus(new RegularStimulus(-250000,5,2000,1))
      {
      }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        AbstractCardiacCell* p_cell;

        /*
         * HOW_TO_TAG Cardiac/Problem definition
         * Use a CVODE adaptor solver in a tissue simulation
         */

        // Here we add a new Cvode Adaptor solver for every node.
        boost::shared_ptr<CvodeAdaptor> p_cvode_solver(new CvodeAdaptor());
        p_cvode_solver->SetMinimalReset(true);
        p_cvode_solver->SetTolerances(1e-4,1e-6);// NB These defaulted to different values in AbstractCvodeSystem.
        if (pNode->GetPoint()[0] == 0.0)
        {
            p_cell = new CellShannon2004FromCellML(p_cvode_solver,
                                                   mpStimulus);
        }
        else
        {
            p_cell = new CellShannon2004FromCellML(p_cvode_solver,
                                                   this->mpZeroStimulus);
        }

        return p_cell;
    }
};

/**
 *
 * Cell Factory defining Shannon cells for 1d chain.
 *
 * In this version each node has its own native Cvode solver
 */
class ShannonCvodeNativeCellFactory : public AbstractCardiacCellFactory<1>
{
private:

    boost::shared_ptr<RegularStimulus> mpStimulus;

public:
    ShannonCvodeNativeCellFactory()
    : AbstractCardiacCellFactory<1>(),
      mpStimulus(new RegularStimulus(-250000,5,2000,1))
      {
      }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        AbstractCvodeCell* p_cell;

        /*
         * HOW_TO_TAG Cardiac/Problem definition
         * Use a native CVODE cell in a tissue simulation
         */

        // Here we add a native Cvode cell (which includes its own solver) at every node.
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;

        if (pNode->GetPoint()[0] == 0.0)
        {
            p_cell = new CellShannon2004FromCellMLCvode(p_empty_solver,
                                                        mpStimulus);
        }
        else
        {
            p_cell = new CellShannon2004FromCellMLCvode(p_empty_solver,
                                                        this->mpZeroStimulus);
        }
        // p_cell->SetMinimalReset(true); // This is now done by AbstractCardiacCellFactory for any AbstractCvodeCell.
        p_cell->SetTolerances(1e-4,1e-6); // NB These defaulted to different values in CvodeAdaptor.
        return p_cell;
    }
};
#endif // CHASTE_CVODE

/**
 * This class tests that a CVODE tissue simulation gets comparable results to those
 * of Forward and Backward Euler.
 *
 * At present FE and BE are more different than FE and CVODE.
 */
class Test1dMonodomainShannonCvodeBenchmarks : public CxxTest::TestSuite
{
private:
    bool CompareBenchmarkResults(const std::vector<double>& rReferenceTrace, const std::vector<double>& rTestTrace, double tol)
    {
        assert(rReferenceTrace.size()==rTestTrace.size());

        // See what the maximum difference in voltages recorded at final time at this node is
        unsigned last_time = rReferenceTrace.size() - 1u;
        double difference = fabs(rReferenceTrace[last_time] - rTestTrace[last_time]);

        // I was looping over all of the time trace here, but conduction velocity must
        // be slightly affected by the different solvers, giving rise to massive
        // (height of upstroke ~135mV) differences, so test at 'plateau' instead.

        // Return a false if it is too big.
        bool result = true;
        if (difference > tol)
        {
            std::cout << "\nDifference in recorded voltages = " << difference << std::endl << std::flush;
            result = false;
        }

        // Also return a false if the cell didn't get stimulated (assuming test only runs to plateau).
        if (rTestTrace[last_time] < 0)
        {
            std::cout << "Cell is not depolarized\n" << std::endl << std::flush;
            result = false;
        }

        return result;
    }

public:

    void TestWithDifferentCellsAndSolvers()
    {
        double duration = 25;
        HeartConfig::Instance()->SetSimulationDuration(duration); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        double pde_time_step = 0.1; //ms
        double printing_time_step = 0.1; //ms

        std::vector<double> fe_node_0;
        std::vector<double> fe_node_100;
        std::vector<double> times;

        {
            HeartConfig::Instance()->SetOutputDirectory("ShannonBenchmark/forward_euler");
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0025,pde_time_step,printing_time_step);
            ShannonCardiacCellFactory<CellShannon2004FromCellML> cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();

            std::stringstream message;
            message << "1. Forward Euler for " << duration << " ms took";

            Timer::Reset();
            monodomain_problem.Solve();
            Timer::Print(message.str());

            Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
            times = data_reader.GetUnlimitedDimensionValues();
            fe_node_0 = data_reader.GetVariableOverTime("V", 0);
            fe_node_100 = data_reader.GetVariableOverTime("V", 100);

        }

        /*
         * To test backward Euler you need to generate an extra Chaste variant of the model from the CellML.
         *
         * i.e. go to heart/src/odes/cellml/Shannon2004-conf.xml and add
         * <arg>--backward-euler</arg>
         *
         * then uncomment the //#include "Shannon2004BackwardEuler.hpp" and the below block of code.
         *
         * Same idea for checking optimised models.
         */

//        {
//            HeartConfig::Instance()->SetOutputDirectory("ShannonBenchmark/backward_euler");
//            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(pde_time_step,pde_time_step,printing_time_step);
//            ShannonCardiacCellFactory<CellShannon2004FromCellMLBackwardEuler> cell_factory;
//            MonodomainProblem<1> monodomain_problem( &cell_factory );
//
//            monodomain_problem.Initialise();
//
//            std::stringstream message;
//            message << "2. Backward Euler for " << duration << " ms took";
//
//            Timer::Reset();
//            monodomain_problem.Solve();
//            Timer::Print(message.str());
//
//            Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
//            std::vector<double> be_node_0 = data_reader.GetVariableOverTime("V", 0);
//            std::vector<double> be_node_100 = data_reader.GetVariableOverTime("V", 100);
//
//            TS_ASSERT(CompareBenchmarkResults(fe_node_0, be_node_0, 1.1));     // More than a milliVolt of difference with B.E.
//            TS_ASSERT(CompareBenchmarkResults(fe_node_100, be_node_100, 1.1)); // More than a milliVolt of difference with B.E.
//        }

#ifdef CHASTE_CVODE
        // First test a CVODE adaptor
        // The ODE system remains the same (std::vectors) and standard vectors are converted into N_Vectors
        // every time CVODE talks to the ODE system.
        {
            HeartConfig::Instance()->SetOutputDirectory("ShannonBenchmark/cvode_adaptor");
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(pde_time_step,pde_time_step,printing_time_step);
            ShannonCvodeAdaptorCellFactory cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();

            std::stringstream message;
            message << "2. CVODE adaptor for " << duration << " ms took";
            Timer::Reset();
            monodomain_problem.Solve();
            Timer::Print(message.str());

            Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
            std::vector<double> cvode_node_0 = data_reader.GetVariableOverTime("V", 0);
            std::vector<double> cvode_node_100 = data_reader.GetVariableOverTime("V", 100);

            TS_ASSERT(CompareBenchmarkResults(fe_node_0, cvode_node_0, 0.4));     // Only 0.2 -- 0.4mV difference with CVODE
            TS_ASSERT(CompareBenchmarkResults(fe_node_100, cvode_node_100, 0.4)); // Only 0.2 -- 0.4mV difference with CVODE
        }


        // Now test a native CVODE cell
        // The ODE system uses N_Vectors, and interacts with CVODE optimally.
        // When Chaste needs to talk to the system, it converts to std::vector.
        {
            HeartConfig::Instance()->SetOutputDirectory("ShannonBenchmark/cvode_native");
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(pde_time_step,pde_time_step,printing_time_step);
            ShannonCvodeNativeCellFactory cell_factory;
            MonodomainProblem<1> monodomain_problem( &cell_factory );

            monodomain_problem.Initialise();

            std::stringstream message;
            message << "3. CVODE native for " << duration << " ms took";

            Timer::Reset();
            monodomain_problem.Solve();
            Timer::Print(message.str());

            Hdf5DataReader data_reader = monodomain_problem.GetDataReader();
            std::vector<double> cvode_node_0 = data_reader.GetVariableOverTime("V", 0);
            std::vector<double> cvode_node_100 = data_reader.GetVariableOverTime("V", 100);

            TS_ASSERT(CompareBenchmarkResults(fe_node_0, cvode_node_0, 0.4));     // Only 0.2 -- 0.4mV difference with CVODE
            TS_ASSERT(CompareBenchmarkResults(fe_node_100, cvode_node_100, 0.4)); // Only 0.2 -- 0.4mV difference with CVODE
        }
#endif // CHASTE_CVODE
    }
};

#endif /*TEST1DMONODOMAINSHANNONCVODEBENCHMARKS_HPP_*/
