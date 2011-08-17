#ifndef TESTDELTANOTCHODESYSTEM_HPP_
#define TESTDELTANOTCHODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <ctime>
#include <vector>
#include <iostream>

#include "OutputFileHandler.hpp"
#include "DeltaNotchOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CvodeAdaptor.hpp"

class TestDeltaNotchOdeSystem : public CxxTest::TestSuite
{
public:

    void TestDeltaNotchOdeSystemSetup() throw(Exception)
    {
    	EXIT_IF_PARALLEL;
#ifdef CHASTE_CVODE
        double mean_delta = 0.5;
        DeltaNotchOdeSystem ode_system(mean_delta);

        double h_value = 0.0001;
        CvodeAdaptor cvode_solver;
        OdeSolution solutions;

        std::vector<double> initial_conditions = ode_system.GetInitialConditions();

        double start_time = std::clock();
        solutions = cvode_solver.Solve(&ode_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        double end_time = std::clock();
        double elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout <<  "1. Cvode Elapsed time = " << elapsed_time << " secs for 100 hours\n";

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time...(going into S phase at 7.8 hours)
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Decent results
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 0.9615, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 0.0107, 1e-4);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2], mean_delta, 1e-4);
#endif //CHASTE_CVODE
    }

    void TestArchiving()
    {
    	EXIT_IF_PARALLEL;
#ifdef CHASTE_CVODE
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_ode.arch";

        {
            DeltaNotchOdeSystem ode_system;

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 0.5, 1e-6);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 3u);
            TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 0.5, 1e-6);

            // Tidy up
            delete p_ode_system;
        }
#endif //CHASTE_CVODE
    }

    void TestSetStateVariables()
    {
    	EXIT_IF_PARALLEL;
#ifdef CHASTE_CVODE

    	std::vector<double> state_vars;
    	state_vars.push_back(0.0);
    	state_vars.push_back(1.0);
    	state_vars.push_back(2.0);
    	DeltaNotchOdeSystem ode_system(0.5,state_vars);

    	TS_ASSERT_EQUALS(ode_system.GetStateVariable(0),0.0);
    	TS_ASSERT_EQUALS(ode_system.GetStateVariable(1),1.0);
    	TS_ASSERT_EQUALS(ode_system.GetStateVariable(2),2.0);

#endif //CHASTE_CVODE
   }
};

#endif /*TESTDELTANOTCHODESYSTEM_HPP_*/
