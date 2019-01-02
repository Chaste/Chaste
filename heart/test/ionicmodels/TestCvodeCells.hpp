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

#ifndef TESTCVODECELLS_HPP_
#define TESTCVODECELLS_HPP_

#include <iostream>
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveLocationInfo.hpp"

#include "AbstractCvodeCell.hpp"
#include "LuoRudy1991Cvode.hpp"
#include "LuoRudy1991.hpp"
#include "Shannon2004.hpp"
#include "Shannon2004Cvode.hpp"
#include "FixedModifier.hpp"
#include "DummyModifier.hpp"
#include "SmartPointers.hpp"

#include "HeartConfig.hpp"
#include "RegularStimulus.hpp"
#include "SimpleStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RunAndCheckIonicModels.hpp"
#include "OdeSystemInformation.hpp"
#include "VectorHelperFunctions.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#ifdef CHASTE_CVODE

class ExceptionalCell : public AbstractCvodeCell
{
private:
    bool mNice;
public :
    ExceptionalCell(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                    boost::shared_ptr<AbstractStimulusFunction> pStimulus)
        : AbstractCvodeCell(pOdeSolver, 2, 0, pStimulus)
    {
        mNice = false;
        mpSystemInfo = OdeSystemInformation<ExceptionalCell>::Instance();
        Init();
        // For coverage purposes, call Init twice
        Init();
    }

    void BeNice()
    {
        mNice = true;
    }

    void EvaluateYDerivatives(double time, N_Vector y, N_Vector ydot)
    {
        NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
        NV_Ith_S(ydot, 1) = -NV_Ith_S(y, 0);
        if (!mNice)
        {
            EXCEPTION(DumpState("I'm feeling nasty!"));
        }
    }

    double GetIIonic(const std::vector<double>* pStateVariables=NULL)
    {
        return 0.0;
    }
};

template<>
void OdeSystemInformation<ExceptionalCell>::Initialise(void)
{
    // State variables
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    this->mInitialised = true;
}
#endif // CHASTE_CVODE

class TestCvodeCells : public CxxTest::TestSuite
{
public:

    void TestLuoRudyCvodeCell()
    {
#ifdef CHASTE_CVODE
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<ZeroStimulus> p_zero_stimulus(new ZeroStimulus()); // For coverage of SetStimulusFunction()
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude,
                                                                          duration,
                                                                          period,
                                                                          start));
        double start_time = 0.0;
        double end_time = 1000.0; //One second in milliseconds
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Make a model that uses Cvode directly:
        CellLuoRudy1991FromCellMLCvode lr91_cvode_system(p_solver, p_zero_stimulus);

        // Cover swapping to a proper stimulus
        lr91_cvode_system.SetStimulusFunction(p_stimulus);
        // More "coverage" - CVODE cells don't have a solver
        TS_ASSERT(!lr91_cvode_system.GetSolver());

        boost::shared_ptr<AbstractStimulusFunction> p_abs_stim = lr91_cvode_system.GetStimulusFunction();
        double period_back = boost::static_pointer_cast<RegularStimulus>(p_abs_stim)->GetPeriod();
        TS_ASSERT_DELTA(period_back,period,1e-7);
        TS_ASSERT_EQUALS(p_abs_stim,p_stimulus);

        TS_ASSERT_EQUALS(lr91_cvode_system.GetVoltageIndex(), 0u);
        TS_ASSERT_EQUALS(lr91_cvode_system.GetMaxSteps(), 0); // 0 means 'UNSET' and Cvode uses the default.

        // 'Traditional' Chaste cell model for comparison of results:
        CellLuoRudy1991FromCellML lr91_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        double max_timestep = 1.0;
        double sampling_time = max_timestep;

        TS_ASSERT_DELTA(lr91_cvode_system.GetTimestep(), DOUBLE_UNSET, 1e-9);
        OdeSolution solution_cvode = lr91_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time);
        // NB the timestep is not automatically set and recorded if we call Solve directly
        // as this methods is on AbstractCvodeSystem, not AbstractCvodeCell which does store it!
        TS_ASSERT_DELTA(lr91_cvode_system.GetTimestep(), DOUBLE_UNSET, 1e-9);

        OdeSolution solution_chaste = lr91_ode_system.Compute(start_time, end_time);

        TS_ASSERT_DELTA(lr91_ode_system.GetIIonic(), lr91_cvode_system.GetIIonic(), 1e-3);

        unsigned step_per_row_chaste = 100u;
        bool clean_dir = false;
        solution_cvode.WriteToFile("TestCvodeCells","lr91_cvode","ms",1,clean_dir);
        solution_chaste.WriteToFile("TestCvodeCells","lr91_chaste","ms",step_per_row_chaste,clean_dir);

        double tolerance = 2e-1;
        bool voltage_only = false;
        CompareCellModelResults("lr91_cvode", "lr91_chaste",tolerance, voltage_only, "TestCvodeCells");

        // Clamping
        lr91_cvode_system.SetVoltageDerivativeToZero();
        solution_cvode = lr91_cvode_system.Solve(end_time, end_time+100.0, max_timestep, sampling_time);
        std::vector<double> voltages = solution_cvode.GetVariableAtIndex(lr91_cvode_system.GetVoltageIndex());
        for (unsigned i=0; i<voltages.size(); i++)
        {
            TS_ASSERT_EQUALS(voltages[i], lr91_cvode_system.GetVoltage());
        }

        lr91_cvode_system.SetVoltageDerivativeToZero(false);

        // Check Compute* & SolveAndUpdateState methods from AbstractCardiacCellInterface
        lr91_cvode_system.ResetToInitialConditions();
        HeartConfig::Instance()->SetPrintingTimeStep(sampling_time);
        OdeSolution solution_cvode_2 = lr91_cvode_system.Compute(start_time, end_time);
        solution_cvode_2.WriteToFile("TestCvodeCells","lr91_cvode_2","ms",1,clean_dir);
        CompareCellModelResults("lr91_cvode_2", "lr91_cvode", 1e-8, voltage_only, "TestCvodeCells");
        lr91_cvode_system.SetMaxSteps(10000); // Needed since we're not sampling

        lr91_cvode_system.ResetToInitialConditions();
        TS_ASSERT_EQUALS(lr91_cvode_system.GetMaxSteps(), 10000);
        lr91_cvode_system.SetMaxTimestep(DOUBLE_UNSET); // Use default (set from HeartConfig)
        lr91_cvode_system.SolveAndUpdateState(start_time, end_time);
        // The max time step, which was unset, should now be set to printing time step.
        TS_ASSERT_DELTA(lr91_cvode_system.GetTimestep(), HeartConfig::Instance()->GetPrintingTimeStep(), 1e-9);
        TS_ASSERT_DELTA(lr91_cvode_system.GetVoltage(), solution_cvode_2.rGetSolutions().back()[lr91_cvode_system.GetVoltageIndex()], 1e-4);
        // Note: adaptive solve takes different time steps when not sampling => can't use very tight tolerance

        // Reset CVODE cell to initial conditions, and solve without sampling
        lr91_cvode_system.ResetToInitialConditions();
        lr91_cvode_system.Solve(start_time, end_time, max_timestep);
        TS_ASSERT_DELTA(lr91_cvode_system.GetVoltage(), lr91_ode_system.GetVoltage(), 1e-3);

        // Test parameter
        TS_ASSERT_EQUALS(lr91_cvode_system.GetNumberOfParameters(), 3u);
        TS_ASSERT_EQUALS(lr91_cvode_system.rGetParameterNames()[1], "membrane_fast_sodium_current_conductance");
        TS_ASSERT_EQUALS(lr91_cvode_system.rGetParameterUnits()[1], "milliS_per_cm2");
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameterIndex("membrane_fast_sodium_current_conductance"), 1u);
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameterUnits(1u), "milliS_per_cm2");
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameter(1u), 23.0);
        TS_ASSERT_EQUALS(lr91_cvode_system.GetParameter("membrane_fast_sodium_current_conductance"), 23.0);

        // Parameter exceptions
        TS_ASSERT_THROWS_THIS(lr91_cvode_system.GetParameterIndex("b"), "No parameter named 'b'.");
        TS_ASSERT_THROWS_THIS(lr91_cvode_system.GetParameter(3u), "The index passed in must be less than the number of parameters.");
        TS_ASSERT_THROWS_THIS(lr91_cvode_system.GetParameterUnits(3u), "The index passed in must be less than the number of parameters.");


        // Coverage
        boost::shared_ptr<const AbstractOdeSystemInformation> p_sys_info = lr91_cvode_system.GetSystemInformation();
        TS_ASSERT(p_sys_info->rGetStateVariableNames() == lr91_cvode_system.rGetStateVariableNames());
        TS_ASSERT(p_sys_info->rGetStateVariableUnits() == lr91_cvode_system.rGetStateVariableUnits());
        TS_ASSERT_EQUALS(lr91_cvode_system.GetStateVariableIndex("membrane_voltage"),
                         lr91_cvode_system.GetVoltageIndex());
        TS_ASSERT_EQUALS(lr91_cvode_system.GetStateVariable(0),lr91_cvode_system.GetVoltage());
        TS_ASSERT_EQUALS(lr91_cvode_system.GetStateVariable("membrane_voltage"),lr91_cvode_system.GetVoltage());
        TS_ASSERT_EQUALS(lr91_cvode_system.GetStateVariableUnits(0), "millivolt");

        TS_ASSERT_DELTA(lr91_cvode_system.GetRelativeTolerance(), 1e-5, 1e-10);
        TS_ASSERT_DELTA(lr91_cvode_system.GetAbsoluteTolerance(), 1e-7, 1e-10);
        TS_ASSERT_LESS_THAN_EQUALS(lr91_cvode_system.GetLastStepSize(), max_timestep);

        double old_v = lr91_cvode_system.GetVoltage();
        const double new_v = -1000.0;
        lr91_cvode_system.SetVoltage(new_v);
        TS_ASSERT_DELTA(lr91_cvode_system.GetVoltage(), new_v, 1e-12);
        lr91_cvode_system.SetVoltage(old_v);

        lr91_cvode_system.SetStateVariable("membrane_voltage", new_v);
        TS_ASSERT_DELTA(lr91_cvode_system.GetVoltage(), new_v, 1e-12);
        lr91_cvode_system.SetVoltage(old_v);

        std::vector<double> s_state_vars = lr91_cvode_system.GetStdVecStateVariables();
        N_Vector n_state_vars = lr91_cvode_system.GetStateVariables();
        TS_ASSERT_EQUALS(GetVectorSize(n_state_vars), GetVectorSize(s_state_vars));
        TS_ASSERT_EQUALS(GetVectorSize(n_state_vars), lr91_cvode_system.GetNumberOfStateVariables());
        for (unsigned i=0; i<s_state_vars.size(); i++)
        {
            TS_ASSERT_DELTA(GetVectorComponent(n_state_vars,i), s_state_vars[i], 1e-9);
        }
        DeleteVector(n_state_vars); // We need to clean up CVODE vectors.

        s_state_vars[0] = -1000; // Alter one of the variables to check the below method works...
        lr91_cvode_system.SetStateVariables(s_state_vars);
        for (unsigned i=0; i<s_state_vars.size(); i++)
        {
            TS_ASSERT_DELTA(lr91_cvode_system.GetStateVariable(i), s_state_vars[i], 1e-9);
        }
        lr91_cvode_system.SetVoltage(old_v);

        // Cover errors
        std::cout << "Testing Error Handling... don't worry about diasters below if this passes!\n";
        lr91_cvode_system.SetMaxSteps(2);
        TS_ASSERT_THROWS_CONTAINS(lr91_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time),
                              "CVODE failed to solve system: CV_TOO_MUCH_WORK");
        // Kill the cell -- will trigger exp overflow if FPE exceptions are on
        boost::shared_ptr<SimpleStimulus> p_boom_stimulus(new SimpleStimulus(-50000, 2.0, 1.0));
#ifndef TEST_FOR_FPE
        CellLuoRudy1991FromCellMLCvode lr91_boom(p_solver, p_boom_stimulus);
        TS_ASSERT_THROWS_CONTAINS(lr91_boom.Solve(start_time, end_time, max_timestep, sampling_time),
                                  "CVODE failed to solve system: CV_"); // Can be ERR_FAILURE or CONV_FAILURE...
        lr91_boom.ResetToInitialConditions();
        lr91_boom.SetMaxSteps(10000);
        TS_ASSERT_THROWS_CONTAINS(lr91_boom.Solve(start_time, end_time, max_timestep),
                                  "CVODE failed to solve system: CV_"); // Can be ERR_FAILURE or CONV_FAILURE...
#endif // TEST_FOR_FPE
        // Nasty cell
        ExceptionalCell bad_cell(p_solver, p_boom_stimulus);
        TS_ASSERT_THROWS_CONTAINS(bad_cell.Solve(start_time, end_time, max_timestep, sampling_time),
                              "CVODE failed to solve system: CV_RHSFUNC_FAIL");
        bad_cell.ResetToInitialConditions();
        TS_ASSERT_THROWS_CONTAINS(bad_cell.Solve(start_time, end_time, max_timestep),
                              "CVODE failed to solve system: CV_RHSFUNC_FAIL");

        // This should work now that metadata has been added to the LuoRudy1991 cellML.
        TS_ASSERT_EQUALS(lr91_cvode_system.HasCellMLDefaultStimulus(), true);
        boost::shared_ptr<RegularStimulus> p_cellml_stim = lr91_cvode_system.UseCellMLDefaultStimulus();
        TS_ASSERT_DELTA(p_cellml_stim->GetPeriod(), 1000.0, 1e-9);

        {
            lr91_cvode_system.SetMaxSteps(10000);
            double voltage_was = lr91_cvode_system.GetAnyVariable("membrane_voltage");
            lr91_cvode_system.ComputeExceptVoltage(start_time, end_time);
            TS_ASSERT_EQUALS(lr91_cvode_system.GetAnyVariable("membrane_voltage"), voltage_was);
        }
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestShannon2004()
    {
#ifdef CHASTE_CVODE
        // Set stimulus
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude,
                                                                          duration,
                                                                          period,
                                                                          start));

        double ode_time_step = 0.001;
        HeartConfig::Instance()->SetOdeTimeStep(ode_time_step);
        double start_time = 0.0;
        double end_time = 1000.0; //One second in milliseconds
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Make a model that uses Cvode directly:
        CellShannon2004FromCellMLCvode sh04_cvode_system(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(sh04_cvode_system.GetVoltageIndex(), 0u);
        TS_ASSERT_EQUALS(sh04_cvode_system.GetMaxSteps(), 0); // 0 means 'UNSET' and Cvode uses the default.

        // 'Traditional' Chaste cell model for comparison of results:
        CellShannon2004FromCellML sh04_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        double max_timestep = 1.0;
        double sampling_time = 1.0;

        // Take a copy of the state variables, solve once, reset state variables and solve again (should get same answer)
        N_Vector state_vars = sh04_cvode_system.GetStateVariables(); // This returns the initial conditions
        sh04_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time);
        N_Vector state_vars_after_solve = sh04_cvode_system.GetStateVariables(); // This returns the state variables at the end (for coverage)...

        sh04_cvode_system.SetStateVariables(state_vars);
        OdeSolution solution_cvode = sh04_cvode_system.Solve(start_time, end_time, max_timestep, sampling_time);

        // Check we get the right answer both times
        std::vector<double> final_state_vars = solution_cvode.rGetSolutions().back();
        for (unsigned i=0; i<final_state_vars.size(); i++)
        {
            TS_ASSERT_DELTA(final_state_vars[i], NV_Ith_S(state_vars_after_solve, i), 1e-6);
        }

        // Free memory
        DeleteVector(state_vars);
        DeleteVector(state_vars_after_solve);

        /*
         * HOW_TO_TAG Cardiac/Output
         * Calculating and outputting ionic currents ('derived quantities') in a single cell simulation using
         * class:OdeSolution - see also [wiki:ChasteGuides/CodeGenerationFromCellML#Derivedquantities this page].
         */
        // Note that the ODE system must have "derived quantities" defined (and a method for calculating them).
        // This is auto-generated by PyCML if the quantity is tagged appropriately - see CodeGenerationFromCellML wiki page.

        // Solve using Chaste solvers for comparison.
        OdeSolution solution_chaste = sh04_ode_system.Compute(start_time, end_time, sampling_time);

        // These don't seem to need template arguments, I suppose the compiler figures out which to use from the last argument?
        solution_cvode.CalculateDerivedQuantitiesAndParameters(&sh04_cvode_system);
        solution_chaste.CalculateDerivedQuantitiesAndParameters(&sh04_ode_system);

        // Check methods for directly getting any quantity out of an OdeSolution
        std::vector<double> voltages_cvode = solution_cvode.GetAnyVariable("membrane_voltage");
        std::vector<double> voltages_chaste = solution_chaste.GetAnyVariable("membrane_voltage");
        // and a parameter...
        std::vector<double> param_cvode = solution_cvode.GetAnyVariable("membrane_fast_sodium_current_conductance");
        std::vector<double> param_chaste = solution_chaste.GetAnyVariable("membrane_fast_sodium_current_conductance");
        // and a derived quantity
        std::vector<double> derived_quantity_cvode = solution_cvode.GetAnyVariable("membrane_fast_sodium_current");
        std::vector<double> derived_quantity_chaste = solution_chaste.GetAnyVariable("membrane_fast_sodium_current");


        TS_ASSERT_EQUALS(voltages_cvode.size(),voltages_chaste.size());
        TS_ASSERT_EQUALS(param_cvode.size(),param_chaste.size());
        TS_ASSERT_EQUALS(derived_quantity_cvode.size(),derived_quantity_chaste.size());

        double tolerance = 6e-2;
        for (unsigned i=0; i<voltages_cvode.size(); ++i)
        {
            // The tolerances for voltage are adjusted cleverly at upstroke in CompareCellModelResults below...
            TS_ASSERT_DELTA(voltages_cvode[i],voltages_chaste[i], 10*tolerance);
            TS_ASSERT_DELTA(param_cvode[i], param_chaste[i], 1e-9); // These should be very very similar!
            TS_ASSERT_DELTA(param_cvode[i], 16.0, 1e-9);
            TS_ASSERT_DELTA(derived_quantity_chaste[i], derived_quantity_cvode[i], 70*tolerance);
        }

        bool clean_dir = false;
        unsigned precision = 6;
        bool include_derived_quantities = true;

        solution_cvode.WriteToFile( "TestCvodeCells","sh04_cvode", "ms",1,clean_dir,precision,include_derived_quantities);
        solution_chaste.WriteToFile("TestCvodeCells","sh04_chaste","ms",1,clean_dir,precision,include_derived_quantities);

        bool voltage_only = false;
        CompareCellModelResults("sh04_cvode", "sh04_chaste", tolerance, voltage_only, "TestCvodeCells");

        // Coverage of GetIIonic method.
        TS_ASSERT_DELTA(sh04_ode_system.GetIIonic(), sh04_cvode_system.GetIIonic(), 1e-4);
        TS_ASSERT_DELTA(sh04_cvode_system.GetIIonic(), 0.0004, 1e-4);

        // Clamping V
        sh04_cvode_system.SetVoltageDerivativeToZero();
        solution_cvode = sh04_cvode_system.Solve(end_time, end_time+100.0, max_timestep, sampling_time);
        std::vector<double> voltages = solution_cvode.GetVariableAtIndex(sh04_cvode_system.GetVoltageIndex());
        for (unsigned i=0; i<voltages.size(); i++)
        {
            TS_ASSERT_EQUALS(voltages[i], sh04_cvode_system.GetVoltage());
        }

        sh04_cvode_system.SetVoltageDerivativeToZero(false);

        // Clamp V with a modifier
        sh04_cvode_system.ResetToInitialConditions();
        MAKE_PTR_ABS(AbstractModifier, FixedModifier, p_modifier, (sh04_cvode_system.GetVoltage()));
        sh04_cvode_system.SetModifier("membrane_voltage", p_modifier);
        solution_cvode = sh04_cvode_system.Solve(0.0, start+10.0, max_timestep, sampling_time);
        // The behaviour isn't quite the same - dV/dt isn't clamped in this case - but we should
        // see a reduced response.
        voltages = solution_cvode.GetVariableAtIndex(sh04_cvode_system.GetVoltageIndex());
        double test_v = solution_chaste.GetVariableAtIndex(sh04_cvode_system.GetVoltageIndex())[voltages.size()-1];
        for (unsigned i=0; i<voltages.size(); i++)
        {
            TS_ASSERT_LESS_THAN(voltages[i], test_v);
        }

        // Coverage of mSetVoltageDerivativeToZero in non-CVODE class
        sh04_ode_system.ComputeExceptVoltage(end_time,end_time+0.01);

        // Blocking membrane_fast_sodium_current_conductance (a parameter), but using modifiers
        sh04_cvode_system.ResetToInitialConditions();
        ASSIGN_PTR(p_modifier, FixedModifier, (0.0));
        sh04_cvode_system.SetModifier("membrane_fast_sodium_current_conductance", p_modifier);
        ASSIGN_PTR(p_modifier, DummyModifier, ());
        sh04_cvode_system.SetModifier("membrane_voltage", p_modifier);
        OdeSolution solution_block = sh04_cvode_system.Solve(0.0, 100.0, max_timestep, sampling_time);
        solution_block.WriteToFile("TestCvodeCells","sh04_block_modifier","ms",1,clean_dir);
        // Now set the parameter instead
        sh04_cvode_system.ResetToInitialConditions();
        sh04_cvode_system.SetModifier("membrane_fast_sodium_current_conductance", p_modifier);
        sh04_cvode_system.SetParameter("membrane_fast_sodium_current_conductance", 0.0);
        solution_block = sh04_cvode_system.Solve(0.0, 100.0, max_timestep, sampling_time);
        solution_block.WriteToFile("TestCvodeCells","sh04_block_param","ms",1,clean_dir);
        CompareCellModelResults("sh04_block_param", "sh04_block_modifier", 1e-6, voltage_only, "TestCvodeCells");
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }

    void TestArchivingCvodeCells()
    {
#ifdef CHASTE_CVODE
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("lr91_cvode.arch");

        const double magnitude_stimulus = -3;  // uA/cm2
        const double duration_stimulus = 3;  // ms
        const double period_stimulus = 1000; //ms
        const double start_stimulus = 10.0;   // ms

        // Save
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set stimulus
            boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude_stimulus,
                                                                              duration_stimulus,
                                                                              period_stimulus,
                                                                              start_stimulus));
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;

            // Check Standard
            AbstractCardiacCellInterface* const p_luo_rudy_cell = new CellLuoRudy1991FromCellMLCvode(p_solver, p_stimulus);

            output_arch <<  p_luo_rudy_cell;

//            // These results are in the repository and should be replicated after the load below
//            RunOdeSolverWithIonicModel(p_luo_rudy_cell,
//                           50.0,
//                           "LRCVODEAfterArchiveValidData");

            delete p_luo_rudy_cell;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCellInterface* p_luo_rudy_cell;
            input_arch >> p_luo_rudy_cell;

            TS_ASSERT_EQUALS( p_luo_rudy_cell->GetNumberOfStateVariables(), 8U );
            boost::shared_ptr<RegularStimulus> p_stimulus = boost::static_pointer_cast<RegularStimulus>(p_luo_rudy_cell->GetStimulusFunction());

            TS_ASSERT_DELTA(p_stimulus->GetPeriod(),   period_stimulus,   1e-12);
            TS_ASSERT_DELTA(p_stimulus->GetDuration(), duration_stimulus, 1e-12);
            TS_ASSERT_DELTA(p_stimulus->GetMagnitude(),magnitude_stimulus,1e-12);
            TS_ASSERT_DELTA(p_stimulus->GetStartTime(),start_stimulus,    1e-12);

            RunOdeSolverWithIonicModel(p_luo_rudy_cell,
                                       50.0,
                                       "LR91CvodeAfterArchive");

            CheckCellModelResults("LR91CvodeAfterArchive");

            delete p_luo_rudy_cell;
        }
#else
        std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
    }
};


#endif /*TESTCVODECELLS_HPP_*/
