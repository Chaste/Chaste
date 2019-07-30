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


#ifndef _TESTIONICMODELSLONG_HPP_
#define _TESTIONICMODELSLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "CheckpointArchiveTypes.hpp" // Needed to avoid segfault on older Boosts...

#include "RunAndCheckIonicModels.hpp"

#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"

#include "EulerIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "FoxModel2002.hpp"
#include "FoxModel2002BackwardEuler.hpp"
#include "Maleckar2008.hpp"
#include "Mahajan2008.hpp"
#include "TenTusscher2006Epi.hpp"
#include "CellProperties.hpp"
#include "HeartConfig.hpp"

#include "FakePetscSetup.hpp"

// Note: RunOdeSolverWithIonicModel(), CheckCellModelResults(), CompareCellModelResults()
// are defined in RunAndCheckIonicModels.hpp

class TestIonicModelsLong : public CxxTest::TestSuite
{
public:
    void TestOdeSolverForFox2002WithRegularStimulus(void)
    {
        clock_t ck_start, ck_end;
        // Set stimulus
        double magnitude = -80.0;
        double duration  = 1.0 ;  // ms
        double start = 50.0; // ms
        double period = 500; // ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude, duration, period, start));

        double end_time = 1000.0; //One second in milliseconds


        HeartConfig::Instance()->SetOdeTimeStep(0.002); // 0.005 leads to NaNs.

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        CellFoxModel2002FromCellML fox_ode_system(p_solver, p_stimulus);

        // Solve and write to file
        ck_start = clock();
        RunOdeSolverWithIonicModel(&fox_ode_system,
                                   end_time,
                                   "FoxRegularStimLong",
                                   500);
        ck_end = clock();
        double forward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        CheckCellModelResults("FoxRegularStimLong");

        // Solve using Backward Euler
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        CellFoxModel2002FromCellMLBackwardEuler backward_system(p_solver, p_stimulus);
        ck_start = clock();
        RunOdeSolverWithIonicModel(&backward_system,
                                   end_time,
                                   "BackwardFoxRegularStimLong",
                                   100);
        ck_end = clock();
        double backward = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;

        CompareCellModelResults("FoxRegularStimLong", "BackwardFoxRegularStimLong", 0.15);
        // Mainly for coverage, and to test consistency of GetIIonic
        TS_ASSERT_DELTA(fox_ode_system.GetIIonic(),
                        backward_system.GetIIonic(),
                        1e-6);

        std::cout << "Run times:\n\tForward: " << forward
                  << "\n\tBackward: " << backward
                  << std::endl;

    }

    /**
     *  Here we test the scale factors methiods for the mahajan model.
     *  The idea is to set the scale factors for the 3 different cell types (epi, mid and endo)
     *  and check that the rsulting APD makes sense if compared to experiemntal results
     */
    void TestScaleFactorsForMahajanModel(void)
    {
        double end_time=300;
        double time_step=0.01;
        double sampling_time=time_step;

        // Set stimulus
        double magnitude_stimulus = -70.0;   // pA/pF
        double duration_stimulus = 1.0;  // ms
        double start_stimulus = 10.0;   // ms
        double period=1000;//here, this is ms
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(magnitude_stimulus,
                                                                          duration_stimulus,
                                                                          period,
                                                                          start_stimulus));

        boost::shared_ptr<EulerIvpOdeSolver> p_forward_solver(new EulerIvpOdeSolver); //define the solver

        CellMahajan2008FromCellML forward_model(p_forward_solver, p_stimulus);
        boost::shared_ptr<BackwardEulerIvpOdeSolver> p_backward_solver(new BackwardEulerIvpOdeSolver(
                                    forward_model.GetNumberOfStateVariables()));

        CellMahajan2008FromCellML epicardial_model(p_backward_solver, p_stimulus);
        CellMahajan2008FromCellML endocardial_model(p_backward_solver, p_stimulus);
        CellMahajan2008FromCellML midmyocardial_model(p_backward_solver, p_stimulus);

        epicardial_model.SetParameter("ScaleFactorGks",1.0);
        epicardial_model.SetParameter("ScaleFactorIto",1.0);
        epicardial_model.SetParameter("ScaleFactorGkr",1.0);

        midmyocardial_model.SetParameter("ScaleFactorGks",0.09);
        midmyocardial_model.SetParameter("ScaleFactorIto",1.0);
        midmyocardial_model.SetParameter("ScaleFactorGkr",1.0);

        endocardial_model.SetParameter("ScaleFactorGks",0.86);
        endocardial_model.SetParameter("ScaleFactorIto",0.2);
        endocardial_model.SetParameter("ScaleFactorGkr",1.0);

        std::vector<double> state_variables_epi = epicardial_model.GetInitialConditions();
        std::vector<double> state_variables_endo = endocardial_model.GetInitialConditions();
        std::vector<double> state_variables_mid = midmyocardial_model.GetInitialConditions();

        const std::string mahajan_epi_file = "Mahajan_epi";
        const std::string mahajan_mid_file = "Mahajan_mid";
        const std::string mahajan_endo_file = "Mahajan_endo";

        // Solve and write to file

        OdeSolution epi_solution;
        epi_solution = p_backward_solver->Solve(&epicardial_model, state_variables_epi, 0, end_time, time_step, sampling_time);

        epi_solution.WriteToFile("TestIonicModels",
                                 mahajan_epi_file,
                                 "ms",//time units
                                 100,//steps per row
                                 false);/*true cleans the directory*/

        OdeSolution mid_solution;
        mid_solution = p_backward_solver->Solve(&midmyocardial_model, state_variables_mid, 0, end_time, time_step, sampling_time);

        mid_solution.WriteToFile("TestIonicModels",
                                 mahajan_mid_file,
                                 "ms",//time units
                                 100,//steps per row
                                 false);/*true cleans the directory*/

        OdeSolution endo_solution;
        endo_solution = p_backward_solver->Solve(&endocardial_model, state_variables_endo, 0, end_time, time_step, sampling_time);

        endo_solution.WriteToFile("TestIonicModels",
                                  mahajan_endo_file,
                                  "ms",//time units
                                  100,//steps per row
                                  false);/*true cleans the directory*/


        ColumnDataReader data_reader_epi("TestIonicModels", mahajan_epi_file);
        ColumnDataReader data_reader_mid("TestIonicModels", mahajan_mid_file);
        ColumnDataReader data_reader_endo("TestIonicModels", mahajan_endo_file);

        std::vector<double> times = data_reader_epi.GetValues("Time");
        std::vector<double> v_endo = data_reader_endo.GetValues("membrane_voltage");
        std::vector<double> v_epi = data_reader_epi.GetValues("membrane_voltage");
        std::vector<double> v_mid = data_reader_mid.GetValues("membrane_voltage");

        CellProperties  cell_properties_endo(v_endo, times);
        CellProperties  cell_properties_epi(v_epi, times);
        CellProperties  cell_properties_mid(v_mid, times);

        double epi_APD = cell_properties_epi.GetLastActionPotentialDuration(90);
        double endo_APD = cell_properties_endo.GetLastActionPotentialDuration(90);
        double mid_APD = cell_properties_mid.GetLastActionPotentialDuration(90);

        // Check that percentage increase from epi to mid and endo (roughly*) matches results
        // from McIntosh et al. Card Res, 45:397-409. 200 (Figure 1 and 2)
        // *this is cardiac modelling after all...
        TS_ASSERT_DELTA((mid_APD-epi_APD)*100/epi_APD, 48.6, 2); // new values because gtos and gtof were the wrong way round (in Mahajan paper - they copied and pasted from Shannon wrong!)
        TS_ASSERT_DELTA((endo_APD-epi_APD)*100/epi_APD, 15.7, 2); // ""
    }

    /**
     * Here we test that the scale factors for the TT model do what they are expected to
     * We check that they modify APD in a way that is expected.
     */
    void TestScaleFactorsForTT06(void)
    {
        double simulation_end=500;/*end time, in milliseconds for this model*/

        // Set the stimulus, the following values are appropriate for single cell simulations of this model.
        double magnitude = -38.0;   // pA/pF
        double duration = 1.0;  // ms
        double start = 100;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude,
                                                                        duration,
                                                                        start));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);// with Forward Euler, this must be as small as 0.001.


        const std::string control_file = "TT_epi";
        const std::string mid_file = "TT_mid";
        const std::string endo_file = "TT_endo";
        const std::string LQT_file = "TT_LQT";

        CellTenTusscher2006EpiFromCellML TT_model_epi(p_solver, p_stimulus);

        TT_model_epi.SetParameter("ScaleFactorIto", 1.0);
        TT_model_epi.SetParameter("ScaleFactorGkr", 1.0);
        TT_model_epi.SetParameter("ScaleFactorGks", 1.0);
        //run the model
        RunOdeSolverWithIonicModel(&TT_model_epi,
                                   simulation_end,
                                   control_file,
                                   100,
                                   false);

        CellTenTusscher2006EpiFromCellML TT_model_mid(p_solver, p_stimulus);
        TT_model_mid.SetParameter("ScaleFactorIto", 1.0);
        TT_model_mid.SetParameter("ScaleFactorGkr", 1.0);
        TT_model_mid.SetParameter("ScaleFactorGks", 0.25);

        RunOdeSolverWithIonicModel(&TT_model_mid,
                                   simulation_end,
                                   mid_file,
                                   100,
                                   false);

        CellTenTusscher2006EpiFromCellML TT_model_endo(p_solver, p_stimulus);
        TT_model_endo.SetParameter("ScaleFactorIto", 0.165);
        TT_model_endo.SetParameter("ScaleFactorGkr", 1.0);
        TT_model_endo.SetParameter("ScaleFactorGks", 0.66);

        RunOdeSolverWithIonicModel(&TT_model_endo,
                                   simulation_end,
                                   endo_file,
                                   100,
                                   false);

        CellTenTusscher2006EpiFromCellML TT_model_LQT(p_solver, p_stimulus);
        TT_model_LQT.SetParameter("ScaleFactorIto", 1.0);
        TT_model_LQT.SetParameter("ScaleFactorGkr", 0.0);
        TT_model_LQT.SetParameter("ScaleFactorGks", 1.0);

        RunOdeSolverWithIonicModel(&TT_model_LQT,
                                   simulation_end,
                                   LQT_file,
                                   100,
                                   false);

        ColumnDataReader data_reader1("TestIonicModels", control_file);
        std::vector<double> voltages1 = GetVoltages(data_reader1);
        ColumnDataReader data_reader2("TestIonicModels", mid_file);
        std::vector<double> voltages2 = GetVoltages(data_reader2);
        ColumnDataReader data_reader3("TestIonicModels", endo_file);
        std::vector<double> voltages3 = GetVoltages(data_reader3);
        ColumnDataReader data_reader4("TestIonicModels", LQT_file);
        std::vector<double> voltages4 = GetVoltages(data_reader4);

        TS_ASSERT_EQUALS(voltages1.size(), voltages2.size());
        TS_ASSERT_EQUALS(voltages2.size(), voltages3.size());
        TS_ASSERT_EQUALS(voltages3.size(), voltages4.size());

        //create the times vector
        std::vector<double> times;
        double k =0;
        for (unsigned i=0; i<voltages2.size(); i++)
        {
          times.push_back(k);
          k=k+0.1;
        }

        CellProperties  cell_properties_control(voltages1, times);
        CellProperties  cell_properties_mid(voltages2, times);
        CellProperties  cell_properties_endo(voltages3, times);
        CellProperties  cell_properties_LQT(voltages4, times);

        double control_APD = cell_properties_control.GetLastActionPotentialDuration(90);
        double mid_APD = cell_properties_mid.GetLastActionPotentialDuration(90);
        double endo_APD = cell_properties_endo.GetLastActionPotentialDuration(90);
        double LQT_APD = cell_properties_LQT.GetLastActionPotentialDuration(90);

        TS_ASSERT_DELTA(control_APD, 300.4789, 0.1);
        TS_ASSERT_DELTA(mid_APD, 392.1871, 0.1);
        TS_ASSERT_DELTA(endo_APD, 329.2048, 0.1);
        TS_ASSERT_DELTA(LQT_APD , 347.8374, 0.1);
     }

    void TestScaleFactorsMaleckar(void)
    {
        double end_time =500;
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-5.6,
                                                                          6,
                                                                          1000,
                                                                          4.0));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver); //define the solver
        HeartConfig::Instance()->SetOdeTimeStep(0.001);
        CellMaleckar2008FromCellML atrial_ode_system(p_solver, p_stimulus);

        const std::string control_file = "control";
        const std::string first_set_file = "first_scale_factor_set";
        const std::string second_set_file = "second_scale_factor_set";
        const std::string AZD_file = "AZD_scale_factor_set";

        OdeSolution control_solution;
        OdeSolution first_scale_factor_set_solution;
        OdeSolution second_scale_factor_set_solution;
        OdeSolution AZD_scale_factor_set_solution;

        double time_step=0.001;
        double sampling_time=0.001;
        std::vector<double> state_variables= atrial_ode_system.GetInitialConditions();

        //default values
        atrial_ode_system.SetParameter("ScaleFactorGks",1.0);
        atrial_ode_system.SetParameter("ScaleFactorIto",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGkr",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGna",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAch",1e-24);
        atrial_ode_system.SetParameter("ScaleFactorGNaK",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGNaCa",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGKur",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGK1",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGCaL",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAZD",0.0);

        control_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);
        control_solution.WriteToFile("TestIonicModels",
                              control_file,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/

        //now apply the first scale factor set, decreases outward currents
        atrial_ode_system.SetParameter("ScaleFactorGks",0.8);
        atrial_ode_system.SetParameter("ScaleFactorIto",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGkr",0.9);
        atrial_ode_system.SetParameter("ScaleFactorGna",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAch",1e-24);
        atrial_ode_system.SetParameter("ScaleFactorGNaK",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGNaCa",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGKur",0.7);
        atrial_ode_system.SetParameter("ScaleFactorGK1",0.6);
        atrial_ode_system.SetParameter("ScaleFactorGCaL",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAZD",0.0);

        state_variables= atrial_ode_system.GetInitialConditions();
        first_scale_factor_set_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);
        first_scale_factor_set_solution.WriteToFile("TestIonicModels",
                              first_set_file,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/

        //now apply the secondscale factor set, this one increases inward currents
        atrial_ode_system.SetParameter("ScaleFactorGks",1.0);
        atrial_ode_system.SetParameter("ScaleFactorIto",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGkr",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGna",1.5);
        atrial_ode_system.SetParameter("ScaleFactorAch",1e-24);
        atrial_ode_system.SetParameter("ScaleFactorGNaK",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGNaCa",2.0);
        atrial_ode_system.SetParameter("ScaleFactorGKur",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGK1",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGCaL",1.6);
        atrial_ode_system.SetParameter("ScaleFactorAZD",0.0);

        state_variables= atrial_ode_system.GetInitialConditions();
        second_scale_factor_set_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);
        second_scale_factor_set_solution.WriteToFile("TestIonicModels",
                              second_set_file,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/

        //check the AZD scale factor (vs control)
        atrial_ode_system.SetParameter("ScaleFactorGks",1.0);
        atrial_ode_system.SetParameter("ScaleFactorIto",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGkr",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGna",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAch",1e-24);
        atrial_ode_system.SetParameter("ScaleFactorGNaK",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGNaCa",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGKur",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGK1",1.0);
        atrial_ode_system.SetParameter("ScaleFactorGCaL",1.0);
        atrial_ode_system.SetParameter("ScaleFactorAZD",5.0);

        AZD_scale_factor_set_solution = p_solver->Solve(&atrial_ode_system, state_variables, 0, end_time, time_step, sampling_time);
        AZD_scale_factor_set_solution.WriteToFile("TestIonicModels",
                              AZD_file,
                              "ms",//time units
                              100,//steps per row
                              false);/*true cleans the directory*/


        ColumnDataReader data_reader1("TestIonicModels", control_file);
        std::vector<double> voltages1 = GetVoltages(data_reader1);
        ColumnDataReader data_reader2("TestIonicModels", first_set_file);
        std::vector<double> voltages2 = GetVoltages(data_reader2);
        ColumnDataReader data_reader3("TestIonicModels", second_set_file);
        std::vector<double> voltages3 = GetVoltages(data_reader3);
        ColumnDataReader data_reader4("TestIonicModels", AZD_file);
        std::vector<double> voltages4 = GetVoltages(data_reader4);

        TS_ASSERT_EQUALS(voltages1.size(), voltages2.size());
        TS_ASSERT_EQUALS(voltages2.size(), voltages3.size());
        TS_ASSERT_EQUALS(voltages3.size(), voltages4.size());

        //create the times vector
        std::vector<double> times;
        double k =0;
        for (unsigned i=0; i<voltages2.size(); i++)
        {
          times.push_back(k);
          k=k+0.1;
        }

        CellProperties cell_properties_control(voltages1, times);
        CellProperties cell_properties_first(voltages2, times);
        CellProperties cell_properties_second(voltages3, times);
        CellProperties cell_properties_AZD(voltages4, times);

        double control_APD = cell_properties_control.GetLastActionPotentialDuration(90);
        double first_APD = cell_properties_first.GetLastActionPotentialDuration(90);
        double second_APD = cell_properties_second.GetLastActionPotentialDuration(90);
        double AZD_APD = cell_properties_AZD.GetLastActionPotentialDuration(90);

        //test that the aps are actually longer than control (all interventions were meant to have that effect except the last one)
        TS_ASSERT_LESS_THAN(control_APD, first_APD);
        TS_ASSERT_LESS_THAN(control_APD, second_APD);
        TS_ASSERT_LESS_THAN(AZD_APD, control_APD);

        //leave some hardcoded value for testing
        TS_ASSERT_DELTA(control_APD, 206.1278, 0.1);
        TS_ASSERT_DELTA(first_APD, 403.8293, 0.1);
        TS_ASSERT_DELTA(second_APD, 320.0195, 0.1);
        TS_ASSERT_DELTA(AZD_APD , 187.8073, 0.1);
    }
};


#endif //_TESTIONICMODELSLONG_HPP_
