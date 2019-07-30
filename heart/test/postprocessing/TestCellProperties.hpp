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

#ifndef _TESTCELLPROPERTIES_HPP_
#define _TESTCELLPROPERTIES_HPP_

#include <cxxtest/TestSuite.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "CellProperties.hpp"
#include "ColumnDataReader.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "FileFinder.hpp"
#include "LuoRudy1991.hpp"
#include "OdeSolution.hpp"
#include "PetscTools.hpp"
#include "RegularStimulus.hpp"

#include "PetscSetupAndFinalize.hpp"

/* HOW_TO_TAG Cardiac/Post-processing
 * Compute action potential properties (APD50, APD90, max upstroke velocities, etc) given voltage traces.
 */
class TestCellProperties : public CxxTest::TestSuite
{
private:
    void LoadMeshalyzerOutputTraces(const FileFinder& rFileFinder, const unsigned length, std::vector<double>& rTime, std::vector<double>& rVoltage)
    {
        std::ifstream apd_file((rFileFinder.GetAbsolutePath()).c_str());
        TS_ASSERT(apd_file.is_open());

        // Create the vectors to be passed to the CellProperties object
        rVoltage.resize(length);
        rTime.resize(length);
        for (unsigned i = 0; i < length; ++i)
        {
            apd_file >> rTime[i];
            apd_file >> rVoltage[i];
        }
        apd_file.close();
    }

public:
    void TestExceptionalBehaviour(void)
    {
        // Check throws an exception if no data given
        std::vector<double> empty;
        TS_ASSERT_THROWS_THIS(CellProperties cell_props(empty, empty),
                              "Insufficient time steps to calculate physiological properties.");

        //Creating an artificial flat potential profile
        std::vector<double> times;
        std::vector<double> flat_v;
        for (unsigned i = 0; i < 100; i++)
        {
            times.push_back(i);
            flat_v.push_back(-85.0);
        }

        {
            CellProperties cell_properties(flat_v, times);

            //Should throw exceptions because the cached vector of onset times (mOnsets) is empty
            TS_ASSERT_THROWS_THIS(cell_properties.GetLastActionPotentialDuration(90), "AP did not occur, never exceeded threshold voltage.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetAllActionPotentialDurations(90)[0], "AP did not occur, never exceeded threshold voltage.");

            //Should throw exceptions because upstroke was never crossed
            TS_ASSERT_THROWS_THIS(cell_properties.GetTimeAtLastMaxUpstrokeVelocity(), "AP did not occur, never descended past threshold voltage.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetLastMaxUpstrokeVelocity(), "AP did not occur, never descended past threshold voltage.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetLastPeakPotential(), "AP did not occur, never descended past threshold voltage.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetTimeAtLastPeakPotential(), "AP did not occur, never descended past threshold voltage.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetMaxUpstrokeVelocities(), "AP did not occur, never descended past threshold voltage.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetTimesAtMaxUpstrokeVelocity(), "AP did not occur, never descended past threshold voltage.");
        }

        //Now make it cross the threshold so the onset vector isn't empty any longer
        times.push_back(100);
        flat_v.push_back(20.0);

        {
            CellProperties cell_properties(flat_v, times);

            //Now this should throw an exception because the vectors of APs is empty...
            TS_ASSERT_THROWS_THIS(cell_properties.GetLastActionPotentialDuration(90), "No full action potential was recorded");

            //...but we can calculate peak properties for the last AP (though incomplete)
            // Gary: Added some extra methods because I think this is a bit misleading and should throw!
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastMaxUpstrokeVelocity(), 100, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastMaxUpstrokeVelocity(), 105, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastPeakPotential(), 20, 1e-6);
            TS_ASSERT_THROWS_THIS(cell_properties.GetLastCompleteMaxUpstrokeVelocity(), "No MaxUpstrokeVelocity matching a full action potential was recorded.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetTimeAtLastCompleteMaxUpstrokeVelocity(), "No TimeAtMaxUpstrokeVelocity matching a full action potential was recorded.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetLastCompletePeakPotential(), "No peak potential matching a full action potential was recorded.");
            TS_ASSERT_THROWS_THIS(cell_properties.GetTimeAtLastCompletePeakPotential(), "No peak potential matching a full action potential was recorded.");
        }

        // Stay up for a while...
        times.push_back(120);
        flat_v.push_back(20.0);
        // Then repolarise
        times.push_back(121);
        flat_v.push_back(-85.0);

        {
            CellProperties cell_properties(flat_v, times);

            // Now both the "complete" and "last" entries should be the same
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastMaxUpstrokeVelocity(), 100, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastMaxUpstrokeVelocity(), 105, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastPeakPotential(), 20, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastCompleteMaxUpstrokeVelocity(), 100, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastCompleteMaxUpstrokeVelocity(), 105, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastCompletePeakPotential(), 20, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastCompletePeakPotential(), 100, 1e-6);
        }

        // Stay down for a while
        times.push_back(140);
        flat_v.push_back(-85.0);
        // Then go up to a different voltage and stay up
        times.push_back(141);
        flat_v.push_back(0.0);
        times.push_back(145);
        flat_v.push_back(0.0);

        {
            CellProperties cell_properties(flat_v, times);

            // Now both the "complete" and "last" entries should now be different
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastMaxUpstrokeVelocity(), 141, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastMaxUpstrokeVelocity(), 85, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastPeakPotential(), 0, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastCompleteMaxUpstrokeVelocity(), 100, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastCompleteMaxUpstrokeVelocity(), 105, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastCompletePeakPotential(), 20, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastCompletePeakPotential(), 100, 1e-6);
        }

        // Go down again...
        times.push_back(146);
        flat_v.push_back(-85.0);

        {
            CellProperties cell_properties(flat_v, times);

            // The "complete" and "last" entries should now be the same again.
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastMaxUpstrokeVelocity(), 141, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastMaxUpstrokeVelocity(), 85, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastPeakPotential(), 0, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastCompleteMaxUpstrokeVelocity(), 141, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastCompleteMaxUpstrokeVelocity(), 85, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetLastCompletePeakPotential(), 0, 1e-6);
            TS_ASSERT_DELTA(cell_properties.GetTimeAtLastCompletePeakPotential(), 141, 1e-6);
        }

        times.push_back(999);
        TS_ASSERT_THROWS_THIS(CellProperties bad_cell_properties(flat_v, times),
                              "Time and Voltage series should be the same length. Time.size() = 108, Voltage.size() = 107");
    }

    void TestCellPhysiologicalPropertiesForRegularLr91(void)
    {
        EXIT_IF_PARALLEL;

        /*
         * Set stimulus
         */
        double magnitude_of_stimulus = -80.0;
        double duration_of_stimulus = 0.5; // ms
        double period = 1000.0; // 1s
        double when = 100.0;
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(
            magnitude_of_stimulus,
            duration_of_stimulus,
            period,
            when));

        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        /*
         * Solve
         */
        double start_time = 0.0; // ms
        double end_time = 3450; // ms

        CellLuoRudy1991FromCellML lr91_ode_system(p_solver, p_stimulus);

        OdeSolution solution = lr91_ode_system.Compute(start_time, end_time);

        solution.WriteToFile("", __FUNCTION__, "ms");

        // Now calculate the properties
        std::vector<double> voltage = solution.GetVariableAtIndex(lr91_ode_system.GetStateVariableIndex("membrane_voltage"));
        CellProperties cell_props(voltage, solution.rGetTimes()); // Use default threshold
        double timestep = solution.rGetTimes()[1] - solution.rGetTimes()[0];
        unsigned size = cell_props.GetMaxUpstrokeVelocities().size();

        TS_ASSERT_EQUALS(size, 4u);
        TS_ASSERT_DELTA(cell_props.GetMaxUpstrokeVelocities()[size - 1], 418.4795, 0.001);
        TS_ASSERT_DELTA(cell_props.GetCycleLengths()[size - 2], 1000.00, 0.01); //last apd is not finished, get cycle lengths from before
        TS_ASSERT_DELTA(cell_props.GetPeakPotentials()[size - 1], 43.1665, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetRestingPotentials()[size - 1], -84.4395, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetActionPotentialAmplitudes()[size - 1], 127.606, 0.001);
        TS_ASSERT_DELTA(cell_props.GetLastActionPotentialAmplitude(), 127.606, 0.001);
        TS_ASSERT_DELTA(cell_props.GetLastPeakPotential(), 43.1665, 0.0001);
        TS_ASSERT_DELTA(cell_props.GetTimeAtLastPeakPotential(), 3101.15, timestep);
        TS_ASSERT_DELTA(cell_props.GetLastActionPotentialDuration(20), 6.5202, timestep);
        TS_ASSERT_DELTA(cell_props.GetLastActionPotentialDuration(50), 271.1389, timestep);
        TS_ASSERT_DELTA(cell_props.GetLastActionPotentialDuration(90), 362.0155, timestep); // Should use penultimate AP
        TS_ASSERT_DELTA(cell_props.GetTimesAtMaxUpstrokeVelocity()[size - 1], 3100.7300, 0.001);
    }

    /**
     * Further tests in a more tricky case.
     * The APs tested here are more bumpy and less regular.
     */
    void TestTrickyActionPotential()
    {
        std::ifstream apd_file("heart/test/data/sample_APs/TrickyAPD.dat");
        TS_ASSERT(apd_file.is_open());

        // Create the vectors to be passed to the CellProperties object
        std::vector<double> voltages(15001);
        std::vector<double> times(15001);
        for (unsigned i = 0; i < 15001; i++)
        {
            apd_file >> voltages[i];
            times[i] = i;
        }
        apd_file.close();

        CellProperties cell_properties(voltages, times); // Use default threshold
        std::vector<double> apds = cell_properties.GetAllActionPotentialDurations(50);
        unsigned size = apds.size();
        TS_ASSERT_EQUALS(size, 10u);

        //First check that the "GetPropertyAtLastAP" actually returns a value equal to
        // the last element of the "GetAllOfThem" vector
        TS_ASSERT_EQUALS(apds[size - 1], cell_properties.GetLastActionPotentialDuration(50));

        TS_ASSERT_EQUALS(cell_properties.GetTimesAtMaxUpstrokeVelocity()[size - 1],
                         cell_properties.GetTimeAtLastMaxUpstrokeVelocity());

        TS_ASSERT_EQUALS(cell_properties.GetMaxUpstrokeVelocities()[size - 1],
                         cell_properties.GetLastMaxUpstrokeVelocity());

        // Then check against hardcoded values (checked manually from the file)
        double timestep = times[1] - times[0];
        TS_ASSERT_DELTA(apds[0], 185.37, timestep);
        TS_ASSERT_DELTA(apds[1], 186.3, timestep);
        TS_ASSERT_DELTA(apds[2], 185.3, timestep);
        TS_ASSERT_DELTA(apds[3], 184.3, timestep);
        TS_ASSERT_DELTA(apds[4], 183.3, timestep);
        TS_ASSERT_DELTA(apds[5], 183.3, timestep);
        TS_ASSERT_DELTA(apds[6], 180.2, timestep);
        TS_ASSERT_DELTA(apds[7], 177.3, timestep);
        TS_ASSERT_DELTA(apds[8], 175.3, timestep);
        TS_ASSERT_DELTA(apds[size - 1], 175.3, timestep);

        // Check against hardcoded resting values (checked manually from the file)
        std::vector<double> resting_values = cell_properties.GetRestingPotentials();
        size = resting_values.size();
        TS_ASSERT_EQUALS(size, 10u);

        TS_ASSERT_DELTA(resting_values[0], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[1], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[2], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[3], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[4], -84.6, 0.1);
        TS_ASSERT_DELTA(resting_values[5], -84.7, 0.1);
        TS_ASSERT_DELTA(resting_values[6], -85.0, 0.1);
        TS_ASSERT_DELTA(resting_values[7], -85.1, 0.1);
        TS_ASSERT_DELTA(resting_values[8], -85.1, 0.1);
        TS_ASSERT_DELTA(resting_values[9], -85.2, 0.1);

        double last_resting_value = cell_properties.GetLastRestingPotential();
        TS_ASSERT_DELTA(resting_values[9], last_resting_value, 1e-12);

        //This file comes from a frequency drop tissue simulation.
        //First five beats at 500 bcl and last five at 2500 bcl
        std::vector<double> cycle_lengths = cell_properties.GetCycleLengths();
        size = cycle_lengths.size();

        TS_ASSERT_DELTA(cycle_lengths[0], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[1], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[2], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[3], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[4], 500, 1);
        TS_ASSERT_DELTA(cycle_lengths[5], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[6], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[7], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[8], 2500, 1);
        TS_ASSERT_DELTA(cycle_lengths[size - 1], 2500, 1);
    }

    void TestEadDetection()
    {
        //this file contains 4 Aps
        std::ifstream ead_file("heart/test/data/sample_APs/Ead.dat");
        TS_ASSERT(ead_file.is_open());

        // Create the vectors to be passed to the CellProperties object
        std::vector<double> voltages(2001);
        std::vector<double> times(2001);
        for (unsigned i = 0; i < 2001; i++)
        {
            ead_file >> voltages[i];
            times[i] = i;
        }
        ead_file.close();

        double threshold = -40;
        CellProperties cell_properties(voltages, times, threshold);

        //first, we calculate how many full Aps we have here
        std::vector<double> apds = cell_properties.GetAllActionPotentialDurations(90);
        unsigned size = apds.size();
        //there should be 4 aps in this file
        TS_ASSERT_EQUALS(size, 4u);

        std::vector<unsigned> above_threshold_depo = cell_properties.GetNumberOfAboveThresholdDepolarisationsForAllAps();
        //first AP has just a notch
        TS_ASSERT_EQUALS(above_threshold_depo[0], 1u);
        //second AP has monotonous repolarisation
        TS_ASSERT_EQUALS(above_threshold_depo[1], 0u);
        //third AP has a notch plus 2 EADs
        TS_ASSERT_EQUALS(above_threshold_depo[2], 3u);
        //fourth AP has monotonous repolarisation
        TS_ASSERT_EQUALS(above_threshold_depo[3], 0u);

        unsigned number_of_changes_for_last_ap = cell_properties.GetNumberOfAboveThresholdDepolarisationsForLastAp();
        TS_ASSERT_EQUALS(number_of_changes_for_last_ap, above_threshold_depo[size - 1]);
    }

    void TestVeryLongApDetection()
    {
        // This test is added to deal with a problem where you get long action potentials (overlapping the next stimulus)
        // and then, depending on the threshold, these got reported as one or two action potentials.
        // i.e. same single big AP reported twice instead of once!

        // Desired behaviour - it is always one action potential until it falls below APD_{X} threshold.

        {
            FileFinder file_finder("heart/test/data/sample_APs", RelativeTo::ChasteSourceRoot);
            ColumnDataReader reader(file_finder, "long_AP_two_stim");

            std::vector<double> times = reader.GetValues("Time");
            std::vector<double> voltages = reader.GetValues("membrane_voltage");

            CellProperties cell_properties(voltages, times, -50.0);
            std::vector<double> apd_90s = cell_properties.GetAllActionPotentialDurations(90);
            TS_ASSERT_EQUALS(apd_90s.size(), 1u);
            TS_ASSERT_DELTA(apd_90s[0], 1712.3745, 1e-3);

            std::vector<double> apd_30s = cell_properties.GetAllActionPotentialDurations(30);
            TS_ASSERT_EQUALS(apd_30s.size(), 1u);
            TS_ASSERT_DELTA(apd_30s[0], 505.4778, 1e-3);
        }

        {
            // With this one we have a threshold which means it makes sense to report two APD30s but only
            // one APD90.
            FileFinder file_finder("heart/test/data/sample_APs", RelativeTo::ChasteSourceRoot);
            ColumnDataReader reader(file_finder, "long_AP_two_stim_2");

            std::vector<double> times = reader.GetValues("Time");
            std::vector<double> voltages = reader.GetValues("membrane_voltage");

            CellProperties cell_properties(voltages, times, -50.0);
            std::vector<double> apd_90s = cell_properties.GetAllActionPotentialDurations(90);
            TS_ASSERT_EQUALS(apd_90s.size(), 1u);
            TS_ASSERT_DELTA(apd_90s[0], 1453.6436, 1e-3);

            std::vector<double> apd_30s = cell_properties.GetAllActionPotentialDurations(30);
            TS_ASSERT_EQUALS(apd_30s.size(), 2u);
            TS_ASSERT_DELTA(apd_30s[0], 768.903, 1e-3);
            TS_ASSERT_DELTA(apd_30s[1], 194.72, 1e-3);
        }
    }

    void TestActionPotentialCalculations()
    {
        /*
        * In this simulation the stimulus was introduced at t=1ms.
        *
        * The only difference between the two files is that Mahajan2008 starts at t=0
        * and Mahajan 2008Immediate starts at t=1 (I just removed the top line).
        *
        * (They should therefore return the same action potential properties,
        * as it is exactly the same trace in each one).
        */

        FileFinder file_finder("heart/test/data/sample_APs", RelativeTo::ChasteSourceRoot);
        double threshold = -70;
        // We now ignore the threshold when calculating the APD - it is now only used for
        // detecting the upstroke and calculating cycle lengths. We now interpolate back to
        // find the time at which the target voltage is exceeded. So the calculation is now
        // robust to this.

        /*
        * I plotted the graph in gnuplot and my back of the envelope calculations are as follows:
        * Peak voltage = 48mV
        * Start and final voltages = -87mV
        * Difference = 135mV
        * Therefore 50% repolarisation = -19.5mV
        *                      and 90% = -73.5mV
        * These are crossed at 4.4, 259.2ish and 2.96, 305.1ish ms
        * Therefore APD50 and 90 should be about 254.8 and 302.2ms you'd think.
        *
        */
        double target_apd_50 = 254.8;
        double target_apd_90 = 302.16;
        double tolerance = 0.1; //ms

        { // Stimulus applied to Mahajan model after 1 ms

            ColumnDataReader reader(file_finder, "Mahajan2008");
            std::vector<double> times = reader.GetValues("Time");
            std::vector<double> voltages = reader.GetValues("membrane_voltage");

            CellProperties cell_properties(voltages, times, threshold);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(50), target_apd_50, tolerance);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(90), target_apd_90, tolerance);
        }

        { // Stimulus applied to Mahajan model immediately

            ColumnDataReader reader(file_finder, "Mahajan2008Immediate");
            std::vector<double> times = reader.GetValues("Time");
            std::vector<double> voltages = reader.GetValues("membrane_voltage");

            CellProperties cell_properties(voltages, times, threshold);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(50), target_apd_50, tolerance);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(90), target_apd_90, tolerance);
        }

        // Different voltage units for the phenomenological models, and time happens to be scaled strangely too...
        threshold = 0.4;
        tolerance = 0.1; //ms - these weren't exactly the same cells from meshalyzer but should be very close...
        target_apd_50 = 197.902;
        target_apd_90 = 237.746;

        { // Stimulus applied to a phenomenological model after 1ms, stimulated cell APD
            FileFinder finder("heart/test/data/sample_APs/phenomenological_delayed_stim.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT_EQUALS(finder.Exists(), true);
            std::vector<double> voltages;
            std::vector<double> times;
            LoadMeshalyzerOutputTraces(finder, 5001, times, voltages);

            CellProperties cell_properties(voltages, times, threshold);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(50), target_apd_50, tolerance);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(90), target_apd_90, tolerance);
        }

        target_apd_50 = 200.8365; // More of a difference in APD50 because of different height peaks and therefore thresholds.
        target_apd_90 = 237.3066;

        { // Stimulus applied to a phenomenological model immediately, stimulated cell APD
            FileFinder finder("heart/test/data/sample_APs/phenomenological_immediate_stim.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT_EQUALS(finder.Exists(), true);
            std::vector<double> voltages;
            std::vector<double> times;
            LoadMeshalyzerOutputTraces(finder, 5001, times, voltages);

            CellProperties cell_properties(voltages, times, threshold);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(50), target_apd_50, tolerance);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(90), target_apd_90, tolerance);
        }

        /**
         * Now we move to phenomenological APDs from the outside of the tissue
         */
        target_apd_50 = 195.646;
        target_apd_90 = 230.361;

        { // Stimulus applied to a phenomenological model after 1ms, outer cell APD
            FileFinder finder("heart/test/data/sample_APs/phenomenological_delayed_stim_outer.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT_EQUALS(finder.Exists(), true);
            std::vector<double> voltages;
            std::vector<double> times;
            LoadMeshalyzerOutputTraces(finder, 5001, times, voltages);

            CellProperties cell_properties(voltages, times, threshold);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(50), target_apd_50, tolerance);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(90), target_apd_90, tolerance);
        }

        { // Stimulus applied to a phenomenological model immediately, outer cell APD
            FileFinder finder("heart/test/data/sample_APs/phenomenological_immediate_stim_outer.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT_EQUALS(finder.Exists(), true);
            std::vector<double> voltages;
            std::vector<double> times;
            LoadMeshalyzerOutputTraces(finder, 5001, times, voltages);

            CellProperties cell_properties(voltages, times, threshold);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(50), target_apd_50, tolerance);
            TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(90), target_apd_90, tolerance);
        }

        {
            // This is a dataset with a small spike of about -52 mV followed by a full action potential.
            // From a 0.5Hz pacing protocol.
            FileFinder finder("heart/test/data/sample_APs/Tricky_alternans.dat", RelativeTo::ChasteSourceRoot);
            TS_ASSERT_EQUALS(finder.IsFile(), true);
            std::vector<double> voltages;
            std::vector<double> times;
            LoadMeshalyzerOutputTraces(finder, 22674, times, voltages);

            //std::cout << "With -50 threshold:\n";
            threshold = -50.0; // Last one above, first one below
            {
                CellProperties cell_properties(voltages, times, threshold);

                std::vector<double> dV_dt_max = cell_properties.GetMaxUpstrokeVelocities();
                TS_ASSERT_EQUALS(dV_dt_max.size(), 1u);

                // This was failing, prompting ticket #2844
                std::vector<double> apds = cell_properties.GetAllActionPotentialDurations(90);
                TS_ASSERT_EQUALS(apds.size(), 1u);

                TS_ASSERT_DELTA(apds[0], 209.4662, tolerance);

                std::vector<double> peak_Vs = cell_properties.GetPeakPotentials();
                TS_ASSERT_EQUALS(peak_Vs.size(), 1u);
                TS_ASSERT_DELTA(peak_Vs[0], 20.742, tolerance);

                std::vector<double> amplitudes = cell_properties.GetActionPotentialAmplitudes();
                TS_ASSERT_EQUALS(amplitudes.size(), 1u);
                TS_ASSERT_DELTA(amplitudes[0], 106.2128, tolerance);

                std::vector<double> upstroke_times = cell_properties.GetTimesAtMaxUpstrokeVelocity();
                TS_ASSERT_EQUALS(upstroke_times.size(), 1u);
                TS_ASSERT_DELTA(upstroke_times[0], 2012.5, tolerance);

                std::vector<double> resting_potentials = cell_properties.GetRestingPotentials();
                TS_ASSERT_EQUALS(resting_potentials.size(), 1u);
                TS_ASSERT_DELTA(resting_potentials[0], -85.4714, tolerance);
            }

            //std::cout << "With -60 threshold:\n";
            threshold = -60.0; // Both above threshold - all these properties were already calculated OK
            {
                CellProperties cell_properties(voltages, times, threshold);

                std::vector<double> dV_dt_max = cell_properties.GetMaxUpstrokeVelocities();
                TS_ASSERT_EQUALS(dV_dt_max.size(), 2u);

                std::vector<double> apds = cell_properties.GetAllActionPotentialDurations(90);
                TS_ASSERT_EQUALS(apds.size(), 2u);

                TS_ASSERT_DELTA(apds[0], 20.3029, tolerance);
                TS_ASSERT_DELTA(apds[1], 209.4662, tolerance);

                std::vector<double> peak_Vs = cell_properties.GetPeakPotentials();
                TS_ASSERT_EQUALS(peak_Vs.size(), 2u);
                TS_ASSERT_DELTA(peak_Vs[0], -52.03, tolerance);
                TS_ASSERT_DELTA(peak_Vs[1], 20.742, tolerance);

                std::vector<double> peak_V_times = cell_properties.GetTimesAtPeakPotentials();
                TS_ASSERT_EQUALS(peak_V_times.size(), 2u);
                TS_ASSERT_DELTA(peak_V_times[0], 8.5, tolerance);
                TS_ASSERT_DELTA(peak_V_times[1], 2018.5, tolerance);

                TS_ASSERT_DELTA(cell_properties.GetTimeAtLastPeakPotential(), 2018.5, tolerance);

                std::vector<double> amplitudes = cell_properties.GetActionPotentialAmplitudes();
                TS_ASSERT_EQUALS(amplitudes.size(), 2u);
                TS_ASSERT_DELTA(amplitudes[0], 33.39, tolerance);
                TS_ASSERT_DELTA(amplitudes[1], 106.2128, tolerance);

                std::vector<double> upstroke_times = cell_properties.GetTimesAtMaxUpstrokeVelocity();
                TS_ASSERT_EQUALS(upstroke_times.size(), 2u);
                TS_ASSERT_DELTA(upstroke_times[0], 0, tolerance);
                TS_ASSERT_DELTA(upstroke_times[1], 2012.5, tolerance);

                std::vector<double> resting_potentials = cell_properties.GetRestingPotentials();
                TS_ASSERT_EQUALS(resting_potentials.size(), 2u);
                TS_ASSERT_DELTA(resting_potentials[0], -85.4714, tolerance);
                TS_ASSERT_DELTA(resting_potentials[1], -85.4714, tolerance);
            }
        }
    }
};

#endif //_TESTCELLPROPERTIES_HPP_
