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


#ifndef _RUNANDCHECKIONICMODELS_HPP_
#define _RUNANDCHECKIONICMODELS_HPP_

#include <vector>
#include <string>

#include "OdeSolution.hpp"

#include "ColumnDataWriter.hpp"
#include "ColumnDataReader.hpp"

#include "AbstractCardiacCell.hpp"
#include "HeartConfig.hpp"

void RunOdeSolverWithIonicModel(AbstractCardiacCell* pOdeSystem,
                                double endTime,
                                std::string filename,
                                int stepPerRow=100,
                                bool doComputeExceptVoltage=true,
                                bool useSamplingInterval=false);

void CheckCellModelResults(const std::string& rBaseResultsFilename,
                           std::string validResultsBasename = "",
                           double tolerance = 1e-3);

std::vector<double> GetVoltages(ColumnDataReader& rReader);

void CompareCellModelResults(std::string baseResultsFilename1, std::string baseResultsFilename2,
                             double tolerance, bool vOnly=false, std::string folderName="TestIonicModels");


/*
 * Note: we have to have the function definitions here, rather than in a .cpp file, if we're
 * to include them in heart/src, as otherwise linking of the heart library fails because
 * CxxTest routines are not defined.
 */

#include <cmath>

void RunOdeSolverWithIonicModel(AbstractCardiacCell* pOdeSystem,
                                double endTime,
                                std::string filename,
                                int stepPerRow,
                                bool doComputeExceptVoltage,
                                bool useSamplingInterval)
{
    double start_time = 0.0;

    if (doComputeExceptVoltage)
    {
        // Store the current system state
        std::vector<double> state_variables_copy = pOdeSystem->rGetStateVariables();

        // Test ComputeExceptVoltage
        double v_init = pOdeSystem->GetVoltage();
        pOdeSystem->ComputeExceptVoltage(start_time, start_time + 100.0 /*ms*/);
        double v_end = pOdeSystem->GetVoltage();
        TS_ASSERT_DELTA(v_init, v_end, 1e-6);

        // Test SetVoltage
        pOdeSystem->SetVoltage(1e6);
        TS_ASSERT_DELTA(pOdeSystem->GetVoltage(), 1e6, 1e-6);

        // Reset the system
        pOdeSystem->SetStateVariables(state_variables_copy);
    }

    // Solve and write to file
    if (useSamplingInterval)
    {
        OdeSolution solution = pOdeSystem->Compute(start_time, endTime, HeartConfig::Instance()->GetOdeTimeStep() * stepPerRow);
        solution.WriteToFile("TestIonicModels", filename, "ms", 1, false, 4);
    }
    else
    {
        OdeSolution solution = pOdeSystem->Compute(start_time, endTime);
        solution.WriteToFile("TestIonicModels", filename, "ms", stepPerRow, false, 4);
    }
}

std::vector<double> GetVoltages(ColumnDataReader& rReader)
{
    //Rather Ugly, we can't currently guarantee what the name of the voltage column is,
    //hence we try to cover the most common possibilities
    std::vector<double> voltages;
    if (rReader.HasValues("V"))
    {
        voltages = rReader.GetValues("V");
    }
    else if (rReader.HasValues("membrane_voltage"))
    {
        voltages = rReader.GetValues("membrane_voltage");
    }
    else if (rReader.HasValues("membrane__V"))
    {
        voltages = rReader.GetValues("membrane__V");
    }
    else
    {
        EXCEPTION("Model membrane voltage not recognised.");
    }
    return voltages;
}

std::vector<double> GetIntracellularCalcium(ColumnDataReader& rReader)
{
    //Rather Ugly, we can't currently guarantee what the name of the calcium column is,
    //hence we try to cover the most common possibilities
    std::vector<double> cai;
    if (rReader.HasValues("CaI"))
    {
        cai = rReader.GetValues("CaI");
    }
    else if (rReader.HasValues("cytosolic_calcium_concentration"))
    {
        cai = rReader.GetValues("cytosolic_calcium_concentration");
    }
    else if (rReader.HasValues("Cai"))
    {
        cai = rReader.GetValues("Cai");
    }
    else
    {
        EXCEPTION("Model intracellular calcium is not recognised.");
    }
    return cai;
}

std::vector<double> GetHGate(ColumnDataReader& rReader)
{
    //Rather Ugly, we can't currently guarantee what the name of the h gate column is,
    //hence we try to cover the most common possibilities
    std::vector<double> h_values;
    if (rReader.HasValues("h"))
    {
        h_values = rReader.GetValues("h");
    }
    else if (rReader.HasValues("fast_sodium_current_h_gate__h"))
    {
        h_values = rReader.GetValues("fast_sodium_current_h_gate__h");
    }
    else if (rReader.HasValues("membrane_fast_sodium_current_h_gate"))
    {
        h_values = rReader.GetValues("membrane_fast_sodium_current_h_gate");
    }
    else
    {
        EXCEPTION("Model h gate is not recognised.");
    }
    return h_values;
}

/*
 * Check the cell model against a previous version
 * or another source e.g. Alan's COR
 */
void CheckCellModelResults(const std::string& rBaseResultsFilename,
                           std::string validResultsBasename,
                           double tolerance)
{
    // read data entries for the new file
    ColumnDataReader data_reader("TestIonicModels", rBaseResultsFilename);
    std::vector<double> times = data_reader.GetValues("Time");
    std::vector<double> voltages = GetVoltages(data_reader);

    if (validResultsBasename == "")
    {
        validResultsBasename = rBaseResultsFilename;
    }

    ColumnDataReader valid_reader("heart/test/data/ionicmodels", validResultsBasename + "ValidData",
                                  false);
    std::vector<double> valid_times = valid_reader.GetValues("Time");
    std::vector<double> valid_voltages = GetVoltages(valid_reader);

    TS_ASSERT_EQUALS(times.size(), valid_times.size());
    for (unsigned i=0; i<valid_times.size(); i++)
    {
        TS_ASSERT_DELTA(times[i], valid_times[i], 1e-12);
        TS_ASSERT_DELTA(voltages[i], valid_voltages[i], tolerance);
    }
}

void CompareCellModelResults(std::string baseResultsFilename1, std::string baseResultsFilename2,
                             double tolerance, bool vOnly, std::string folderName)
{
    // Compare 2 sets of results, e.g. from 2 different solvers for the same model.
    // If the time series differ, the finer resolution must be given first.
    ColumnDataReader data_reader1(folderName, baseResultsFilename1);
    std::vector<double> times1 = data_reader1.GetValues("Time");
    std::vector<double> voltages1 = GetVoltages(data_reader1);
    std::vector<double> calcium1;
    std::vector<double> h1;

    ColumnDataReader data_reader2(folderName, baseResultsFilename2);
    std::vector<double> times2 = data_reader2.GetValues("Time");
    std::vector<double> voltages2 = GetVoltages(data_reader2);
    std::vector<double> calcium2;
    std::vector<double> h2;

    if (!vOnly)
    {
        calcium1 = GetIntracellularCalcium(data_reader1);
        h1 = GetHGate(data_reader1);
        calcium2 = GetIntracellularCalcium(data_reader2);
        h2 = GetHGate(data_reader2);
    }

    TS_ASSERT(times1.size() >= times2.size());
    double last_v = voltages2[0];
    double tol = tolerance;
    for (unsigned i=0, j=0; i<times2.size(); i++)
    {
        // Find corresponding time index
        while (j<times1.size() && times1[j] < times2[i] - 1e-12)
        {
            j++;
        }

        // Set tolerance higher in upstroke
        if (fabs(voltages2[i] - last_v) > 0.05)
        {
            tol = tolerance * 25;
        }
        else
        {
            tol = tolerance;
        }
        last_v = voltages2[i];

        TS_ASSERT_DELTA(times1[j], times2[i], 1e-12);
        // adjust tol to data
        TS_ASSERT_DELTA(voltages1[j], voltages2[i], tol);
        if (!vOnly)
        {
            TS_ASSERT_DELTA(calcium1[j],  calcium2[i],  tol/100);
            TS_ASSERT_DELTA(h1[j],        h2[i],        tol/10);
        }
    }
}


#endif //_RUNANDCHECKIONICMODELS_HPP_
