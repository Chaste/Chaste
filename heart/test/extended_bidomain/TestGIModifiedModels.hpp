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


#ifndef TESTGIMODIFIEDMODELS_HPP_
#define TESTGIMODIFIEDMODELS_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Exception.hpp"
#include "RunAndCheckIonicModels.hpp"
#include "SimpleStimulus.hpp"
#include "RegularStimulus.hpp"
#include "ZeroStimulus.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta2IvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "HeartConfig.hpp"
#include "ColumnDataReader.hpp"
#include "CellProperties.hpp"
#include "CorriasBuistSMCModified.hpp"
#include "CorriasBuistICCModified.hpp"

// Note: RunOdeSolverWithIonicModel(), CheckCellModelResults(), CompareCellModelResults()
// are defined in RunAndCheckIonicModels.hpp

class TestGIModifiedModels : public CxxTest::TestSuite
{
public:



    void TestICCmodelModified(void) throw (Exception)
    {
        HeartConfig::Instance()->Reset();

        boost::shared_ptr<ZeroStimulus> stimulus(new ZeroStimulus()); //define the stimulus
        boost::shared_ptr<EulerIvpOdeSolver> solver(new EulerIvpOdeSolver); //define the solver

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,0.1);
        CorriasBuistICCModified icc_ode_system(solver, stimulus);
        icc_ode_system.SetIP3Concentration(0.000635);
        icc_ode_system.SetFractionOfVDDRInPU(0.04);

        TS_ASSERT_EQUALS(icc_ode_system.GetCarbonMonoxideScaleFactor(), 1.0); //coverage
        icc_ode_system.SetCarbonMonoxideScaleFactor(1.0);
        icc_ode_system.SetSercaPumpScaleFactor(1.0);

        std::string filebasename = "ICCmodelMartinCode";
        // Solve and write to file
        RunOdeSolverWithIonicModel(&icc_ode_system,
                                   20000,/*end time, in milliseconds*/
                                   filebasename,
                                   1000);

        ////////////////////////////////////////////////////////////////////
        //Compare with valid data. The data are valid because of the following:
        // valid data plotted against the results of Martin's code are shown in projects/alberto/test/data/MartinCodeVsChasteSingleCell.eps
        // The two solutions match very closely, hence the code that generated those files is valid.
        ///////////////////////////////////////////////////////////////////
        double tolerance = 1e-3;

        // read data entries for the new file
        ColumnDataReader data_reader("TestIonicModels", filebasename, true);
        std::vector<double> times = data_reader.GetValues("Time");
        std::vector<double> voltages = data_reader.GetValues("Vm");

        ColumnDataReader valid_reader("heart/test/data/extendedbidomain", filebasename + "ValidData",false);
        std::vector<double> valid_times = valid_reader.GetValues("Time");
        std::vector<double> valid_voltages = valid_reader.GetValues("Vm");

        TS_ASSERT_EQUALS(times.size(), valid_times.size());
        for (unsigned i=0; i<valid_times.size(); i++)
        {
            TS_ASSERT_DELTA(times[i], valid_times[i], 1e-12);
            TS_ASSERT_DELTA(voltages[i], valid_voltages[i], tolerance);
        }
    }

    void TestSMCmodelModified(void) throw (Exception)
    {
        // Set stimulus (no stimulus in this case)
        double magnitude_stimulus = 0.0;   // dimensionless
        double duration_stimulus = 0.05;  // ms
        double start_stimulus = 0.01;   // ms
        double period=1;//
        boost::shared_ptr<RegularStimulus> stimulus(new RegularStimulus(magnitude_stimulus,
                                                                        duration_stimulus,
                                                                        period,
                                                                        start_stimulus));

        boost::shared_ptr<EulerIvpOdeSolver> solver(new EulerIvpOdeSolver); //define the solver

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,1);
        CorriasBuistSMCModified smc_ode_system(solver, stimulus);

        TS_ASSERT_EQUALS(smc_ode_system.GetCarbonMonoxideScaleFactor(), 1.0); //coverage
        smc_ode_system.SetCarbonMonoxideScaleFactor(1.0);//coverage
        // Solve and write to file
        RunOdeSolverWithIonicModel(&smc_ode_system,
                                   60000,/*end time, in milliseconds*/
                                   "SMCmodel_modified",
                                   1000);

        // Calculate APDs and compare with hardcoded values
        ColumnDataReader data_reader1("TestIonicModels", "SMCmodel_modified");
        std::vector<double> voltages = data_reader1.GetValues("Vm_SM");
        //create the times vector
        double k =0;
        std::vector<double> times;
        for (unsigned i=0; i<voltages.size(); i++)
        {
          times.push_back(k);
          k=k+100;
        }

        CellProperties  cell_properties(voltages, times, -45.0);//note the lower threshold for SMC calculation of 'AP');

        TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(50), 8331.3835, 0.1);
        TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(70), 8571.1671, 0.1);
        TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(80), 8722.6543, 0.1);
        TS_ASSERT_DELTA(cell_properties.GetLastActionPotentialDuration(90), 8955.8126, 0.1);

        CorriasBuistSMCModified smc_ode_system_no_stim(solver, stimulus);

        //check the default value of the presence of ICC stimulus
        TS_ASSERT_EQUALS(smc_ode_system_no_stim.GetFakeIccStimulusPresent(), true);
        //Switch off ICC stimulus and check the voltage is almost flat
        smc_ode_system_no_stim.SetFakeIccStimulusPresent(false);
        //check member variable was switched correctly
        TS_ASSERT_EQUALS(smc_ode_system_no_stim.GetFakeIccStimulusPresent(), false);
        // Solve and write to file
        RunOdeSolverWithIonicModel(&smc_ode_system_no_stim,
                                   60000,/*end time, in milliseconds*/
                                   "SMCmodel_modified_nostim",
                                   1000);

        // Calculate APDs and compare with hardcoded values
        ColumnDataReader data_reader2("TestIonicModels", "SMCmodel_modified_nostim");
        std::vector<double> voltages_flat = data_reader2.GetValues("Vm_SM");

        CellProperties  cell_properties_2(voltages_flat, times, -45.0);

        //it should be flat now, no AP and cell properties should throw
        TS_ASSERT_THROWS_THIS(cell_properties_2.GetLastActionPotentialDuration(90),
                            "AP did not occur, never exceeded threshold voltage.");

     }

    void TestArchiving(void) throw(Exception)
    {
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename =  ArchiveLocationInfo::GetProcessUniqueFilePath("GI.arch");

        // Save
        {

            boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(0.0,1.0,0.5));
            boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
            double time_step = 0.01;
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(time_step, time_step, time_step);

            // icc and smc
            AbstractCardiacCell* const p_smc = new CorriasBuistSMCModified(p_solver, p_stimulus);
            AbstractCardiacCell* const p_icc = new CorriasBuistICCModified(p_solver, p_stimulus);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch <<  p_smc;
            output_arch <<  p_icc;

            delete p_smc;
            delete p_icc;
        }
        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCell* p_smc;
            AbstractCardiacCell* p_icc;
            input_arch >> p_smc;
            input_arch >> p_icc;

            TS_ASSERT_EQUALS( p_smc->GetNumberOfStateVariables(), 14U );
            TS_ASSERT_EQUALS( p_icc->GetNumberOfStateVariables(), 18U );

            delete p_smc;
            delete p_icc;
        }
     }

};
#endif //TESTGIMODIFIEDMODELS_HPP_
