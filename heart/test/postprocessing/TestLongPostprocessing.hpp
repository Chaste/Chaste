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

#ifndef TESTLONGPOSTPROCESSING_HPP_
#define TESTLONGPOSTPROCESSING_HPP_


#include <cxxtest/TestSuite.h>
#include <ctime>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>

#include "CardiacSimulationArchiver.hpp"

#include "AbstractCardiacCell.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "RegularStimulus.hpp"
#include "LuoRudy1991.hpp"
#include "HeartConfig.hpp"
#include "MonodomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"

#include "PetscSetupAndFinalize.hpp"


template <unsigned DIM>
class PointStimulusCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<RegularStimulus> mpStimulus;
    double mAreaOrVolume;

public:
    PointStimulusCellFactory(const double& rStimMag, const double& rStimDuration, const double& rPacingCycleLength, const double& rAreaOrVolume )
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new RegularStimulus(rStimMag, rStimDuration, rPacingCycleLength, 0)),// Default stimulus introduced at t=0.
          mAreaOrVolume(rAreaOrVolume)
    {
        assert(DIM==1 || DIM==2);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];

        if (x < 0.03)
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

class TestLongPostprocessing : public CxxTest::TestSuite
{

public:


    void Test2DSimulations() throw(Exception)
    {
        double conductivity_scale = 1;
        double h = 0.01; // cm
        double ode_time_step = 0.005; //ms
        double pde_time_step = 0.01; //ms
        unsigned num_stims = 1;

        TetrahedralMesh<2,2> mesh;
        unsigned num_elem_x = (unsigned)(0.5/h);  // num elements to make 5mm
        unsigned num_elem_y = (unsigned)(0.5/h);  // num elements to make 5mm
        //unsigned num_elem_z = (unsigned)(0.15/h);// Num elements to make 0.3cm
        double pacing_cycle_length = 350;
        double stim_mag = -500000;
        double stim_dur = 3;
        double area = 0.005;

        mesh.ConstructRectangularMesh(num_elem_x, num_elem_y);
        mesh.Scale(h,h); // Get mesh into units of cm.

        std::string archive_dir_base("LongPostprocessing_archives/archive");
        std::string archive_dir_current;

        // Setup
        HeartConfig::Instance()->SetSimulationDuration(pacing_cycle_length); //ms
        HeartConfig::Instance()->SetOutputDirectory("LongPostprocessing");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        // These lines make postprocessing fast or slow.
        //HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_time_step, pde_time_step, 1);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_time_step, pde_time_step, 0.01);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.4*conductivity_scale*1.171, 1.4*conductivity_scale*1.171));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400.0); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2

        std::vector<std::pair<double,double> > apds_requested;
        apds_requested.push_back(std::pair<double, double>(90,-30)); //repolarisation percentage and threshold
        HeartConfig::Instance()->SetApdMaps(apds_requested);
//        std::vector<double> excitation_threshold;
//        excitation_threshold.push_back(-30.0);
//        HeartConfig::Instance()->SetUpstrokeTimeMaps(excitation_threshold);
//        HeartConfig::Instance()->SetMaxUpstrokeVelocityMaps(excitation_threshold);

        for (unsigned stim_counter=0; stim_counter < num_stims; stim_counter++ )
        {
            // Load problem
            MonodomainProblem<2> *p_monodomain_problem;
            if (stim_counter==0)
            {
                PointStimulusCellFactory<2> cell_factory(stim_mag, stim_dur, pacing_cycle_length, area);
                p_monodomain_problem = new MonodomainProblem<2>( &cell_factory );
                p_monodomain_problem->SetMesh(&mesh);
                p_monodomain_problem->Initialise();
            }
            else
            {
                p_monodomain_problem = CardiacSimulationArchiver<MonodomainProblem<2> >::Load(archive_dir_current);
            }

            HeartConfig::Instance()->SetSimulationDuration((double) (stim_counter+1)*pacing_cycle_length); //ms

            // set new directories to work from
            std::stringstream stringoutput;
            stringoutput << stim_counter;
            std::string stim_counter_string = stringoutput.str();

            archive_dir_current = archive_dir_base + "_" + stim_counter_string;
            HeartConfig::Instance()->SetOutputFilenamePrefix("results_" + stim_counter_string);

            // Solve problem (this does the postprocessing too when HeartConfig options are set).
            p_monodomain_problem->Solve();

            HeartEventHandler::Headings();
            HeartEventHandler::Report();

            // Save problem to archive
            CardiacSimulationArchiver<MonodomainProblem<2> >::Save(*p_monodomain_problem, archive_dir_current, false);
            std::cout << "Archived to " << archive_dir_current << "\n" << std::flush;

            if (PetscTools::AmMaster()) // Best to have only the master processor messing with files
            {
                // Copy the postprocessing results into the archive folders so they aren't wiped.
                std::vector<std::string> files;
                files.push_back("Apd_90_-30_Map");
//                files.push_back("MaxUpstrokeVelocityMap_-30");
//                files.push_back("UpstrokeTimeMap_-30");
                for (unsigned i=0; i<files.size(); i++)
                {
                    std::string command = "mv " +  OutputFileHandler::GetChasteTestOutputDirectory() + HeartConfig::Instance()->GetOutputDirectory() + "/output/" + files[i] + ".dat "
                                          + OutputFileHandler::GetChasteTestOutputDirectory() + archive_dir_current + "/" + files[i] + ".dat";

                    // If we are in an PetscTools::AmMaster() block we need to use this instead of EXPECT0.
                    ABORT_IF_NON0(system, command);
                }
            }
            PetscTools::Barrier();
        }// close for loop
    }//close void Test2dSimulations

};

#endif /* TESTLONGPOSTPROCESSING_HPP_ */

