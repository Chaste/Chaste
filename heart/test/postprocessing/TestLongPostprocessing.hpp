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

#ifndef TESTLONGPOSTPROCESSING_HPP_
#define TESTLONGPOSTPROCESSING_HPP_


#include <cxxtest/TestSuite.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>

#include "CardiacSimulationArchiver.hpp"

#include "AbstractCardiacCell.hpp"
#include "AbstractCardiacCellWithModifiers.hpp"
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        double x = pNode->rGetLocation()[0];

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


    void Test2DSimulations()
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
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_time_step, pde_time_step, 10); // Leads to 10MB VTK file
        //HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(ode_time_step, pde_time_step, 0.01); // Leads to 1GB VTK file

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.4*conductivity_scale*1.171, 1.4*conductivity_scale*1.171));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400.0); // 1/cm
        HeartConfig::Instance()->SetCapacitance(1.0); // uF/cm^2
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer();
#ifdef CHASTE_VTK
        HeartConfig::Instance()->SetVisualizeWithVtk();
#endif

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
            OutputFileHandler archive_directory(archive_dir_current, true); // Clean a folder for new results
            HeartConfig::Instance()->SetOutputFilenamePrefix("results_" + stim_counter_string);

            // Solve problem (this does the postprocessing too when HeartConfig options are set).
            p_monodomain_problem->Solve();

            HeartEventHandler::Headings();
            HeartEventHandler::Report();

            // Save problem to archive
            CardiacSimulationArchiver<MonodomainProblem<2> >::Save(*p_monodomain_problem, archive_dir_current, false);
            std::cout << "Archived to " << archive_dir_current << "\n" << std::flush;

            // Copy the postprocessing results into the archive folders so they aren't wiped.
            std::vector<std::string> files;
            files.push_back("Apd_90_minus_30_Map");
//                files.push_back("MaxUpstrokeVelocityMap_-30");
//                files.push_back("UpstrokeTimeMap_-30");

            for (unsigned i=0; i<files.size(); i++)
            {
                FileFinder file_to_copy(HeartConfig::Instance()->GetOutputDirectory() + "/output/" + files[i] + ".dat", RelativeTo::ChasteTestOutput);
                TS_ASSERT(file_to_copy.IsFile());
                archive_directory.CopyFileTo(file_to_copy);
            }
        }// close for loop
    }//close void Test2dSimulations
};

#endif /* TESTLONGPOSTPROCESSING_HPP_ */
