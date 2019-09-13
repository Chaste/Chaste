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

#ifndef TESTCVODEWITHJACOBIAN_HPP_
#define TESTCVODEWITHJACOBIAN_HPP_

#include <cxxtest/TestSuite.h>
#include "Timer.hpp"

#include "CvodeAdaptor.hpp"
#include "DiFrancescoNoble1985.hpp"
#include "DiFrancescoNoble1985Cvode.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "FaberRudy2000.hpp"
#include "FaberRudy2000Cvode.hpp"
#include "FoxModel2002.hpp"
#include "FoxModel2002Cvode.hpp"
#include "HodgkinHuxley1952.hpp"
#include "HodgkinHuxley1952Cvode.hpp"
#include "LuoRudy1991.hpp"
#include "LuoRudy1991Cvode.hpp"
#include "Mahajan2008.hpp"
#include "Mahajan2008Cvode.hpp"
#include "Maleckar2008.hpp"
#include "Maleckar2008Cvode.hpp"
#include "NobleVargheseKohlNoble1998a.hpp"
#include "NobleVargheseKohlNoble1998aCvode.hpp"
#include "RegularStimulus.hpp"
#include "Shannon2004.hpp"
#include "Shannon2004Cvode.hpp"
#include "TenTusscher2006Epi.hpp"
#include "TenTusscher2006EpiCvode.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCvodeWithJacobian : public CxxTest::TestSuite
{
public:
    void TestTimingsWithAndWithoutJacobian()
    {
#ifdef CHASTE_CVODE
        // Set up a default solver and a stimulus
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new RegularStimulus(-25, 5, 1000, 1));
        double simulation_duration = 10000; // This increased to 1e6 to give the timings shown on #1795.
        double max_time_step = 1.0;

        for (unsigned i = 0; i < 10; i++)
        {
            boost::shared_ptr<CvodeAdaptor> p_cvode_adaptor(new CvodeAdaptor()); // moving this outside the loop causes seg-faults!
            p_cvode_adaptor->SetTolerances(1e-5, 1e-7); // Match AbstractCvodeCell
            boost::shared_ptr<AbstractCardiacCell> p_cell_cvode_adaptor;
            boost::shared_ptr<AbstractCvodeCell> p_cell_cvode;
            switch (i)
            {
                case 0:
                {
                    // Fox model
                    p_cell_cvode_adaptor.reset(new CellFoxModel2002FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellFoxModel2002FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 1:
                {
                    // Shannon model
                    p_cell_cvode_adaptor.reset(new CellShannon2004FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellShannon2004FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 2:
                {
                    // HH1952 model
                    p_cell_cvode_adaptor.reset(new CellHodgkinHuxley1952FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellHodgkinHuxley1952FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 3:
                {
                    // Luo Rudy model
                    p_cell_cvode_adaptor.reset(new CellLuoRudy1991FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellLuoRudy1991FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 4:
                {
                    // Mahajan model
                    p_cell_cvode_adaptor.reset(new CellMahajan2008FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellMahajan2008FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 5:
                {
                    // Maleckar model
                    p_cell_cvode_adaptor.reset(new CellMaleckar2008FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellMaleckar2008FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 6:
                {
                    // FaberRudy2000 model
                    p_cell_cvode_adaptor.reset(new CellFaberRudy2000FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellFaberRudy2000FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 7:
                {
                    // DiFrancescoNoble1985 model
                    p_cell_cvode_adaptor.reset(new CellDiFrancescoNoble1985FromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellDiFrancescoNoble1985FromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 8:
                {
                    // NobleVargheseKohlNoble1998a model
                    p_cell_cvode_adaptor.reset(new CellNobleVargheseKohlNoble1998aFromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellNobleVargheseKohlNoble1998aFromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                case 9:
                {
                    // TenTusscher2006Epi model
                    p_cell_cvode_adaptor.reset(new CellTenTusscher2006EpiFromCellML(p_cvode_adaptor, p_stimulus));
                    p_cell_cvode.reset(new CellTenTusscher2006EpiFromCellMLCvode(p_solver, p_stimulus));
                    break;
                }
                default:
                    EXCEPTION("No model");
            }

            std::cout << p_cell_cvode_adaptor->GetSystemName() << ": " << p_cell_cvode_adaptor->GetNumberOfStateVariables() << " ODEs.\n";

            // A standard CVODE adaptor solve
            p_cvode_adaptor->SetMaxSteps(0xffffffff);
            p_cell_cvode_adaptor->SetTimestep(max_time_step);

            Timer::Reset();
            p_cell_cvode_adaptor->SolveAndUpdateState(0, simulation_duration);
            Timer::Print(" 1. CVODE adaptor (numeric Jacobian)");

            // A standard native CVODE solve
            p_cell_cvode->SetMaxSteps(0xffffffff);
            p_cell_cvode->SetMaxTimestep(max_time_step);
            p_cell_cvode->ForceUseOfNumericalJacobian(true);
            TS_ASSERT(!p_cell_cvode->GetUseAnalyticJacobian());

            Timer::Reset();
            p_cell_cvode->SolveAndUpdateState(0, simulation_duration);
            Timer::Print(" 2. CVODE native with Numeric Jacobian");

            // A jacobian native CVODE solve
            p_cell_cvode->ResetToInitialConditions();
            p_cell_cvode->ForceUseOfNumericalJacobian(false);
            TS_ASSERT(p_cell_cvode->GetUseAnalyticJacobian());

            Timer::Reset();
            p_cell_cvode->SolveAndUpdateState(0, simulation_duration);
            Timer::Print(" 3. CVODE native with Analytic Jacobian");
        }
#else
        std::cout << "CVODE is not installed or Chaste hostconfig is not using it." << std::endl;
#endif
    }
};

#endif // TESTCVODEWITHJACOBIAN_HPP_
