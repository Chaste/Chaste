/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "CvodeAdaptor.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RegularStimulus.hpp"
#include "HodgkinHuxley1952.hpp"
#include "HodgkinHuxley1952Cvode.hpp"
#include "HH1952WithJacobian.hpp"
#include "HH1952WithJacobianCvode.hpp"

class TestCvodeWithJacobian : public CxxTest::TestSuite
{
public:
    void TestTimingsWithAndWithoutJacobian() throw (Exception)
    {
#ifdef CHASTE_CVODE
        // Set up a default solver and a stimulus
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver());
        boost::shared_ptr<CvodeAdaptor> p_cvode_adaptor(new CvodeAdaptor());
        boost::shared_ptr<CvodeAdaptor> p_cvode_adaptor_jacobian(new CvodeAdaptor());
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new RegularStimulus(-25,5,1000,1));
        double simulation_duration = 10000;

        boost::shared_ptr<AbstractCardiacCell> hh_1952_cvode_adaptor(new CellHodgkinHuxley1952FromCellML(p_cvode_adaptor,p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> hh_1952_cvode(new CellHodgkinHuxley1952FromCellMLCvode(p_solver,p_stimulus));
        boost::shared_ptr<AbstractCardiacCell> hh_1952_cvode_adaptor_jacobian(new HH1952WithJacobian(p_cvode_adaptor_jacobian,p_stimulus));
        boost::shared_ptr<AbstractCvodeCell> hh_1952_cvode_jacobian(new HH1952WithJacobianCvode(p_solver,p_stimulus));


        double start_time, end_time, elapsed_time = 0.0;

        // A standard CVODE adaptor solve
        p_cvode_adaptor->SetMaxSteps(0xffffffff);
        start_time = (double) std::clock();
        hh_1952_cvode_adaptor->SolveAndUpdateState(0, simulation_duration);
        end_time = (double) std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout << "1. CVODE adaptor, elapsed time = " << elapsed_time << " secs for " << simulation_duration << " ms\n";

        // A standard native CVODE solve
        hh_1952_cvode->SetMaxSteps(0xffffffff);
        start_time = (double) std::clock();
        hh_1952_cvode->SolveAndUpdateState(0, simulation_duration);
        end_time = (double) std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout << "2. CVODE native, elapsed time = " << elapsed_time << " secs for " << simulation_duration << " ms\n";

/// \todo # Not yet implemented (it would work, but wouldn't use a Jacobian).
//        // A jacobian CVODE adaptor solve
//        p_cvode_adaptor_jacobian->SetMaxSteps(1e10);
//        start_time = std::clock();
//        hh_1952_cvode_adaptor_jacobian->SolveAndUpdateState(0, simulation_duration);
//        end_time = std::clock();
//        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
//        std::cout << "3. CVODE adaptor with Jacobian, elapsed time = " << elapsed_time << " secs for " << simulation_duration << " ms\n";

        // A jacobian native CVODE solve
        hh_1952_cvode_jacobian->SetMaxSteps(0xffffffff);
        start_time = (double) std::clock();
        hh_1952_cvode_jacobian->SolveAndUpdateState(0, simulation_duration);
        end_time = (double) std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout << "4. CVODE native with Jacobian, elapsed time = " << elapsed_time << " secs for " << simulation_duration << " ms\n";
#else
        std::cout << "CVODE is not installed or Chaste hostconfig is not using it." << std::endl;
#endif
    }
};

#endif // TESTCVODEWITHJACOBIAN_HPP_
