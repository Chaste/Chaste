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


#ifndef TESTHODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_
#define TESTHODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>
//#include <iostream>
#include <vector>

#include "HodgkinHuxley1952.hpp"
#include "SimpleStimulus.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestHodgkinHuxleySquidAxon1952OriginalOdeSystem: public CxxTest::TestSuite
{
public:

    void TestHHModelAtSingularities()
    {
        /*
        * Set stimulus
        */
        double magnitude_stimulus = 0.0;  // uA/cm2
        double duration_stimulus = 0.;  // ms
        double start_stimulus = 0.0;   // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(
                magnitude_stimulus,
                duration_stimulus,
                start_stimulus));

        boost::shared_ptr<AbstractIvpOdeSolver> p_solver; // We don't actually need a solver

        CellHodgkinHuxley1952FromCellML hh52_ode_system(p_solver, p_stimulus);

        double v_singularity[2];
        v_singularity[0]=-65;
        v_singularity[1]=-50;

        for (int i=0; i<2; i++)
        {

            std::vector<double> yleft;

            //mVariableNames.push_back("V");
            //mVariableUnits.push_back("mV");
            yleft.push_back(v_singularity[i]+0.1);

            //mVariableNames.push_back("n");
            //mVariableUnits.push_back("");
            yleft.push_back(0.325);

            //mVariableNames.push_back("h");
            //mVariableUnits.push_back("");
            yleft.push_back(0.6);

            //mVariableNames.push_back("m");
            //mVariableUnits.push_back("");
            yleft.push_back(0.05);


            std::vector<double> rhsleft(yleft.size());
            hh52_ode_system.EvaluateYDerivatives (0.0, yleft, rhsleft);

            std::vector<double> yright;

            //mVariableNames.push_back("V");
            //mVariableUnits.push_back("mV");
            yright.push_back(v_singularity[i]-0.1);

            //mVariableNames.push_back("n");
            //mVariableUnits.push_back("");
            yright.push_back(0.325);

            //mVariableNames.push_back("h");
            //mVariableUnits.push_back("");
            yright.push_back(0.6);

            //mVariableNames.push_back("m");
            //mVariableUnits.push_back("");
            yright.push_back(0.05);

            std::vector<double> rhsright(yright.size());
            hh52_ode_system.EvaluateYDerivatives (0.0, yright, rhsright);


            std::vector<double> y_at_singularity;

            //mVariableNames.push_back("V");
            //mVariableUnits.push_back("mV");
            y_at_singularity.push_back(v_singularity[i]);

            //mVariableNames.push_back("n");
            //mVariableUnits.push_back("");
            y_at_singularity.push_back(0.325);

            //mVariableNames.push_back("h");
            //mVariableUnits.push_back("");
            y_at_singularity.push_back(0.6);

            //mVariableNames.push_back("m");
            //mVariableUnits.push_back("");
            y_at_singularity.push_back(0.05);

            std::vector<double> rhs_at_singularity(y_at_singularity.size());
            hh52_ode_system.EvaluateYDerivatives (0.0, y_at_singularity, rhs_at_singularity);

            for (int j=0; j<4; j++)
            {
                //std::cout << j << "\t" << rhsright[j] << "\t" << rhsleft[j] << "\t" << rhs_at_singularity[j] << std::endl;
                TS_ASSERT_DELTA((rhsright[j]+rhsleft[j])/2, rhs_at_singularity[j], 0.1);
            }
        }
    }
};

#endif /*TESTHODGKINHUXLEYSQUIDAXON1952ORIGINALODESYSTEM_HPP_*/
