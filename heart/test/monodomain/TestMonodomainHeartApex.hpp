/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef _TESTMONODOMAINHEARTAPEX_HPP_
#define _TESTMONODOMAINHEARTAPEX_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "SimpleStimulus.hpp"

class PointStimulusHeartCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
public:
    PointStimulusHeartCellFactory()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-1000*1000, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
    }

    void FinaliseCellCreation(std::vector<AbstractCardiacCell* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        /*
         * Here's the list of stimulated cells from the original mesh file with tetgen numbering:
37483   1.95075 0.02458 0.007709
37498   1.974   0.0669055       0.0212167
37776   1.92132 -0.0185282      0.0264612
37778   1.90362 -0.0457586      0.0653502
38007   1.93232 -0.0006106      0.0718023
38331   1.9199  -0.0110326      -0.0303119
38586   1.90586 -0.0489975      -0.0275832
38587   1.90054 -0.0704444      0.0187846
39311   1.95105 0.0306952       -0.0127931
39313   1.97209 0.072277        -0.0302588
39642   1.93571 0.0272909       -0.0672191
40587   1.95069 0.0286633       -0.0049338
40589   1.97168 0.0738751       -0.0122153
63884   1.93084 0       0
         * Here's the list of stimulated cells from the new apex mesh file with tetgen numbering:
875     1.950749999999999984    0.024580000000000001154 0.00770899999999999970680
890     1.9739999999999999769   0.066905500000000006633 0.02121670000000000153 0
1168    1.9213199999999999168   -0.018528200000000001613        0.026461200000000000693 0
1170    1.9036200000000000898   -0.045758600000000003438        0.06535019999999999718  0
1361    1.9323200000000000376   -0.00061059999999999998825      0.071802299999999999458 0
1660    1.91989999999999994     -0.011032600000000000046        -0.030311899999999999261        0
1795    1.9058600000000001096   -0.048997499999999999387        -0.027583199999999998692        0
1796    1.9005399999999998961   -0.070444400000000004236        0.018784599999999998521 0
2192    1.951049999999999951    0.030695199999999998874 -0.0127930999999999999580
2194    1.9720899999999998986   0.072276999999999994029 -0.0302587999999999988920
2452    1.9357100000000000417   0.027290899999999999881 -0.0672191000000000038470
3257    1.950690000000000035    0.028663299999999999196 -0.0049338000000000003381       0
3259    1.9716800000000000992   0.073875099999999999101 -0.0122153000000000001160
5109    1.9308399999999998897   0       0       0


 Failure occurs on node
  5045    1.4275700000000000056   0.24721699999999999231  0.05276329999999999909 0
[     */

        int stimulated_cells[] = {
                                    875,
                                    890,
                                    1168,
                                    1170,
                                    1361,
                                    1660,
                                    1795,
                                    1796,
                                    2192,
                                    2194,
                                    2452,
                                    3257,
                                    3259,
                                    5109
                                 };


        for (unsigned i=0; i<14; i++)
        {
            int global_index = stimulated_cells[i];
            if ((global_index>=(int)lo) && (global_index<(int)hi))
            {
                int local_index = global_index - lo;
                (*pCellsDistributed)[ local_index ]->SetStimulusFunction(mpStimulus);
            }
        }
    }
};

class TestMonodomainHeart : public CxxTest::TestSuite
{

public:
    void TestMonodomainDg0Heart() throw(Exception)
    {
        ///////////////////////////////////////////////////////////////////////
        // Solve
        ///////////////////////////////////////////////////////////////////////
        double end_time = 10.0;        // ms

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetPrintingTimeStep(end_time/100);
        HeartConfig::Instance()->SetPdeTimeStep(0.01);
        HeartConfig::Instance()->SetOdeTimeStep(0.01/3.0);
        HeartConfig::Instance()->SetSimulationDuration(end_time);
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/HeartApex"); // note that this is the full heart mesh (not fifthheart)
        HeartConfig::Instance()->SetOutputDirectory("MonoDg0HeartApex");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_HeartApex");

        PointStimulusHeartCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetWriteInfo();

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        ///////////////////////////////////////////////////////////////////////
        ///todo: now reread the data and check verify that one of the stimulated
        ///nodes was actually stimulated, and that the propagation spread to
        ///a nearby node
        ///////////////////////////////////////////////////////////////////////
        /*
         *
        ColumnDataReader data_reader("MonoDg0Heart","MonodomainLR91_Heart");

        // get the voltage values at stimulated node
        std::vector<double> voltage_values_at_node_37483 = data_reader.GetValues("V", 37484-1);
        // get the voltage values at a nearby unstimulated node
        std::vector<double> voltage_values_at_node_500 = data_reader.GetValues("V", 501-1);
        bool stimulated_node_was_excited = false;
        bool unstimulated_node_was_excited = false;

        for (unsigned i=0; i<voltage_values_at_node_37483.size(); i++)
        {
            if (voltage_values_at_node_37483[i] > 0)
            {
                stimulated_node_was_excited = true;
            }
            if (voltage_values_at_node_500[i] > 0)
            {
                unstimulated_node_was_excited = true;
            }
        }
        TS_ASSERT(stimulated_node_was_excited);
        TS_ASSERT(unstimulated_node_was_excited);
        *
        */
    }
};

#endif //_TESTMONODOMAINHEARTAPEX_HPP_
