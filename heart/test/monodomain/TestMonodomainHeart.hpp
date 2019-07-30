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


#ifndef _TESTMONODOMAINHEART_HPP_
#define _TESTMONODOMAINHEART_HPP_

#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include <petscvec.h>
#include <vector>

#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "Hdf5DataReader.hpp"
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
    }

    void FinaliseCellCreation(std::vector<AbstractCardiacCellInterface* >* pCellsDistributed, unsigned lo, unsigned hi)
    {
        /* Here's the list of stimulated cells from the mesh file with tetgen numbering:
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
        */
        int stimulated_cells[] = {
                                     37484-1,
                                     37499-1,
                                     37777-1,
                                     37779-1,
                                     38008-1,
                                     38332-1,
                                     38587-1,
                                     38588-1,
                                     39312-1,
                                     39314-1,
                                     39643-1,
                                     40588-1,
                                     40590-1,
                                     63885-1
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
    void TestMonodomainDg0Heart()
    {
        ///////////////////////////////////////////////////////////////////////
        // Solve
        ///////////////////////////////////////////////////////////////////////
        double pde_time_step = 0.01;  // ms
        double end_time = 100;        // ms

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetPrintingTimeStep(end_time/100);
        HeartConfig::Instance()->SetPdeTimeStep(pde_time_step);
        HeartConfig::Instance()->SetOdeTimeStep(pde_time_step/4.0);
        HeartConfig::Instance()->SetSimulationDuration(end_time); //ms
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/UCSD_heart"); // note that this is the full heart mesh (not fifthheart)
        HeartConfig::Instance()->SetOutputDirectory("MonoDg0Heart");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_Heart");

        PointStimulusHeartCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        monodomain_problem.SetWriteInfo();

        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        ///////////////////////////////////////////////////////////////////////
        // now reread the data and check verify that one of the stimulated
        // nodes was actually stimulated, and that the propagation spread to
        // a nearby node
        ///////////////////////////////////////////////////////////////////////
        Hdf5DataReader data_reader = monodomain_problem.GetDataReader();

        // get the voltage values at stimulated node
        std::vector<double> voltage_values_at_node_37483 = data_reader.GetVariableOverTime("V", 37484-1);
        // get the voltage values at a nearby unstimulated node
        std::vector<double> voltage_values_at_node_500 = data_reader.GetVariableOverTime("V", 501-1);
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
    }
};

#endif //_TESTMONODOMAINHEART_HPP_
