/*

Copyright (c) 2005-2026, University of Oxford.
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


#ifndef _TESTBIDOMAINHEART_HPP_
#define _TESTBIDOMAINHEART_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>


#include "BidomainProblem.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "TrianglesMeshReader.hpp"
#include "SimpleStimulus.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class PointStimulusHeartCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
public:
    PointStimulusHeartCellFactory()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-1000.0*500, 0.5))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        // Stimulate the apex
        if (pNode->rGetLocation()[0] > 0.94)
        {
            return new CellLuoRudy1991FromCellML(mpSolver,mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver,mpZeroStimulus);
        }
    }
};

class TestBidomainHeart : public CxxTest::TestSuite
{
public:

    void TestBidomainDg0Heart()
    {
        PointStimulusHeartCellFactory cell_factory;
        BidomainProblem<3> bidomain_problem(&cell_factory);

        //The conductivities were in the Metis test (not the plain test)
        bidomain_problem.SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        bidomain_problem.SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        bidomain_problem.SetOdePdeAndPrintingTimeSteps(0.0025, 0.005, 0.1);
        bidomain_problem.SetSimulationDuration(100.0);  //ms
        bidomain_problem.SetKspSolverType("symmlq");
        bidomain_problem.SetKspPreconditionerType("bjacobi");
        bidomain_problem.SetOutputFilenamePrefix("BidomainLR91HalfHeart");
        bidomain_problem.SetMeshFileName("heart/test/data/scaled_UCSD_heart");
        bidomain_problem.SetOutputDirectory("BiDg0Heart");

        PetscTools::SetOption("-options_table", "");

        bidomain_problem.SetWriteInfo();

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};

#endif //_TESTBIDOMAINHEART_HPP_
