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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        // Stimulate the apex
        if (GetMesh()->GetNode(node)->rGetLocation()[0] > 0.94)
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
private:
    void SetParameters()
    {
        HeartConfig::Instance()->Reset();
        //The conductivities were in the Metis test (not the plain test)
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.0025, 0.005, 0.1);
        HeartConfig::Instance()->SetSimulationDuration(100.0);  //ms

        HeartConfig::Instance()->SetKSPSolver("symmlq");
        HeartConfig::Instance()->SetKSPPreconditioner("bjacobi");
        HeartConfig::Instance()->SetOutputFilenamePrefix("BidomainLR91HalfHeart");
        PetscOptionsSetValue("-options_table", "");
    }
public:

    void TestBidomainDg0Heart() throw (Exception)
    {
        SetParameters();

        HeartConfig::Instance()->SetMeshFileName("heart/test/data/halfheart");
        HeartConfig::Instance()->SetOutputDirectory("BiDg0Heart");

        PointStimulusHeartCellFactory cell_factory;
        BidomainProblem<3> bidomain_problem(&cell_factory);

        bidomain_problem.SetWriteInfo();

        bidomain_problem.Initialise();
        bidomain_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};

#endif //_TESTBIDOMAINHEART_HPP_
