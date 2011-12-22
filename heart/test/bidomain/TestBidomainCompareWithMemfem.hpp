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


#ifndef TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_
#define TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <iostream>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "SimpleStimulus.hpp"

class BidomainPointStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
public:
    BidomainPointStimulusCellFactory()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-1000.0*1000, 1))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        if (node==19)
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};

class TestBidomainCompareWithMemfem :  public CxxTest::TestSuite
{
public:

    void TestBidomainCompareWithMemfemBasic()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.19, 0.19, 1.79));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(2.36, 2.36, 6.25));
        HeartConfig::Instance()->SetOdeTimeStep(0.001);
        HeartConfig::Instance()->SetSimulationDuration(50.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("heart/test/data/memfem_mesh/simple");
        HeartConfig::Instance()->SetOutputDirectory("Bidomain3d_CompareWithMemfem");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain3d");

        BidomainPointStimulusCellFactory bidomain_cell_factory;

        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        // set the back face (nodes 468-506) to have phi_e fixed to zero
        std::vector<unsigned> fixed_nodes;
        for (unsigned i=468;i<507;i++)
        {
            fixed_nodes.push_back(i);
        }
        bidomain_problem.SetFixedExtracellularPotentialNodes(fixed_nodes);

        bidomain_problem.SetWriteInfo();

        bidomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1500); //    1/cm

        try
        {
            TS_FAIL("Doesn't yet agree with Memfem");
         /*
         * We don't test anything, since we haven't managed to get memfem to agree
         * with chaste - probably because we couldn't find any identical cell models
         * to those we have, and issues setting identical stimuli.  (But we suspect
         * memfem and chaste won't agree anyway, and given all our other tests we
         * should probably assume that it's memfem that incorrect? dunno).
         */
            //bidomain_problem.Solve();
        }
        catch (Exception e)
        {
            std::cout << e.GetMessage() << "\n";
        }


    }
};


#endif /*TESTBIDOMAINCOMPAREWITHMEMFEM_HPP_*/
