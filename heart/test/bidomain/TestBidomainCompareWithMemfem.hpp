/*

Copyright (c) 2005-2012, University of Oxford.
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
