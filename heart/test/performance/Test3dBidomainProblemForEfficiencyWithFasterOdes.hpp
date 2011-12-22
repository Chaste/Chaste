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


#ifndef TEST3DBIDOMAINFOREFFICIENCYWITHFASTERODES_HPP_
#define TEST3DBIDOMAINFOREFFICIENCYWITHFASTERODES_HPP_




#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "BidomainProblem.hpp"
#include <petscvec.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "RegularStimulus.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "SimpleStimulus.hpp"

class BidomainFaceStimulusCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    boost::shared_ptr<RegularStimulus> mpRegStimulus;

public:
    //Pdetime step is (by default) 0.01
    //Odetime step set below to 0.01 as backward Euler should be stable
    BidomainFaceStimulusCellFactory()
        : AbstractCardiacCellFactory<3>(),
          mpStimulus(new SimpleStimulus(-900.0*1000, 0.5)),
          mpRegStimulus(new RegularStimulus(-900.0*1000, 0.5, 100.0, 0.0))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned node)
    {
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        if (GetMesh()->GetNode(node)->GetPoint()[0] == 0.0)
        {
            //std::cout << node+1 << "\n";
            return new CellLuoRudy1991FromCellMLBackwardEuler(p_solver, mpRegStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellMLBackwardEuler(p_solver, mpZeroStimulus);
        }
    }
};

class Test3dBidomainProblemForEfficiencyWithFasterOdes :  public CxxTest::TestSuite
{
public:

    void TestBidomain3d() throw (Exception)
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 1.75, 1.75));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 7.0, 7.0));
        HeartConfig::Instance()->SetSimulationDuration(150.0);
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/3D_0_to_.5mm_1889_elements_irregular");

        BidomainFaceStimulusCellFactory bidomain_cell_factory;

        BidomainProblem<3> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.PrintOutput(false);

        PetscOptionsSetValue("-log_summary", "");

        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        Vec voltage=bidomain_problem.GetSolution();
        ReplicatableVector voltage_replicated;
        voltage_replicated.ReplicatePetscVector(voltage);

        /*
         * Test the top right node against the right one in the 1D case,
         * comparing voltage, and then test all the nodes on the right hand
         * face of the cube against the top right one, comparing voltage.
         */
        bool need_initialisation = true;
        double probe_voltage=-9999.;

        need_initialisation = true;

        // Test the RHF of the mesh
        for (unsigned i = 0; i < bidomain_problem.rGetMesh().GetNumNodes(); i++)
        {
            if (bidomain_problem.rGetMesh().GetNode(i)->GetPoint()[0] == 0.05)
            {
                // x = 0 is where the stimulus has been applied
                // x = 0.05cm is the other end of the mesh and where we want to
                //       to test the value of the nodes

                if (need_initialisation)
                {
                    probe_voltage = voltage_replicated[2*i];
                    need_initialisation = false;
                }
                else
                {
                    // the voltage at the end face varies a little because
                    // of drift due to the orientation of the tets in the mesh,
                    // hence the tolerance of 0.02
                    TS_ASSERT_DELTA(voltage_replicated[2*i], probe_voltage, 0.02);
                }

                // Check against hard coded value
                // For 50 ms test TS_ASSERT_DELTA(voltage_replicated[2*i],  7.3, 0.2);
                // For 150 ms test
                TS_ASSERT_DELTA(voltage_replicated[2*i], -1.735, 0.1);
            }
        }
    }
};


#endif /*TEST3DBIDOMAINFOREFFICIENCYWITHFASTERODES_HPP_*/
