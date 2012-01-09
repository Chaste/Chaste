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

#ifndef TESTOPERATORSPLITTINGMONODOMAINSOLVER_HPP_
#define TESTOPERATORSPLITTINGMONODOMAINSOLVER_HPP_


#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>
#include "MonodomainProblem.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "AbstractCardiacCellFactory.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "MonodomainSolver.hpp"
#include "TenTusscher2006Epi.hpp"
#include "Mahajan2008.hpp"

// stimulate a block of cells (an interval in 1d, a block in a corner in 2d)
template<unsigned DIM>
class BlockCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BlockCellFactory()
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(-1000000.0, 0.5))
    {
        assert(DIM<3);
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];

        if (fabs(x)<0.02+1e-6)
        {
            return new CellLuoRudy1991FromCellMLBackwardEuler(this->mpSolver, this->mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellMLBackwardEuler(this->mpSolver, this->mpZeroStimulus);
        }
    }
};


class TestOperatorSplittingMonodomainSolver : public CxxTest::TestSuite
{
public:

    // The operator splitting and normal methods should agree closely with very small dt and h, but this takes
    // too long to run in the continuous build (see instead TestOperatorSplittingMonodomainSolverLong)
    //
    // Here we run on a fine (as opposed to v fine) mesh and with a normal dt, and check that the solutions
    // are near.
    void TestComparisonOnNormalMeshes() throw(Exception)
    {
        ReplicatableVector final_voltage_normal;
        ReplicatableVector final_voltage_operator_splitting;

        HeartConfig::Instance()->SetSimulationDuration(4.0); //ms
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.1);

        double h = 0.01;

        // Normal
        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, 1.0);
            HeartConfig::Instance()->SetOutputDirectory("MonodomainCompareWithOperatorSplitting_normal");
            BlockCellFactory<1> cell_factory;

            MonodomainProblem<1> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
            monodomain_problem.Solve();

            final_voltage_normal.ReplicatePetscVector(monodomain_problem.GetSolution());
        }

        // Operator splitting
        {
            TetrahedralMesh<1,1> mesh;
            mesh.ConstructRegularSlabMesh(h, 1.0);
            HeartConfig::Instance()->SetOutputDirectory("MonodomainCompareWithOperatorSplitting_splitting");
            BlockCellFactory<1> cell_factory;

            HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplitting();

            MonodomainProblem<1> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(&mesh);
            monodomain_problem.Initialise();
            monodomain_problem.Solve();

            final_voltage_operator_splitting.ReplicatePetscVector(monodomain_problem.GetSolution());
        }

        // hardcoded value to check nothing has changed
        TS_ASSERT_DELTA(final_voltage_operator_splitting[30], 5.0577, 1e-3);

        bool some_node_depolarised = false;
        assert(final_voltage_normal.GetSize()==final_voltage_operator_splitting.GetSize());
        for(unsigned j=0; j<final_voltage_normal.GetSize(); j++)
        {
            // this tolerance means the wavefronts are not on top of each other, but not too far
            // separated (as otherwise max difference between the voltages across space would be
            // greater than 80).
            double tol=25;

            TS_ASSERT_DELTA(final_voltage_normal[j], final_voltage_operator_splitting[j], tol);

            if(final_voltage_normal[j]>-80)
            {
                // shouldn't be exactly equal, as long as away from resting potential
                TS_ASSERT_DIFFERS(final_voltage_normal[j], final_voltage_operator_splitting[j]);
            }

            if(final_voltage_normal[j]>0.0)
            {
                some_node_depolarised = true;
            }
        }
        assert(some_node_depolarised);
    }
};

#endif /* TESTOPERATORSPLITTINGMONODOMAINSOLVER_HPP_ */
