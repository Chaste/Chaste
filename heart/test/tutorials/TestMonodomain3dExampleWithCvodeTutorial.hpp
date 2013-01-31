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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTMONODOMAIN3DEXAMPLEWITHCVODETUTORIAL_HPP_
#define TESTMONODOMAIN3DEXAMPLEWITHCVODETUTORIAL_HPP_

/*
 * = 3D monodomain example using CVODE for ODE solution =
 *
 * This tutorial is based on UserTutorials/Monodomain3dExample except this time we will
 * use CVODE solvers. To highlight the changes needed to run with CVODE we omit the usual
 * explanations of the rest of the code - see UserTutorials/Monodomain3dExample for these.
 *
 * First include the headers
 */
#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "LuoRudy1991.hpp"
/* Chaste actually has two ways of using CVODE for solution of cardiac action potential model ODEs.
 *  1. via a `CvodeAdaptor` solver - this would work on the usual cell model that was included above.
 *  1. via an `AbstractCvodeCell` instead of an `AbstractCardiacCell` - this class uses native CVODE vectors.
 *
 * In order to generate CVODE cells please see ChasteGuides/CodeGenerationFromCellML.
 * Note that recent improvements (becoming available from release 3.2) mean that
 * maple can be used to generate an '''analytic jacobian''' which is then made available in the
 * native AbstractCvodeCell, and this will provide a speed up of between 5-30% (depending on the size of
 * the ODE system).
 *
 * So here we do the #include to import the native CVODE version of the cell model.
 */
#include "LuoRudy1991Cvode.hpp"
/* then include the rest of the headers as usual */
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"


class BenchmarkCellFactory : public AbstractCardiacCellFactory<3> // <3> here
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BenchmarkCellFactory()
        : AbstractCardiacCellFactory<3>(), // <3> here as well!
          mpStimulus(new SimpleStimulus(-100000.0, 2))
    {
    }

    /*
     * The following method definition changes to return an `AbstractCvodeCell`
     * instead of an `AbstractCardiacCell`.
     */
    AbstractCvodeCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        /*
         * Purely in order to maintain a consistent interface,
         * an AbstractCvodeCell expects an AbstractIvpOdeSolver in its
         * constructor, but it is not used (CVODE is instead). So an empty
         * pointer can be passed.
         */
        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;

        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];
        double z = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[2];

        if ( (x<0.1+1e-6) && (y<0.1+1e-6) && (z<0.1+1e-6) )
        {
            return new CellLuoRudy1991FromCellMLCvode(p_empty_solver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellMLCvode(p_empty_solver, mpZeroStimulus);
        }
    }
};

/*
 * The rest of the test is almost identical to the non-CVODE cell case,
 * - just note the #ifdef tag and the comment about ODE timesteps.
 */
class TestMonodomain3dExample : public CxxTest::TestSuite
{
public:
    void TestMonodomain3d() throw(Exception)
    {
        /*
         * Since CVODE is an optional extra dependency for Chaste - albeit now
         * one that is highly recommended - see InstallGuide.
         *
         * EMPTYLINE
         *
         * We guard any code that relies upon it with the following `#ifdef`.
         * This CHASTE_CVODE flag is set automatically if your hostconfig file
         * (in python/hostconfig) sets `use_cvode` and calls `DetermineCvodeVersion(<path to CVODE includes>)`.
         * See the end of the file python/hostconfig/default.py for an example of this.
         *
         */
#ifdef CHASTE_CVODE
        TetrahedralMesh<3,3> mesh;
        double h=0.02;
        mesh.ConstructRegularSlabMesh(h, 0.8 /*length*/, 0.3 /*width*/, 0.3 /*depth*/);

        HeartConfig::Instance()->SetSimulationDuration(5); //ms
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3dExample");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        /*
         * Note - when using CVODE in cardiac tissue simulations the ODE timestep
         * should be set to the same as the PDE timestep.
         *
         * CVODE will take as many adaptive internal timesteps as it requires each time
         * it is called, so we should just call it once per PDE timestep - i.e. set the
         * ODE and PDE timesteps to be the same.
         *
         * EMPTYLINE
         *
         * ''''NB'''': CVODE will only give you a big speedup when the ODE/PDE timestep is larger than
         * a typical Forward Euler timestep would be for that model. But it doesn't
         * seem to be any slower than Forward Euler, even at this PDE resolution.
         *
         * A convergence analysis should be performed to ensure that the PDE is being solved
         * accurately before reducing the step just to get faster ODE solution!
         */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);

        /*
         * The rest of the code is unchanged.
         */
        BenchmarkCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem( &cell_factory );

        monodomain_problem.SetMesh(&mesh);

        bool partial_output = false;
        if(partial_output)
        {
            std::vector<unsigned> nodes_to_be_output;
            nodes_to_be_output.push_back(0);
            nodes_to_be_output.push_back((unsigned)round( (mesh.GetNumNodes()-1)/2 ));
            nodes_to_be_output.push_back(mesh.GetNumNodes()-1);
            monodomain_problem.SetOutputNodes(nodes_to_be_output);
        }

        monodomain_problem.SetWriteInfo();


        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        ReplicatableVector voltage(monodomain_problem.GetSolution());
        /*
         * '''NB''': CVODE gives a more accurate ODE solution than Forward Euler, so this
         * result has been tweaked from previous tutorial (34.9032mV previously).
         */
        TS_ASSERT_DELTA(voltage[0], 34.8548, 1e-2);
        /*
         * Here we add a visual warning in case CVODE is not installed and/or set up.
         *
         * Otherwise the test passes, as it is not doing anything!
         */
#else
        std::cout << "CVODE is not installed, or CHASTE is not configured to use it, check your hostconfig settings." << std::endl;
#endif // CHASTE_CVODE
    }
};
#endif /*TESTMONODOMAIN3DEXAMPLEWITHCVODETUTORIAL_HPP_*/
