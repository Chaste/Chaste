/*

Copyright (c) 2005-2014, University of Oxford.
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
#ifndef TESTMONODOMAIN3DRABBITHEARTTUTORIAL_HPP_
#define TESTMONODOMAIN3DRABBITHEARTTUTORIAL_HPP_

/*
 * = 3D monodomain rabbit heart example =
 *
 * In this tutorial we show how to run a 3D simulation using the monodomain equation.
 * To go from monodomain to bidomain or vice versa is trivial, and for 2d to 3d is
 * also very easy. This tutorial runs a simulation on a whole heart mesh. Note that this
 * mesh is far too coarse for converged simulations, but provides a useful example.
 *
 * First include the headers, `MonodomainProblem` this time.
 */
#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "GenericMeshReader.hpp"
#include "DistributedTetrahedralMesh.hpp"

/* Here we define a cell factory that gives stimuli to all cells
 * below height z = 0.042... this corresponds to the apex of the heart.
 */
class RabbitHeartCellFactory : public AbstractCardiacCellFactory<3> // <3> here
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    RabbitHeartCellFactory()
        : AbstractCardiacCellFactory<3>(), // <3> here as well!
          mpStimulus(new SimpleStimulus(-80000.0, 2))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double z = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[2];

        if (z <= 0.04248)
        {
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellMLBackwardEuler(mpSolver, mpZeroStimulus);
        }
    }
};

/* Now define the test */
class TestMonodomain3dRabbitHeartTutorial : public CxxTest::TestSuite
{
public:
    void TestMonodomain3dRabbitHeart() throw(Exception)
    {
        /*
         * HOW_TO_TAG Cardiac/Problem definition
         * Read in a mesh from file, via `HeartConfig`.
         */
        HeartConfig::Instance()->SetMeshFileName("apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um",
                                                 cp::media_type::Axisymmetric);

        /*
         * Specify the conductivity vector to use in the simulation. Since this is going to be
         * a monodomain simulation, we only specify intra-cellular conductivities.
         * And additionally, because this is an Axi-symmetric mesh then we must specify
         * conductivity along sheet and normal directions to be the same (an exception will be
         * thrown if not).
         *
         * The 3rd entry would be different for an orthotropic mesh with fibre, sheet and
         * normal directions there would just be one entry for no fibre directions.
         */
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.19, 0.19));

        /*
         * Set the simulation duration, and output directory etc.
         *
         * We have set the simulation duration to be very short here so this tutorial runs
         * fast, increase it to see decent propagation of the wavefront.
         */
        HeartConfig::Instance()->SetSimulationDuration(1); //ms
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3dRabbitHeart");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        /* The ODE and PDE timesteps should be refined when using this code for real
         * scientific simulations. The values here are sufficient to ensure stability
         * in this case, but not converged numerical behaviour.
        */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.02, 0.1, 0.5);

        /*
         * Here we create an instance of our cell factory, which will tell the `MonodomainProblem`
         * class which action potential models to use at which nodes.
         */
        RabbitHeartCellFactory cell_factory;

        /* Now we declare the problem class, `MonodomainProblem<3>` instead of `BidomainProblem<2>`.
         * The interface for both is the same.
         */
        MonodomainProblem<3> monodomain_problem( &cell_factory );

        /* `SetWriteInfo` is a useful method that means that the min/max voltage is
         * printed as the simulation runs (useful for verifying that cells are stimulated
         * and the wave propagating, for example) (although note scons does buffer output
         * before printing to screen) */
        monodomain_problem.SetWriteInfo();

        /* the problem must be initialised to set up the mesh, cell models and PDE. */
        monodomain_problem.Initialise();
        /*
         * Call Solve to run the simulation.
         */
        monodomain_problem.Solve();
    }
};

/*
 * '''Note''' if you were doing a 'real' scientific simulation you would want to use a higher
 * resolution mesh. A version of this can be found on the [http://www.cs.ox.ac.uk/chaste/download.html Chaste download website]
 *
 * Navigate to the "Data" tab, and download either
 *  * [source:/data/public/OxfordRabbitHeart/OxfordRabbitHeart_binary.tgz OxfordRabbitHeart_binary.tgz]  - 599MB, or
 *  * [source:/data/public/OxfordRabbitHeart/OxfordRabbitHeartWithBath_binary.tgz OxfordRabbitHeartWithBath_binary.tgz]  - 846MB.
 *
 * These will probably require HPC resources, and finer ODE and PDE time steps than we used here.
 */

#endif /*TESTMONODOMAIN3DRABBITHEARTTUTORIAL_HPP_*/
