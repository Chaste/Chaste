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
 * This tutorial runs a simulation on a whole rabbit heart mesh. Note that this
 * mesh is far too coarse for converged simulations, but provides a useful example.
 *
 * First include the headers, `MonodomainProblem` this time.
 */
#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "GenericMeshReader.hpp"
#include "SimpleStimulus.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

/** TEMPORARY FOR DEBUGGING */
///\todo #2739
#include <unistd.h>
#include <sys/resource.h>
//#include "Debug.hpp"

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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        double z = pNode->rGetLocation()[2];

        if (z <= 0.05)
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
    double GetMemoryUsage()
    {
       struct rusage rusage;
       getrusage( RUSAGE_SELF, &rusage );

       return (double)(rusage.ru_maxrss)/(1024);// Convert KB to MB
    }

    void TestMonodomain3dRabbitHeart()
    {
        /*
         * HOW_TO_TAG Cardiac/Problem definition
         * Read in a mesh from file, via `HeartConfig`.
         */
        HeartConfig::Instance()->SetMeshFileName("apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um",
                                                 cp::media_type::Axisymmetric);

//        HeartConfig::Instance()->SetMeshFileName("OxfordRabbitHeart_ascii",
//                                                         cp::media_type::Axisymmetric);


        /*
         * Specify the conductivity vector to use in the simulation. Since this is going to be
         * a monodomain simulation, we only specify intra-cellular conductivities.
         * Additionally, because this is an Axi-symmetric mesh then we must specify
         * conductivity along the sheet and normal directions to be the same (an exception will be
         * thrown if not).
         *
         * The 3rd entry would be different for an orthotropic mesh with fibre, sheet and
         * normal directions. For a simulation without fibre directions, there should be one value.
         */
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.19, 0.19));

        /*
         * Set the simulation duration, output directory, filename and VTK visualization.
         *
         * We have set the simulation duration to be very short here so this tutorial runs
         * quickly, increase it to see decent propagation of the wavefront.
         */
        HeartConfig::Instance()->SetSimulationDuration(2); //ms
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3dRabbitHeart");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        /* The ODE and PDE timesteps should be refined when using this code for real
         * scientific simulations. The values here are sufficient to ensure stability
         * in this case, but not sufficient for converged numerical behaviour.
        */
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.02, 0.1, 0.2);

        /*
         * Here we create an instance of our cell factory, which will tell the `MonodomainProblem`
         * class which action potential models to use at which nodes. The rest of the problem is set up
         * identically to [wiki:UserTutorials/Monodomain3dExample Monodomain3dExample].
         */
        RabbitHeartCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem( &cell_factory );
        monodomain_problem.SetWriteInfo();
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        /* We can access nodes in the mesh using a `NodeIterator`. Here, we check that each node
         * has not been assigned to bath, and throw an error if it has. This is not a particularly useful test,
         * but it does demonstrate the principle.
         * */

        AbstractTetrahedralMesh<3,3>* p_mesh = &(monodomain_problem.rGetMesh());

        /** \todo #2739
        double before_init = GetMemoryUsage();
        monodomain_problem.Initialise();
        double after_init = GetMemoryUsage();
        PRINT_VARIABLE(after_init - before_init);
        */
        for (AbstractTetrahedralMesh<3,3>::NodeIterator iter = p_mesh->GetNodeIteratorBegin();
                         iter != p_mesh->GetNodeIteratorEnd(); ++iter)
        {
            TS_ASSERT_EQUALS(HeartRegionCode::IsRegionBath( iter->GetRegion() ), false);
        }
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
 *
 * To visualize these results, see ChasteGuides/VisualisationGuides/ParaviewForCardiac. The colour axes in Paraview
 * may need to be rescaled in order to see the voltage changes.
 */

#endif /*TESTMONODOMAIN3DRABBITHEARTTUTORIAL_HPP_*/
