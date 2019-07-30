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
#ifndef TESTDYNAMICVENTILATIONTUTORIAL_HPP_
#define TESTDYNAMICVENTILATIONTUTORIAL_HPP_

/* HOW_TO_TAG Lung/Simulation
 * Simulate ventilation on a full lung geometry with acinar dynamics over a breathing cycle.
 */

/*
 * = An example showing how to simulate ventilation over a full breathing cycle. =
 *
 * In this tutorial we demonstrate the use of !DynamicVentilationProblem to simulate breathing over a
 * single breathing cycle. Acini are represented by linear balloon models and the simulation is driven
 * by an oscillating pleural pressure. Simulation results are written to a VTK unstructured grid file
 * for visualisation.
 *
 * NB. UMFPACK or KLU is required for this test to execute. The test has a reasonably long run time (~6 minutes on a quad core i7).
 * Progress can be followed by watching the file $CHASTE_TEST_OUTPUT/TestDynamicVentilationTutorial/progress_status.txt
 */

/* The usual headers are included */
#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"

/* !DynamicVentilationProblem does most of the work in calculating a ventilation distribution. */
#include "DynamicVentilationProblem.hpp"

/* A number of acinar models could be used. Here we include the simplest possible model: a linear elastic balloon. */
#include "SimpleBalloonAcinarUnit.hpp"

/* Note that this tutorial only works with UMFPACK or KLU -- we need to warn the user if it's not installed. */
#include "Warnings.hpp"

/* !DynamicVentilationProblem uses the Petsc solver library. This setups up Petsc ready for use. */
#include "PetscSetupAndFinalize.hpp"

/*
 * == Acinar Unit Factory ==
 *
 * !DynamicVentilationProblem couples a ventilation model on an airway tree to acinar dynamics models. The
 * distal ends of the airway tree have an acinar model associated with them. These models are specified
 * by subclassing !AbstractAcinarUnitFactory as shown here.
 */
class SimpleAcinarUnitFactory : public AbstractAcinarUnitFactory
{
public:
    SimpleAcinarUnitFactory(double acinarCompliance,
                            double frequency = 0.5) : mAcinarCompliance(acinarCompliance),
                                                      mFrequency(frequency)
    {}

    /*
     * !DynamicVentilationProblem calls this method once for each terminal node.
     * The user must create, configure and return an Acinar model. The created
     * model is subsequently used in the ventilation simulation.
     */
    virtual AbstractAcinarUnit* CreateAcinarUnitForNode(Node<3>* pNode)
    {
        /*
         * Here we use the simplest possible acinar model: a linear elastic balloon.
         */
        SimpleBalloonAcinarUnit* p_acinus = new SimpleBalloonAcinarUnit;

        /*
         * The acinar model can be configured in different ways.
         * Here we simply set a constant spatially homogeneous compliance and
         * a spatially homogeneous initial volume. However,
         * things like gravitational gradients can easily be configured here by
         * altering the compliance (or other parameters) based on the nodal
         * location (obtained from pNode->rGetLocation()).
         */
        p_acinus->SetCompliance(mAcinarCompliance);
        p_acinus->SetUndeformedVolume(mAcinarCompliance*500);


        return p_acinus;
    }

    /*
     * !DynamicVentilationProblems are driven by a change in Pleural pressure. !DynamicVentilationProblem
     * calls this method to determine what the current pleural pressure is. Here
     * we specify a spatially homogeneous oscillating pressure in the tidal breathing range.
     */
    virtual double GetPleuralPressureForNode(double time, Node<3>* pNode)
    {
        return -750 - 250*sin(2*M_PI*mFrequency*(time - 0.25));
    }

private:
    double mAcinarCompliance;
    double mFrequency;
};


/* Define the test */
class TestDynamicVentilationTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestSimulateTidalBreathing()
    {
        /*
         * !DynamicVentilationProblem is not (yet) parallel.
         */
        EXIT_IF_PARALLEL;

        /*
         * == IMPORTANT ==
         * See the note below about use of direct solvers. This tutorial cannot be run without a direct solver.
         */
#if defined(LUNG_USE_UMFPACK) || defined(LUNG_USE_KLU)

        /* First we need to create an acinar unit factory object from the class we specified earlier.
         * An acinar compliance is specified in Pa/m^3^. The given value is roughly equal to a lung compliance
         * of 0.1 cmH2O/L (assuming there are 30000 acini).
         */
        double acinar_compliance = 0.1/98.0665/1e3/30000;
        SimpleAcinarUnitFactory factory(acinar_compliance);

        /*
         * We now create a !DynamicVentilationProblem object that does most of the work in simulating ventilation.
         * The factory we just created is passed to the constructor along with the location of an airways mesh.
         */
        DynamicVentilationProblem problem(&factory, "lung/test/data/simplified_airways", 0u);

        /* We assign a zero pressure boundary condition at the entrance to the trachea. */
        problem.rGetMatrixVentilationProblem().SetOutflowPressure(0.0);

        /* The mesh we are using is specified in millimetres rather than in metres. */
        problem.rGetMatrixVentilationProblem().SetMeshInMilliMetres();

        /* The mesh we are using specifies airway radii on the edges rather than the nodes. */
        problem.rGetMatrixVentilationProblem().SetRadiusOnEdge();

        /* Tell the solver to use a more accurate Pedley based dynamic resistance scheme. */
        problem.rGetMatrixVentilationProblem().SetDynamicResistance();

        /* Here we tell the solver the time step size to use. The given value will typically
         * result in a stable numerical scheme. However, users should assess the time step
         * size required to achieve a suitable level of numerical convergence when setting
         * up their own simulations.
         */
        problem.SetTimeStep(0.02);

        /* Tell the solver where to write its output to.
         * The solver will also write out a progress_status.txt file to
         * this directory to allow the user to monitor progress.
         * We specify output in VTK format for
         * easy visualisation. */
        problem.SetOutputDirectory("TestDynamicVentilationTutorial");
        problem.SetOutputFilenamePrefix("tidal_breathing");
        problem.SetWriteVtkOutput();

        /* Tell the solver how often to write output. Here we ask for output every 5 time steps. */
        problem.SetSamplingTimeStepMultiple(5u);

        /* Specify when to end the simulation (in seconds) and solve. Note that
         * after solving it is possible to set a new end time and solve again if needed.
         *
         * For execution speed reasons we only simulate an inspiration. Try setting a later
         * end time to see a full breathing cycle.
         */
        problem.SetEndTime(1.0);
        problem.Solve();

        /*
         * It is now possible to analyse the data produced by the ventilation simulation. Typically
         * this will be done using the output written to disk and an external program. This simulation
         * will output a file in $CHASTE_TEST_OUTPUT/TestDynamicVentilationTutorial/tidal_breathing.vtu for easy visualisation.
         * For demonstration purposes we perform a simple check of the final lung volume here.
         */
        std::map<unsigned, AbstractAcinarUnit*>& acinar_map = problem.rGetAcinarUnitMap();
        double lung_volume = 0;

        for (std::map<unsigned, AbstractAcinarUnit*>::iterator iter = acinar_map.begin();
             iter != acinar_map.end();
             ++iter)
        {
            lung_volume += iter->second->GetVolume();
        }

        std::cout << "The total lung volume at the end of the simulation is " << lung_volume*1e3 << " L. " << std::endl;
#else
        WARNING("Not compiled with a direct solver. Test not executed.");
#endif
    }

    /*
     *
     * == IMPORTANT: Using UMFPACK/KLU ==
     *
     * Ventilation problems lead to very badly conditioned matrices. Iterative solvers such as GMRES can stall on these
     * matrices. When running problems on large airway trees it is vital that to change the linear solver to a direct
     * solver such as UMFPACK or KLU. UMFPACK and KLU are not pre-requisites for installing Chaste, hence this is not (currently)
     * the default linear solver for ventilation problems.
     *
     * ''UMFPACK or KLU should be considered pre-requisites for large ventilation problems''
     *
     * To use UMFPACK or KLU, you need to have PETSc installed with UMFPACK/KLU.
     *
     * To switch on UMFPACK or KLU on within chaste, set "ccflags='-DLUNG_USE_UMFPACK'" or
     * "ccflags='-DLUNG_USE_KLU'" in your local.py or
     * open the file `lung/src/ventilation/MatrixVentilationProblem.hpp` and uncomment the line
     * #define LUNG_USE_UMFPACK near the top of the file.
     *
     */
};

#endif /*TESTDYNAMICVENTILATIONTUTORIAL_HPP_*/


