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
#ifndef TESTSTATICVENTILATIONTUTORIAL_HPP_
#define TESTSTATICVENTILATIONTUTORIAL_HPP_

/* HOW_TO_TAG Lung/Simulation
 * Calculate ventilation distribution in an airway tree for a given flow rate at the trachea
 */

/*
 * = An example showing how to calculate ventilation distribution in an airway tree for a given flow rate at the trachea =
 *
 * In this tutorial we demonstrate the use of !MatrixVentilationProblem to calculate airflow distribution in an airway
 * tree model. Homogeneous pressure boundary conditions are used at the terminals of the tree and a flow boundary condition
 * is used at the trachea. We demonstrate how to calculate the total bronchial pressure drop at different flow rates.
 *
 */

/* The usual headers are included */
#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"

/* !MatrixVentilationProblem does most of the work in calculating a ventilation distribution. */
#include "MatrixVentilationProblem.hpp"

/* Note that this tutorial only works with UMFPACK or KLU -- we need to warn the user if it's not installed */
#include "Warnings.hpp"

/* !MatrixVentilationProblem uses the Petsc solver library. This setups up Petsc ready for use. */
#include "PetscSetupAndFinalize.hpp"

/* Define the test */
class TestStaticVentilationTutorial : public CxxTest::TestSuite
{
public: // Tests should be public!

    void TestCalculatePressureDrop()
    {
        EXIT_IF_PARALLEL;

        /* First we setup a !MatrixVentilationProblem and tell it to load an airway centerline mesh.
         *
         * Note that the airway centerline  mesh cannot have intermediate nodes/elements within an airway.
         * Typically intermediate nodes are found in imaging derived airways centerlines. These can be removed
         * using the !AirwayRemesher class.
         *
         * == IMPORTANT ==
         * See the note below about use of UMFPACK or KLU. If UMFPACK or KLU are not available we use a different (not realistic) airways
         * mesh.
         */
#if defined(LUNG_USE_UMFPACK) || defined(LUNG_USE_KLU)
        MatrixVentilationProblem problem("lung/test/data/simplified_airways", 0u);
#else
        WARNING("Not compiled with UMFPACK or KLU.  Using non-realistic airway tree.");
        MatrixVentilationProblem problem("mesh/test/data/y_branch_3d_mesh", 0u);
#endif

        /* Matrix ventilation problem uses SI units but the mesh is specified in mm. This method allows the solver
         * to handle this discrepancy. */
        problem.SetMeshInMilliMetres();

        /* Airway meshes can have radii defined either on nodes (default) or on edges. Here we tell the solver
         * that the mesh has radii defined on edges.
         */
        problem.SetRadiusOnEdge();

        /* The ventilation solver can use different flow models. By default it uses the simplest Poiseuille flow
         * model. Whilst this is quick to calculate the results are not accurate for high flow rates. Here we
         * tell the solver to use a more accurate flow model based on Pedley 1970. Whilst this is more accurate
         * simulations take longer!
         */
        problem.SetDynamicResistance();

        /* We define some tracheal flow rates to test. The values are specified in m^3/s
         * These give a range from 10 L/min to 100 L/min.
         */
        std::vector<double> flows;
        flows.push_back(0.00017);
        flows.push_back(0.00083);
        flows.push_back(0.00167);
        flows.push_back(0.003);

        /* Loop over the tracheal flow rates and solve */
        for (unsigned i = 0; i < flows.size(); ++i)
        {
            /* This sets the boundary conditions: flow at the trachea and homogeneous zero pressure at the terminal airways */
            problem.SetOutflowFlux(flows[i]);
            problem.SetConstantInflowPressures(0.0);

            /* Calculates flow on the tree */
            problem.Solve();

            /* Here we obtain vectors containing the calculated pressures and fluxes at all nodes and elements in the mesh.
             */
            std::vector<double> flux, pressure;
            problem.GetSolutionAsFluxesAndPressures(flux, pressure);

            /* The vectors can be used to analyse the flow distribution in detail.
             * However, for this example we just output the total bronchial pressure drop between the trachea and the terminal
             * airways.
             */
            std::cout << "Total bronchial pressure drop for a tracheal flow rate of " << flows[i] << " m^3/s is " << pressure[0] << " Pa.\n";


           // cout << pressure[1] << endl;


        }
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

#endif /*TESTSTATICVENTILATIONTUTORIAL_HPP_*/
