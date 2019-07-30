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
#ifndef TESTSOLVINGELASTICITYPROBLEMSTUTORIAL_HPP_
#define TESTSOLVINGELASTICITYPROBLEMSTUTORIAL_HPP_

/*
 * = Solving solid mechanics problems =
 *
 * In this tutorial we show how Chaste can be used to solve solid mechanics problems.
 * We assume the reader has some familiarity with solid mechanics problems (the
 * equations of nonlinear elasticity are given in the PDF on equations
 * and finite element implementations (see ChasteGuides --> Miscellaneous information)). It is also best
 * to have had a look at the solving linear PDEs tutorials.
 *
 * In brief, there several facets to solid mechanics models:
 *  * Time-dependent problems versus static problems
 *  * Linear elasticity versus nonlinear elasticity
 *  * Compressible versus incompressible materials
 *  * The type of material behaviour (elastic, visco-elastic, etc..)
 *  * Specification of geometry, material law, body force, displacement boundary conditions, and traction boundary conditions
 *
 * The solvers currently implemented are STATIC (time-independent) and use NONLINEAR ELASTICITY. There are solvers for
 * both INCOMPRESSIBLE and COMPRESSIBLE deformations. The material behaviour is
 * assumed to be ELASTIC (stress is just a function of strain, not strain-rate etc), and in particular HYPER-ELASTIC
 * (stress is a function of strain via a 'strain energy function', for which stress is obtained by differentiating the
 * strain energy function with respect to strain).
 *
 * To solve a mechanics problem we need to
 *  * Choose the solver (compressible or incompressible)
 *  * Specify the geometry (ie the mesh)
 *  * Specify the material law (ie the strain-energy function) (incompressible or compressible as appropriate).
 *  * Specify the BODY FORCE -- this is a force density acting throughout the body (eg. acceleration due to gravity),
 *  and also the mass density.
 *  * Specify some DISPLACEMENT BOUNDARY CONDITIONS -- some part of the boundary must have the displacement specified on it
 *  * Specify TRACTION BOUNDARY CONDITIONS (if non-zero) on the rest of the boundary -- tractions are forces per unit area applied
 *  the rest of the surface of the deformable object.
 *
 *  '''VERY IMPORTANT NOTE:''' For incompressible problems, make sure you read the comment about HYPRE below before going to
 *  3D or refining meshes in these tests.
 *
 * '''Another note:''' mechanics problems are not currently implemented to scale in parallel yet.
 *
 * As always we include this first class as a test suite */
#include <cxxtest/TestSuite.h>
/* On some systems there is a clash between Boost Ublas includes and PETSc.  This can be
 * resolved by making sure that Chaste's interface to the Boost libraries are included
 * as early as possible.
 */
#include "UblasCustomFunctions.hpp"
/* The incompressible solver is called `IncompressibleNonlinearElasticitySolver` */
#include "IncompressibleNonlinearElasticitySolver.hpp"
/* The simplest incompressible material law is the Mooney-Rivlin material law (of which
 * Neo-Hookean laws are a subset)  */
#include "MooneyRivlinMaterialLaw.hpp"
/* Another incompressible material law */
#include "ExponentialMaterialLaw.hpp"
/* This is a useful helper class */
#include "NonlinearElasticityTools.hpp"
/* For visualising results in Paraview */
#include "VtkNonlinearElasticitySolutionWriter.hpp"
/* As before: !PetscSetupAndFinalize.hpp must be included in every test that uses PETSc. Note that it
 * cannot be included in the source code. */
#include "PetscSetupAndFinalize.hpp"


/*
 * HOW_TO_TAG Continuum mechanics
 * Solve nonlinear elasticity problems
 */

/*
 *
 * == Simple incompressible deformation: 2D shape hanging under gravity ==
 *
 */
class TestSolvingElasticityProblemsTutorial : public CxxTest::TestSuite
{
public:
    /* In the first test we use INCOMPRESSIBLE nonlinear elasticity. For incompressible elasticity, there
     * is a constraint on the deformation, which results in a pressure field (a Lagrange multiplier)
     * which is solved for together with the deformation.
     *
     * All the mechanics solvers solve for the deformation using the finite element method with QUADRATIC
     * basis functions for the deformation. This necessitates the use of a `QuadraticMesh` - such meshes have
     * extra nodes that aren't vertices of elements, in this case midway along each edge. (The displacement
     * is solved for at ''each node'' in the mesh (including internal [non-vertex] nodes), whereas the pressure
     * is only solved for at each vertex - in FEM terms, quadratic interpolation for displacement, linear
     * interpolation for pressure, which is required for stability. The pressure at internal nodes is computed
     * by linear interpolation).
     *
     */
    void TestSimpleIncompressibleProblem()
    {
        /* First, define the geometry. This should be specified using the `QuadraticMesh` class, which inherits from `TetrahedralMesh`
         * and has mostly the same interface. Here we define a 0.8 by 1 rectangle, with elements 0.1 wide.
         * (`QuadraticMesh`s can also be read in using `TrianglesMeshReader`; see next tutorial/rest of code base for examples of this).
         */
        QuadraticMesh<2> mesh;
        mesh.ConstructRegularSlabMesh(0.1 /*stepsize*/, 0.8 /*width*/, 1.0 /*height*/);

        /* We use a Mooney-Rivlin material law, which applies to isotropic materials and has two parameters.
         * Restricted to 2D however, it only has one parameter, which can be thought of as the total
         * stiffness. We declare a Mooney-Rivlin law, setting the parameter to 1.
         */
        MooneyRivlinMaterialLaw<2> law(1.0);

        /* Next, the body force density. In realistic problems this will either be
         * acceleration due to gravity (ie b=(0,-9.81)) or zero if the effect of gravity can be neglected.
         * In this problem we apply a gravity-like downward force.
         */
        c_vector<double,2> body_force;
        body_force(0) =  0.0;
        body_force(1) = -2.0;


        /* Two types of boundary condition are required: displacement and traction. As with the other PDE solvers,
         * the displacement (Dirichlet) boundary conditions are specified at nodes, whereas traction (Neumann) boundary
         * conditions are specified on boundary elements.
         *
         * In this test we apply displacement boundary conditions on one surface of the mesh, the upper (Y=1.0) surface.
         * We are going to specify zero-displacement at these nodes.
         * We do not specify any traction boundary conditions, which means that (effectively) zero-traction boundary
         * conditions (ie zero pressures) are applied on the three other surfaces.
         *
         * We need to get a `std::vector` of all the node indices that we want to fix. The `NonlinearElasticityTools`
         * has a static method for helping do this: the following gets all the nodes for which Y=1.0. The second
         * argument (the '1') indicates Y . (So, for example, `GetNodesByComponentValue(mesh, 0, 10)` would get the nodes on X=10).
         */
        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh, 1, 1.0);

        /*
         * Before creating the solver we create a `SolidMechanicsProblemDefinition` object,  which contains
         * everything that defines the problem: mesh, material law, body force,
         * the fixed nodes and their locations, any traction boundary conditions, and the density
         * (which multiplies the body force, otherwise isn't used).
         */
        SolidMechanicsProblemDefinition<2> problem_defn(mesh);

        /* Set the material problem on the problem definition object, saying that the problem, and
         * the material law, is incompressible. All material law files can be found in
         * `continuum_mechanics/src/problem/material_laws`. */
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);

        /* Set the fixed nodes, choosing zero displacement for these nodes (see later for how
         * to provide locations for the fixed nodes). */
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        /* Set the body force and the density. (Note that the second line isn't technically
         * needed, as internally the density is initialised to 1)
         */
        problem_defn.SetBodyForce(body_force);
        problem_defn.SetDensity(1.0);

        /* Now we create the (incompressible) solver, passing in the mesh, problem definition
         * and output directory
         */
        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "SimpleIncompressibleElasticityTutorial");

        /* .. and to compute the solution, just call `Solve()` */
        solver.Solve();


        /* '''Visualisation'''. Go to the folder `SimpleIncompressibleElasticityTutorial` in your test-output directory.
         * There should be 2 files, initial.nodes and solution.nodes. These are the original nodal positions and the deformed
         * positions. Each file has two columns, the x and y locations of each node. To visualise the solution in say
         * Matlab or Octave, you could do: `x=load('solution.nodes'); plot(x(:,1),x(:,2),'k*')`. For Cmgui output, see below.
         *
         * To get the actual solution from the solver, use these two methods. Note that the first
         * gets the deformed position (ie the new location, not the displacement). They are both of size
         * num_total_nodes.
         */
        std::vector<c_vector<double,2> >& r_deformed_positions = solver.rGetDeformedPosition();
        std::vector<double>& r_pressures = solver.rGetPressures();
        /* Let us obtain the values of the new position, and the pressure, at the bottom right corner node. */
        unsigned node_index = 8;
        assert( fabs(mesh.GetNode(node_index)->rGetLocation()[0] - 0.8) < 1e-6); // check that X=0.8, ie that we have the correct node,
        assert( fabs(mesh.GetNode(node_index)->rGetLocation()[1] - 0.0) < 1e-6); // check that Y=0.0, ie that we have the correct node,
        std::cout << "New position: " << r_deformed_positions[node_index](0) << " " << r_deformed_positions[node_index](1) << "\n";
        std::cout << "Pressure: " << r_pressures[node_index] << "\n";

        /* HOW_TO_TAG Continuum mechanics
         * Visualise nonlinear elasticity problems solutions, including visualisng strains
         */

        /* One visualiser is Cmgui. This method can be used to convert all the output files to Cmgui format.
         * They are placed in `[OUTPUT_DIRECTORY]/cmgui`. A script is created to easily load the data: in a
         * terminal cd to this directory and call `cmgui LoadSolutions.com`. (In this directory, the initial position is given by
         * solution_0.exnode, the deformed by solution_1.exnode).
         */
        solver.CreateCmguiOutput();

        /* The recommended visualiser is Paraview, for which Chaste must be installed with VTK. With paraview, strains (and in the future
         * stresses) can be visualised on the undeformed/deformed geometry). We can create VTK output using
         * the `VtkNonlinearElasticitySolutionWriter` class. The undeformed geometry, solution displacement, and pressure (if incompressible
         * problem) are written to file, and below we also choose to write the deformation tensor C for each element.
         */
        VtkNonlinearElasticitySolutionWriter<2> vtk_writer(solver);
        vtk_writer.SetWriteElementWiseStrains(DEFORMATION_TENSOR_C); // other options are DEFORMATION_GRADIENT_F and LAGRANGE_STRAIN_E
        vtk_writer.Write();


        /* These are just to check that nothing has been accidentally changed in this test.
         * Newton's method (with damping) was used to solve the nonlinear problem, and we check that
         * 4 iterations were needed to converge.
         */
        TS_ASSERT_DELTA(r_deformed_positions[node_index](0),  0.7980, 1e-3);
        TS_ASSERT_DELTA(r_deformed_positions[node_index](1), -0.1129, 1e-3);
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 4u);
    }
    /* ''Exercise'': convert to a compressible solver and compare the resultant deformations.
     * The next tutorial describes how to solve for a compressible deformation,
     * but the changes are essentially trivial: `IncompressibleNonlinearElasticitySolver` needs to be changed to
     * `CompressibleNonlinearElasticitySolver`, the line `problem_defn.SetMaterialLaw(..)` needs changing, and
     * the material law itself should be of type `AbstractCompressibleMaterialLaw`. An example is
     * `CompressibleMooneyRivlinMaterialLaw`. Also `solver.rGetPressures()` doesn't exist (or make sense)
     * when the solver is an `CompressibleNonlinearElasticitySolver`.
     */

    /* == Incompressible deformation: 2D shape hanging under gravity with a balancing traction ==
     *
     * We now repeat the above test but include a traction on the bottom surface (Y=0). We apply this
     * in the inward direction so that is counters (somewhat) the effect of gravity. We also show how stresses
     * and strains can be written to file.
     */
    void TestIncompressibleProblemWithTractions()
    {
        /* All of this is exactly as above */
        QuadraticMesh<2> mesh;
        mesh.ConstructRegularSlabMesh(0.1 /*stepsize*/, 0.8 /*width*/, 1.0 /*height*/);

        MooneyRivlinMaterialLaw<2> law(1.0);

        c_vector<double,2> body_force;
        body_force(0) =  0.0;
        body_force(1) = -2.0;

        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mesh, 1, 1.0);

        /* Now the traction boundary conditions. We need to collect all the boundary elements on the surface which we want to
         * apply non-zero tractions, put them in a `std::vector`, and create a corresponding `std::vector` of the tractions
         * for each of the boundary elements. Note that the each traction is a 2D vector with dimensions of pressure.
         *
         * First, declare the data structures:
         */
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        /* Create a constant traction */
        c_vector<double,2> traction;
        traction(0) = 0;
        traction(1) = 1.0; // this choice of sign corresponds to an inward force (if applied to the bottom surface)
        /* Loop over boundary elements */
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            /* If the centre of the element has Y value of 0.0, it is on the surface we need */
            if (fabs((*iter)->CalculateCentroid()[1] - 0.0) < 1e-6)
            {
                /* Put the boundary element and the constant traction into the stores. */
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        /* A quick check */
        assert(boundary_elems.size() == 8u);

        /* Now create the problem definition object, setting the material law, fixed nodes and body force as
         * before (this time not calling `SetDensity()`, so using the default density of 1.0,
         * and also calling a method for setting tractions, which takes in the boundary elements
         * and tractions for each of those elements.
         */
        SolidMechanicsProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetBodyForce(body_force);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);

        /* Create solver as before */
        IncompressibleNonlinearElasticitySolver<2> solver(mesh,
                                                          problem_defn,
                                                          "IncompressibleElasticityWithTractionsTutorial");

        /* In this test we also output the stress and strain. For the former, we have to tell the solver to store
         * the stresses that are computed during the solve.
         */
        solver.SetComputeAverageStressPerElementDuringSolve();


        /* Call `Solve()` */
        solver.Solve();

        /* If VTK output is written (discussed above) strains can be visualised. Alternatively, we can create text files
         * for strains and stresses by doing the following.
         *
         * Write the final deformation gradients to file. The i-th line of this file provides the deformation gradient F,
         * written as 'F(0,0) F(0,1) F(1,0) F(1,1)', evaluated at the centroid of the i-th element. The first variable
         * can also be DEFORMATION_TENSOR_C or LAGRANGE_STRAIN_E to write C or E. The second parameter is the file name.
         */
        solver.WriteCurrentStrains(DEFORMATION_GRADIENT_F,"deformation_grad");

        /* Since we called `SetComputeAverageStressPerElementDuringSolve`, we can write the stresses to file too. However,
         * note that for each element this is not the stress evaluated at the centroid, but the mean average of the stresses
         * evaluated at the quadrature points - for technical cardiac electromechanics reasons, it is difficult to
         * define the stress at non-quadrature points.
         */
        solver.WriteCurrentAverageElementStresses("2nd_PK_stress");


        /* Another quick check */
        TS_ASSERT_EQUALS(solver.GetNumNewtonIterations(), 4u);

        /* Visualise as before by going to the output directory and doing
         * `x=load('solution.nodes'); plot(x(:,1),x(:,2),'m*')` in Matlab/octave, or by using Cmgui.
         * The effect of the traction should be clear (especially when compared to
         * the results of the first test).
         *
         * Create Cmgui output
         */
        solver.CreateCmguiOutput();

        /* This is just to check that nothing has been accidentally changed in this test */
        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[8](0), 0.8561, 1e-3);
        TS_ASSERT_DELTA(solver.rGetDeformedPosition()[8](1), 0.0310, 1e-3);
    }
};
    /* More examples are given in the next tutorial
     *
     * == IMPORTANT: Using HYPRE ==
     *
     * Mechanics solves being nonlinear are expensive, so it is recommended you also use `CMAKE_BUILD_TYPE=Release` when running with CMake
     * on larger problems.
     *
     * When running '''incompressible''' problems in 3D, or with more elements, it is vital to also change the linear solver to use HYPRE, an algebraic multigrid
     * solver. Without HYRPE, the linear solve (i) may become very very slow; or (ii) may not converge, in which case the nonlinear
     * solve will (probably) not converge. HYPRE is (currently) not a pre-requisite for installing Chaste, hence this is not (currently)
     * the default linear solver for incompressible mechanics problems, although this will change in the future.
     *
     * ''HYPRE should be considered a pre-requisite for large incompressible mechanics problems.''
     *
     * To use HYPRE, you need to have PETSc installed with HYPRE. However, if you followed installation
     * instructions for Chaste 2.1 or later, you probably do already have PETSc installed with HYPRE.
     *
     * To switch on HYPRE, open the file `continuum_mechanics/src/solver/AbstractNonlinearElasticitySolver` and uncomment the line
     * #define MECH_USE_HYPRE
     * near the top of the file (currently: line 57).
     *
     * Note: PETSc unfortunately doesn't quit if you try to use HYPRE without it being installed, but it spew lots of error messages.
     *
     */

#endif /*TESTSOLVINGELASTICITYPROBLEMSTUTORIAL_HPP_*/
