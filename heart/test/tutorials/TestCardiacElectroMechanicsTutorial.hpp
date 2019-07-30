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
#ifndef TESTCARDIACELECTROMECHANICSTUTORIAL_HPP_
#define TESTCARDIACELECTROMECHANICSTUTORIAL_HPP_

/*
 * = Cardiac Electro-mechanical Problems =
 *
 * == Introduction ==
 *
 * The tutorial explains how electro-mechanics problems can be solved in Chaste. The reader should certainly read
 * the electro-physiological tutorials before this tutorial, and really they should have also had a look at
 * the tutorial(s) on solving general solid mechanics problems.
 *
 * The equations of cardiac electro-mechanics are written down in Section 4.2 of the PDF on equations and
 * finite element implementations in ChasteGuides -> Miscellaneous information. '''Note:''' By default we do
 * not solve these full equations: the mechanics information is not coupled back to electrics, ie by default
 * the conductivities do not depend on deformation, and cell models do not get affected by stretch.
 * This has to be switched on if required, as will be described further below.
 *
 * Before going to the code, we list the sub-models/parameters that need to be set, or can be varied,
 * in electro-mechanical problems. The last five of the below are mechanics-specific.
 *  * The geometry (see note 1 below)
 *  * The region electrically stimulated
 *  * The cell model
 *  * Electro-physiological parameters (conductivity, capacitance, surface-area-to-volume ratio)
 *  * Electro-physiological timesteps: ode and pde (but not printing timestep) (see note 2 below)
 *  * Fibre directions (and maybe sheet/normal directions) (see note 3 below)
 *  * The part of the boundary that has displacement boundary conditions
 *  * Any pressure or traction boundary conditions
 *  * The contraction model [the model which takes in electrical variables (voltage or calcium typically), and
 *  returns cellular active tension]
 *  * Whether the tissue should be treated as compressible or incompressible. (Although likely technically
 *  incompressible at appropriate scales, cardiac tissue is often treated as compressible due to blood
 *  squeezed out of the coronary vessels during contraction).
 *  * The material law [the strain-energy function]
 *  * Mechanics timesteps: mechanics update timestep, contraction model ode timestep. (see note 4 below)
 *
 * Notes:
 *  * ''Meshes:'' Two meshes for the geometry are required, one for the electrics solve and one for the mechanics.
 * The mechanics mesh would ideally be coarser but any two meshes are technically possible. The meshes should
 * ideally both cover exactly the same geometry (ie either mesh being contained in the other), but the meshes
 * not completely overlapping is allowed - some extrapolation of quantities will then occur.
 *  * ''The electro-physiology printing timestep:'' This is not used in electro-mechanics problems; output is
 * instead written after every mechanics solve, so effectively the mechanics update timestep is equal to
 * the printing timestep.
 *  * ''Fibres:'' In electro-physiological simulations the fibre direction is in the X-direction
 * by default, but if isotropic conductivities are used the fibre direction won't be used. In mechanics
 * solves, the fibres will always be used as it determines the direction of contraction. It defaults to the
 * X-direction, so this is the direction the tissue will contract, unless a fibre file is given.
 * If the material law is transversely isotropic, the problem is independent of sheet & normal directions.
 * If the material law is anisotropic, the problem is dependent of sheet & normal directions.
 *  * ''Timesteps:'' Should-divide rules are: (a) ode_timestep should-divide pde_timestep should-divide
 *  mechanics_update_timestep and (b) contraction_model_ode_timestep should-divide mechanics_update_timestep.
 *
 * '''Another important note:''' mechanics problems are not currently implemented to scale in parallel yet. This
 * is work in progress.
 *
 * The basic includes are */
#include <cxxtest/TestSuite.h>
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "LuoRudy1991.hpp"
/* Some other includes that are used */
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "FileComparison.hpp"
#include "FileFinder.hpp"

/*
 * HOW_TO_TAG Cardiac/Electro-mechanics
 * Run basic electro-mechanics simulations; specify different models, boundary conditions, fibres
 */

/*
 * == IMPORTANT: using HYPRE ==
 *
 * Mechanics solves being nonlinear are expensive, so it is recommended you also use a `Release` build type for `cmake`
 * (or `build=GccOpt_ndebug` when running on the old `scons` build system)
 * on larger problems. Also:
 *
 * Mechanics solves involve solving a nonlinear system, which is broken down into a sequence of linear solves.
 * When running '''incompressible''' problems '''in 3D, or with more elements than in the first test below''',
 * it is vital to change the linear solver to use HYPRE, an algebraic multigrid solver.
 * Without HYRPE, the linear solve (i) may become very very slow; or
 * (ii) may not converge, in which case the nonlinear solve will (probably) not converge. See the comments on using
 * HYPRE in the first solid mechanics tutorial.
 *
 * == Simple 2d test ==
 *
 * This test shows how to use the `CardiacElectroMechProbRegularGeom` class, which
 * inherits from a more general class, `CardiacElectroMechanicsProblem`, and
 * sets up a square or cubic geometry for you. Using
 * `CardiacElectroMechProbRegularGeom` is not really recommended, as the functionality
 * it allows is very limited - it is better to use `CardiacElectroMechanicsProblem`, which
 * is shown in the following tests. We use `CardiacElectroMechProbRegularGeom`
 * in this first tutorial just to illustrate a simulation with a few lines (four!) of code.
 */
class TestCardiacElectroMechanicsTutorial : public CxxTest::TestSuite
{
public:
    void TestCardiacElectroMechanicsExample()
    {
        /* All electro-mechanics problems require a cell factory as normal. This particular
         * factory stimulates the LHS side (X=0) surface. */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000);

        /* Electro-physiology parameters, such as the cell-model ODE timestep, the monodomain PDE timestep,
         * the conductivities, capacitance etc, are set using `HeartConfig` as in electro-physiological
         * (ie not mechanical) simulations. We use the defaults for all of these. The one variable that
         * has to be set on `HeartConfig` is the end time of the simulation.
         */
        HeartConfig::Instance()->SetSimulationDuration(40.0);

        /* The main solver class for electro-mechanics, equivalent to `MonodomainProblem` or `BidomainProblem`,
         * is `CardiacElectroMechanicsProblem`. We will show how to use this class in later tests. The
         * subclass of `CardiacElectroMechanicsProblem` called `CardiacElectroMechProbRegularGeom`
         * can be used to quickly set up simulations on a square geometry. It is only present for
         * convenience, and doesn't allow for much flexibility or configurability. The constructor of this
         * class takes in whether an incompressible or compressible problem should be solved,
         * information about the geometry to be created, and the some information about the
         * mechanics: which contraction model to use, what ODE timestep to use with it, and how often
         * to solve the mechanics. In this class the equation that describes the electrics is assumed to be
         * the monodomain one.
         */
        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     0.1,  // width of square (cm)
                                                     5,    // Number mechanics elements in each direction
                                                     10,   // Number electrics elements in each direction
                                                     &cell_factory,
                                                     KERCHOFFS2003,  // The contraction model (see below)
                                                     1.0,  // mechanics solve timestep
                                                     0.01, // contraction model ode timestep
                                                     "TestCardiacElectroMechanicsExample" /* output directory */);
        /* The contraction model chosen above is 'KERCHOFFS2003' (Kerchoffs, Journal of Engineering Mathematics, 2003). Other possibilities
         * are 'NHS' (Niederer, Hunter, Smith, 2006), and 'NASH2004' (Nash, Progress in Biophysics and Molecular Biology, 2004).
         *
         * Two meshes are created, one with five elements in each direction for the mechanics (so 5*5*2 triangles in total),
         * and a finer one for the electrics.
         *
         * This leaves the material law, fibres direction and fixed nodes from the list above: the material
         * law is the default incompressible material law (pole-zero), the fibre direction is by default
         * the X-direction, and the fixed nodes are automatically set be those satisfying X=0, ie
         * the left-hand edge. No surface tractions are set. To do something more general, `CardiacElectroMechanicsProblem`
         * must be used.
         *
         * All we now have to do is call Solve.
         */
        problem.Solve();

        /* Go to the output directory. There should be log file (which, note, can be used to watch progress
         * during a simulation), and a directory for the electrics output and the mechanics output. The electrics
         * directory is not the same as when running an electrics solve: the basic HDF5 data is there but
         * there is no Meshalyzer output, and there is always Cmgui output, including the electrics solution downsampled
         * onto the mechanics mesh.
         * The deformation output directory contains the deformed solution each timestep in several simple
         * MATLAB-readable files, and a cmgui output directory. The latter has a script for automatically loading
         * all the results.
         *
         * HOW_TO_TAG Cardiac/Electro-mechanics
         * Visualise results in Cmgui (very brief description)
         *
         * Visualise the results by calling `cmgui LoadSolutions.com` in the directory
         * `TestCardiacElectroMechanicsExample/deformation/cmgui` . The electrics data can be visualised on the
         * deforming mesh by using the Scene and Spectrum Editors.  See also ChasteGuides/UsingCmgui.
         *
         * (See cmgui website for information on how
         * to use cmgui, but very briefly: graphics -> scene editor -> select surfaces -> add, then check 'Data'. Then
         * graphics -> Spectrum editor -> min=-90, max=50.).
         *
         * The undeformed shape is visualied as as the t=-1 configuration. Here this is the same as the t=0 configuration,
         * but in some cases that might not be the case - see eg a later tutorials involving internal pressures and an
         * annulus.
         *
         * To observe the tissue relaxing you can re-run the simulation with an end time of more than 350ms.
         *
         */
    }

    /* == Same simulation, this time using `CardiacElectroMechanicsProblem` ==
     *
     * Let us repeat the above test using `CardiacElectroMechanicsProblem`. */
    void TestCardiacElectroMechanicsExampleAgain()
    {
        /* This lines is as above */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000);

        /* Create two meshes, one for the electrics, one for the mechanics, covering the same
         * region, with different mesh resolutions. The first mesh should be a `TetrahedralMesh`,
         * (as used in monodomain/bidomain), the second should be a `QuadraticMesh` (as used
         * in mechanics problems).
         */
        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01/*stepsize*/, 0.1/*length*/, 0.1/*width*/, 0.1/*depth*/);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1, 0.1 /*as above with a different stepsize*/);

        /* Set the end time as above */
        HeartConfig::Instance()->SetSimulationDuration(40.0);

        /* In the solid mechanics tutorials, you can see how to use the class `SolidMechanicsProblemDefinition`
         * to set up a mechanics problem to be solved. The class allows you to specify things like: material law,
         * fixed nodes, traction boundary conditions, gravity, and so on. For electro-mechanics problems, we use
         * the  class `ElectroMechanicsProblemDefinition`, which inherits from `SolidMechanicsProblemDefinition`
         * (and therefore has the same functionality), as well as a few electro-mechanics specific methods.
         *
         * We choose to fix the nodes on X=0. For this the `NonlinearElasticityTools` class
         * is helpful. The static method called below returns all nodes for which the X value
         * (indicated by the '0' ('0' for X, '1' for Y, '2' for Z)) is equal to 0.0.
         */
        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.0); // all the X=0.0 nodes

        /* Now we create the problem definition class, tell it about the fixed nodes, the contraction model to be used,
         * that we want to use the default cardiac material law, and the mechanics
         * solve timestep (how often the mechanics is solved). An error would occur if we failed to provide
         * information about any of these. Optional other things that could have been set are gravity, tractions,
         * and whether to use mechano-electric feedback. The material law used below is the default incompressible material law, the Pole-Zero law. This is
         * defined in `continuum_mechanics/src/problem/material_laws/NashHunterPoleZeroLaw`. (All material
         * laws are in this folder). Note that the parameters values in this law are such that the
         * material law is transversely isotropic, so sheet and normal directions do not matter.
         */
        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,0.01/*contraction model ODE timestep*/);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        /* Now create the problem class, passing in the compressibility type (COMPRESSIBLE or INCOMPRESSIBLE),
         * the type of electrics propagation equation (MONODOMAIN in this case),the meshes, the cell factory, and the problem_definition class,
         * and call solve. The first template parameter (2) is the dimension of the space, the second one is the number of unknowns
         * in the electrics problem (1 for MONODOMAIN, 2 for BIDOMAIN)*/
        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestCardiacElectroMechanicsExample2");

        problem.Solve();
        /* Visualise as above.
         *
         * Some comments: to use compressibility instead of incompressibility, just change the two
         * 'INCOMPRESSIBLE's to 'COMPRESSIBLE'. (Note that this leads to a completely different type
         * of problem, and a completely different type of solver - the incompressible problem involves
         * solving for displacement and pressure, and mixed formulations and saddle-point problems,
         * the compressible problem does not. Compressible problems  are far less computationally-demanding).
         * A compressible example is given later in this tutorial.
         *
         * The default incompressible material law is the pole-zero law, and the default
         * compressible material law is an exponential law. To pass in your own choice of
         * material law, call `SetMaterialLaw()`, as in a normal solid mechanics simulation. For example:
         */
        CompressibleMooneyRivlinMaterialLaw<2> law(2.0,1.0); // random (non-cardiac) material law
        problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
        /* As mentioned above, by default the deformation does '''not''' couple back to the electrics.
         * The stretch is not passed to the cell model to allow for stretch-activated channels (M.E.F.),
         * and the deformation is not used in altering the conductivity tensor (the latter simplifications has
         * little effect in
         * in simple propagation problems - see "A numerical method for cardiac mechano-electric simulations",
         * Annals of Biomedical Engineering). To set the solver to use either of these, do, for example */
        problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, true /*deformation affects cell models*/);
        /* before calling `problem.Solve()`. Deformation affecting cell models is described in more detail
         * later in this tutorial. For deformation affecting conductivity, note that the electrics solve will
         * slow down, since the linear system matrix now varies with time (as conductivities depend
         * on deformation), and has to be recomputed after every mechanics update. The set-up cost in this
         * case currently requires optimisation, also.
         *
         * Finally, `SetNoElectricsOutput` is a method that is sometimes useful with a fine electrics mesh. */
        problem.SetNoElectricsOutput();

        /* The final position of the nodes can be obtained as follows (same interface in described in the solid mechanics tutorials). */
        TS_ASSERT_DELTA(problem.rGetDeformedPosition()[5](0), 0.090464, 1e-4);

        /* Ignore these tests, they are they to check nothing has changed in this tutorial */
        FileFinder finder1("TestCardiacElectroMechanicsExample/deformation/solution_40.nodes", RelativeTo::ChasteTestOutput);
        FileFinder finder2("TestCardiacElectroMechanicsExample2/deformation/solution_40.nodes", RelativeTo::ChasteTestOutput);
        FileComparison comparer(finder1,finder2);
        TS_ASSERT(comparer.CompareFiles());
    }

    /* == Twisting cube: 3d example with varying fibre directions ==
     *
     * The third test is a longer running 3d test - the 'dont' in the name of the test
     * means it isn't run automatically. To run, remove the 'dont'. It is worth running
     * with `build=GccOpt_ndebug`; and see the comments about HYPRE above if you change
     * this to an incompressible solve.
     *
     * This test shows how to do 3d simulations (trivial changes), and how to pass in
     * fibre directions for the mechanics mesh. It also uses a compressible law.
     */
    void dontTestTwistingCube()
    {
        /* Cell factory as normal */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 3> cell_factory(-1000*1000);

        /* Set up two meshes of 1mm by 1mm by 1mm, one a `TetrahedralMesh`
         * for the electrics solve, one a (coarser) `QuadraticMesh` for the mechanics
         * solve. */
        TetrahedralMesh<3,3> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.01/*stepsize*/, 0.1/*length*/, 0.1/*width*/, 0.1/*depth*/);

        QuadraticMesh<3> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1, 0.1 /*as above with a different stepsize*/);

        /* Collect the nodes on Z=0 */
        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<3>::GetNodesByComponentValue(mechanics_mesh, 2, 0.0);

        /* Set the simulation end time as before */
        HeartConfig::Instance()->SetSimulationDuration(50.0);

        /* Create the problem definition object as before (except now the template parameter is 3). */
        ElectroMechanicsProblemDefinition<3> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(COMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        /* The default fibre direction is the X-direction (and the default sheet plane is the XY plane). Now we show
         * how this can be changed.
         *
         * Fibre files should be .ortho files, not .axi file (see file formats documentation if you haven't come across these files,
         * basically .axi files specify the fibre directions; .ortho the fibre sheet and normal directions).
         * For mechanics problems, the .ortho file
         * can be used to either define the fibre information PER-ELEMENT or PER-QUADRATURE-POINT (ie all the quadrature points
         * in all the elements). The latter provides a higher resolution description of fibres.
         *
         * In this tutorial, we will generate both types of fibre files, using our own choice of fibre fibre for the cubic tissue.
         * To generate a fibre file prescribing fibres that depend on the X-coordinate, one fibre definition per element,
         * we can do:
         */
        OutputFileHandler handler("TutorialFibreFiles");
        out_stream p_file = handler.OpenOutputFile("5by5by5_fibres.ortho");

        *p_file << mechanics_mesh.GetNumElements() << "\n"; // first line is number of entries
        for (unsigned i=0; i<mechanics_mesh.GetNumElements(); i++)
        {
            double X = mechanics_mesh.GetElement(i)->CalculateCentroid()(0);
            double theta = M_PI/3 - 10*X*2*M_PI/3; // 60 degrees when X=0, -60 when X=0.1;
            *p_file <<  "0 " << cos(theta)  << " " << sin(theta)  // first three entries are fibre direction
                    << " 0 " << -sin(theta) << " " << cos(theta)  // next three are sheet direction
                    << " 1 0 0\n";                                // then normal to sheet direction
        }
        p_file->close();
        /* This will generate a file, TutorialFibreFiles/5by5by5_fibres.ortho. Note that out_streams are essentially
         * pointers to a C++ ofstream.
         *
         * More advanced: we can also generate the same type of file, but where there is one line for each quadrature point.
         * By default there are, per element, 3 quadrature points in each direction, so in this 3D problem there are
         * (3^3^)*num_elem quadrature points. Here's how we can obtain their positions, and set-up the analogous
         * fibre file, which we name similarly to the above but change the extension. */
        out_stream p_file2 = handler.OpenOutputFile("5by5by5_fibres.orthoquad");

        // Mechanics deformation solvers use 3rd order quadrature rules
        GaussianQuadratureRule<3> quad_rule(3);
        QuadraturePointsGroup<3> quad_points(mechanics_mesh, quad_rule);

        *p_file2 << quad_points.Size() << "\n";
        for (unsigned i=0; i<quad_points.Size(); i++)
        {
            double X = quad_points.rGet(i)(0);
            double theta = M_PI/3 - 10*X*2*M_PI/3;
            *p_file2 <<  "0 " << cos(theta)  << " " << sin(theta)
                     << " 0 " << -sin(theta) << " " << cos(theta)
                     << " 1 0 0\n";
        }
        p_file2->close();

        /*
         * We use the `FileFinder` class to identify locations of files.
         * `OutputFileHandler` has a handy method called `FindFile()` which returns a `FileFinder` to a file in the folder it points to.
         */
        FileFinder finder = handler.FindFile("5by5by5_fibres.orthoquad");
        problem_defn.SetVariableFibreSheetDirectionsFile(finder, true);

        /* Create the problem object */
        CardiacElectroMechanicsProblem<3,1> problem(COMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestCardiacElectroMech3dTwistingCube");

        /* Now call `Solve`. This will take a while to run, so watch progress using the log file to estimate when
         * it will finish. `build=GccOpt_ndebug` will speed this up by a factor of about 5. Visualise in Cmgui as usual.
         */
        problem.Solve();
    }
};

#endif /*TESTCARDIACELECTROMECHANICSTUTORIAL_HPP_*/

