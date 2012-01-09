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

/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */
#ifndef TESTBIDOMAINWITHBATHANDFIBRESTUTORIAL_HPP_
#define TESTBIDOMAINWITHBATHANDFIBRESTUTORIAL_HPP_
/*
 * = Running a bidomain simulation with a bath and fibres =
 *
 * In this tutorial we run a bidomain simulation with both a bath and fibres
 *
 * EMPTYLINE
 *
 * We include the same headers as in the previous fibre tutorial
 */
#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"

/* Define the test class as before
 */
class TestBidomainWithBathAndFibresTutorial : public CxxTest::TestSuite
{
public:
    void TestSimulation() throw(Exception)
    {
        HeartConfig::Instance()->SetSimulationDuration(5.0);  //ms
        HeartConfig::Instance()->SetOutputDirectory("BidomainTutorialWithBathAndFibres");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        /* Bath problems seem to require decreased ODE timesteps. We use the
         * Backward Euler version of the Luo-Rudy model (see below) instead to
         * improve code performance.
         */
        HeartConfig::Instance()->SetOdeTimeStep(0.01);  //ms

        /* Use the {{{PlaneStimulusCellFactory}}} to define a set of Luo-Rudy cells, in this
         * case with a Backward Euler solver. We pass the stimulus magnitude as 0.0
         * as we don't want any stimulated cells.
         */
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellMLBackwardEuler,2> cell_factory(0.0);

        /*
         * Note that in the previous bath example, a mesh was read in and elements where then set to be
         * bath elements in the test. With fibres as well, in a bath simulation, it is better to read in a
         * mesh that has all the information: this mesh has bath elements defined as an extra column in the
         * .ele file, and a .ortho file which defines the fibre direction for each element. Note that the
         * .ortho file should include fibre information for bath elements as well, but they won't be used
         * in the simulation. (The fibres read here are the same 'kinked' fibres as in the previous fibre
         * tutorial).
         */
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements_bath_sides", cp::media_type::Orthotropic);

        /* Set anistropic conductivities.
         */
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.175));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(7.0, 0.7));

        /* and now we define the electrodes.. */
        double magnitude = -9.0e3; // uA/cm^2
        double start_time = 0.0;
        double duration = 2; //ms
        HeartConfig::Instance()->SetElectrodeParameters(false, 0, magnitude, start_time, duration);

        /* Now create the problem class, using the cell factory and passing
         * in `true` as the second argument to indicate we are solving a bath
         * problem, and solve.
         */
        BidomainProblem<2> bidomain_problem( &cell_factory, true );
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
};

#endif /*TESTBIDOMAINWITHBATHANDFIBRESTUTORIAL_HPP_*/
