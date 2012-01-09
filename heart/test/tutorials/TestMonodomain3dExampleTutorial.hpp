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
#ifndef TESTMONODOMAIN3DEXAMPLETUTORIAL_HPP_
#define TESTMONODOMAIN3DEXAMPLETUTORIAL_HPP_

/*
 * = 3D monodomain example =
 *
 * In this tutorial we show how to run a 3D simulation using the monodomain equation.
 * To go from monodomain to bidomain or vice versa is trivial, and for 2d to 3d is
 * also very easy.
 *
 * EMPTYLINE
 *
 * First include the headers, `MonodomainProblem` this time.
 */
#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "TetrahedralMesh.hpp"

/* Here we define a cell factory that gives stimuli to cells in the block
 * 0<x<0.1, 0<y<0.1, 0<z<0.1. Note that it inherits from `AbstractCardiacCellFactory<3>`
 * this time (not `<2>`).
 */
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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];
        double z = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[2];

        if ( (x<0.1+1e-6) && (y<0.1+1e-6) && (z<0.1+1e-6) )
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }
    }
};

/* Now define the test */
class TestMonodomain3dExample : public CxxTest::TestSuite
{
public:
    void TestMonodomain3d() throw(Exception)
    {
        /* We will auto-generate a mesh this time, and pass it in, rather than
         * provide a mesh file name. This is how to generate a cuboid mesh with
         * a given spatial stepsize h */
        TetrahedralMesh<3,3> mesh;
        double h=0.02;
        mesh.ConstructRegularSlabMesh(h, 0.8 /*length*/, 0.3 /*width*/, 0.3 /*depth*/);
        /* (In 2D the call is identical, but without the depth parameter).
         *
         * EMPTYLINE
         *
         * Set the simulation duration, etc, and create an instance of the cell factory.
         * One thing that should be noted for monodomain problems, the ''intracellular
         * conductivity'' is used as the monodomain effective conductivity (not a
         * harmonic mean of intra and extracellular conductivities). So if you want to
         * alter the monodomain conductivity call
         * `HeartConfig::Instance()->SetIntracellularConductivities`
         */
        HeartConfig::Instance()->SetSimulationDuration(5); //ms
        HeartConfig::Instance()->SetOutputDirectory("Monodomain3dExample");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, 0.01, 0.1);

        BenchmarkCellFactory cell_factory;

        /* Now we declare the problem class, `MonodomainProblem<3>` instead of `BidomainProblem<2>`.
         * The interface for both is the same.
         */
        MonodomainProblem<3> monodomain_problem( &cell_factory );

        /* If a mesh-file-name hasn't been set using `HeartConfig`, we have to pass in
         * a mesh using the `SetMesh` method (must be called before `Initialise`). */
        monodomain_problem.SetMesh(&mesh);

        /* By default data for all nodes is output, but for big simulations, sometimes this
         * might not be required, and the action potential only at certain nodes required.
         * The following code shows how to output the results at the first, middle and last
         * nodes, for example. (The output is written to the HDF5 file; no meshalyzer output
         * will be made. HDF5 files can be read using Matlab). We are not using this in this
         * simulation however (hence the boolean being set to false).
         */
        bool partial_output = false;
        if(partial_output)
        {
            std::vector<unsigned> nodes_to_be_output;
            nodes_to_be_output.push_back(0);
            nodes_to_be_output.push_back((unsigned)round( (mesh.GetNumNodes()-1)/2 ));
            nodes_to_be_output.push_back(mesh.GetNumNodes()-1);
            monodomain_problem.SetOutputNodes(nodes_to_be_output);
        }

        /* `SetWriteInfo` is a useful method that means that the min/max voltage is
         * printed as the simulation runs (useful for verifying that cells are stimulated
         * and the wave propagating, for example) (although note scons does buffer output
         * before printing to screen) */
        monodomain_problem.SetWriteInfo();

        /* Finally, call `Initialise` and `Solve` as before */
        monodomain_problem.Initialise();
        monodomain_problem.Solve();

        /* This part is just to check nothing has accidentally been changed in this example */
        ReplicatableVector voltage(monodomain_problem.GetSolution());
        TS_ASSERT_DELTA(voltage[0], 35.0939, 1e-2);
    }
};
#endif /*TESTMONODOMAIN3DEXAMPLETUTORIAL_HPP_*/
