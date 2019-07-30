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

#ifndef TESTCARDIACELECTROMECHANICSFURTHERFUNCTIONALITY_HPP_
#define TESTCARDIACELECTROMECHANICSFURTHERFUNCTIONALITY_HPP_

#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include <petscvec.h>
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "LuoRudy1991.hpp"
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "Hdf5DataReader.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"

class TestCardiacElectroMechanicsFurtherFunctionality : public CxxTest::TestSuite
{
public:
    void TestDeterminingWatchedNodes()
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

        HeartConfig::Instance()->SetSimulationDuration(1.0);

        CardiacElectroMechProbRegularGeom<2> problem(INCOMPRESSIBLE,
                                                     1.0, /* width (cm) */
                                                     1,   /* mech elem each dir */
                                                     96,  /* elec elem each dir */
                                                     &cell_factory,
                                                     NHS,
                                                     1.0,  /* mechanics solve timestep */
                                                     0.01, /* contraction model ode timestep */
                                                     "");

        c_vector<double,2> pos;
        pos(0) = 1.0;
        pos(1) = 0.0;
        problem.SetWatchedPosition(pos);
        problem.SetNoElectricsOutput();
        problem.Initialise();

        // have checked these hardcoded values correspond to the nodes
        // at (1,0);
        TS_ASSERT_EQUALS(problem.mWatchedElectricsNodeIndex, 96u);
        TS_ASSERT_EQUALS(problem.mWatchedMechanicsNodeIndex, 1u);

        //// would like to do the following....
        //CardiacElectroMechanicsProblem<2> problem2(&cell_factory,
        //                                           1, 10, 100, 0.01,
        //                                           "");
        //pos(1) = 1.1;
        //problem2.SetWatchedPosition(pos);
        //TS_ASSERT_THROWS_THIS(problem2.Initialise(), "");
        //// ... but the exception causes a segmentation fault and had to be replaced
        //// with an assert(0);
    }

    // Sets up a short simulation on a square with zero stimulus, but a model with stretch activated channels.
    // Hacks the mechanics initial condition to correspond to some stretch, which should create a bit of
    // SAC activity and increased voltage
    void TestWithMechanoElectricFeedback()
    {
        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory(0.0);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1);

        QuadraticMesh<2> mechanics_mesh(0.1, 0.1, 0.1);

        HeartConfig::Instance()->SetSimulationDuration(1.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(NASH2004,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        std::vector<unsigned> fixed_nodes;
        std::vector<c_vector<double,2> > locations;
        for ( TetrahedralMesh<2,2>::BoundaryNodeIterator iter = mechanics_mesh.GetBoundaryNodeIteratorBegin();
              iter != mechanics_mesh.GetBoundaryNodeIteratorEnd();
              ++iter)
        {
            double X = (*iter)->rGetLocation()[0];
            double Y = (*iter)->rGetLocation()[1];

            fixed_nodes.push_back((*iter)->GetIndex());

            c_vector<double,2> new_position;
            new_position(0) = 1.2*X;
            new_position(1) = Y/1.2;
            locations.push_back(new_position);
        }
        problem_defn.SetFixedNodes(fixed_nodes, locations);

        // Set deformation not affecting conductivity, but affecting cell models.
        problem_defn.SetDeformationAffectsElectrophysiology(false,true);

        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestNobleSacActivatedByStretchTissue");

        problem.Initialise();

        // hack into the mechanics solver and set up the current solution so that it corresponds to
        // the square of tissue being stretched - see comment below
        for (unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
            // the 3 here is the problem dim (DIM+1 as incompressible)
            problem.mpMechanicsSolver->rGetCurrentSolution()[3*i]   = X*0.2;
            problem.mpMechanicsSolver->rGetCurrentSolution()[3*i+1] = Y*(1.0/1.2 - 1);
        }

        // we are going to get the modified conductivity tensor directly, without (initially) calling solve,
        // so need to do the following, which is normally done inside the Solve
        problem.mpCardiacMechSolver->ComputeDeformationGradientAndStretchInEachElement(problem.mDeformationGradientsForEachMechanicsElement, problem.mStretchesForEachMechanicsElement);

        // just get the default conductivity so don't to hardcode it (1.75 at the moment)
        c_vector<double, 2> conductivities;
        HeartConfig::Instance()->GetIntracellularConductivities(conductivities);
        assert((conductivities(0)-conductivities(1))<1e-8);
        double default_conductivity = conductivities(0);

        // test directly that the conductivity hasn't been modified
        for (unsigned i=0; i<electrics_mesh.GetNumElements(); i++)
        {
            const c_matrix<double,2,2>& r_tensor = problem.mpElectricsProblem->GetTissue()->rGetIntracellularConductivityTensor(i);
            TS_ASSERT_DELTA(r_tensor(0,0), default_conductivity, 1e-9);
            TS_ASSERT_DELTA(r_tensor(0,1), 0.0,                  1e-9);
            TS_ASSERT_DELTA(r_tensor(1,0), 0.0,                  1e-9);
            TS_ASSERT_DELTA(r_tensor(1,1), default_conductivity, 1e-9);
        }


        // Now we call Solve. The solver will, before the first timestep, first the equilibrium solution
        // given BCs and loading, which obviously won't require any newton steps as the current solution
        // was manually set to the right thing above. Then it will solve one timestep.
        problem.Solve();

        // Get the voltage at the start and end of the simulation, check the stretch was passed down to the
        // cell model and caused increased voltage

        Hdf5DataReader reader("TestNobleSacActivatedByStretchTissue/electrics", "voltage");
        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
        TS_ASSERT_EQUALS(num_timesteps, 2u);

        Vec start_voltage = PetscTools::CreateVec(36);
        Vec end_voltage = PetscTools::CreateVec(36);
        reader.GetVariableOverNodes(start_voltage, "V", 0);
        reader.GetVariableOverNodes(end_voltage, "V", 1);
        ReplicatableVector start_voltage_repl(start_voltage);
        ReplicatableVector end_voltage_repl(end_voltage);

        for (unsigned i=0; i<start_voltage_repl.GetSize(); i++)
        {
            TS_ASSERT_LESS_THAN(start_voltage_repl[i], -90.0);
            TS_ASSERT_LESS_THAN(-90, end_voltage_repl[i]);
        }

        PetscTools::Destroy(start_voltage);
        PetscTools::Destroy(end_voltage);
    }

    // Similar to first part of above test, except the deformation is not homogeneous, and
    // also here we say that the deformation DOES affect conductivity.
    void TestWithMefAndAlteredConductivitesHeterogeneousStretch()
    {
        // irrelevant, not going to call solve
        PlaneStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory(0.0);

        // set up two meshes of 1mm by 1mm by 1mm
        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.025, 0.1, 0.1); // the h should be such that there are an even number of elements per dim, so
                                                                  // electrical element centroids don't lie on mechanics element boundaries
                                                                  // (for the below test to not have to worry about boundary cases).

        QuadraticMesh<2> mechanics_mesh(0.1, 0.1, 0.1); // 2 elements

        // irrelevant, not going to call solve
        std::vector<unsigned> fixed_nodes
          = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh,0,0);

        HeartConfig::Instance()->SetSimulationDuration(1.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(NASH2004,1.0);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);
        problem_defn.SetDeformationAffectsElectrophysiology(true,true);


        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "");


        problem.Initialise();

        // hack into the mechanics solver and set up the current solution so that it corresponds to
        // the some stretch in the upper element
        for (unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
        {
            double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
            double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
            // the 3 here is problem dim (DIM+1 as incompressible)
            problem.mpMechanicsSolver->rGetCurrentSolution()[3*i]   = X*Y*0.1; // displacement is zero except for (1,1) node
            problem.mpMechanicsSolver->rGetCurrentSolution()[3*i+1] = X*Y*0.1; // displacement is zero except for (1,1) node
        }

        // we are going to get the modified conductivity tensor directly, without (initially) calling solve,
        // so need to do the following, which is normally done inside the Solve
        problem.mpCardiacMechSolver->ComputeDeformationGradientAndStretchInEachElement(problem.mDeformationGradientsForEachMechanicsElement, problem.mStretchesForEachMechanicsElement);

        c_matrix<double,2,2> F;
        F(0,0) = 1.01;
        F(1,0) = 0.01;
        F(0,1) = 0.01;
        F(1,1) = 1.01;
        c_matrix<double,2,2> inverse_C = prod(Inverse(F), trans(Inverse(F)));

        // just get the default conductivity so don't to hardcode it (1.75 at the moment)
        c_vector<double, 2> conductivities;
        HeartConfig::Instance()->GetIntracellularConductivities(conductivities);
        assert((conductivities(0)-conductivities(1))<1e-8);
        double default_conductivity = conductivities(0);

        // test directly that the conductivity is being computed using the deformation
        for (unsigned i=0; i<electrics_mesh.GetNumElements(); i++)
        {
            // sigma = F^{-1} sigma_undef F^{-T},
            const c_matrix<double,2,2>& r_tensor = problem.mpElectricsProblem->GetTissue()->rGetIntracellularConductivityTensor(i);
            c_vector<double,2> centroid = electrics_mesh.GetElement(i)->CalculateCentroid();
            if (centroid(0)+centroid(1)<0.1)
            {
                // in mechanics element with corners (0,0),(0,1),(1,0) -- F=I here
                TS_ASSERT_DELTA(r_tensor(0,0), default_conductivity, 1e-9);
                TS_ASSERT_DELTA(r_tensor(0,1), 0.0,                  1e-9);
                TS_ASSERT_DELTA(r_tensor(1,0), 0.0,                  1e-9);
                TS_ASSERT_DELTA(r_tensor(1,1), default_conductivity, 1e-9);
            }
            else
            {
                // in mechanics element with corners (0,1),(1,0),(1,1) -- F \ne I here.
                // sigma_def = F^{-1} sigma_undef F^{-T}, but this is 1.75 C^{-1} since sigma_undef = 1.75*I
                TS_ASSERT_DELTA(r_tensor(0,0), default_conductivity*inverse_C(0,0), 1e-9);
                TS_ASSERT_DELTA(r_tensor(0,1), default_conductivity*inverse_C(0,1), 1e-9);
                TS_ASSERT_DELTA(r_tensor(1,0), default_conductivity*inverse_C(1,0), 1e-9);
                TS_ASSERT_DELTA(r_tensor(1,1), default_conductivity*inverse_C(1,1), 1e-9);
            }
        }
    }

    // Run a where the domain in long and thin, and held squashed in the X-direction, by Dirichlet
    // boundary conditions on every node. Run with and without deformation affecting the
    // conductivity - in the latter as the conductivity will be increased in the X-direction, the
    // wave should travel a little bit faster.
    void TestDeformationAffectingConductivity()
    {
        unsigned num_stimulated_nodes[2];

        for (unsigned sim=0; sim<2; sim++)
        {
            PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-1000*1000);

            double tissue_init_width = 0.5;

            TetrahedralMesh<2,2> electrics_mesh;
            electrics_mesh.ConstructRegularSlabMesh(0.0125, tissue_init_width, 0.025);

            QuadraticMesh<2> mechanics_mesh;
            mechanics_mesh.ConstructRegularSlabMesh(0.025, tissue_init_width, 0.025);

            // Fix every single node, to the deformation x=alpha X, y = Y/alpha
            std::vector<unsigned> fixed_nodes;
            std::vector<c_vector<double,2> > fixed_node_locations;

            for (unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
            {
                double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
                double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];

                c_vector<double,2> new_position;
                new_position(0) = 0.8*X;
                new_position(1) = Y/0.8;
                fixed_nodes.push_back(i);
                fixed_node_locations.push_back(new_position);
            }

            ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
            problem_defn.SetContractionModel(KERCHOFFS2003,1.0);
            problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
            problem_defn.SetFixedNodes(fixed_nodes, fixed_node_locations);
            problem_defn.SetMechanicsSolveTimestep(0.5);

            std::string dir;
            if (sim==0)
            {
                problem_defn.SetDeformationAffectsElectrophysiology(false,false);
                dir = "TestCardiacEmDeformationNotAffectingConductivity";
            }
            else
            {
                problem_defn.SetDeformationAffectsElectrophysiology(true,false);
                dir = "TestCardiacEmDeformationAffectingConductivity";
            }

            // Small enough end-time so that the wavefront doesn't reach the other side..
            HeartConfig::Instance()->SetSimulationDuration(5.0);

            CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                      MONODOMAIN,
                                                      &electrics_mesh,
                                                      &mechanics_mesh,
                                                      &cell_factory,
                                                      &problem_defn,
                                                      dir);


            problem.Initialise();

            // hack into the mechanics solver and set up the solution to be true solution
            // because the solution is quite far from the initial guess.
            for (unsigned i=0; i<mechanics_mesh.GetNumNodes(); i++)
            {
                double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
                double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
                problem.mpMechanicsSolver->rGetCurrentSolution()[3*i]   = 0.8*X - X;
                problem.mpMechanicsSolver->rGetCurrentSolution()[3*i+1] = Y/0.8 - Y;
            }

            problem.Solve();

            // count the number of nodes that were stimulated at the last timestep
            Hdf5DataReader reader(dir+"/electrics","voltage");
            unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
            Vec voltage = PetscTools::CreateVec(electrics_mesh.GetNumNodes());
            reader.GetVariableOverNodes(voltage, "V", num_timesteps-1);

            ReplicatableVector voltage_repl(voltage);

            num_stimulated_nodes[sim] = 0;
            for (unsigned i=0; i<voltage_repl.GetSize(); i++)
            {
                if (voltage_repl[i]>0.0)
                {
                    num_stimulated_nodes[sim]++;
                }
            }
            PetscTools::Destroy(voltage);
        }

        // check the number of stimulated nodes is greater in the second case,
        // where conductivity is affected by the deformation
        TS_ASSERT_EQUALS(num_stimulated_nodes[0], 90u);
        TS_ASSERT_EQUALS(num_stimulated_nodes[1], 108u);

    }


    void TestCardiacElectroMechanicsHeterogeneousMaterialLaws()
    {
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000);

        TetrahedralMesh<2,2> electrics_mesh;
        electrics_mesh.ConstructRegularSlabMesh(0.02/*stepsize*/, 0.1/*length*/, 0.1/*width*/);

        QuadraticMesh<2> mechanics_mesh;
        mechanics_mesh.ConstructRegularSlabMesh(0.02, 0.1, 0.1 /*as above with a different stepsize*/);

        std::vector<unsigned> fixed_nodes
            = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.0);


        HeartConfig::Instance()->SetSimulationDuration(20.0);

        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(KERCHOFFS2003,0.01);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(1.0);

        // Create two laws to have stiff and soft tissue
        std::vector<AbstractMaterialLaw<2>*> laws;
        CompressibleMooneyRivlinMaterialLaw<2> stiff_law(1.0/*effects stiffness*/, 1.0/*affects amount of compressibility*/);
        CompressibleMooneyRivlinMaterialLaw<2> soft_law(1.0/5.0/*effects stiffness*/, 1.0/*affects amount of compressibility*/);
        for (TetrahedralMesh<2,2>::ElementIterator iter = mechanics_mesh.GetElementIteratorBegin();
             iter != mechanics_mesh.GetElementIteratorEnd();
             ++iter)
        {
            if (((iter)->CalculateCentroid()[1] >= 0.04)
                 && ((iter)->CalculateCentroid()[1] <= 0.06))
            {
                laws.push_back(&soft_law);
            }
            else
            {
                laws.push_back(&stiff_law);
            }
        }

        problem_defn.SetMaterialLaw(COMPRESSIBLE,laws);


        CardiacElectroMechanicsProblem<2,1> problem(COMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &electrics_mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "TestCardiacElectroMechanicsHeterogeneousMaterialLaws" /* output directory */);


        problem.Solve();

        // test by checking the length of the tissue against hardcoded value
        std::vector<c_vector<double,2> >& r_deformed_position = problem.rGetDeformedPosition();

        // check node 5 starts at (1,0)
        assert(fabs(mechanics_mesh.GetNode(5)->rGetLocation()[0] - 0.1)<1e-6);
        assert(fabs(mechanics_mesh.GetNode(5)->rGetLocation()[1]      )<1e-6);

        // Visualised solution to check heterogeneous stiffnesses are taken into account.
        // Here we just have a hardcoded test to check nothing has changed.
        // The effect of the weak region is small but noticeable, compared to a simulation
        // with stiff law everywhere - the weak region contracts a tiny bit more.
        TS_ASSERT_DELTA(r_deformed_position[5](0),  0.0917, 1e-4);
        TS_ASSERT_DELTA(r_deformed_position[5](1), -0.0002, 1e-4);
    }
};


#endif // TESTCARDIACELECTROMECHANICSFURTHERFUNCTIONALITY_HPP_
